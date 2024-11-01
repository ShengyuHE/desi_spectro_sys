import numpy as np
from astropy.table import Table
import fitsio
import sys
import os
import argparse
sys.path.append('/global/homes/s/shengyu/project_rc/main/desihub/')
import LSS.common_tools as common
from LSS.globals import main

parser = argparse.ArgumentParser()
parser.add_argument("--par", help="run different random number in parallel?",default='n')
parser.add_argument("--tracer", help="tracer type to be selected", choices=['ELG_LOPnotqso','LRG','QSO'])
parser.add_argument("--datadir", help="base directory for input",default=None)
parser.add_argument("--addcatas",help="apply catastrophics in the clean mock redshift 'Z', 'realistic' is the observed pattern, 'failures' is the upper limit 1%, slitless is the 5% assumed catastrophics",nargs='*',type=str,choices=['realistic','failures','slitless'],default=None)
parser.add_argument('--remove_zerror', help='the suffix of redshift column without the redsihft error', type=str, default=None)
args = parser.parse_args()

if args.remove_zerror == "None":
    args.remove_zerror = None

c = 299792 # speed of light in km/s
fn_data = args.tracer+'_{}GC_clustering.dat.fits'
fn_rand = args.tracer+'_{}GC_{}_clustering.ran.fits'
fn_nz   = args.tracer+'_{}GC_nz{}.txt'

if args.tracer == 'QSO':
    P0 = 6000
if args.tracer[:3] == 'LRG':
    P0 = 10000
if args.tracer[:3] == 'ELG':
    P0 = 4000
if args.tracer[:3] == 'BGS':
    P0 = 7000

if args.addcatas is not None:
    types = args.addcatas
if args.remove_zerror is not None:
    types = args.remove_zerror

## the catalogues are huge, so memory might be a problem
## this code changes NZ(z) to NZ(z_catas)
## TODO: NZ(z_catas) -> NZ_catas(z_catas)
def catasfun(par):
    catalog_fn,sourcedir = par[0],par[1]
    catalog=Table(fitsio.read(catalog_fn))
    if 'FKP_failures' in catalog.colnames:
        print('catastrophics FKP ready:',catalog_fn)
    else:
        import time
        T0=time.time()
        for catas_type in types:
            catalog[f'FKP_{catas_type}'] = catalog['WEIGHT_FKP']*1
            # change WEIGHT_FKP of catastrophics
            # change WEIGHT_FKP accordingly
            nz        = np.loadtxt(sourcedir+fn_nz.format(GC,f'_{catas_type}'))
            ## find the nz of the true redshift
            ind_rawNZ = np.argmin(abs(catalog['Z']-nz[:,0][:,np.newaxis]),axis=0)
            ## caluclate the completeness rescaling of nz for FKP weight
            norm = catalog['NX']/nz[ind_rawNZ,3]
            dv   = (catalog[f'Z_{catas_type}']-catalog['Z'])/(1+catalog['Z'])*c
            dz   = (catalog[f'Z_{catas_type}']-catalog['Z'])
            # note that FKP changes are fewer than z changes, this is because
            ## 1. 1% of the extra catas has dv<1000km/s for z=1.32 catas(negligible)
            ## 2. 99% of the extra catas was not from 0.8<Z_RAW<1.6=> need to find the norm of the closest NX for FKP calculation as they had norm==0
            tmp      = np.argsort(catalog,order=['RA', 'DEC'])
            catalog  = catalog[tmp]
            norm     = norm[tmp]
            dv       = dv[tmp]
            NX       = catalog['NX']*1
            norm[norm==0] = np.nan
            print('there are {} samples to find new FKP'.format(np.sum((dv!=0)&(np.isnan(norm)))))
            for ID in np.where((dv!=0)&(np.isnan(norm)))[0]:
                if (2<ID)&(ID<len(catalog)-2):
                    norm[ID] = np.nanmedian(norm[[ID-2,ID-1,ID+1,ID+2]])
                elif ID<2:
                    norm[ID] = np.nanmedian(norm[[ID+1,ID+2]])
                elif ID>len(catalog)-2:
                    norm[ID] = np.nanmedian(norm[[ID-2,ID-1]])

                ## update NX for norm ==0
                ind_newNZ = np.argmin(abs(catalog[f'Z_{catas_type}'][ID]-nz[:,0]))
                NX[ID] = norm[ID]*nz[ind_newNZ,3]
            ## select all catastrophics
            sel = dv!=0
            ind_newNZ = np.argmin(abs(catalog[f'Z_{catas_type}'][sel]-nz[:,0][:,np.newaxis]),axis=0)
            ## update NX and WEIGHT_FKP columns for all catastrophics
            NX[sel]  = norm[sel]*nz[ind_newNZ,3]
            catalog[f'FKP_{catas_type}'][sel] = 1 / (NX[sel]*P0+1)
            catalog[f'FKP_{catas_type}'][np.isnan(catalog[f'FKP_{catas_type}'])] = 1
            # find the nz of the catastrophics redshift
            print('implement catastrophophics took time:{:.2f}s'.format(time.time()-T0))
        catalog.write(catalog_fn,overwrite=True)
        print('catastrophics FKP corrected')

#catalog_fns = [datadir.format(mockid,mockid)+fn_rand.format(GC,ranid) for GC in ['N','S'] for mockid in range(Nrand) for ranid in range(18)]+[datadir.format(mockid,mockid)+fn_data.format(GC) for GC in ['N','S'] for mockid in range(Nrand) ]
#sourcedirs  = [datadir.format(mockid,mockid) for GC in ['N','S'] for mockid in range(Nrand) for ranid in range(18)]+[datadir.format(mockid,mockid) for GC in ['N','S'] for mockid in range(Nrand)]
Nrand = 0
catalog_fns = [args.datadir.format(mockid,mockid)+fn_rand.format(GC,ranid) for GC in ['N','S'] for mockid in range(Nrand,25) for ranid in range(18)]+[datadir.format(mockid,mockid)+fn_data.format(GC) for GC in ['N','S'] for mockid in range(Nrand,25) ]
sourcedirs  = [args.datadir.format(mockid,mockid) for GC in ['N','S'] for mockid in range(Nrand,25) for ranid in range(18)]+[datadir.format(mockid,mockid) for GC in ['N','S'] for mockid in range(Nrand,25)]

if args.par == 'n':
    for catalog_fn,sourcedir in zip(catalog_fns,sourcedirs):
        catasfun([catalog_fn,sourcedir])
if args.par == 'y':
    from multiprocessing import Pool
    with Pool(processes = 64) as pool:
        pool.map(catasfun, zip(catalog_fns,sourcedirs))

