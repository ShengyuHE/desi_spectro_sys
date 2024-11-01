# check the MPI
from mpi4py import MPI
mpicomm = MPI.COMM_WORLD
mpiroot = 0

import os
import numpy as np
import sys
import argparse
import fitsio
from astropy.table import Table

c = 299792 # speed of light in km/s
P0_values = {'ELG_LOPnotqso': 4000, 'LRG': 10000, 'QSO': 6000, 'BGS': 7000}
NRAN_values = {'ELG_LOPnotqso':10, 'LRG':8, 'QSO':2}

def FKPupdate(catalog_fn, nz_fn, catas_type):
    # change WEIGHT_FKP of catastrophics
    import time 
    T0=time.time()
    catalog=Table(fitsio.read(catalog_fn))
    try:
        catalog[f'Z_{catas_type}']
    except KeyError:
        raise ValueError(f"Invalid Zcatas type: '{catas_type}'.")
    
    # if f'FKP_{catas_type}' in catalog.colnames:
    #     print(f'{catas_type} catastrophics FKP ready:', catalog_fn)
    if 1:
        catalog[f'FKP_{catas_type}'] = catalog['WEIGHT_FKP']*1
        # find the nz of the true redshift
        nz = np.loadtxt(nz_fn)
        ind_rawNZ = np.argmin(abs(catalog['Z']-nz[:,0][:,np.newaxis]),axis=0)
        # caluclate the completeness rescaling of nz for FKP weight
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
            # update NX for norm ==0
            ind_newNZ = np.argmin(abs(catalog[f'Z_{catas_type}'][ID]-nz[:,0]))
            NX[ID] = norm[ID]*nz[ind_newNZ,3]
        #select all catastrophics
        sel = dv!=0
        ind_newNZ = np.argmin(abs(catalog[f'Z_{catas_type}'][sel]-nz[:,0][:,np.newaxis]),axis=0)
        # update NX and WEIGHT_FKP columns for all catastrophics
        NX[sel]  = norm[sel]*nz[ind_newNZ,3]
        catalog[f'FKP_{catas_type}'][sel] = 1 / (NX[sel]*P0+1)
        catalog[f'FKP_{catas_type}'][np.isnan(catalog[f'FKP_{catas_type}'])] = 1
        # find the nz of the catastrophics redshift
        print('implement {} catastrophophics took time: {:.2f}s'.format(catas_type, time.time()-T0))
        catalog.write(catalog_fn, overwrite=True)
        print(f'{catas_type} catastrophics FKP corrected')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--targDir", help="base directory for input",default=None)
    parser.add_argument("--tracer", help="tracer type to be selected", choices=['ELG_LOPnotqso','LRG','QSO'],default=None)
    parser.add_argument("--region", help="region to be selected", choices=['NGC','SGC','NGC SGC'], default ='NGC SGC')
    parser.add_argument("--mock_type", help="mocks type", choices=['dat','ran','all'], default='all')
    parser.add_argument("--addcatas",help="apply catastrophics in the clean mock redshift 'Z', 'realistic' is the observed pattern, 'failures' is the upper limit 1%, slitless is the 5% assumed catastrophics",nargs='*',type=str,choices=['realistic','failures','slitless'], default=['realistic','failures'])
    parser.add_argument('--remove_zerror', help='the suffix of redshift column without the redsihft error', type=str, default=None)
    parser.add_argument("--parallel", help="run different random number in parallel?",default='n')
    args = parser.parse_args()

    if args.remove_zerror == "None":
        args.remove_zerror = None

    filedir = args.targDir
    tracer = args.tracer
    mock_types = ['dat','ran'] if args.mock_type == 'all' else [args.mock_type]
    regions = ['NGC', 'SGC'] if args.region == 'NGC SGC' else [args.region]

    P0 = P0_values.get(tracer, None)
    NRAN = NRAN_values.get(tracer, None)

    if args.addcatas is not None:
        types = args.addcatas
    if args.remove_zerror is not None:
        types = args.remove_zerror

    for region in regions:
        for mock_type in mock_types:
            if mock_type == 'dat':
                catalog_fn = filedir+ f'/{tracer}_{region}_clustering.{mock_type}.fits'
                for catas_type in types:
                    print(f'\nupdate FKP_{catas_type} weight for {tracer} data mocks')
                    nz_fn = filedir+f'/{tracer}_{region}_nz_{catas_type}.txt'
                    FKPupdate(catalog_fn, nz_fn, catas_type)
            elif mock_type == 'ran':
                for num in range(NRAN):
                    catalog_fn = filedir+ f'/{tracer}_{region}_{num}_clustering.{mock_type}.fits'
                    for catas_type in types:
                        print(f'\nupdate FKP_{catas_type} weight for {tracer} random {num} mocks')
                        nz_fn = filedir+f'/{tracer}_{region}_nz_{catas_type}.txt'
                        FKPupdate(catalog_fn, nz_fn, catas_type)
