import os, sys
import fitsio
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table
import matplotlib.pyplot as plt
import random
import h5py

from Y3_dv import vsmear

from pypower import CatalogFFTPower,mpi
from pycorr import TwoPointCorrelationFunction

from triumvirate.catalogue import ParticleCatalogue
from triumvirate.threept import compute_bispec_in_gpp_box, compute_3pcf_in_gpp_box
from triumvirate.parameters import fetch_paramset_template, ParameterSet

pt_type = sys.argv[1].split(',') #2pt,3pt

mockname = f'{os.environ["SCRATCH"]}/SHAM/catalog/UNIT_SHAM/DESI/SHAM_DESI_LRG_z0.8z1.1_at0.51450_seed15000.hdf5'
f = h5py.File(mockname,"r")
DATA = np.zeros((len(f["gal"]['X'][:]),3))
i = 0
for key in f['gal'].keys():
    if (key == 'X')|(key == 'Y')|(key == 'Z_RSD_true'):
        DATA[:,i] = f[f'gal/{key}'][...]
        i = i+1

def statistics_2pt(catalogue,source,output):
    # compute mps
    fn_mps = source+output.format('mps')+'.npy'
    if not os.path.exists(fn_mps):
        result_mps = TwoPointCorrelationFunction('smu', smuedges, data_positions1=catalogue,engine='corrfunc', boxsize=boxsize, los='z',position_type='xyz', mpicomm=mpicomm, mpiroot=mpiroot)
        result_mps.save(fn_mps)
    else:
        result_mps = TwoPointCorrelationFunction.load(fn_mps)
    
    # compute mps log scales
    fn_mpslog = source+output.format('mpslog')+'.npy'
    if not os.path.exists(fn_mpslog):
        result_mps = TwoPointCorrelationFunction('smu', slogmuedges, data_positions1=catalogue,engine='corrfunc', boxsize=boxsize, los='z',position_type='xyz', mpicomm=mpicomm, mpiroot=mpiroot)
        result_mps.save(fn_mpslog)
    else:
        result_mps = TwoPointCorrelationFunction.load(fn_mpslog)

    # compute pk
    fn_pk = source+output.format('pk')+'.npy'
    if not os.path.exists(fn_pk):
        result_pk = CatalogFFTPower(data_positions1=catalogue,edges=kedges, ells=ells, interlacing=3, boxsize=boxsize, nmesh=512, resampler='tsc',los='z', position_type='xyz', mpiroot=mpiroot,mpicomm=mpicomm)
        result_pk.save(fn_pk)
    else:
        result_pk = CatalogFFTPower.load(fn_pk)
    #return result_mps,result_pk

def statistics_3pt(catalogue,source,output):
    CATALOGUE = ParticleCatalogue(catalogue[0],catalogue[1],catalogue[2], nz=len(DATA)/boxsize**3)
    
    # compute bispectrum
    fn_bispec = source+output.format('bk000_diag')+'.npz'
    if not os.path.exists(fn_bispec):
        bispec_param.update(directories={'measurements':source},tags={'output':output[2:]})
        result_bispec = compute_bispec_in_gpp_box(CATALOGUE,paramset=bispec_param,save='.npz')

    # compute 3pt
    fn_3pt = source+output.format('zeta000_diag')+'.npz'
    if not os.path.exists(fn_3pt):
        threept_param.update(directories={'measurements':source},tags={'output':output[2:]})
        result_threept = compute_3pcf_in_gpp_box(CATALOGUE,paramset=threept_param,save='.npz')

outputdir= f'{os.environ["SCRATCH"]}/DESI_spectroscopic_systematics/dv-test/'
outputfn = '{}_{}_z{}-{}_vsmear_dv-{}'
a_t      = '0.49220'
z        = 1/float(a_t)-1
Om       = 0.3089
Ode      = 1-Om
H        = 100*np.sqrt(Om*(1+z)**3+Ode)
zmin     = 0.8
zmax     = 1.1
ells     = (0, 2)
boxsize  = 1000

if '2pt' in pt_type:
    mpicomm = mpi.COMM_WORLD
    mpiroot = 0
    if mpicomm.rank == mpiroot:
        print('recommended setting: srun -N 2 -n 16 -c 16')
    # Pk settings
    kedges   = np.arange(0.,0.4001,0.001)
    # mps settings
    smuedges  = (np.linspace(0., 200, 201), np.linspace(-1., 1., 201))
    slogmuedges= (np.geomspace(0.01, 100., 49), np.linspace(-1., 1., 201))
    statistics_2pt(DATA.T,outputdir,outputfn.format({},'standard','min','zmax','dvmode'))
if '3pt' in pt_type:
    os.environ['OMP_NUM_THREADS'] = '128'
    print('recommended setting: srun -N 1 -n 1 -c 128')
    # bispec settings
    bispec_dict = fetch_paramset_template('dict')
    bispec_dict.update({
        'degrees': {'ELL': 0, 'ell1': 0, 'ell2': 0},
        'boxsize': {'x': boxsize, 'y': boxsize, 'z': boxsize},
        'statistic_type': 'bispec',
        'ngrid': {'x': 512, 'y': 512, 'z': 512},
        'num_bins': 10,
        'range': [0.005, 0.105],
    })
    bispec_param = ParameterSet(param_dict=bispec_dict)

    # 3pt settings
    bispec_dict.update({
        'statistic_type': '3pcf',
        'range': [0,50],
    })
    threept_param = ParameterSet(param_dict=bispec_dict)

    statistics_3pt(DATA.T,outputdir,outputfn.format({},'standard','min','zmax','dvmode'))

for tracer in ['BGS','LRG','ELG']:
    for dvmode in ['model','obs']:
        for dzmode in ['dz-LSS', 'dz-0.1']:
            zmode = 'LSS' if dzmode == 'dz-LSS' else 'fine'
                # AR per-tracer setting
            if tracer[:3] == "BGS":
                if dzmode.find('0.1')!=-1:
                    zmins  = np.arange(0.1,0.31,0.1)
                    zmaxs   = np.arange(0.2,0.41,0.1)
                elif dzmode.find('LSS')!=-1:
                    zmins, zmaxs = [0.1], [0.4]
            elif tracer == 'LRG': 
                if dzmode.find('0.1')!=-1:
                    zmins  = np.arange(0.4,1.01,0.1)
                    zmaxs   = np.arange(0.5,1.11,0.1)
                elif dzmode.find('LSS')!=-1:
                    zmins  = [0.4,0.6,0.8]
                    zmaxs  = [0.6,0.8,1.1]
            elif tracer[:3] == 'ELG':
                if dzmode.find('0.1')!=-1:
                    zmins  = np.arange(0.8,1.51,0.1)
                    zmaxs   = np.arange(0.9,1.61,0.1)
                elif dzmode.find('LSS')!=-1:
                    zmins  = [0.8,1.1,1.3]
                    zmaxs  = [1.1,1.3,1.6]
            for zmin,zmax in zip(zmins,zmaxs):
                # apply smearing
                Z_RSD_vsmear = (DATA[:,2]+vsmear(tracer,zmin,zmax,len(DATA),zmode=zmode,dvmode=dvmode)*(1+z)/H)%boxsize
                if '2pt' in pt_type:
                    statistics_2pt(np.array([DATA[:,0],DATA[:,1],Z_RSD_vsmear]),outputdir,outputfn.format('{}',tracer,f'{zmin:.1f}',f'{zmax:.1f}',dvmode))
                if '3pt' in pt_type:
                    #statistics_3pt(np.array([DATA[:,0],DATA[:,1],Z_RSD_vsmear]),outputdir,outputfn.format('{}',tracer,f'{zmin:.1f}',f'{zmax:.1f}',dvmode))
                    statistics_3pt(np.array([DATA[:,0],DATA[:,1],Z_RSD_vsmear]),outputdir,outputfn.format('{}',tracer,int(zmin/0.1),int(zmax/0.1),dvmode))
