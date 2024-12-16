import os
import sys
import argparse
import fitsio
import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d

from pypower import CatalogFFTPower,mpi, setup_logging
from pycorr import TwoPointCorrelationFunction, setup_logging
# from triumvirate.catalogue import ParticleCatalogue
# from triumvirate.threept import compute_bispec_in_gpp_box, compute_3pcf_in_gpp_box
# from triumvirate.parameters import fetch_paramset_template, ParameterSet

setup_logging()

sys.path.append('../')
from helper import REDSHIFT_VSMEAR, REDSHIFT_CUBICBOX, EDGES

# from mpi4py import MPI

# mpicomm = mpi.COMM_WORLD
# mpiroot = 0

boxsize = 2000.
kedges   = np.arange(0.,0.4001,0.001); ells = (0, 2)
smuedges  = (np.linspace(0., 200, 201), np.linspace(-1., 1., 201))
slogmuedges= (np.geomspace(0.01, 100., 49), np.linspace(-1., 1., 201))

def statistics_2pt(position, fn):
    # compute mps
    fn_mps = fn.format('xipoles')
    if not os.path.exists(fn_mps):
        result_mps = TwoPointCorrelationFunction('smu', smuedges, data_positions1=position, engine='corrfunc', 
                                                 boxsize=boxsize, los='z', position_type='xyz',
                                                 gpu=True, nthreads = 4)
                                                #  mpiroot=mpiroot, mpicomm=mpicomm)
        result_mps.save(fn_mps)
    else:
        result_mps = TwoPointCorrelationFunction.load(fn_mps)
    # compute mps log scales
    fn_mpslog = fn.format('mpslog')
    if not os.path.exists(fn_mpslog):
        result_mps = TwoPointCorrelationFunction('smu', slogmuedges, data_positions1=position, engine='corrfunc', 
                                                 boxsize=boxsize, los='z', position_type='xyz',
                                                 gpu=True, nthreads = 4)
                                                #  mpiroot=mpiroot, mpicomm=mpicomm)
        result_mps.save(fn_mpslog)
    else:
        result_mps = TwoPointCorrelationFunction.load(fn_mpslog)
    # compute pk
    fn_pk = fn.format('pkpoles')
    if not os.path.exists(fn_pk):
        result_pk = CatalogFFTPower(data_positions1=position, edges=kedges, ells=ells, interlacing=3, 
                                    boxsize=boxsize, nmesh=512, resampler='tsc',los='z', position_type='xyz',)
                                    # mpiroot=mpiroot, mpicomm=mpicomm)
        result_pk.save(fn_pk)
    else:
        result_pk = CatalogFFTPower.load(fn_pk)
    return result_mps,result_pk

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument("--nthreads", type = int, default = 4)
    # parser.add_argument("--datDir", help="base directory for void data catalogs", default=None)
    # parser.add_argument("--outputDir", help="base directory for void random catalogs", default=None)
    parser.add_argument("--mockid", type=str, default="0-24", help="Mock ID range or list")
    parser.add_argument("--tracers", help="tracer type to be selected", type = str, choices=['LRG','ELG','QSO'], default=['LRG','ELG','QSO'], nargs = '+')
    args = parser.parse_args()

    # Convert mockid string input to a list
    if '-' in args.mockid:
        start, end = map(int, args.mockid.split('-'))
        mockids = list(range(start, end + 1))
    else:
        mockids = list(map(int, args.mockid.split(',')))
    for tracer in args.tracers:
        for (zmin, zmax) in REDSHIFT_VSMEAR[tracer]:
            for mock_id in mockids:
                mock_id03 =  f"{mock_id:03}"
                basedir = f'/pscratch/sd/s/shengyu/galaxies/catalogs/cosmosim/AbacusHOD_mocks_v1/CubicBox/{tracer}/obs_z{zmin:.1f}-{zmax:.1f}/AbacusSummit_base_c000_ph{mock_id03}'
                catalog_fn = basedir+f'/catalog_rsd_xi2d_{tracer}_z{zmin:.1f}-{zmax:.1f}_velbias_B_s_mockcov.fits'
                catalog=Table(fitsio.read(catalog_fn))
                for sysmodel in ['standard', 'dv-obs']:
                    if sysmodel == 'standard':
                        positions = [catalog['x'],catalog['y'], catalog['z']]
                        fn = basedir+f'/mpspk/{{}}_{tracer}_z{zmin:.1f}-{zmax:.1f}_{sysmodel}.npy'
                        statistics_2pt(positions, fn)

                    if sysmodel == 'dv-obs':
                        positions = [catalog['x'],catalog['y'], catalog['z_dv']]
                        fn = basedir+f'/mpspk/{{}}_{tracer}_z{zmin:.1f}-{zmax:.1f}_{sysmodel}.npy'
                        statistics_2pt(positions, fn)