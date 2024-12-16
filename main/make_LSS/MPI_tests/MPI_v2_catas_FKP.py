import os
import glob
import time
import fitsio
import argparse
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from scipy.spatial import cKDTree

c = 299792 # speed of light in km/s
P0_values = {'ELG_LOPnotqso': 4000, 'LRG': 10000, 'QSO': 6000, 'BGS': 7000}
NRAN_values = {'ELG_LOPnotqso':10, 'LRG':8, 'QSO':4}

def NXnorm(catalog_fn, nz_fn):
    #Calculate the NXnorm for the catalog, norm = NX/n(z)
    catalog = Table(fitsio.read(catalog_fn)) 
    nz      = np.loadtxt(nz_fn)
    tree    = cKDTree(nz[:, 0].reshape(-1, 1))
    _, ind_rawNZ = tree.query(catalog['Z'].reshape(-1, 1))
    norm         = catalog['NX']/nz[ind_rawNZ,3]
    return norm

def FKPupdate(catalog_fn, nz_fn, catas_type, norm, mpicomm, mpiroot=0):
    T0 = time.time()
    rank = mpicomm.rank
    size = mpicomm.size
    print('rank:', rank)
    print('size:', size)

    n_galaxies = fitsio.FITS(catalog_fn)[1].get_nrows()

    # Check for Z_{catas_type} existence
    if rank == mpiroot:
        catalog_columns = fitsio.FITS(catalog_fn)[1].get_colnames()
        if f'Z_{catas_type}' not in catalog_columns:
            raise ValueError(f"Invalid Zcatas type: '{catas_type}'.")
    else:
        catalog_columns = None
    catalog_columns = mpicomm.bcast(catalog_columns, root=mpiroot)

    # Calculate the indices for each rank
    # counts = n_galaxies // size
    # remainder = n_galaxies % size
    # start = rank * counts + min(rank, remainder)
    # end = start + counts + (1 if rank < remainder else 0)

    indices = np.array_split(np.arange(n_galaxies), size)
    local_indices = indices[rank]

    # Each rank reads its own chunk of the catalog
    catalog = Table(fitsio.read(catalog_fn, rows=local_indices))
    print('indices', local_indices)
    print('rank:', rank, len(catalog))
    norm_chunk = norm[local_indices]

    # Add 'ROW_INDEX' to keep track of original order
    # catalog['ROW_INDEX'] = np.arange(start, end)
    catalog[f'FKP_{catas_type}'] = catalog['WEIGHT_FKP'].copy()

    # Load the nz_catas and create the cKDTree
    nz = np.loadtxt(nz_fn)
    tree_nz = cKDTree(nz[:, 0].reshape(-1, 1))

    # Calculate dv and dz for FKP weight
    dv = (catalog[f'Z_{catas_type}'] - catalog['Z']) / (1 + catalog['Z']) * c
    dz = catalog[f'Z_{catas_type}'] - catalog['Z']

    # Sort local_catalog if necessary
    tmp = np.argsort(catalog['RA'], kind='mergesort')
    catalog = catalog[tmp]
    norm_chunk = norm_chunk[tmp]
    dv = dv[tmp]

    NX = catalog['NX'].copy()
    norm_chunk[norm_chunk == 0] = np.nan

    # Handle NaN values in norm_chunk
    for idx in np.where((dv != 0) & (np.isnan(norm_chunk)))[0]:
        if 2 < idx < len(catalog) - 2:
            norm_chunk[idx] = np.nanmedian(norm_chunk[[idx - 2, idx - 1, idx + 1, idx + 2]])
        elif idx <= 2:
            norm_chunk[idx] = np.nanmedian(norm_chunk[[idx + 1, idx + 2]])
        elif idx >= len(catalog) - 2:
            norm_chunk[idx] = np.nanmedian(norm_chunk[[idx - 2, idx - 1]])
        # Update NX for norm == 0
        ind_newNZ = np.argmin(abs(catalog[f'Z_{catas_type}'][idx] - nz[:, 0]))
        NX[idx] = norm_chunk[idx] * nz[ind_newNZ, 3]

    # Update NX and FKP weights
    sel = dv != 0
    _, ind_newNZ = tree_nz.query(catalog[f'Z_{catas_type}'][sel].reshape(-1, 1))
    NX[sel] = norm_chunk[sel] * nz[ind_newNZ, 3]
    catalog[f'FKP_{catas_type}'][sel] = 1 / (NX[sel] * P0 + 1)
    # Replace NaNs in FKP weights with 1
    nan_indices = np.isnan(catalog[f'FKP_{catas_type}'])
    catalog[f'FKP_{catas_type}'][nan_indices] = 1

    mpicomm.Barrier()
    # Gather the updated local_catalogs back to the root rank
    all_catalogs = mpicomm.gather(catalog, root=mpiroot)

    if mpicomm.rank == mpiroot:
        # Concatenate all updated data
        final_catalog = vstack(all_catalogs)
        # final_catalog.sort('ROW_INDEX')
        # final_catalog.remove_column('ROW_INDEX')
        # Write the updated catalog to file
        final_catalog.write(catalog_fn, overwrite=True)
        print(f"[FKP weight]: Processed in {time.time() - T0:.2f}s", flush=True)
    
    return 0


def FKPupdate_old(catalog_fn, nz_fn, catas_type, norm, mpicomm, mpiroot=0):
    T0 = time.time()
    catalog=Table(fitsio.read(catalog_fn))
    # Check for Z_{catas_type} existence
    if f'Z_{catas_type}' not in catalog.colnames:
        raise ValueError(f"Invalid Zcatas type: '{catas_type}'.")

    # Add 'ROW_INDEX' to keep track of original order
    catalog['ROW_INDEX'] = np.arange(len(catalog))
    
    # Distribute galaxies among ranks
    n_galaxies = len(catalog)
    indices = np.arange(n_galaxies)
    indices_per_rank = np.array_split(indices, mpicomm.size)
    local_indices = indices_per_rank[mpicomm.rank]

    # Each rank processes its subset of galaxies
    local_catalog = catalog[local_indices]
    local_norm = norm[local_indices]

    local_catalog[f'FKP_{catas_type}'] = local_catalog['WEIGHT_FKP'].copy()

    # Load the nz_catas and create the cKDTree
    nz = np.loadtxt(nz_fn)
    tree = cKDTree(nz[:, 0].reshape(-1, 1))

    # caluclate the completeness rescaling of nz for FKP weight
    dv = (local_catalog[f'Z_{catas_type}'] - local_catalog['Z']) / (1 + local_catalog['Z']) * c
    dz = local_catalog[f'Z_{catas_type}'] - local_catalog['Z']
    tmp      = np.argsort(local_catalog['RA'], kind='mergesort')
    # tmp      = np.argsort(catalog,order=['RA', 'DEC'])
    local_catalog  = local_catalog[tmp]
    local_norm     = local_norm[tmp]
    dv       = dv[tmp]
    NX       = local_catalog['NX'].copy()
    local_norm[local_norm==0] = np.nan
    # print('there are {} samples to find new FKP'.format(np.sum((dv!=0)&(np.isnan(norm)))), flush=True)
    for ID in np.where((dv!=0)&(np.isnan(local_norm)))[0]:
        if (2<ID)&(ID<len(local_catalog)-2):
            local_norm[ID] = np.nanmedian(local_norm[[ID-2,ID-1,ID+1,ID+2]])
        elif ID<2:
            local_norm[ID] = np.nanmedian(local_norm[[ID+1,ID+2]])
        elif ID>len(local_catalog)-2:
            local_norm[ID] = np.nanmedian(local_norm[[ID-2,ID-1]])
        # update NX for norm ==0
        ind_newNZ = np.argmin(abs(local_catalog[f'Z_{catas_type}'][ID]-nz[:,0]))
        NX[ID] = local_norm[ID]*nz[ind_newNZ,3]
    # update NX and WEIGHT_FKP columns for all catastrophics
    sel = dv != 0
    _, ind_newNZ = tree.query(local_catalog[f'Z_{catas_type}'][sel].reshape(-1, 1))
    NX[sel] = local_norm[sel] * nz[ind_newNZ, 3]
    local_catalog[f'FKP_{catas_type}'][sel] = 1 / (NX[sel] * P0 + 1)
    local_catalog[f'FKP_{catas_type}'][np.isnan(local_catalog[f'FKP_{catas_type}'])] = 1

    # Gather the updated local_catalogs back to the root rank
    mpicomm.Barrier()
    all_local_catalogs = mpicomm.gather(local_catalog, root=mpiroot)

    if mpicomm.rank == mpiroot:
        updated_catalog = vstack(all_local_catalogs)
        updated_catalog.sort('ROW_INDEX')
        updated_catalog.remove_column('ROW_INDEX')
        # Write the updated catalog to file
        updated_catalog.write(catalog_fn, overwrite=True)
        print('[FKP weight]: implement {} catastrophophics took time: {:.2f}s'.format(catas_type, time.time()-T0), flush=True)
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--targDir", help="base directory for input",default=None)
    parser.add_argument("--tracer", help="tracer type to be selected", choices=['ELG_LOPnotqso','LRG','QSO'],default=None)
    parser.add_argument("--region", help="region to be selected", choices=['NGC','SGC','NGC SGC'], default ='NGC SGC')
    parser.add_argument("--mock_type", help="mocks type", choices=['dat','ran','all'], default='all')
    parser.add_argument("--addcatas",help="apply catastrophics in the clean mock redshift 'Z', 'realistic' is the observed pattern, 'failures' is the upper limit 1%, slitless is the 5% assumed catastrophics",nargs='*',type=str,choices=['realistic','failures','slitless'], default=['realistic','failures'])
    parser.add_argument('--remove_zerror', help='the suffix of redshift column without the redsihft error', type=str, default=None)
    parser.add_argument("--parallel", help="MPI parallelization",choices=['y','n'], default='y')
    args = parser.parse_args()

    if args.remove_zerror == "None":
        args.remove_zerror = None

    T1 = time.time()

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

    if args.parallel == 'y':
        from mpi4py import MPI
        mpicomm = MPI.COMM_WORLD
        mpiroot = 0
        rank = mpicomm.Get_rank()
        size = mpicomm.Get_size()
        if rank == mpiroot:
            print("[FKP weight] Enable MPI parallization")
    else:
        rank = 0
        size = 1

    for region in regions:
        nz_fiducial_fn = filedir+f'/{tracer}_{region}_nz.txt'
        for mock_type in mock_types:
            if mock_type == 'dat':
                catalog_fn = filedir+ f'/{tracer}_{region}_clustering.{mock_type}.fits'
                norm = NXnorm(catalog_fn, nz_fiducial_fn)
                for catas_type in types:
                    nz_fn = filedir+f'/{tracer}_{region}_nz_{catas_type}.txt'
                    FKPupdate(catalog_fn, nz_fn, catas_type, norm, mpicomm, mpiroot)
                    if args.parallel == 'y':
                        # MPI initialization...
                        mpicomm.Barrier()
                    if rank == mpiroot:
                        print(f'[FKP weight]: {tracer}_{region} dat mocks FKP_{catas_type} weight updated', flush =True)
            elif mock_type == 'ran':
                for num in range(4):
                    catalog_fn = filedir+ f'/{tracer}_{region}_{num}_clustering.{mock_type}.fits'
                    norm = NXnorm(catalog_fn, nz_fiducial_fn)
                    for catas_type in types:
                        nz_fn = filedir+f'/{tracer}_{region}_nz_{catas_type}.txt'
                        FKPupdate(catalog_fn, nz_fn, catas_type, norm, mpicomm, mpiroot)
                        if args.parallel == 'y':
                        # MPI initialization...
                            mpicomm.Barrier()
                        if rank == mpiroot:
                            print(f'[FKP weight]: {tracer}_{region} random {num} mocks FKP_{catas_type} weight updated', flush =True)

    if rank == 0:
        print('[FKP weight]: Implement all catastrophophics takes time: {:.2f}s'.format(time.time()-T1), flush=True)