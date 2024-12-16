import os
import fitsio
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d

from desitarget.targetmask import desi_mask, bgs_mask
from desitarget.targetmask import zwarn_mask as zmtl_zwarn_mask

def vsmear_modelling(tracer,zmin,zmax,dvfn='./'):
    """
    vsmear_modelling function:

    This function provide the probability distribution function (PDF) and the cumulative probability function (CDF) \
    of the redshift difference Delta_v from repeated observations of DESI Y3 tracers (e.g., LRG, ELG, QSO, BGS) \
    at a given redshift range [zmin, zmax]. This is to prepare for applying them to the mocks.

    Parameters:
    ----------
    - tracer : str
        The tracer type. Valid values are 'LRG', 'ELG', 'QSO', and 'BGS' or their subsamples.
    
    - zmin : float
        The minimum redshift of the range to be considered for selection.
    
    - zmax : float
        The maximum redshift of the range to be considered for selection.
    
    - dvfn : str, optional
        The directory path to store the output files (default is './').

    Returns:
    --------
    This function does not return any values but saves output to files in the specified directory (dvfn).
        - The PDF and CDF of DV for the specified redshift range are saved as `.npz` files.


    Example usage:
    --------------
    vsmear_modelling('LRG', 0.4, 0.6, dvfn='./output')
    """
    # Validate the tracer input
    if tracer[:3] in ['LRG', 'ELG', 'QSO']:
        repeatdir= '/global/cfs/projectdirs/desi/users/jiaxiyu/repeated_observations/EDR_vs_Y3/LSS-scripts_repeats/main-repeats-kibo-dark-pairs.fits'
        mask, mask_key = desi_mask, "DESI_TARGET"
        effkey, effmin, effmax, effxlim = (
            "TSNR2_LRG",
            0.85 * 1000,
            1.5 * 1000,
            (500, 1500),
        )
    elif tracer[:3] == 'BGS':
        repeatdir= '/global/cfs/projectdirs/desi/users/jiaxiyu/repeated_observations/EDR_vs_Y3/LSS-scripts_repeats/main-repeats-kibo-bright-pairs.fits'
        mask, mask_key = bgs_mask, "BGS_TARGET"
        effkey, effmin, effmax, effxlim = "TSNR2_BGS", 0.85 * 180, 1.5 * 180, (0, 500)
    else:
        raise ValueError(f"Invalid tracer: {tracer[:3]}. Must be in/be a subsample of ['BGS', 'LRG', 'ELG', 'QSO'].")
    # Create the output directory if not existed
    if not os.path.exists(dvfn):
        os.system(f'mkdir -p {dvfn}')

    # Read the observed pairs 
    d     = Table.read(repeatdir)

    # Set parameters for the selection and calculation process
    catasmin, catasmax, catasbin = -3, 6, 0.2
    goodkey = f"GOOD_{tracer[:3]}"
    if tracer == "QSO":
        catasmin = -2
    
    # Redrock redshift selections:
    if tracer != 'QSO':
        # efftime_spec calculation for selections
        snr2time = d.meta["{}SNR2T".format(effkey.split("_")[1])]
        efftime0s = snr2time * d["{}_0".format(effkey)]
        efftime1s = snr2time * d["{}_1".format(effkey)]
        # zmtl_zwarn_mask nodata + bad selections
        nodata0 = (d["ZMTL_ZWARN_0"] & zmtl_zwarn_mask["NODATA"]) > 0
        nodata1 = (d["ZMTL_ZWARN_1"] & zmtl_zwarn_mask["NODATA"]) > 0
        badqa0 = (d["ZMTL_ZWARN_0"] & zmtl_zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA")) > 0
        badqa1 = (d["ZMTL_ZWARN_1"] & zmtl_zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA")) > 0
        # Apply the selection criteria to clean the data
        sel = (d[mask_key] & mask[tracer]) > 0
        sel &= (d["COADD_FIBERSTATUS_0"] == 0) & (d["COADD_FIBERSTATUS_1"] == 0)
        sel &= (~nodata0) & (~nodata1)
        sel &= (~badqa0) & (~badqa1)
        sel &= (efftime0s > effmin) & (efftime1s > effmin)
        sel &= (efftime0s < effmax) & (efftime1s < effmax)
        sel &= (d["{}_0".format(goodkey)]) & (d["{}_1".format(goodkey)])
        sel &= (d['SURVEY_0']=='main')&(d['SURVEY_1']=='main')
        if tracer[:3] == 'ELG':
            sel &= (~d["GOOD_QSO_0"]) & (~d["GOOD_QSO_1"])
    else:
        qsofn    = repeatdir[:-5]+'_QSO'+repeatdir[-5:]
        d        = Table(fitsio.read(qsofn))
        d        = d[~np.isnan(d['DV'])]
        sel      = np.full(len(d),True)
        
    # cut on redshift range 
    selz = ((zmin<d["Z_0"])&(d["Z_0"]<zmax))|((zmin<d["Z_1"])&(d["Z_1"]<zmax))
    dv_final = np.log10(abs(d['DV'][sel&selz]))  
    
    # provide the Delta_velocity distributions
    dens,bins = np.histogram(dv_final,bins=np.arange(catasmin,catasmax,catasbin),density=True)

    ## keep none-zero elements
    sel_clean = dens>0
    vmid      = (bins[1:]+bins[:-1])/2
    vmid      = vmid[sel_clean]
    ## save the observed PDF and CDF
    cdffn_data= f'{dvfn}/{tracer[:3]}_z{zmin:.1f}-{zmax:.1f}_CDF'
    if not os.path.exists(cdffn_data+'.npz'):       
        vbin_fine = 0.005
        dens_fine,bins_fine=np.histogram(dv_final,bins=np.arange(catasmin,catasmax,vbin_fine),density=True)
        cdf_data     = np.cumsum(dens_fine) * vbin_fine 
        np.savez(cdffn_data, vbin=(bins_fine[1:]+bins_fine[:-1])/2, pdf=dens_fine, cdf=cdf_data)
                    
    # interpolation for the observed Delta_v distribution
    vnewbin = 0.005
    vnew    = np.arange(vmid[0]-catasbin/2,vmid[-1]+catasbin/2+0.01,vnewbin)
    vnewmid = (vnew[1:]+vnew[:-1])/2
    kernel  = 0.3
    ## compute the modelled PDF and CDF and save them
    cdffn= f'{dvfn}/{tracer[:3]}_z{zmin:.1f}-{zmax:.1f}_kernel{kernel}_CDF'
    if not os.path.exists(cdffn+'.npz'):       
        interp  = gaussian_kde(dv_final, bw_method=kernel)
        pdf     = interp(vnewmid)
        cdf     = np.cumsum(pdf) * vnewbin  
        np.savez(cdffn, vbin=vnewmid, pdf=pdf, cdf=cdf)
        

def vsmear(tracer,zmin,zmax,Ngal,zmode='LSS',dvmode='obs',seed=42):
    """
    vsmear function:

    This function generate an array to smear the velocity along the line of sight of a tracer (e.g., LRG, ELG, QSO, BGS) at a given redshift range [zmin,zmax]. \
    This is to model the effect of redshift uncertainties and catastrophics on realistic galayx mocks. 
    The distribution of the array is based on cumulative distribution functions (CDF) derived from pre-existing data generated by Y3_redshift_systematics.py

    Parameters:
    -----------
    - tracer : str
        The tracer type. Valid values are 'LRG', 'ELG', 'QSO', and 'BGS' or their subsamples.
    
    - zmin, zmax : float, float
        The redshift lower and upper limits. 
      - For 'BGS', valid [zmin,zmax] are [0.1, 0.4] for zmode=='LSS', or are chosen from np.arange(0.1,0.41,0.1) for zmode=='fine'.
      - For 'LRG', valid [zmin,zmax] are [(0.4, 0.6), (0.6, 0.8), (0.8, 1.1)] for zmode=='LSS', or are chosen from np.arange(0.4,1.11,0.1) for zmode=='fine'.
      - For 'ELG', valid [zmin,zmax] are [(0.8, 1.1), (1.1, 1.3), (1.3, 1.6)] for zmode=='LSS' or are chosen from np.arange(0.8,1.61,0.1) for zmode=='fine'.
      - For 'QSO', valid [zmin,zmax] are [(0.8, 1.1), (1.1, 1.4), (1.4, 1.7), (1.7, 2.1), (1.1, 1.3), (1.3, 1.7)] for zmode=='LSS' or are chosen from np.arange(0.8,2.11,0.1) for zmode=='fine'.

    - Ngal : int
        The number of galaxies of your mock. 

    - zmode : str, optional (default 'LSS')
        The redshift mode to use for selecting redshift ranges. Options are 'LSS' (the zbinning DESI LSS analysis used) and 'fine' (dz=0.1 binning).

    - dvmode : str, optional (default 'model')
        The mode for selecting the velocity dispersion distribution. Options are 'obs' (observational) and 'model' (theoretical model).

    Returns:
    --------
    - dv : numpy.ndarray
        An array of velocity dispersions (Delta_v) generated based on the specified tracer and redshift range.

    Example usage:
    --------------
    dv = vsmear('LRG', 0.4, 0.6, Ngal=100000, zmode='LSS', dvmode='obs')
    """
    np.random.seed(seed)
    if tracer[:3] not in ['BGS', 'LRG', 'ELG', 'QSO']:
        raise ValueError(f"Invalid tracer: {tracer[:3]}. Must be one of ['BGS', 'LRG', 'ELG', 'QSO'].")
    
    # Define redshift ranges in LSS mode and in fine-binning mode
    #repeatdir = '/global/homes/j/jiaxi/DESI_spectroscopic_systematics/Y3'
    repeatdir = '/global/cfs/projectdirs/desi/users/jiaxiyu/repeated_observations/EDR_vs_Y3/LSS-scripts_repeats'

    # Extract tracer key and range info 
    tracer_key = tracer[:3]

    # check if the CDF exists
    if   dvmode == 'obs':
        fn_cdf = f'{repeatdir}/{tracer_key}_z{zmin:.1f}-{zmax:.1f}_CDF.npz'
    elif dvmode == 'model':
        fn_cdf = f'{repeatdir}/{tracer_key}_z{zmin:.1f}-{zmax:.1f}_kernel0.3_CDF.npz'  
    print(f'load {fn_cdf} to add redshift uncertainties and catastrophics')

    if not os.path.exists(fn_cdf):
        raise ValueError(f"No prepared file: {fn_cdf}. \
        Please use function Y3_redshift_systematics.vsmear_modelling to prepare them!")
    
    # Load the saved data
    data     = np.load(fn_cdf)
    xgrid_raw= data["vbin"]
    cdf_raw  = data["cdf"]

    # inverse the CDF with interpolation 
    ## Remove duplicate values in the CDF
    cdf, ind = np.unique(cdf_raw, return_index=True)
    x_grid   = xgrid_raw[ind]
    inv_cdf = interp1d(cdf/cdf[-1], x_grid, bounds_error=False, fill_value=(x_grid[0], x_grid[-1]), kind='cubic')

    # Generate the exponents of Delta_v using the inverse CDF
    random_uniform_samples = np.random.uniform(0, 1, int(Ngal/2))  # Uniform samples [0, 1]
    exponent = inv_cdf(random_uniform_samples)  # Transform to match the redshift error distributions
    dv       = np.append(10**exponent,-10**exponent) # get the actual Delta_v
    np.random.shuffle(dv) # shuffle dv to get random positive and negative dv
    if Ngal%2 ==1:
        dv = np.append(np.zeros(1),dv)
    return dv
