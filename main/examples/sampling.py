'''
## MCMC Sampling

This notebook performs a posterior analysis of cosmological parameters using MCMC sampling:
1. Load power spectrum from ELG_LOPnotqso and covariance matrices from EZmocks.
2. Run MCMC chains to explore the posterior distributions.
3. Visualize the posterior distributions and compare them to true values.

The focus is on examining the resulting posterior distributions and evaluating their consistency with expected cosmological values.
'''

import os
import glob
import numpy as np
import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from getdist import plots

# from mockfactory import Catalog
from cosmoprimo.fiducial import DESI
from desilike.theories.galaxy_clustering import FOLPSTracerPowerSpectrumMultipoles, FOLPSRCTracerPowerSpectrumMultipoles
from desilike.observables.galaxy_clustering import TracerPowerSpectrumMultipolesObservable
from desilike.emulators import EmulatedCalculator, Emulator, TaylorEmulatorEngine
from desilike.likelihoods import ObservablesGaussianLikelihood
from desilike.samplers.emcee import EmceeSampler
from desilike.samples import plotting, Chain
from desilike import setup_logging
setup_logging()  # for logging messages

# vairables for pk calculation
redshift = 1.0
tracer = "ELG_LOPnotqso" 
ran_mock_num = "10" # ELG:10; LRG:8; QSO:4
region = "NGC"  # NGC or SGC
mocks_fn = os.path.join(os.environ['SCRATCH'], "mocks/")
pk_dir = os.path.join(os.environ['SCRATCH'], "data/pk/")
# pk_dir = os.path.join(os.environ['HOME'], "project_rc/main/data/pk/")
output_dir = '/global/homes/s/shengyu/project_rc/main/results'
cov_dir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/desipipe/v1/ffa/2pt/'

# load the pk file, the pk calculated from the pkrun.py
pk_fn = pk_dir+f'pkpoles_{tracer}_{region}_0.8_1.1_default_FKP_lin_thetacut0.05.npy'
wmatrix_fn = pk_dir+f'wmatrix_smooth_{tracer}_{region}_0.8_1.1_default_FKP_lin_thetacut0.05.npy'

# load the covariance matrix from the Ezmocks
cov_fns = []
cov_fn = f'pkpoles_ELG_LOP_{region}_z0.8-1.1_default_FKP_lin_nran10_cellsize6_boxsize9000.npy'
for i in range(100, 200):
    mock_dir = f'mock{i}/pk/'
    cov_fns.extend(glob.glob(os.path.join(cov_dir,  mock_dir, cov_fn), recursive=True))


def set_true_values(params, catalogue='fiducial'):
    # true values for the QUIJOTE simulation (not true for the Abacus)
    update_values = {
        'fiducial': {'h': 0.6711, 'omega_cdm': 0.1209, 'Omega_cdm': 0.2685, 'omega_b':0.02207,'logA': 3.0631, 'm_ncdm': 0.0, 'n_s':0.9624, 'w0_fld':-1.0, 'fc':0.01},
    }
    if catalogue in update_values:
        truth_values = update_values[catalogue]
    return [truth_values[param] for param in params if param in truth_values]


def plot_covergence_walk(fn, params):
    ndim            = len(params)
    chain = Chain.load(fn).remove_burnin()[::]
    chain_samples   = dict(zip(chain.basenames(), chain.data))
    samples         = np.array([chain_samples[p] for p in params])
    medians         = np.array(chain.median(params=params))
    true_values     = set_true_values(params)
    fig, ax = plt.subplots(ndim, sharex=True, figsize=(16, 2 * ndim))
    for i in range(nwalkers):
        for j in range(ndim):
            ax[j].plot(samples[j, :, i], c = 'green', lw=0.3)
            ax[j].set_ylabel(params[j], fontsize=15)
            ax[j].grid(True)
            ax[j].axhline(medians[j], c='blue', lw=1.2)
            ax[j].axhline(true_values[j], c='red', lw=1.2)
    fig.savefig(output_dir+'./plots/covergence.png',dpi=300)
    print('plot the covergence walk of parameters')


def plot_corner(fn, params):
    burnin      = 0.50
    slice_step  = 60
    chain = Chain.load(fn).remove_burnin(burnin)[::slice_step]
    
    g = plots.get_subplot_plotter()
    g.settings.fig_width_inch= 8
    g.settings.legend_fontsize = 20
    g.settings.axes_labelsize = 20
    g.settings.axes_fontsize = 16
    g.settings.figure_legend_frame = False
    plotting.plot_triangle(chain, title_limit=0, filled = True, params = params,
                            #    legend_labels = [r'$FOLPS+f_c$', r'$FOLPS$'], legend_loc= 'upper right',

                                # contour_ls = lss, contour_lws = lws, contour_colors = colors, 
                                # param_limits=param_limits, 
                                smoothed=True, show=False, g=g)
    true_values     = set_true_values(params)
    for i in range(len(true_values)):
        for j in range(i+1):
            g.subplots[i,j].axvline(true_values[j], c = 'k', ls = ':', lw = 1.2)
            if i != j:
                g.subplots[i,j].axhline(true_values[i], c = 'k', ls = ':', lw = 1.2)
    g.export(output_dir+'./plots/corner.png',dpi=300)
    
if __name__ == "__main__":
    # set the k bins
    kmin     = 0.005
    kmax     = 0.205
    binning  = 0.005
    k_ev     = np.arange(kmin, kmax, binning)
    klen     = len(k_ev)
    klim     = {ell*2: (kmin,kmax,binning) for ell in range(2)}

    cosmology = 'LCDM' #LCDM, nuCDM, wCDM
    theory_model = 'FOLPS'  # TNS, FOLPS, FOLPSRC

    emulator_fn = f'./results/emulators/emulator_{cosmology}_{theory_model}.npy'
    chain_fn = f'./results/samples/{cosmology}/chain_{tracer}_{region}_z{redshift}.npy'
    
    # MCMC sampling
    nwalkers = 120 
    interations = 1501 # save every 300 iterations
    
    if theory_model == 'FOLPS':
        theory_el = FOLPSTracerPowerSpectrumMultipoles(pt=EmulatedCalculator.load(emulator_fn))
    if theory_model == 'FOLPSRC':
        theory_el = FOLPSRCTracerPowerSpectrumMultipoles(pt=EmulatedCalculator.load(emulator_fn))  
    observable = TracerPowerSpectrumMultipolesObservable(data= pk_fn, 
                                                        klim=klim, 
                                                        covariance = cov_fns,
                                                        theory = theory_el,
                                                        kin=np.arange(0.001, 0.35, 0.002),
                                                        wmatrix=wmatrix_fn)                    
    likelihood = ObservablesGaussianLikelihood(observable)
    likelihood()
    sampler = EmceeSampler(likelihood, seed=42, nwalkers=nwalkers, save_fn =chain_fn)
    sampler.run(check={'max_eigen_gr': 0.05}, max_iterations = interations) # save every 300 iterations
    print("Sampling finished")