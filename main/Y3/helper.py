import numpy as np

REDSHIFT_VSMEAR = dict(LRG = [(0.4, 0.6), (0.6, 0.8), (0.8, 1.1)], 
                    ELG = [(0.8, 1.1), (1.1, 1.3), (1.3, 1.6)],
                    QSO = [(0.8, 1.1), (1.1, 1.4), (1.4, 1.7), (1.7, 2.1)])

REDSHIFT_CUBICBOX = dict(LRG = [0.500, 0.800, 0.800],
                         ELG= [0.800, 1.100, 1.100],
                         QSO = [0.800, 1.100, 1.100, 1.100])

TRACER_CUBICBOX = dict(LRG = 'lrg', ELG= 'elg', QSO = 'qso')

EDGES = dict(pk=np.arange(0.,0.4001,0.001),
    xi=(np.linspace(0., 200, 201), np.linspace(-1., 1., 201)),
    mpslog=(np.geomspace(0.01, 100., 49), np.linspace(-1., 1., 201)))

def bins(statistics):
    if statistics == 'xipoles':
        binning   = 4
        rmin      = 20
        rmax      = 200
        lenr      = 45
        return (rmin, rmax, binning, lenr)
    elif statistics == 'pkpoles':
        kmin     = 0.0
        kmax     = 0.401
        binning  = 0.005
        lenk = 80
        return (kmin, kmax, binning, lenk)
    elif statistics == 'mpslog':
        smin      = 0.20
        smax      = 30
        return (smin, smax, binning, 0)
    elif statistics == 'wp':
        pmin = 0.20
        pmax = 30
        return (pmin, pmax, 0, 0)