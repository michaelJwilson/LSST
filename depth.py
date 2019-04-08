import numpy             as      np
import pylab             as      pl
import matplotlib.pyplot as      plt

from   utils             import  latexify
from   collections       import  OrderedDict


latexify(columns=1, equal=True, fontsize=12, ggplot=True, usetex=True)

Coadd_m5             = {}
colors               = plt.rcParams['axes.prop_cycle'].by_key()['color']

single_m5            = OrderedDict()

single_m5['u']       = 23.98
single_m5['g']       = 24.91
single_m5['r']       = 24.42
single_m5['i']       = 23.97 
single_m5['z']       = 23.38 
single_m5['y']       = 22.47

kraken_visPerFilter  = {'u': 64,    'g':    90, 'r': 206,   'i':   204, 'z': 186,   'y':   188}
kraken_diffPerFilter = {'u': 2.258, 'g': 2.442, 'r': 2.892, 'i': 2.887, 'z': 2.837, 'y': 2.843}

##  Single visit overwrite, to match draft.  Fiducial 10 year depths:
FidDepths            = {'u': 25.3,  'g': 26.84, 'r': 27.04, 'i': 26.35, 'z': 25.22, 'y': 24.47}

##  Single visits    = Final - kraken_diffPerFilter = Final - 2.5 * np.log10(np.sqrt(kraken_visPerFilter[f]))  
SigDepths            = OrderedDict()

SigDepths['u']       = 23.042
SigDepths['g']       = 24.398
SigDepths['r']       = 24.148
SigDepths['i']       = 23.463
SigDepths['z']       = 22.383
SigDepths['y']       = 21.627

colors               = [colors[1], colors[5], colors[0], colors[2], colors[3], colors[4]]

for ii, f in enumerate(SigDepths):
  Ys      = np.arange(0., 10.1, 0.1)
  results = []

  for Y in Ys:
    Coadd_m5[f] = SigDepths[f] + 2.5 * np.log10(np.sqrt( Y * kraken_visPerFilter[f] / 10. ))

    results.append(Coadd_m5[f])

  results = np.array(results)  

  pl.plot(2023 + Ys, results, c=colors[ii], label=r'$%s$' % f)


pl.xlim(2022.5, 2033.5)
pl.legend(ncol=3, loc=4, frameon=False)

pl.xlabel('Year')
pl.ylabel('Depth')

plt.tight_layout()

pl.savefig('plots/LSST_timeline.pdf')
