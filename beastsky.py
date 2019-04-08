import sqlite3
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

from   astropy             import units as u
from   astropy.coordinates import SkyCoord
from   plot_map            import radec2project


def remap_long(llong):
  if llong.value > 180.:
    return llong - 360. * u.deg

  else:
    return llong
  
print('\n\nWelcome to BEAST sky.\n\n')

##  https://github.com/LSSTScienceCollaborations/survey_strategy_wp/blob/master/tools/WFD%2Bbigarea.ipynb
##  http://opsim.lsst.org/runs/reference_run/minion_1016/summary_minion_1016.pdf
totalNvis  = 2447931  ## with snaps, (2x15s exposures/visit).
totalNvis *=    1.08  ## without snaps, (1x30s/visit).  Assuming 8% increase in visits. 

##  90% assumed to be available. 
percentTotal =         0.90
totalNvis   *= percentTotal

print("The number of visits available is %.2fM." % (totalNvis/1000000))

conn         =   sqlite3.connect('/global/cscratch1/sd/mjwilson/LSST/minion_1016_sqlite.db')
c            =   conn.cursor()

##  Field properties;  set of 5280 hexagons and 12 pentagons inscribed in circular fields having a 3.5-degree diameter                                                                                                                                                                                   
c.execute('SELECT {coi1},{coi2},{coi3},{coi4} FROM {tn}'.\
          format(coi1='fieldID', coi2='fieldFov', coi3='fieldRA', coi4='fieldDec', tn='Field'))

##  Max field is 5292.                                                                                                                                                                                                                                                                                    
_FieldID, _FieldFOV, _FieldRa, _FieldDec = np.split(np.array(c.fetchall()), 4, 1)

##  Closing the connection to the database file.                                                                                                                                                                                                                                                           
conn.close()

bigsky  = []
wfd     = []
bigNwfd = []

##  And plot.                                                                                                                                                                                                                                                                                            
fig     =  plt.figure(figsize=(8, 8))
ax      =  plt.subplot(111, projection="aitoff")

for i, id in enumerate(_FieldID):
  c     = SkyCoord(ra=_FieldRa[i] * u.degree, dec= _FieldDec[i] * u.degree, frame='icrs')                                                                                                                                                                                                             
  cg    = c.galactic                                                                                                                                                                                                                                                                                  
  l, b  = remap_long(cg.l), cg.b  
    
  if (_FieldDec[i] >= -72.25) & (_FieldDec[i] <= 12.4) & ((cg.b.value >= 15.) | (cg.b.value <= -15.)):
    wfd.append(id)

    ##  x, y = radec2project(_FieldRa[i], _FieldDec[i])
    ##  cax  = ax.scatter(x, y, c='r', rasterized=True, alpha=0.5)

  elif (_FieldDec[i] >= -90.) & (_FieldDec[i] <= 32.):
    bigsky.append(id)

    ##  x, y = radec2project(_FieldRa[i], _FieldDec[i])
    ##  cax  = ax.scatter(x, y, c='b', rasterized=True, alpha=0.5)

  else:
    pass

## 
for id in bigsky:
  if id not in wfd:
    bigNwfd.append(id)

##  Assuming 825 visits per field to WFD. 
wfdvisits   = len(wfd) * 825

##  Remaining footprint for Big Sky. 
availNvis   = totalNvis - wfdvisits

#  Exact number of visits per field in Big not WFD. 
visPerField = availNvis / len(bigNwfd)

#  Round the number of visits per field to an integer
visPerField = int(round(visPerField))

##  Additional visits to WFD.
propTotal   = visPerField * len(bigNwfd)

##  Current expected performance
single_m5 = {'u': 23.98, 'g': 24.91, 'r': 24.42, 'i': 23.97, 'z': 23.38, 'y': 22.47}

##  Relative distribution of filters for Big Sky Not WFD.
fractionsPerFilter = {'u': 0.2, 'g': 0.2, 'r': 0.20, 'i': 0.10, 'z': 0.10, 'y': 0.10}

visPerFilter = {}
fieldTotal   = 0

print('\nTotal visits to each field in BEAST sky:\n')

for f in fractionsPerFilter:
    visPerFilter[f] = int(round(fractionsPerFilter[f] * visPerField))
    fieldTotal += visPerFilter[f]

    print('Visits in %s: %d' % (f, visPerFilter[f]))

##  Jiggle by hand based on output. 
print('\nTotal per field: %d (compared to potential %d per field previously calculated)' %(fieldTotal, visPerField))

coadd_m5 = {}

print('\nCoadded depths achieved:\n')

for f in visPerFilter:
    coadd_m5[f] = single_m5[f] + 2.5 * np.log10(np.sqrt(visPerFilter[f]))
    print("Coadded depth in %s: %.2f" % (f, coadd_m5[f]))

# Using estimates from kraken_2026:
kraken_single_m5    = {'u': 23.78, 'g': 24.81, 'r': 24.35, 'i': 23.92, 'z': 23.34, 'y': 22.45}
kraken_visPerFilter = {'u': 64, 'g': 90, 'r': 206, 'i': 204, 'z': 186, 'y': 188}
opsimCoaddM5        = {'u': 25.65, 'g': 27.15, 'r': 27.20, 'i': 26.62, 'z': 25.72, 'y': 24.91}
kraken_coadd_m5 = {}

offset = {}

print('\nOpsim depth corrections\n')

for f in visPerFilter:
    kraken_coadd_m5[f] = kraken_single_m5[f] + 2.5 * np.log10(np.sqrt(kraken_visPerFilter[f]))
    offset[f]          = kraken_coadd_m5[f] - opsimCoaddM5[f]

    print("Coadded depth in %s: %.2f - correction is %.2f" % (f, kraken_coadd_m5[f], offset[f]))

saved_offsets = {'u': 0.39, 'g': 0.10, 'r': 0.04, 'i': 0.19, 'z': 0.46, 'y': 0.38}

print('\nCorrected depths in Beast Sky (not WFD).\n')

# Means we might achieve these coadded values:
for f in visPerFilter:
    coadd_m5[f] = single_m5[f] + 2.5 * np.log10(np.sqrt(visPerFilter[f])) - saved_offsets[f]

    print("Coadded depth in %s: %.2f" % (f, coadd_m5[f]))

##  And plot.
for i, id in enumerate(bigsky):
  id    = np.int(id)

  c     = SkyCoord(ra=_FieldRa[id][0] * u.degree, dec= _FieldDec[id][0] * u.degree, frame='icrs')
  cg    = c.galactic
  l, b  = remap_long(cg.l), cg.b

  x, y  = radec2project(_FieldRa[id][0], _FieldDec[id][0])
  cax   =    ax.scatter(x, y, c=coadd_m5['u'], cmap=plt.cm.coolwarm, rasterized=True, vmin=22., vmax=26., s=6)

plt.grid(True)
pl.tight_layout()
pl.colorbar(cax, shrink=0.4)
pl.show(block=True)

print('\n\nDone.\n\n')
