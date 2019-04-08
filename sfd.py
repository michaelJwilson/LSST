import  json
import  sfdmap
import  numpy               as  np
import  pylab               as  pl
import  matplotlib.pyplot   as  plt

from    astropy             import units as u
from    astropy.coordinates import SkyCoord
from    plot_map            import radec2project
from    numpy.random        import shuffle


sfd   = sfdmap.SFDMap('/global/homes/m/mjwilson/SFD/', scaling=1.0)

#c     = SkyCoord(ra=_FieldRa[i] * u.degree, dec= _FieldDec[i] * u.degree, frame='icrs')
#cg    = c.galactic
#l, b  = cg.l, cg.b

##  And plot.                                                                                                                                             
fig   =  plt.figure(figsize=(8, 8))
ax    =  plt.subplot(111, projection="aitoff")

ras   =  np.linspace(  0., 360., 5000)
decs  =  np.linspace(-90.,  90., 5000)

np.random.shuffle(ras)
np.random.shuffle(decs)

##  Get E(B-V) value at RA [deg], Dec [deg], J2000.
EBV   = np.array([sfd.ebv(ra, dec) for ra, dec in zip(ras, decs)]) 
cax   = plt.scatter(np.radians(ras) - np.pi, np.radians(decs), c=EBV, cmap=plt.cm.coolwarm, rasterized=True, vmin=0.0, vmax=0.5)
    
plt.grid(True)
pl.tight_layout()
  
pl.colorbar(cax, shrink=0.4)
pl.show()
