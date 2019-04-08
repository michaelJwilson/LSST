import  sys 
import  json
import  pickle
import  sfdmap
import  sqlite3
import  numpy                           as  np
import  pylab                           as  pl
import  matplotlib.pyplot               as  plt

from    astropy                         import units as u
from    astropy.coordinates             import SkyCoord
from    schlafly                        import get_SchlaflyExt
from    astropy.coordinates.tests.utils import randomly_sample_sphere


def radec2project(ra, dec):
    return (np.radians(ra) - np.pi, np.radians(dec))

if __name__ == '__main__':
  if len(sys.argv) > 1:
    ##  python run.py YEAR filter test PROPID                                                                                     
    YEAR    =  np.int(sys.argv[1])
    filter  =  '%s' % sys.argv[2]
    test    =  np.int(sys.argv[3])
    PROPID  =  np.int(sys.argv[4])

    if PROPID == 99:
       PROPID  = 'ALL' 

  else:
    filter  =    'i'
    PROPID  =  'ALL'
    YEAR    =     1
    test    =     0

  toplot = 'CoAdd5sig_corr'  ## ['EBV', 'CoAdd5sig_corr', 'CoAdd5sig', 'A' + filter, 'EqualDepthCoadd']

  with open('json/lsst%s_YR%d_PROPID_%s.json' % (filter, YEAR, PROPID)) as json_file:  
    hitmap  = json.load(json_file)

  for id in hitmap:
      ##  Correct for galactic extinction. 
      hitmap[id]['CoAdd5sig_corr']  = hitmap[id]['CoAdd5sig']       - hitmap[id]['A' + filter]
      hitmap[id]['EqualDepthCoadd'] = hitmap[id]['EqualDepthCoadd'] - hitmap[id]['A' + filter]
      
  ##  Current expected single-visit performance.                                                                                                    
  ##  single_m5 = {'u': 23.98, 'g': 24.91, 'r': 24.42, 'i': 23.97, 'z': 23.38, 'y': 22.47}

  ##  Tile cut on i band depth.
  if filter  == 'i':
    ilims     = {1: 24.5, 5: 25.25, 10: 26.0}  

    ilim      =  ilims[YEAR]
    inLSST    =           []
    ninLSST   =           []  ## NOT inlsst.                                                                                                                

    field_ras =    []
    field_dcs =    []

    for field in hitmap.keys():
      if (hitmap[field]['CoAdd5sig_corr'] > ilim) & (hitmap[field]['EBV'] <= 0.2):
        field_ras.append(hitmap[field]['RA'])
        field_dcs.append(hitmap[field]['DEC'])

        inLSST.append(field)

      else:
        ninLSST.append(field)

    pickle.dump(inLSST,    open( "pickle/inlsst_PROPID_%s_YR%d.p"    % (PROPID, YEAR), "wb" ))
    pickle.dump(ninLSST,   open( "pickle/ninlsst_PROPID_%s_YR%d.p"   % (PROPID, YEAR), "wb" ))
    pickle.dump(field_ras, open( "pickle/field_ras_PROPID_%s_YR%d.p" % (PROPID, YEAR), "wb" ))
    pickle.dump(field_dcs, open( "pickle/field_dcs_PROPID_%s_YR%d.p" % (PROPID, YEAR), "wb" ))

  else:
    ##  Load the i-band limited depth tiles. 
    inLSST    = pickle.load(open( "pickle/inlsst_PROPID_%s_YR%d.p"     % (PROPID, YEAR), "rb" ))
    ninLSST   = pickle.load(open( "pickle/ninlsst_PROPID_%s_YR%d.p"    % (PROPID, YEAR), "rb" ))

  notin = []

  ##  Check that the i-band limited tiles are available for this filter. 
  for i, field in enumerate(inLSST):
    if field not in hitmap.keys():
      notin.append(field)

  ##  And remove from available tiles. 
  for x in notin:
    inLSST.remove(x)
    ninLSST.append(x)

  if filter is not 'i':
    ##  Save for this filter. 
    pickle.dump(inLSST,   open( "pickle/inlsst%s_PROPID_%s_YR%d.p"     % (filter, PROPID, YEAR), "wb" ))
    pickle.dump(ninLSST,  open( "pickle/ninlsst%s_PROPID_%s_YR%d.p"    % (filter, PROPID, YEAR), "wb" ))

  ##  Plot the galactic extinction values. 
  if toplot in ['A' + filter]:
      label =  toplot
      med   =     0.0

      ##  Bottom out colourscale at EBV=0.2
      vmin  =  get_SchlaflyExt(0.2, filter, Rv=None)
      vmax  =   10.0

  ##  Plot the galactic colour excess / relative extinction.
  elif toplot == 'EBV':
      label = toplot
      med   =    0.0

      vmin  =    0.2
      vmax  =    0.5

  else:
      ##  Plot the coadded or N-visit single visit depths.   
      med   = np.array([hitmap[field][toplot] for field in inLSST])
      med   = np.median(med)

      label = r'$%s$ - %.2lf' % (filter, med)
    
      vmin  = -0.5
      vmax  =  0.5

  ##  Area calc.                                                                                                                                      
    
  Area                       =         0.0
  nrand                      =  np.int(5e5)

  ##  LSST FOV radius in deg. 
  LSST_rad                   = 3.5 / 2.

  field_cat                  = SkyCoord(ra=[hitmap[field]['RA'] * u.degree for field in inLSST],\
                                        dec=[hitmap[field]['DEC'] * u.degree for field in inLSST])

  rand_ra, rand_dec, _       = randomly_sample_sphere(nrand)

  rand_ra                    =  rand_ra.to(u.deg)
  rand_dec                   = rand_dec.to(u.deg)

  rand_cat                   = SkyCoord(ra = rand_ra, dec = rand_dec)

  idxc, idxcatalog, d2d, _   = rand_cat.search_around_sky(field_cat, LSST_rad * u.deg)

  nrand_inlsst               = len(rand_cat[idxcatalog])

  ##  Can check by plotting random distribution.                                                                                                       
  ##  
  ##  for i, x in enumerate(rand_cat[idxcatalog]):                                                                                                       
  ##    x, y    = radec2project(x.ra.value, x.dec.value)                                                                                                     
  ##    ax.scatter(x, y, c='c', marker='^', rasterized=True, s=10, alpha=0.4)                                                                            

  DegOnSky                   = 41253.
  Area                       = (1. * nrand_inlsst / nrand) * DegOnSky

  print('For Y%d in %s and %s, the area is %.2lf deg2, with a median depth of %.2lf' % (YEAR, filter, PROPID, Area, med))

  ##  Plot
  fig       =  plt.figure(figsize=(8, 8))
  ax        =  plt.subplot(111, projection="aitoff")

  for field in inLSST:
    x, y    = radec2project(hitmap[field]['RA'], hitmap[field]['DEC'])
    cax     = ax.scatter(x, y, c=np.array(hitmap[field][toplot]) - med, cmap=plt.cm.coolwarm, rasterized=True, vmin=-0.4, vmax=0.4, s=6)

  '''
  for field in ninLSST:
    x, y    = radec2project(hitmap[field]['RA'], hitmap[field]['DEC'])
    ax.scatter(x, y, c='k', marker='o', rasterized=True, s=6, alpha=0.2)  
  '''

  plt.grid(True)
  pl.title(r'           Y%d (%s) Area:  %.1lf deg$^2$' % (YEAR, PROPID, Area))
  pl.tight_layout()
  pl.colorbar(cax, shrink=0.4, label=label)
    
  pl.show()
  ##  pl.savefig('plots/%sDepthMap_YR%s_PROPID_%s.pdf' % (filter, YEAR, PROPID))

  print('\n\nDone.\n\n')


