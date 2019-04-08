import  os 
import  numpy  as  np
import  pylab  as  pl


test = 0

for YEAR in [1]:                                  ##  [1, 5, 10]
  ##  i first as defines detection band.
  for filter in ['i', 'u', 'g', 'r', 'z', 'y']:   ##  ['i', 'u', 'g', 'r', 'z', 'y']:
    for PROPID in [54]:                           ##  54:  WFD, 55:  NES, 99:  ALL.
      ##  
      ##  os.system('python run.py %d %s %d %d'  % (YEAR, filter, test, PROPID))
      os.system('python plot_map.py %d %s %d %d' % (YEAR, filter, test, PROPID))

      break
