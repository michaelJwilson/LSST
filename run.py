import  sys
import  json
import  sfdmap
import  sqlite3
import  numpy               as  np
import  pylab               as  pl
import  matplotlib.pyplot   as  plt

from    astropy             import units as u
from    astropy.coordinates import SkyCoord
from    schlafly            import get_SchlaflyExt


def radec2project(ra, dec):
    return (np.radians(ra) - np.pi, np.radians(dec))

##  https://arxiv.org/pdf/1809.01669.pdf
##  Appendix C1. 

##  https://www.lsst.org/scientists/simulations/opsim/opsim-v335-benchmark-surveys
##  Knutt Olsen LSST notebooks.

##  http://opsim.lsst.org/runs/reference_run/minion_1016/summary_minion_1016.pdf


verbose       =   True

##  
if len(sys.argv) > 1:
  ##  python run.py YEAR filter test PROPID
  YEAR        =  np.int(sys.argv[1])
  filter      =  '"%s"' % sys.argv[2]
  test        =  np.int(sys.argv[3])
  PROPID      =  np.int(sys.argv[4])

else: 
  YEAR        =      10
  filter      =   '"r"'  
  test        =   True    ##  Row limit for test runs;  Limits the number of visits considered.
  PROPID      =     99    ##  WFD: 54;  North ecliptic spur:  55 (griz only); 99 mapped to ALL.

print('Solving for: YEAR %d, filter %s, test %d, PROPID %d' % (YEAR, filter, test, PROPID))

##  http://ops2.lsst.org/docs/current/architecture.html                                                                                                   
conn          =   sqlite3.connect('/global/cscratch1/sd/mjwilson/LSST/minion_1016_sqlite.db')
c             =   conn.cursor()

if verbose:
  for Table in ['ObsHistory', 'Field', 'Proposal', 'ObsHistory_Proposal', 'Proposal_Field']:
    c.execute('PRAGMA TABLE_INFO({})'.format(Table))

    print([tup[1] for tup in c.fetchall()])

'''                                                                                                                                                         
http://opsim.lsst.org/runs/reference_run/minion_1016/summary_minion_1016.pdf
    -- Total number of visits is 2,447,931, with 85.1% spent on the main wide-fast-deep (WFD).  2447931 * 0.851 = 2083758. 
    -- The median single-visit depths for WFD fields are 23.14, 24.47, 24.16, 23.40, 22.23, 21.57 in the ugrizy bands4

##  Pg. 27 of https://github.com/LSSTScienceCollaborations/ObservingStrategy/blob/pdf/whitepaper/LSST_Observing_Strategy_White_Paper.pdf                   

(52, u'conf/survey/GalacticPlaneProp.conf', u'WL', 4367689744, u'minion', 1016)                 --                                                          
(53, u'conf/survey/SouthCelestialPole-18.conf', u'WL', 4367691152, u'minion', 1016)             --                                                          
(54, u'conf/survey/Universal-18-0824B.conf', u'WLTSS', 4367690832, u'minion', 1016)             --  WFD.                                                    
(55, u'conf/survey/NorthEclipticSpur-18c.conf', u'WLTSS', 4367690768, u'minion', 1016)          --                                                          
(56, u'conf/survey/DDcosmology1.conf', u'WLTSS', 4367691088, u'minion', 1016)                   --  Deep Drilling fields                                     
'''

##  Field properties;  set of 5280 hexagons and 12 pentagons inscribed in circular fields having a 3.5-degree diameter                                                                                                                  
c.execute('SELECT {coi1},{coi2},{coi3},{coi4} FROM {tn}'.\
          format(coi1='fieldID', coi2='fieldFov', coi3='fieldRA', coi4='fieldDec', tn='Field'))

##  Max field is 5292.                                                                                                                                                                                                                    
_FieldID, _FieldFOV, _FieldRa, _FieldDec = np.split(np.array(c.fetchall()), 4, 1)

##  (1, 5292)
##  print(_FieldID.min(), _FieldID.max())

if PROPID in [54, 55]:
    ##  A many-to-many relationship table that stores the fields fieldID from the Field table that maps to the field centers as specified by proposal propID.
    c.execute('SELECT {coi1},{coi2} FROM {tn} WHERE {cn}={pn}'.format(coi1='Field_fieldID', coi2='proposal_field_id', tn='Proposal_Field', cn='Proposal_propID', pn=PROPID))

    ##  WFD:  Field_fieldID, proposal_field_id.
    WFD_FieldID, WFD_PropFieldID = np.split(np.array(c.fetchall()), 2, 1)

    ## (2293, 2776, 39847, 42986).
    ## print(len(WFD_FieldID), WFD_FieldID.max(), WFD_PropFieldID.min(), WFD_PropFieldID.max())

    ##  This table maps visits to a field to the proposal or proposals that requested it.
    c.execute('SELECT {coi1},{coi2} FROM {tn} WHERE {cn} = {pn}'.\
              format(coi1='ObsHistory_obsHistID', coi2='obsHistory_propID', tn='ObsHistory_Proposal', cn='Proposal_propID', pn=PROPID))

    _ObsHistID, _PropIDs = np.split(np.array(c.fetchall()), 2, 1)

    ##  Note:  first number:  visits in the WFD program. 
    ##  WFD (54):  (2083758, 2447931, 21064332, 29833372)
    ##  NES (55):  ( 158912, 1225657, 21075714, 27870824)
    ##  print(PROPID, len(_ObsHistID), _ObsHistID.max(), _PropIDs.min(), _PropIDs.max())

    ##  This table keeps a record of each visit made by the telescope during a simulated survey.
    ##  Note:  all visits (2 x 15s exposures) where filter is filter.
    if test:
        c.execute('SELECT {coi1},{coi2},{coi3},{coi4},{coi5} FROM {tn} WHERE {cn}={dn} LIMIT {ln}'.\
                  format(coi1='Field_fieldID', coi2='expMJD', coi3='visitExpTime', coi4='fiveSigmaDepth', coi5='obsHistID', tn='ObsHistory', cn='filter', dn=filter, ln=str(1e4)))

    else:
        c.execute('SELECT {coi1},{coi2},{coi3},{coi4},{coi5} FROM {tn} WHERE {cn}={dn}'.\
                  format(coi1='Field_fieldID', coi2='expMJD', coi3='visitExpTime', coi4='fiveSigmaDepth', coi5='obsHistID', tn='ObsHistory', cn='filter', dn=filter))

else:
    PROPID = 'ALL'

    ##  ** IN ALL PROGRAMS ** 
    ##  A many-to-many relationship table that stores the fields fieldID from the Field table that maps to the field centers as specified by proposal propID.                                                                                                                    
    c.execute('SELECT {coi1},{coi2} FROM {tn}'.format(coi1='Field_fieldID', coi2='proposal_field_id', tn='Proposal_Field', cn='Proposal_propID'))

    ##  WFD:  Field_fieldID, proposal_field_id.                                                                                                                                                                                                                                  
    WFD_FieldID, WFD_PropFieldID = np.split(np.array(c.fetchall()), 2, 1)

    ## (2293, 2776, 39847, 42986).                                                                                                                                                                                                                                                   ## print(len(WFD_FieldID), WFD_FieldID.max(), WFD_PropFieldID.min(), WFD_PropFieldID.max())                                                                                                                                                                                  
    ##  This table maps visits to a field to the proposal or proposals that requested it.                                                                                                                                                                                        
    c.execute('SELECT {coi1},{coi2} FROM {tn}'.\
              format(coi1='ObsHistory_obsHistID', coi2='obsHistory_propID', tn='ObsHistory_Proposal', cn='Proposal_propID'))

    _ObsHistID, _PropIDs = np.split(np.array(c.fetchall()), 2, 1)

    ##  Note:  first number:  visits in the WFD program.                                                                                                                                                                                                                         
    ## (2083758, 2447931, 21064332, 29833372)                                                                                                                                                                                                                                    
    ##  print(len(_ObsHistID), _ObsHistID.max(), _ObsHistPropIDs.min(), _ObsHistPropIDs.max(), _PropIDs.min(), _PropIDs.max())                                                                                                                                                    
    ##  This table keeps a record of each visit made by the telescope during a simulated survey.                                                                                                                                                                                 
    ##  Note:  all visits (2 x 15s exposures) where filter is filter.                                                                                                                                                                                                             
    if test:
        ##  Return properties for each visit if in requested filter.  Test:  limited to XX rows. 
        c.execute('SELECT {coi1},{coi2},{coi3},{coi4},{coi5} FROM {tn} WHERE {cn}={dn} LIMIT {ln}'.\
                  format(coi1='Field_fieldID', coi2='expMJD', coi3='visitExpTime', coi4='fiveSigmaDepth', coi5='obsHistID', tn='ObsHistory', cn='filter', dn=filter, ln=str(1e4)))

    else:
        c.execute('SELECT {coi1},{coi2},{coi3},{coi4},{coi5} FROM {tn} WHERE {cn}={dn}'.\
                  format(coi1='Field_fieldID', coi2='expMJD', coi3='visitExpTime', coi4='fiveSigmaDepth', coi5='obsHistID', tn='ObsHistory', cn='filter', dn=filter))


##  FieldIDs cut on filter. 
FieldIDs, MJDs, ExpTimes, FiveSigs, ObsHistID = np.split(np.array(c.fetchall()), 5, 1)
ObsHistID                                     = ObsHistID.astype(np.int)

assert np.all(np.diff(MJDs)       >= 0.)
assert np.all(np.diff(ObsHistID)  >= 0.)
assert np.all(np.diff(_ObsHistID) >= 0.)

##  Closing the connection to the database file                                                                                                                                     
conn.close()
  
##  (1.0, 3793.0)
##  print(FieldIDs.min(), FieldIDs.max())

##  (2063, 94527), cut on filter; 2063, 2064, 2065 ...
##  print(ObsHistID.min(), ObsHistID.max()

##  Cut to the ObvsIDs available up to a given year in the requested proposals (for programs).  Y1: 10% etc.
_ObsHistID   = _ObsHistID[:np.round(YEAR * 0.1 * len(_ObsHistID)).astype(np.int)]
_PropIDs     =   _PropIDs[:np.round(YEAR * 0.1 * len(_PropIDs)).astype(np.int)]

WFDFIELDIDS  = []

for i, x in enumerate(ObsHistID):
    
  if x in _ObsHistID:  ## In the requested programs, e.g. WFD. 
    WFDFIELDIDS.append([x[0], FieldIDs[i][0], FiveSigs[i][0]])

    completeness  = 100. * i / len(ObsHistID)

    if completeness > 5.:
        ##  break
        pass

  print('Percentage complete: %.2lf' % (100. * i / len(ObsHistID)))

##  
WFDFIELDIDS    = np.array(WFDFIELDIDS)
uWFDIDs, uVis  = np.unique(WFDFIELDIDS[:,1], return_counts=True)  ##  Unique WFD field ids and number of visits. 

##  Visits per WFD field.                                                                                                                                                                           
hitmap         = dict(zip(uWFDIDs, uVis))

for id in hitmap:
  hitmap[id]   = {'Visits': hitmap[id], 'CoAdd5sig': 0.0}

##  Current expected per visit performance.                                                                                                                                                        
single_m5      = {'u': 23.14, 'g': 24.47, 'r': 24.16, 'i': 23.40, 'z': 22.23, 'y': 21.57}
                              
##  By default, a scaling of 0.86 is applied to the map values to reflect the recalibration by Schlafly & Finkbeiner (2011)                                                                         
##  https://github.com/kbarbary/sfdmap.                                                                                                                                                              
sfd            = sfdmap.SFDMap('/global/homes/m/mjwilson/SFD/')
          
for i, x in enumerate(WFDFIELDIDS[:,1]):                                                                                                                         
    x        = np.int(x)

    ##  RA, DEC [degrees].
    EBV      = sfd.ebv(_FieldRa[x][0], _FieldDec[x][0])
    Ab       = get_SchlaflyExt(EBV, filter.replace('"', ''), Rv=None)
    
    ##  c    = SkyCoord(ra=_FieldRa[i] * u.degree, dec= _FieldDec[i] * u.degree, frame='icrs')
    ##  cg   = c.galactic
    ##  l,b  = cg.l, cg.b

    ##  Coadded 5 sig. depth:  One. Visit. 5. Sigma. Depth. + 2.5 * np.log10(np.sqrt(Vists. Per. Filter.))
    ##  Coadded 5 sig. depth:  sigm_COADD ** 2 = [\sum_visit (1 / visit depth) ** 2] ^-1.
    hitmap[x] = {'CoAdd5sig': hitmap[x]['CoAdd5sig'] + 10.**(-0.8*(single_m5[filter.replace('"','')] - FiveSigs[i][0])),\
                 'Visits': int(hitmap[x]['Visits']),\
                 'RA': _FieldRa[x][0],\
                 'DEC': _FieldDec[x][0],\
                 'FOV': _FieldFOV[x][0],\
                 'EBV': EBV,\
                 'A%s' % filter.replace('"',''): Ab,\
                 'FieldID': int(x),\
                 'SingleVisitDepth': single_m5[filter.replace('"','')],\
                 'EqualDepthCoadd': single_m5[filter.replace('"','')] + 2.5 * np.log10(np.sqrt(hitmap[x]['Visits']))\
                 }    
    
for id in hitmap:
  hitmap[id]['CoAdd5sig'] = single_m5[filter.replace('"','')] + 2.5 * np.log10(hitmap[id]['CoAdd5sig']) / 2.
    
print(hitmap)

with open('json/lsst%s_YR%d_PROPID_%s.json' % (filter.replace('"', ''), YEAR, str(PROPID)), 'w') as outfile:
  json.dump(hitmap, outfile)

print('\n\nDone.\n\n')
