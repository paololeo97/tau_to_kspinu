 ###########################################################################
# Selection of 3%1 prong taupair decays with loose cuts
#
#                 l(h) nu nu
#                  |
#        e+ e- -> tau+ tau-
#                       |
#                     h h h nu
#
# Contributors:
# Michel Hernandez-Villanueva, Petar Rados, Ami Rostomyan, Francesco Tenchini
#
# Use cases: tau analysis
#
# last modified: February 2020 
###########################################################################

# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import basf2 as b2
import modularAnalysis as ma
import variables.collections as vc
import variables.utils as vu
from variables import variables
import sys

######################################################
# Select some options...
######################################################

### flags unique to mc
sys_track_eff = [False, 0.0069]
correct_pi0_eff = [True, 0.914, 0.020]

### flags unique to data
sys_track_scale = [False, 1.00056]
sys_photon_energy = False

applyMomentumSF = 'sf'  # sf down or up. Will be set to False for MC

######################################################
# Specify path to data sample and database
######################################################
main = b2.create_path()

arg_dataORmc = str(sys.argv[1])
arg_input = str(sys.argv[2])
arg_output = str(sys.argv[3])

if arg_dataORmc == 'data':
    b2.B2INFO("CHECK: analysing data, is this correct?")
    ma.inputMdst('default', '%s/*.root' % arg_input, path=main)
elif arg_dataORmc == 'mc':
    b2.B2INFO("CHECK: analysing MC, is this correct?")
    ma.inputMdst('default', arg_input, path=main)
else:
    b2.B2INFO("WARNING: unkown input...")
    
######################################################
# create and fill the ParticleLists
######################################################

ma.fillParticleList('e-:all', '', path=main)
ma.fillParticleList('mu-:all', '', path=main)
ma.fillParticleList('pi+:all', '', path=main)

ma.fillParticleList('gamma:all', '', path=main)

if sys_photon_energy:
    ma.correctEnergyBias(inputListNames=['gamma:all'], tableName="PhotonEnergyBiasCorrection_MC13a_Sep2021" , path=main)

##################################################################################
# pi0 reconstruction with loose photons 
##################################################################################

gammaDetectorLocation = {
    "FWD" : "clusterReg == 1",
    "BRL" : "clusterReg == 2",
    "BWD" : "clusterReg == 3"
}

gammaForPi0lists = []
for g in gammaDetectorLocation.keys():
    
    gammaForPi0Cuts = gammaDetectorLocation[g]
    gammaForPi0Cuts += ' and E > 0.1'
    gammaForPi0Cuts += ' and abs(clusterTiming) < 200'
    gammaForPi0Cuts += ' and -0.8660 < cosTheta < 0.9563'
    gammaForPi0Cuts += ' and clusterNHits > 1.5'
    gammaForPi0Cuts += ' and (minC2TDist > 40 or E > 0.4)'
    
    gammaForPi0 = 'gamma:looseForPi0{}'.format(g)
    gammaForPi0lists.append(gammaForPi0)
    
    ma.cutAndCopyLists(gammaForPi0, 'gamma:all', gammaForPi0Cuts, path=main)

# cos of opening angle between the photons
variables.addAlias('cosAngle2Photons', 
'formula((daughter(0, px) * daughter(1, px) + daughter(0, py) * daughter(1, py) + daughter(0, pz) * daughter(1, pz) ) / daughter(0, p) / daughter(1, p) )')

# other variables used for cuts
variables.addAlias('leadingclusterE', 'formula(max(daughter(0, clusterE ),daughter(1, clusterE)))' )
variables.addAlias('subleadingclusterE', 'formula(min(daughter(0, clusterE ),daughter(1, clusterE)))')

# Determine pi0 reco for individual Detector Parts (values from Zuzana Master thesis)
Pi0CutLabel = ["leadingclusterE","subleadingclusterE","cosAngle2Photons","p"]
Pi0CutValue = {
"FWD,FWD" : [0.5625, 0.1625, 0.9458, 0.9444],
"BRL,BRL" : [0.4125, 0.0625, 0.8875, 0.6333],
"BWD,BWD" : [0.4125, 0.1125, 0.8708, 0.6111],
"BRL,FWD" : [0.3625, 0.0875, 0.8875, 0.5889],
"BRL,BWD" : [0.3625, 0.0875, 0.8875, 0.5889]
}

Pi0lists = []
for cut in Pi0CutValue.keys():
    gammalists = cut.split(",")
    CurrentPi0List = 'pi0:fromLooseGammas{}{}'.format(gammalists[0],gammalists[1])
    Pi0lists.append(CurrentPi0List)
    Pi0Cut = '0.115 < M < 0.152'
    for i,c in enumerate(Pi0CutLabel):
        Pi0Cut += ' and {} > {}'.format(c,Pi0CutValue[cut][i])
    ma.reconstructDecay('{} -> gamma:looseForPi0{} gamma:looseForPi0{}'.format(CurrentPi0List,gammalists[0],gammalists[1]), Pi0Cut, path=main)


ma.copyLists('pi0:fromLooseGammas', Pi0lists, path=main)
variables.addAlias('nPi0', 'countInList(pi0:fromLooseGammas)')

ma.cutAndCopyLists('gamma:pi0', gammaForPi0lists, 'isDescendantOfList(pi0:fromLooseGammas) == 1', path=main)
variables.addAlias('nPhotonsFromPi0', 'countInList(gamma:pi0)')

variables.addAlias('nPhotonsAll', 'countInList(gamma:all)')

if correct_pi0_eff[0]:
    pi0_ineff_cut = 'random < '+str(correct_pi0_eff[1])
    ma.cutAndCopyLists('pi0:fromIneffCorr', 'pi0:fromLooseGammas', pi0_ineff_cut, path=main)
    variables.addAlias('nPi0_Ineff', 'countInList(pi0:fromIneffCorr)')
    ####
    pi0_ineff_cut_up = 'random < '+str(correct_pi0_eff[1] + correct_pi0_eff[2])
    ma.cutAndCopyLists('pi0:fromIneffCorrUp', 'pi0:fromLooseGammas', pi0_ineff_cut_up, path=main)
    variables.addAlias('nPi0_IneffUp', 'countInList(pi0:fromIneffCorrUp)')
    ####
    pi0_ineff_cut_dn = 'random < '+str(correct_pi0_eff[1] - correct_pi0_eff[2])
    ma.cutAndCopyLists('pi0:fromIneffCorrDn', 'pi0:fromLooseGammas', pi0_ineff_cut_dn, path=main)
    variables.addAlias('nPi0_IneffDn', 'countInList(pi0:fromIneffCorrDn)')
else:
    ma.copyList('pi0:fromIneffCorr', 'pi0:fromLooseGammas', path=main)
    variables.addAlias('nPi0_Ineff', 'countInList(pi0:fromIneffCorr)')
    ####
    ma.copyList('pi0:fromIneffCorrUp', 'pi0:fromLooseGammas', path=main)
    variables.addAlias('nPi0_IneffUp', 'countInList(pi0:fromIneffCorrUp)')
    ####
    ma.copyList('pi0:fromIneffCorrDn', 'pi0:fromLooseGammas', path=main)
    variables.addAlias('nPi0_IneffDn', 'countInList(pi0:fromIneffCorrDn)')
    
############################################################################################################
# track (not used for reconstruction of K_0S) and photon (not used for reconstruction of pi0) cuts 
############################################################################################################
trackCuts = '-3.0 < dz < 3.0'
trackCuts += ' and dr < 1.0'

variables.addAlias('EoverP', 'formula( ifNANgiveX( clusterE, -1 )/p )')

# TODO pid selections (to develop later, if the data and MC PIDs would agree )
eIDCuts  = 'EoverP > 0.8'
piIDCuts = 'EoverP < 0.8'
muIDCuts = piIDCuts
    
ma.cutAndCopyLists('pi+:tracks', 'pi+:all', trackCuts, path=main)
ma.cutAndCopyLists('pi+:pid', 'pi+:tracks', piIDCuts, path=main)
ma.cutAndCopyLists('e-:pid', 'e-:all', trackCuts+' and '+eIDCuts, path=main)

gammaCuts = 'E > 0.1'
gammaCuts += ' and -0.8660 < cosTheta < 0.9563'
gammaCuts += ' and clusterNHits > 1.5'
gammaCuts += ' and isDaughterOfList(pi0:fromLooseGammas) == 0'
gammaCuts += ' and abs(clusterTiming) < 200'
gammaCuts += ' and (minC2TDist > 40 or E > 0.4)'

ma.cutAndCopyLists('gamma:notPi0', 'gamma:all', gammaCuts, path=main)

if arg_dataORmc == 'mc':
    #can I use fillParticleList to add gen particles?
    main.add_module('ParticleLoader', decayStringsWithCuts=[('pi0:gen', '')], useMCParticles=True )

######################################################
# tracking related corrections and systematics
######################################################

if sys_track_eff[0]:
    ma.trackingEfficiency(inputListNames=['e-:pid', 'pi+:pid'], fraction=sys_track_eff[1], path=main)

if sys_track_scale[0]:
    ma.scaleTrackMomenta(inputListNames=['e-:pid', 'pi+:pid'], scale=sys_track_scale[1], path=main)    
    
######################################################
# event based cut - 4 tracks in event
######################################################
variables.addAlias('nGoodPhotons', 'countInList(gamma:notPi0)')

variables.addAlias('nGoodTracks', 'countInList(pi+:tracks)')
ma.applyEventCuts('nGoodTracks == 4', path=main)

#######################################################
# EventShape and EventKinamatics modules
#######################################################
ma.buildEventShape(['pi+:tracks','gamma:pi0', 'gamma:notPi0'], path=main)
ma.buildEventKinematics(['pi+:tracks','gamma:pi0', 'gamma:notPi0'], path=main)

######################################################
# Signal and tag sides
#######################################################
# -- 3 prong
ma.reconstructDecay('tau+:3prong -> pi+:pid pi-:pid pi+:pid', '', path=main) 

# ---- pseudomass of 3 prong
variables.addAlias('Mmin_nom', 'formula( [ daughter(0,M) * daughter(0,M) + 2. * \
                                       ( 5.289699 - useCMSFrame(daughter(0,E)) ) * \
                       ( useCMSFrame(daughter(0,E)) - useCMSFrame(daughter(0,p)) ) \
                     ]^0.5 )' ) 

variables.addAlias('Mmin', 'formula( [ daughter(0,M) * daughter(0,M) + 2. * \
                                       ( Ecms*0.5 - useCMSFrame(daughter(0,E)) ) * \
                       ( useCMSFrame(daughter(0,E)) - useCMSFrame(daughter(0,p)) ) \
                     ]^0.5 )' )

variables.addAlias('mcMmin', 'formula( [ mcDaughter(0,M) * mcDaughter(0,M) + 2. * \
                                       ( Ecms*0.5 - useCMSFrame(mcDaughter(0,E)) ) * \
                       ( useCMSFrame(mcDaughter(0,E)) - useCMSFrame(mcDaughter(0,p)) ) \
                     ]^0.5 )' )

# ---- invariant masses for 3prong decay
variables.addAlias('invM12',  'daughterInvariantMass(0, 1)')
variables.addAlias('invM23',  'daughterInvariantMass(1, 2)')


# -- 1 prong
ma.reconstructDecay('tau-:e -> e-:pid', '', path=main, dmID=11)
ma.reconstructDecay('tau-:pimu -> pi-:pid', '', path=main, dmID=2113)
variables.addAlias('dmID_1prong', 'daughter(1, extraInfo(decayModeID))') #reconstructed 1-prong decay mode
ma.copyLists('tau-:1prong',['tau-:e','tau-:pimu'], path=main)

# -- gamma* -> tau+ tau-
ma.reconstructDecay('vpho -> tau+:3prong tau-:1prong', '', path=main)

######################################################
# pions and lepton on the opposide sides
######################################################
variables.addAlias('prod1',
                   'formula(daughter(0, daughter(0, cosToThrustOfEvent))*daughter(1, daughter(0,cosToThrustOfEvent)))')
variables.addAlias('prod2',
                   'formula(daughter(0, daughter(1, cosToThrustOfEvent))*daughter(1, daughter(0,cosToThrustOfEvent)))')
variables.addAlias('prod3',
                   'formula(daughter(0, daughter(2, cosToThrustOfEvent))*daughter(1, daughter(0,cosToThrustOfEvent)))')

ma.applyCuts('vpho', 'prod1 < 0 and prod2 < 0 and prod3 < 0', path=main)

############################################################################################
# number of photons and pi0 on 3-prong and 1-prong
# put requirements on the number of pi0s on 1-prong to choose the decay mode
############################################################################################
ma.copyList('gamma:1prong', 'gamma:notPi0', path=main)
ma.copyList('gamma:3prong', 'gamma:notPi0', path=main)

ma.copyList('pi0:1prong', 'pi0:fromLooseGammas', path=main)
ma.copyList('pi0:3prong', 'pi0:fromLooseGammas', path=main)

# 1-prong on positive or negative side of thrust axis
variables.addAlias('1prongInPosThrust', 'countInList(vpho, daughter(1, daughter(0,cosToThrustOfEvent)) > 0)')
variables.addAlias('1prongInNegThrust', 'countInList(vpho, daughter(1, daughter(0,cosToThrustOfEvent)) < 0)')

positiveThrust = b2.create_path()
negativeThrust = b2.create_path()

ma.applyCuts('gamma:1prong', 'cosToThrustOfEvent > 0', path=positiveThrust)
ma.applyCuts('gamma:3prong', 'cosToThrustOfEvent < 0', path=positiveThrust)

ma.applyCuts('gamma:1prong', 'cosToThrustOfEvent < 0', path=negativeThrust)
ma.applyCuts('gamma:3prong', 'cosToThrustOfEvent > 0', path=negativeThrust)

ma.applyCuts('pi0:1prong', 'cosToThrustOfEvent > 0', path=positiveThrust)
ma.applyCuts('pi0:3prong', 'cosToThrustOfEvent < 0', path=positiveThrust)

ma.applyCuts('pi0:1prong', 'cosToThrustOfEvent < 0', path=negativeThrust)
ma.applyCuts('pi0:3prong', 'cosToThrustOfEvent > 0', path=negativeThrust)

# take different paths if 1-prong in cosToThrustOfEvent > or < 0
sigThrustModule = main.add_module('VariableToReturnValue', variable='1prongInPosThrust')
sigThrustModule.if_value('> 0', positiveThrust, b2.AfterConditionPath.CONTINUE)
sigThrustModule = main.add_module('VariableToReturnValue', variable='1prongInNegThrust')
sigThrustModule.if_value('> 0', negativeThrust, b2.AfterConditionPath.CONTINUE)

# the number of photons and pi0s in 3prong and 1prong hemispheres
variables.addAlias('nPhotons_1prong', 'nParticlesInList(gamma:1prong)')
variables.addAlias('nPhotons_3prong', 'nParticlesInList(gamma:3prong)')

variables.addAlias('nPi0s_1prong', 'nParticlesInList(pi0:1prong)')
variables.addAlias('nPi0s_3prong', 'nParticlesInList(pi0:3prong)')

######################################################
# kinematics of 3-prong and 1-prong
######################################################
variables.addAlias('photonE_3prong', 'totalEnergyOfParticlesInList(gamma:3prong)')
variables.addAlias('photonE_1prong', 'totalEnergyOfParticlesInList(gamma:1prong)')

variables.addAlias('photonECMS_3prong', 'useCMSFrame(totalEnergyOfParticlesInList(gamma:3prong))')
variables.addAlias('photonECMS_1prong', 'useCMSFrame(totalEnergyOfParticlesInList(gamma:1prong))')

# select the tracks that enter the reconstruction of tau+ and tau- and vpho
ma.cutAndCopyList('tau+:3prongFromVpho', 'tau+:3prong', 'isDaughterOfList(vpho)', path=main)
ma.cutAndCopyList('tau-:1prongFromVpho', 'tau-:1prong', 'isDaughterOfList(vpho)', path=main)

# energy 3prong
variables.addAlias('ECMS_3prong', 'useCMSFrame(totalEnergyOfParticlesInList(tau+:3prongFromVpho))')
variables.addAlias('ECMS_3prong_pi0s',
                   'formula(useCMSFrame(totalEnergyOfParticlesInList(tau+:3prongFromVpho)) + useCMSFrame(totalEnergyOfParticlesInList(pi0:3prong)))')
variables.addAlias('ECMS_3prong_photons',
                   'formula(useCMSFrame(totalEnergyOfParticlesInList(tau+:3prongFromVpho)) + useCMSFrame(totalEnergyOfParticlesInList(gamma:3prong)))')
variables.addAlias('ECMS_3prong_photons_pi0s',
                   'formula(useCMSFrame(totalEnergyOfParticlesInList(tau+:3prongFromVpho)) + useCMSFrame(totalEnergyOfParticlesInList(gamma:3prong)) + useCMSFrame(totalEnergyOfParticlesInList(pi0:3prong)) )')

# energy 1prong
variables.addAlias('ECMS_1prong', 'useCMSFrame(totalEnergyOfParticlesInList(tau-:1prongFromVpho))')
variables.addAlias('ECMS_1prong_photons',
                   'formula(useCMSFrame(totalEnergyOfParticlesInList(tau-:1prongFromVpho)) + useCMSFrame(totalEnergyOfParticlesInList(gamma:1prong)))')

variables.addAlias('M_3prong', 'invMassInLists(tau+:3prongFromVpho)')
variables.addAlias('M_3prong_pi0s', 'invMassInLists(tau+:3prongFromVpho, pi0:3prong)')
variables.addAlias('M_3prong_photons', 'invMassInLists(tau+:3prongFromVpho, gamma:3prong)')
variables.addAlias('M_3prong_photons_pi0s', 'invMassInLists(tau+:3prongFromVpho, gamma:3prong, pi0:3prong)')

variables.addAlias('M_1prong', 'invMassInLists(tau-:1prongFromVpho)')
variables.addAlias('M_1prong_photons', 'invMassInLists(tau-:1prongFromVpho, gamma:1prong)')


if arg_dataORmc == 'mc':
    variables.addAlias('nGenPi0', 'countInList(pi0:gen)')

######################################################
# perform MC matching for MC samples
######################################################
if arg_dataORmc == 'mc':
    ma.matchMCTruth('vpho', path=main)
    ma.labelTauPairMC(path=main)

#####################################################
# select the variables to be stored in the ntuple
#####################################################
# CMS variables
variables.addAlias('E_CMS','useCMSFrame(E)')
variables.addAlias('p_CMS','useCMSFrame(p)')
variables.addAlias('px_CMS','useCMSFrame(px)')
variables.addAlias('py_CMS','useCMSFrame(py)')
variables.addAlias('pz_CMS','useCMSFrame(pz)')
variables.addAlias('pt_CMS','useCMSFrame(pt)')
variables.addAlias('theta_CMS','useCMSFrame(theta)')
variables.addAlias('phi_CMS','useCMSFrame(phi)')
# expert PID
variables.addAlias('muid_KLM','pidProbabilityExpert(13, KLM)')
variables.addAlias('eid_KLM','pidProbabilityExpert(11, KLM)')
variables.addAlias('piid_KLM','pidProbabilityExpert(211, KLM)')
variables.addAlias('muid_ECL','pidProbabilityExpert(13, ECL)')
variables.addAlias('eid_ECL','pidProbabilityExpert(11, ECL)')
variables.addAlias('piid_ECL','pidProbabilityExpert(211, ECL)')
variables.addAlias('pimuid_KLM','pidDeltaLogLikelihoodValueExpert(211, 13, KLM)')
# ancestors
variables.addAlias('fromTau', 'hasAncestor(15)')
variables.addAlias('fromKS', 'hasAncestor(310)')
## ancestors
variables.addAlias('fromTau', 'hasAncestor(15)')
variables.addAlias('fromKS', 'hasAncestor(310)')
variables.addAlias("genMotherPDG0", "genMotherPDG(0)")
variables.addAlias("genMotherPDG1", "genMotherPDG(1)")
variables.addAlias("genMotherPDG2", "genMotherPDG(2)")
variables.addAlias("genMotherPDG3", "genMotherPDG(3)")
variables.addAlias("genMotherPDG4", "genMotherPDG(4)")
ancestorVariables  =  ['fromTau', 'fromKS', 'genMotherPDG0', 'genMotherPDG1', 'genMotherPDG2', 'genMotherPDG3']

# -- event based variables
eventVariables = ['dmID_1prong',
                  'nGoodTracks',
                  'nGoodPhotons', 'nPhotons_3prong', 'nPhotons_1prong', 'nPhotonsFromPi0',
                  'nPi0', 'nPi0s_3prong', 'nPi0s_1prong',
                  #'nK0S'
                  ]
eventVariables += ['thrust',
                   'visibleEnergyOfEventCMS',
                   'missingMomentumOfEvent', 'missingMomentumOfEvent_theta',
                   'missingMomentumOfEventCMS', 'missingMomentumOfEventCMS_theta',
                   'missingMass2OfEvent',
                   'Mmin', 'Mmin_nom',
                  ]
eventVariables += ['photonE_3prong', 'photonECMS_3prong',
                   'photonE_1prong', 'photonECMS_1prong',
                   'ECMS_3prong', 'ECMS_3prong_photons', 'ECMS_3prong_photons_pi0s', 'ECMS_3prong_pi0s',
                   'ECMS_1prong', 'ECMS_1prong_photons', 
                   'M_3prong', 'M_3prong_photons', 'M_3prong_photons_pi0s', 'M_3prong_pi0s',
                   'M_1prong', 'M_1prong_photons'
                   ]

eventVariables += ['Ecms', 'beamE', 'beamPx', 'beamPy', 'beamPz']

commonVariables = vc.kinematics 
commonVariables += ['theta', 'cosTheta', 'phi']
commonVariables += ['E_CMS', 'p_CMS', 'px_CMS', 'py_CMS', 'pz_CMS', 'pt_CMS', 'theta_CMS', 'phi_CMS']
commonVariables += ['charge', 'cosToThrustOfEvent']

# -- tau candidate variables
tauVariables =  vc.inv_mass + vc.vertex + ['chiProb']

tau3prongVariables = ['invM12', 'invM23']

# -- track level variables
trackVariables = ['clusterE', 'EoverP']
trackVariables += vc.pid + ['muid_KLM', 'eid_KLM', 'piid_KLM', 'muid_ECL', 'eid_ECL', 'piid_ECL', 'pimuid_KLM']
trackVariables += vc.track_hits
trackVariables += ['dz', 'dr']

# -- MC specific info
if arg_dataORmc == 'mc':
  # -- event variables
  eventVariables += ['tauPlusMCMode', 'tauMinusMCMode', 'tauPlusMCProng', 'tauMinusMCProng', 'mcMmin', 'nGenPi0'] 
  # -- common variables   
  commonVariables += vc.mc_variables + vc.mc_truth + ancestorVariables
  # -- tau vertex truth
  tauVariables += vc.mc_vertex
  tauVariables += vc.mc_variables

# -- data specific info
if arg_dataORmc == 'data':
# -- L1 
    variables.addAlias('psnm_atleastone', 'L1Trigger')
    eventVariables += ['psnm_atleastone']
    n_trigs = 135
    for ix in range(n_trigs):
      variables.addAlias('psnm_%i' % ix, 'L1PSNMBit(%i)' % ix)
      eventVariables += ['psnm_%i' % ix]
      variables.addAlias('ftdl_%i' % ix, 'L1FTDLBit(%i)' % ix)
      eventVariables += ['ftdl_%i' % ix]
    n_triginput = 146
    for iy in range(n_triginput):
        variables.addAlias('triginput_%i' % iy, 'L1InputBit(%i)' % iy)
        eventVariables += ['triginput_%i' % iy]            
        
vphoVariableList = vu.create_aliases_for_selected(list_of_variables=eventVariables,
                                                  decay_string='^vpho') + \
                   vu.create_aliases_for_selected(list_of_variables=commonVariables + tauVariables,
                                                  decay_string='vpho -> ^tau+ ^tau-', 
                                                  prefix=['tau_3prong','tau_1prong']) + \
                   vu.create_aliases_for_selected(list_of_variables=tau3prongVariables,
                                                  decay_string='vpho -> ^tau+ tau-', 
                                                  prefix=['tau_3prong']) + \
                   vu.create_aliases_for_selected(list_of_variables=commonVariables + trackVariables,
                                                  decay_string='vpho -> [tau+ -> ^pi+ ^pi- ^pi+] [tau- -> ^pi-]',
                                                  prefix=['track1_3prong','track2_3prong','track3_3prong','track_1prong']) 



######################################################
# Write flat ntuples
######################################################

ma.variablesToNtuple(decayString='vpho',
                     variables=vphoVariableList,
                     filename=arg_outfile,
                     treename='tau3x1',
                     path=main)


# Process the events
b2.print_path(main)
b2.process(main)

#from basf2.utils import pretty_print_module
#for m in main.modules():
#    pretty_print_module(m, m.name())
print(b2.statistics)
