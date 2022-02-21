import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process("cscRechitCluster",Run2_2018)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

process.GlobalTag.globaltag = '102X_dataRun2_v12'
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/eos/uscms/store/user/kkwok/llp/ggH_HToSSTo4Tau_MH-125_TuneCP5_13TeV-powheg-pythia8_GEN-SIM-RECO.root' # choose your input file here
        #'file:/eos/uscms/store/user/lpclonglived/HLT/ggH_HToSSTo4Tau_MH-125_TuneCP5_13TeV-powheg-pythia8/0A05C63C-7C96-934A-ABA9-F4F21B1BA39E.root' # choose your input file here
        "file:/eos/cms/store/group/phys_susy/razor/run2/RAW/CSCShowerTriggerSkim/sixie/PassL1MuShowerNominal_2021213/merged/pickEvents_L1CSCSkim_Merged_0.root",
        "file:/eos/cms/store/group/phys_susy/razor/run2/RAW/CSCShowerTriggerSkim/sixie/PassL1MuShowerNominal_2021213/merged/pickEvents_L1CSCSkim_Merged_1.root",
        "file:/eos/cms/store/group/phys_susy/razor/run2/RAW/CSCShowerTriggerSkim/sixie/PassL1MuShowerNominal_2021213/merged/pickEvents_L1CSCSkim_Merged_2.root",
        "file:/eos/cms/store/group/phys_susy/razor/run2/RAW/CSCShowerTriggerSkim/sixie/PassL1MuShowerNominal_2021213/merged/pickEvents_L1CSCSkim_Merged_3.root",
        "file:/eos/cms/store/group/phys_susy/razor/run2/RAW/CSCShowerTriggerSkim/sixie/PassL1MuShowerNominal_2021213/merged/pickEvents_L1CSCSkim_Merged_4.root",
        "file:/eos/cms/store/group/phys_susy/razor/run2/RAW/CSCShowerTriggerSkim/sixie/PassL1MuShowerNominal_2021213/merged/pickEvents_L1CSCSkim_Merged_5.root",
        "file:/eos/cms/store/group/phys_susy/razor/run2/RAW/CSCShowerTriggerSkim/sixie/PassL1MuShowerNominal_2021213/merged/pickEvents_L1CSCSkim_Merged_6.root",
        "file:/eos/cms/store/group/phys_susy/razor/run2/RAW/CSCShowerTriggerSkim/sixie/PassL1MuShowerNominal_2021213/merged/pickEvents_L1CSCSkim_Merged_7.root",
        "file:/eos/cms/store/group/phys_susy/razor/run2/RAW/CSCShowerTriggerSkim/sixie/PassL1MuShowerNominal_2021213/merged/pickEvents_L1CSCSkim_Merged_8.root",
        "file:/eos/cms/store/group/phys_susy/razor/run2/RAW/CSCShowerTriggerSkim/sixie/PassL1MuShowerNominal_2021213/merged/pickEvents_L1CSCSkim_Merged_9.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_aa.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ab.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ac.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ad.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ae.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_af.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ag.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ah.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ai.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_aj.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ak.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_al.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_am.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_an.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ao.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ap.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_aq.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ar.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_as.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_at.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_au.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_av.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_aw.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ax.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ay.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_az.root",
        #"file:/eos/uscms/store/user/kkwok/llp/zerobias18/Events_ComparatorSkim_ba.root",

    )
)

runOnRaw=True

if runOnRaw:
    process.dummy1 = cms.ESSource("EmptyESSource",
                                      recordName = cms.string("CSCIndexerRecord"),
                                      firstValid = cms.vuint32(1),
                                      iovIsRunNotTime = cms.bool(True)
                                  )
    
    process.dummy2 = cms.ESSource("EmptyESSource",
                                      recordName = cms.string("CSCChannelMapperRecord"),
                                      firstValid = cms.vuint32(1),
                                      iovIsRunNotTime = cms.bool(True)
                                  )
    
    process.load('RecoMuon.MuonRechitClusterProducer.MuonSystemRawToReco_cff')
    process.CSCIndexerESProducer = cms.ESProducer("CSCIndexerESProducer", AlgoName = cms.string("CSCIndexerStartup") )
    process.CSCChannelMapperESProducer = cms.ESProducer("CSCChannelMapperESProducer", AlgoName = cms.string("CSCChannelMapperStartup") )

    process.unpack = cms.Path(process.muonSystemClusterSelSeq)
    process.reco = cms.Path(process.muonSystemClusterSelSeq*process.csc2DRecHits)


process.load('RecoMuon.MuonRechitClusterProducer.rechitCluster_cff')


if runOnRaw:
   #outf = "test_output_raw.root"
    outf = "~/eos/HLT/zerobias/output_data_123X_Feb14.root"
else:
    outf = "test_output_reco.root"

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(outf), # choose your output file
    outputCommands = cms.untracked.vstring(
              'drop *',
              'keep *_gen*_*_*',
              'keep *_csc2DRecHits_*_*',
              'keep *_DTRecHits_*_*',
              'keep *_*_*_cscRechit*',
              'keep *_*_*_cscJet*',
    )
)

process.producer = cms.Path(
                            #process.ca4cscRechitClusters*process.ca3cscRechitClusters*process.ca2cscRechitClusters*
                            #process.ca4dtRechitClusters*process.ca3dtRechitClusters*process.ca2dtRechitClusters
                            process.ca4CSCrechitClusters*
                            process.ca4DTrechitClusters
                        )
process.end = cms.EndPath(process.out)

if  runOnRaw:
    process.schedule = cms.Schedule( process.reco, process.producer, process.end)
else:
    process.schedule = cms.Schedule( process.producer, process.end)
