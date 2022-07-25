

#cmsRun runCSCTriggerPrimitiveProducer_cfg.py maxEvents=-1 unpack=true l1=true run3=true dqm=true saveEdmOutput=true preTriggerAnalysis=true  inputFiles=file:~/eos/HLT/Commissioning2022/877caf18-4328-415e-bb74-b40adf51f9b4.root runNumber=351405
#cmsRun runCSCTriggerPrimitiveProducer_cfg.py maxEvents=-1 unpack=true l1=true run3=true dqm=true saveEdmOutput=true preTriggerAnalysis=true  inputFiles=file:~/eos/HLT/Commissioning2022/877caf18-4328-415e-bb74-b40adf51f9b4.root runNumber=351405
#cmsRun runCSCTriggerPrimitiveProducer_cfg.py maxEvents=-1 unpack=true l1=true run3=true dqm=true saveEdmOutput=true preTriggerAnalysis=true  inputFiles=file:/eos/cms/tier0/store/data/Run2022A/ZeroBias7/RAW/v1/000/353/689/00000/a2a1e5ca-c7ba-479b-b148-5d4d05620fae.root runNumber=353689

#cmsRun runCSCTriggerPrimitiveProducer_cfg.py maxEvents=-1 unpack=true l1=true run3=true dqm=true saveEdmOutput=true preTriggerAnalysis=true  inputFiles=file:/eos/cms/store/express/Run2022A/ExpressCosmics/FEVT/Express-v1/000/353/685/00000/556e513e-da56-4c82-9e25-42e68ec9c63a.root runNumber=353685

#cmsRun runCSCTriggerPrimitiveProducer_cfg.py maxEvents=-1 unpack=true l1=true run3=true dqm=true saveEdmOutput=true preTriggerAnalysis=true  inputFiles=list_351405.txt runNumber=351405
#cmsRun runCSCL1TDQMClient_cfg.py maxEvents=10000
#cmsRun runCSCTriggerPrimitiveAnalyzer_cfg.py dataVsEmulation=true dataVsEmulationFile=file:DQM_V0001_R000350619__Global__CMSSW_X_Y_Z__RECO.root runNumber=350619 maxEvents=10000
## Physics stream, no L1 seed
#cmsRun runCSCTriggerPrimitiveProducer_cfg.py maxEvents=10000 unpack=true l1=true run3=true dqm=true saveEdmOutput=true preTriggerAnalysis=true  inputFiles=file:/eos/cms/store/express/Run2022A/ExpressPhysics/FEVT/Express-v1/000/353/018/00000/bb9a164a-c7e8-4f54-9fe5-5248e458e561.root runNumber=353018
## ZeroBias stream
#cmsRun runCSCTriggerPrimitiveProducer_cfg.py maxEvents=-1 unpack=true l1=true run3=true dqm=true saveEdmOutput=true preTriggerAnalysis=true  inputFiles=file:/eos/cms/store/data/Run2022A/ZeroBias/RAW/v1/000/353/018/00000/3c488acf-ea53-4f64-8770-23fc2f9e8ce6.root runNumber=353018

#cmsRun runCSCShowerAnalyzer_cfg.py inputFiles=file:./lcts2.root

#cmsRun runCSCShowerAnalyzer_cfg.py inputFiles=list_351405.txt runNumber=351407 maxEvent=10
#cmsRun runCSCShowerAnalyzer_cfg.py inputFiles=file:~/eos/HLT/Commissioning2022/877caf18-4328-415e-bb74-b40adf51f9b4.root runNumber=351407 maxEvents=100

#cmsRun runCSCShowerAnalyzer_cfg.py inputFiles=file:/eos/cms/store/express/Run2022B/ExpressPhysics/FEVT/Express-v1/000/355/404/00000/4062441b-32e0-4d73-81f6-ae95bd78dc5d.root runNumber=355404 maxEvents=-1
#cmsRun runCSCShowerAnalyzer_cfg.py inputFiles=file:~/eos/HLT/Commissioning2022/877caf18-4328-415e-bb74-b40adf51f9b4.root runNumber=351407 maxEvents=-1


#cmsRun runCSCShowerAnalyzer_cfg.py inputFiles=file:/eos/cms/tier0/store/data/Run2022C/ZeroBias/RAW/v1/000/355/921/00000/b766054e-a1e4-4af4-b30b-f0f34573b7b5.root runNumber=355921 maxEvents=100


#python3 printcount.py -i plots_353689.root -o run353689/
#python3 printcount.py -i plots_353689.root -o run355100/

#python3 printcount.py -i /eos/cms/store/user/kakwok/HLT/Commissioning2022/run355100/plots.root -o run355100/
#python3 printcount.py -i /eos/cms/store/user/kakwok/HLT/Commissioning2022/run355404/plots.root -o run355404/
python3 printcount.py -i /eos/cms/store/user/kakwok/HLT/Commissioning2022/run355921/plots.root -o run355921/
#python3 printcount.py -i /eos/cms/store/user/kakwok/HLT/Commissioning2022/run355404_zb/plots.root -o run355404_zb/

#python3 analyze_data.py -i ~/eos/HLT/Commissioning2022/run355100/output_reco.root -o analyze_355100.root --isData
#python3 analyze_data.py -i ~/eos/HLT/Commissioning2022/run355404/output_reco.root -o analyze_355404.root --isData
#python3 analyze_data.py -i ~/eos/HLT/Commissioning2022/run355100/output_reco.root -o analyze_355100.root --isData
#python3 analyze_data.py -i ~/eos/HLT/Commissioning2022/run355404_zb/output_reco.root -o analyze_355404_zb.root --isData
#python3 analyze_data.py -i ~/eos/HLT/Commissioning2022/run355558/output_reco.root -o analyze_355558.root --isData

#python3 analyze_data.py -i ~/eos/HLT/Commissioning2022/run355769/output_reco.root -o analyze_355769.root --isData

#python3 analyze_data.py -i ./output_reco.root -o analyze_351407.root --isData
#hadd -f hmt_ntuple_351407.root  analyze_351407.root ./plots_351407.root 


#hadd hmt_ntuple_355558.root  analyze_355558.root ~/eos/HLT/Commissioning2022/run355558/plots.root
#hadd hmt_ntuple_355100.root  analyze_355100.root ~/eos/HLT/Commissioning2022/run355100/plots.root
#hadd -f hmt_ntuple_355404.root  analyze_355404.root ~/eos/HLT/Commissioning2022/run355404/plots.root
#hadd hmt_ntuple_355404_zb.root  analyze_355404_zb.root ~/eos/HLT/Commissioning2022/run355404_zb/plots.root
#hadd -f hmt_ntuple_355769.root  analyze_355769.root ~/eos/HLT/Commissioning2022/run355769/plots.root
