import FWCore.ParameterSet.Config as cms




from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register ('useMiniAOD',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "use miniAOD rather than AOD")
options.parseArguments()
useMiniAOD=options.useMiniAOD
process = cms.Process("GenAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")


process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(100)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            #'file:root://xrootd-cms.infn.it///store/user/mcampana/RunIISummer20UL18_SingleLQ_ueLQue_M1000_Lambda0p1_POWHEG-HerwigV7/RunIISummer20UL18_MiniAOD/211116_105025/0000/MiniAOD_1.root'
             			open('MiniAOD.list').readlines()
			   )
                            )



from LeptoquarkAnalysis.MiniAODAnalyzer.tools import setupVIDForHEEPV70
setupVIDForHEEPV70(process,useMiniAOD=useMiniAOD)


#from LeptoquarkAnalysis.MiniAODAnalyzer.tools import addHEEPV70ElesMiniAOD
#addHEEPV70ElesMiniAOD(process,useStdName=True)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple_lepto.root")
)


process.GenAnalysis = cms.EDAnalyzer('MiniAODAnalyzer',
                                     generatorInfo= cms.InputTag("generator"),
                                     genparticles    = cms.untracked.InputTag("prunedGenParticles", "", "PAT"),
                                     muons = cms.untracked.InputTag("slimmedMuons","","PAT"),
                                     electrons = cms.untracked.InputTag("slimmedElectrons","","PAT"),
                                       vid=cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
                                       vidBitmap=cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70Bitmap")
                              )

process.p = cms.Path(
	process.egmGsfElectronIDSequence*
	process.GenAnalysis)

