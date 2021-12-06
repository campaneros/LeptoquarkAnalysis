// -*- C++ -*-
//
// Package:    LeptoquarkAnalysis/MiniAODAnalyzer
// Class:      MiniAODAnalyzer
//
/**\class MiniAODAnalyzer MiniAODAnalyzer.cc LeptoquarkAnalysis/MiniAODAnalyzer/plugins/MiniAODAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mattia Campana
//         Created:  Fri, 03 Dec 2021 13:52:30 GMT
//
//



// system include files
#include <memory>
//
// // user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
//#include <SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h>

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEgamma/ElectronIdentification/interface/ElectronIDAlgo.h"



#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

// root include files
#include "TTree.h"
#include "TLorentzVector.h"
#include <TStyle.h>
#include <TCanvas.h>



#include "LeptoquarkAnalysis/MiniAODAnalyzer/interface/CutNrs.h"
#include "LeptoquarkAnalysis/MiniAODAnalyzer/interface/VIDCutCodes.h"




namespace{
  struct NrPassFail {
    NrPassFail():nrPass(0),nrFail(0){}
    mutable std::atomic<int> nrPass;
    mutable std::atomic<int> nrFail;
  };
}



struct tree_struc_{
  int run;
  int lumi;
  long unsigned int event;
  int nEle;
  int nMu;
  float TQ_genMass;
  float TQ_recoMass;
  std::vector<float>            genLep_pt;
  std::vector<float>            genLep_eta;
  std::vector<float>            genLep_mass;
  std::vector<float>            genLep_phi;
  std::vector<int>            genLep_pdgId;
  std::vector<float>            genMom_pt;
  std::vector<float>            genMom_eta;
  std::vector<float>            genMom_mass;
  std::vector<float>            genMom_phi;
  std::vector<int>            genMom_pdgId;

  std::vector<float>            recoMu_pt;
  std::vector<float>            recoMu_eta;
  std::vector<float>            recoMu_mass;
  std::vector<float>            recoMu_phi;
  std::vector<float>            recoMu_dR;

  std::vector<float>            recoEle_pt;
  std::vector<float>            recoEle_eta;
  std::vector<float>            recoEle_mass;
  std::vector<float>            recoEle_phi;
  std::vector<float>            recoEle_dR;


};




//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class MiniAODAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MiniAODAnalyzer(const edm::ParameterSet&);
      ~MiniAODAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void initTreeStructure();
      void clearVectors(); 

      // ----------member data ---------------------------
      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
        edm::Service<TFileService> fs;
  edm::EDGetTokenT<GenEventInfoProduct>genInfoToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> >genparticleToken_;
  edm::EDGetTokenT<std::vector<pat::Muon> > muonToken_;
  edm::EDGetTokenT<std::vector<pat::Electron> > eleToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > vidPassToken_;
  edm::EDGetTokenT<edm::ValueMap<unsigned int> > vidBitmapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> >  vidResultToken_;


 //setup tree;
  TTree* tree;
    tree_struc_ tree_;


   TH1F* HEEP_pt;  
   TH1F* Electron;
   TH1F* Ratio;
     

   TH1F* HEEP_single;  
   TH1F* HEEP;  
   TCanvas* c;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MiniAODAnalyzer::MiniAODAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
     genInfoToken_= consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generatorInfo"));
  genparticleToken_    = consumes<std::vector<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genparticles"));
  muonToken_=consumes<std::vector<pat::Muon>>(iConfig.getUntrackedParameter<edm::InputTag>("muons"));
  eleToken_=consumes<std::vector<pat::Electron>>(iConfig.getUntrackedParameter<edm::InputTag>("electrons"));
  vidPassToken_=consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("vid"));
  vidBitmapToken_=consumes<edm::ValueMap<unsigned int> >(iConfig.getParameter<edm::InputTag>("vidBitmap"));
  vidResultToken_=consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("vid"));
}




MiniAODAnalyzer::~MiniAODAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MiniAODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(genInfoToken_, genInfo);  

  edm::Handle<std::vector<reco::GenParticle> > genparticles;
  iEvent.getByToken(genparticleToken_, genparticles);

  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByToken(eleToken_, electrons);
  edm::Handle<edm::ValueMap<bool> > vidPass;
  edm::Handle<edm::ValueMap<unsigned int> > vidBitmap;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > vidResult;

  iEvent.getByToken(vidPassToken_,vidPass);
  iEvent.getByToken(vidBitmapToken_,vidBitmap);
  iEvent.getByToken(vidResultToken_,vidResult);


  std::vector<float>genLep_pt;
  std::vector<float>genLep_mass;
  std::vector<float>genLep_eta;
  std::vector<float>genLep_phi;
  std::vector<int>genLep_pdgId;

  std::vector<float>genMom_pt;
  std::vector<float>genMom_mass;
  std::vector<float>genMom_eta;
  std::vector<float>genMom_phi;
  std::vector<int>genMom_pdgId;


  std::vector<float>recoEle_pt;
  std::vector<float>recoEle_mass;
  std::vector<float>recoEle_eta;
  std::vector<float>recoEle_phi;
  std::vector<float>recoEle_dR;
  std::vector<float>recoEle_index;


  unsigned long int event = iEvent.id().event();   
  int run                 = iEvent.id().run();
  int lumi                = iEvent.luminosityBlock();
  int nEle=0;
  int nMu=0;




 for (const auto & genpar_iter : *genparticles){

 if((abs(genpar_iter.pdgId())!=13 && abs(genpar_iter.pdgId())!=11) && genpar_iter.status()!=1)continue;
 if(genpar_iter.mother(0)->pdgId()!=9911561) continue; 

      genLep_pt.push_back(genpar_iter.pt());
      genLep_mass.push_back(genpar_iter.mass());
      genLep_eta.push_back(genpar_iter.eta());
      genLep_phi.push_back(genpar_iter.phi());
      genLep_pdgId.push_back(genpar_iter.pdgId());
      

      genMom_pt.push_back(genpar_iter.mother(0)->pt());
      genMom_mass.push_back(genpar_iter.mother(0)->mass());
      genMom_eta.push_back(genpar_iter.mother(0)->eta());
      genMom_phi.push_back(genpar_iter.mother(0)->phi());
      genMom_pdgId.push_back(genpar_iter.mother(0)->pdgId());
      
      if(abs(genpar_iter.pdgId())==11)nEle++;
      if(abs(genpar_iter.pdgId())==13)nMu++;

}

  for(unsigned int i=0;i<genLep_pdgId.size();i++){
    if(abs(genLep_pdgId[i])!=11)continue;
    float dR=999.;
    float dRMin=999.;
    float recopt;
    float recoeta;
    float recophi;
    float recomass;
    int index;   
   // int heepIDBits;
  // bool HEEP_ind= false;
  // bool heepID;
   bool passHEEPV70=false;
   size_t count =0;
   unsigned int heepV70Bitmap;

    //for (size_t eleNr=0;eleNr<electrons->size();eleNr++){
  //	edm::Ptr<reco::GsfElectron> elePtr(electrons,eleNr);
//	std::cout<<"eleNr   "<<eleNr<<std::endl; 
//	count++	
//}
      for (const auto & el_iter: *electrons){
   // for (size_t eleNr=0;eleNr<electrons->size();eleNr++){
     // const auto & el_iter : 
       edm::Ptr<reco::GsfElectron> elePtr(electrons,count);      
      float eta=el_iter.eta();
      float phi=el_iter.phi();
      
      
      dR= deltaR(eta,phi,genLep_eta[i],genLep_phi[i]);
      if(dR<dRMin){
	recopt=el_iter.pt();
	recoeta=el_iter.eta();
	recophi=el_iter.phi();
	recomass=el_iter.mass();
	dRMin=dR;
        bool heepID=(*vidPass)[elePtr];
	passHEEPV70=heepID;
        heepV70Bitmap = (*vidBitmap)[elePtr];	
        //heepID=el_iter.userInt("heepElectronID-HEEPV70");
        //HEEP_ind = heepID;
      }
  // std::cout<<"count   "<<count<<std::endl;
   count++;  
    } 
    if (dRMin<0.2){
      recoEle_pt.push_back(recopt);
      recoEle_mass.push_back(recomass);
      recoEle_eta.push_back(recoeta);
      recoEle_phi.push_back(recophi);
      recoEle_dR.push_back(dRMin);
      Electron->Fill(recopt);
      if (passHEEPV70) HEEP_pt->Fill(recopt);
       int bincounter=1;
    using HEEPV70 = VIDCutCodes<cutnrs::HEEPV70>;
 bool pass_singleET = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET});
 bool pass_singleETA = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ETA});
 bool pass_singleDETAINSEED = HEEPV70::pass(heepV70Bitmap,{HEEPV70::DETAINSEED});
 bool pass_singleDPHIIN = HEEPV70::pass(heepV70Bitmap,{HEEPV70::DPHIIN});
 bool pass_singleSIGMAIETAIETA = HEEPV70::pass(heepV70Bitmap,{HEEPV70::SIGMAIETAIETA});
 bool pass_singleE2X5OVER5X5 = HEEPV70::pass(heepV70Bitmap,{HEEPV70::E2X5OVER5X5});
 bool pass_singleHADEM = HEEPV70::pass(heepV70Bitmap,{HEEPV70::HADEM});
 bool pass_singleTRKISO = HEEPV70::pass(heepV70Bitmap,{HEEPV70::TRKISO});
 bool pass_singleEMHADD1ISO = HEEPV70::pass(heepV70Bitmap,{HEEPV70::EMHADD1ISO});
 bool pass_singleDXY = HEEPV70::pass(heepV70Bitmap,{HEEPV70::DXY});
 bool pass_singleMISSHITS = HEEPV70::pass(heepV70Bitmap,{HEEPV70::MISSHITS});
 bool pass_singleECALDRIVEN = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ECALDRIVEN});

 bool passETA = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET,HEEPV70::ETA});
 bool passDETAINSEED = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET,HEEPV70::ETA,HEEPV70::DETAINSEED});
 bool passDPHIIN = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET,HEEPV70::ETA,HEEPV70::DETAINSEED});
 bool passSIGMAIETAIETA = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET,HEEPV70::ETA,HEEPV70::DETAINSEED,HEEPV70::DPHIIN,HEEPV70::SIGMAIETAIETA});
 bool passE2X5OVER5X5 = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET,HEEPV70::ETA,HEEPV70::DETAINSEED,HEEPV70::DPHIIN,HEEPV70::SIGMAIETAIETA,HEEPV70::E2X5OVER5X5});
 bool passHADEM = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET,HEEPV70::ETA,HEEPV70::DETAINSEED,HEEPV70::DPHIIN,HEEPV70::SIGMAIETAIETA,HEEPV70::E2X5OVER5X5,HEEPV70::HADEM});
 bool passTRKISO = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET,HEEPV70::ETA,HEEPV70::DETAINSEED,HEEPV70::DPHIIN,HEEPV70::SIGMAIETAIETA,HEEPV70::E2X5OVER5X5,HEEPV70::HADEM,HEEPV70::TRKISO});
 bool passEMHADD1ISO = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET,HEEPV70::ETA,HEEPV70::DETAINSEED,HEEPV70::DPHIIN,HEEPV70::SIGMAIETAIETA,HEEPV70::E2X5OVER5X5,HEEPV70::HADEM,HEEPV70::TRKISO,HEEPV70::EMHADD1ISO});
 bool passDXY = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET,HEEPV70::ETA,HEEPV70::DETAINSEED,HEEPV70::DPHIIN,HEEPV70::SIGMAIETAIETA,HEEPV70::E2X5OVER5X5,HEEPV70::HADEM,HEEPV70::TRKISO,HEEPV70::EMHADD1ISO,HEEPV70::DXY});
 bool passMISSHITS = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET,HEEPV70::ETA,HEEPV70::DETAINSEED,HEEPV70::DPHIIN,HEEPV70::SIGMAIETAIETA,HEEPV70::E2X5OVER5X5,HEEPV70::HADEM,HEEPV70::TRKISO,HEEPV70::EMHADD1ISO,HEEPV70::DXY,HEEPV70::MISSHITS});
 bool passECALDRIVEN = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET,HEEPV70::ETA,HEEPV70::DETAINSEED,HEEPV70::DPHIIN,HEEPV70::SIGMAIETAIETA,HEEPV70::E2X5OVER5X5,HEEPV70::HADEM,HEEPV70::TRKISO,HEEPV70::EMHADD1ISO,HEEPV70::DXY,HEEPV70::MISSHITS,HEEPV70::ECALDRIVEN});
 bool passHEEP = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET,HEEPV70::ETA,HEEPV70::DETAINSEED,HEEPV70::DPHIIN,HEEPV70::SIGMAIETAIETA,HEEPV70::E2X5OVER5X5,HEEPV70::HADEM,HEEPV70::TRKISO,HEEPV70::EMHADD1ISO,HEEPV70::DXY,HEEPV70::MISSHITS,HEEPV70::ECALDRIVEN});


        HEEP->Fill(bincounter-1);
        HEEP_single->Fill(bincounter-1);

HEEP->GetXaxis()->SetBinLabel(bincounter,"NoCuts");
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"NoCuts");
                ++bincounter;

          if (pass_singleET)
                {
        HEEP->Fill(bincounter-1);
        HEEP_single->Fill(bincounter-1);
        HEEP->GetXaxis()->SetBinLabel(bincounter,"E_{T}");
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"E_{T}");
                }
        ++bincounter;

                  if (pass_singleETA)
                {
        HEEP_single->Fill(bincounter-1);
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"#eta");
                }
                  if (passETA)
                {
        HEEP->Fill(bincounter-1);
        HEEP->GetXaxis()->SetBinLabel(bincounter,"#eta");
                }
        ++bincounter;


                  if (pass_singleDETAINSEED)
                {
        HEEP_single->Fill(bincounter-1);
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"#eta in seed");
                }
                 if (passDETAINSEED)
                {
        HEEP->Fill(bincounter-1);
        HEEP->GetXaxis()->SetBinLabel(bincounter,"#eta in seed");
                }
        ++bincounter;


                  if (pass_singleDPHIIN)
                {
        HEEP_single->Fill(bincounter-1);
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"#Delta#phi");
                }
                  if (passDPHIIN)
                {
        HEEP->Fill(bincounter-1);
        HEEP->GetXaxis()->SetBinLabel(bincounter,"#Delta#phi");
                }
        ++bincounter;


                  if (pass_singleSIGMAIETAIETA)
                {
        HEEP_single->Fill(bincounter-1);
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"#sigma_{i#eta}");
                }
                  if (passSIGMAIETAIETA)
                {
        HEEP->Fill(bincounter-1);
        HEEP->GetXaxis()->SetBinLabel(bincounter,"#sigma_{i#eta}");
                }
        ++bincounter;


                  if (pass_singleE2X5OVER5X5)
                {
        HEEP_single->Fill(bincounter-1);
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"E^{2x5}/E^{5x5}");
                }
                  if (passE2X5OVER5X5)
                {
        HEEP->Fill(bincounter-1);
        HEEP->GetXaxis()->SetBinLabel(bincounter,"E^{2x5}/E^{5x5}");
                }
        ++bincounter;
                 if (pass_singleHADEM)
                {
        HEEP_single->Fill(bincounter-1);
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"HADEM");
                }
                  if (passHADEM)
                {
        HEEP->Fill(bincounter-1);
        HEEP->GetXaxis()->SetBinLabel(bincounter,"HADEM");
                }
        ++bincounter;

        if (pass_singleTRKISO)
                {
        HEEP_single->Fill(bincounter-1);
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"TRKISO");
                }
                  if (passTRKISO)
                {
        HEEP->Fill(bincounter-1);
        HEEP->GetXaxis()->SetBinLabel(bincounter,"TRKISO");
                }
        ++bincounter;


        if (pass_singleEMHADD1ISO)
                {
        HEEP_single->Fill(bincounter-1);
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"EMHADD1ISO");
                }
                  if (passEMHADD1ISO)
                {
        HEEP->Fill(bincounter-1);
        HEEP->GetXaxis()->SetBinLabel(bincounter,"EMHADD1ISO");
                }
        ++bincounter;

        if (pass_singleDXY)
                {
        HEEP_single->Fill(bincounter-1);
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"DXY");
                }
                  if (passDXY)
                {
        HEEP->Fill(bincounter-1);
        HEEP->GetXaxis()->SetBinLabel(bincounter,"DXY");
                }
        ++bincounter;


        if (pass_singleMISSHITS)
                {
        HEEP_single->Fill(bincounter-1);
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"MISSHITS");
                }
                  if (passMISSHITS)
                {
        HEEP->Fill(bincounter-1);
        HEEP->GetXaxis()->SetBinLabel(bincounter,"MISSHITS");
                }
        ++bincounter;

        if (pass_singleECALDRIVEN)
                {
        HEEP_single->Fill(bincounter-1);
        HEEP_single->GetXaxis()->SetBinLabel(bincounter,"ECALDRIVEN");
                }
                  if (passECALDRIVEN)
                {
        HEEP->Fill(bincounter-1);
        HEEP->GetXaxis()->SetBinLabel(bincounter,"ECALDRIVEN");
                }
        ++bincounter;
 
    }  
 }




initTreeStructure();

clearVectors();


   for(unsigned int i=0;i<genLep_pt.size();i++){
      tree_.genLep_pt.push_back(genLep_pt[i]);
      tree_.genLep_mass.push_back(genLep_mass[i]);
      tree_.genLep_eta.push_back(genLep_eta[i]);
      tree_.genLep_phi.push_back(genLep_phi[i]);
      tree_.genLep_pdgId.push_back(genLep_pdgId[i]);

      tree_.genMom_pt.push_back(genMom_pt[i]);
      tree_.genMom_mass.push_back(genMom_mass[i]);
      tree_.genMom_eta.push_back(genMom_eta[i]);
      tree_.genMom_phi.push_back(genMom_phi[i]);
      tree_.genMom_pdgId.push_back(genMom_pdgId[i]);
}





tree->Fill();
gStyle->SetOptStat(0);
//HEEP->Sumw2();
//Electron->Sumw2();
HEEP->Draw();
c->SaveAs("LQ_HEEP_bits.png");
HEEP_single->Draw();
c->SaveAs("LQ_HEEP_singlebits.png");
Ratio->Divide(HEEP_pt,Electron,1.,1.,"B");
Ratio->Draw("PE");
Ratio->SetMarkerSize(2.);
Ratio->SetTitle("LQ_HEEP_efficiency");
Ratio->SetAxisRange(0.,1.1, "Y");
Ratio->GetYaxis()->SetTitle("Efficiency");
Ratio->GetYaxis()->SetTitleSize(20);
Ratio->GetYaxis()->SetNdivisions(520);
Ratio->GetYaxis()->SetTitleFont(43);
Ratio->GetYaxis()->SetTitleOffset(1.55);
Ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
Ratio->GetYaxis()->SetLabelSize(15);


Ratio->GetXaxis()->SetTitleSize(20);
Ratio->GetXaxis()->SetTitleFont(43);
Ratio->GetXaxis()->SetTitle("pt (GeV)");
Ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
Ratio->GetXaxis()->SetLabelSize(15);
c->SaveAs("LQ_HEEP_efficiency.png");

}


// ------------ method called once each job just before starting event loop  ------------
void
MiniAODAnalyzer::beginJob()
{


  double bin[26]={0., 35.,100.,200.,300.,400.,500.,600.,700.,800.,900.,1000.,1100.,1200.,1300.,1400.,1500.,1600.,1700.,1800.,1900.,2000.,2250.,2500.,2750.,3000.};
  HEEP_pt = fs->make<TH1F>("LQ_HEEP_pt","LQ_HEEP_pt", 25,bin);
  Electron = fs->make<TH1F>("Electron_Reco","Electron_Reco", 25,bin);
  Ratio = fs->make<TH1F>("Ratio","Ratio", 25,bin);
  c=fs->make<TCanvas>("canvas","canvas",900,600);
  HEEP = fs->make<TH1F>("LQ_HEEP cuts","LQ_HEEP cuts", 13,0,13);
  HEEP_single = fs->make<TH1F>("LQ_HEEP singlecuts","LQ_HEEP singleuts", 13,0,13);


  tree = fs->make<TTree>("tree","tree");
  tree->Branch("run", &tree_.run, "run/I");
  tree->Branch("event", &tree_.event, "event/I");
  tree->Branch("lumi", &tree_.lumi, "lumi/I");
  tree->Branch("nEle", &tree_.nEle, "nEle/I");
  tree->Branch("nMu", &tree_.nMu, "nMu/I");
  tree->Branch("TQ_genMass", &tree_.TQ_genMass, "TQ_genMass/F");
  tree->Branch("TQ_recoMass", &tree_.TQ_recoMass, "TQ_recoMass/F");
  tree->Branch("genLep_pt",     &tree_.genLep_pt);
  tree->Branch("genLep_mass",      &tree_.genLep_mass);
  tree->Branch("genLep_eta",&tree_.genLep_eta);
  tree->Branch("genLep_phi",&tree_.genLep_phi);
  tree->Branch("genLep_pdgId",&tree_.genLep_pdgId);
  tree->Branch("genMom_pt",     &tree_.genMom_pt);
  tree->Branch("genMom_mass",      &tree_.genMom_mass);
  tree->Branch("genMom_eta",&tree_.genMom_eta);
  tree->Branch("genMom_phi",&tree_.genMom_phi);
  tree->Branch("genMom_pdgId",&tree_.genMom_pdgId);
  tree->Branch("recoMu_pt",     &tree_.recoMu_pt);
  tree->Branch("recoMu_mass",      &tree_.recoMu_mass);
  tree->Branch("recoMu_eta",&tree_.recoMu_eta);
  tree->Branch("recoMu_phi",&tree_.recoMu_phi);
  tree->Branch("recoMu_dR",&tree_.recoMu_dR);
  tree->Branch("recoEle_pt",     &tree_.recoEle_pt);
  tree->Branch("recoEle_mass",      &tree_.recoEle_mass);
  tree->Branch("recoEle_eta",&tree_.recoEle_eta);
  tree->Branch("recoEle_phi",&tree_.recoEle_phi);
  tree->Branch("recoEle_dR",&tree_.recoEle_dR);

}

// ------------ method called once each job just after ending the event loop  ------------
void
MiniAODAnalyzer::endJob()
{


}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}


void MiniAODAnalyzer::initTreeStructure()
{
  
  tree_.run= -500.;
  tree_.lumi= -500.;
  tree_.event= 0.;


}

void MiniAODAnalyzer::clearVectors()
{

  tree_.genLep_pt.clear();
  tree_.genLep_mass.clear();
  tree_.genLep_eta.clear();
  tree_.genLep_phi.clear();
  tree_.genLep_pdgId.clear();
  tree_.genMom_pt.clear();
  tree_.genMom_mass.clear();
  tree_.genMom_eta.clear();
  tree_.genMom_phi.clear();
  tree_.genMom_pdgId.clear();

  tree_.recoMu_pt.clear();
  tree_.recoMu_mass.clear();
  tree_.recoMu_eta.clear();
  tree_.recoMu_phi.clear();
  tree_.recoMu_dR.clear();


  tree_.recoEle_pt.clear();
  tree_.recoEle_mass.clear();
  tree_.recoEle_eta.clear();
  tree_.recoEle_phi.clear();
  tree_.recoEle_dR.clear();


}


//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODAnalyzer);
