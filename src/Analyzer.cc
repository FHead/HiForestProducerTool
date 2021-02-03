#include <memory>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"

// Tracking
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

// Vertex 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"

// beamspot
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// Trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// ROOT
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

class Analyzer : public edm::EDAnalyzer
{
public:
   explicit Analyzer(const edm::ParameterSet &);
   ~Analyzer();

   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
   virtual void beginJob();
   virtual void analyze(const edm::Event &, const edm::EventSetup &);
   virtual void endJob();

   virtual void beginRun(const edm::Run &, const edm::EventSetup &);
   virtual void endRun(const edm::Run &, const edm::EventSetup &);
   virtual void beginLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup &);
   virtual void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup &);

   // user routines (detailed description given with the method implementations)
   bool FillEvent(const edm::Event &iEvent);
   bool FillTrack(const edm::Handle<reco::TrackCollection> &Tracks, const reco::VertexCollection::const_iterator &PV);
   bool FillMuon(const edm::Handle<reco::TrackCollection> &muons, const reco::VertexCollection::const_iterator &pv);
   bool FillPrimaryVertex(const edm::Handle<reco::VertexCollection> &Vertex);
   const reco::Candidate* GetFinalState(const reco::Candidate* particle, const int id);
   void FillFourMomentum(const reco::Candidate* particle, float* p);
   void InitBranchVars();

   // input tags
   edm::InputTag InputTagTracks;
   edm::InputTag InputTagMuons;
   edm::InputTag InputTagElectrons;
   edm::InputTag InputTagBtags;
   edm::InputTag InputTagPrimaryVertex;

   // general flags and variables
   int FlagMC;
   int FlagRECO;
   int FlagGEN;
   int NEvents;
   int NEventsSelected;
   int SignLeptonP;
   int SignLeptonM;

   // storage
   TFile* File;
   TTree* Tree;

   // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   // >>>>>>>>>>>>>>>> event variables >>>>>>>>>>>>>>>>>>>>>>>
   // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   // (their description given when tree branches are created)
   // event
   int RunNumber;
   int EventNumber;
   // tracks
   static const int MaxNTrack = 10000;
   int NTrack;
   float TrackPt[MaxNTrack];
   float TrackEta[MaxNTrack];
   float TrackPhi[MaxNTrack];
   int TrackCharge[MaxNTrack];
   int TrackHitsValid[MaxNTrack];
   int TrackHitsPixel[MaxNTrack];
   float TrackDistPV0[MaxNTrack];
   float TrackDistPVz[MaxNTrack];
   float TrackChi2NDOF[MaxNTrack];
   // muons
   static const int MaxNmu = 10;
   int NMu;
   int NMu0;
   float MuPt[MaxNmu];
   float MuEta[MaxNmu];
   float MuPhi[MaxNmu];
   float MuC[MaxNmu];
   float MuIso03[MaxNmu];
   float MuIso04[MaxNmu];
   int MuHitsValid[MaxNmu];
   int MuHitsPixel[MaxNmu];
   float MuDistPV0[MaxNmu];
   float MuDistPVz[MaxNmu];
   float MuTrackChi2NDOF[MaxNmu];
   // primary vertex
   int NPV;
   int PVNDOF;
   float PVZ;
   float PVRho;
};

Analyzer::Analyzer(const edm::ParameterSet &iConfig)
{
   // for proper log files writing (immediate output)
   setbuf(stdout, NULL);

   // input tags
   InputTagMuons = edm::InputTag("globalMuons");
   InputTagTracks = edm::InputTag("generalTracks");
   InputTagPrimaryVertex = edm::InputTag("offlinePrimaryVertices");
   // vertex input tag used for pp collisions = "offlinePrimaryVertices"
   // PbPb might be "hiSelectedVertex"

   // read configuration parameters
   FlagMC = 0;   // iConfig.getParameter<int>("mc"); // true for MC, false for data
   FlagRECO = 1; // iConfig.getParameter<int>("reco"); // if true, RECO level processed
   FlagGEN = 0;  // iConfig.getParameter<int>("gen"); // if true, generator level processed (works only for MC)
   NEvents = 0;  // number of processed events
   NEventsSelected = 0; // number of selected events
   edm::Service<TFileService> FileService;
   Tree = FileService->make<TTree>("EventObjects", "EventObjects"); //make output tree

   // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   // >>>>>>> tree branches >>>>>>>>>>>>
   // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   //
   // event
   Tree->Branch("RunNumber",  &RunNumber, "RunNumber/I"); // run number
   Tree->Branch("EventNumber",  &EventNumber, "EventNumber/I"); // event number

   if(FlagRECO)
   {
      // tracks
      Tree->Branch("NTrack",  &NTrack, "NTrack/I");   // number of tracks
      Tree->Branch("TrackPt", TrackPt, "TrackPt[NTrack]/F");
      Tree->Branch("TrackEta", TrackEta, "TrackEta[NTrack]/F");
      Tree->Branch("TrackPhi", TrackPhi, "TrackPhi[NTrack]/F");
      Tree->Branch("TrackCharge", TrackCharge, "TrackCharge[NTrack]/I");
      Tree->Branch("TrackHitsValid", TrackHitsValid, "TrackHitsValid[NTrack]/I");
      Tree->Branch("TrackHitsPixel", TrackHitsPixel, "TrackHitsPixel[NTrack]/I");
      Tree->Branch("TrackDistPV0", TrackDistPV0, "TrackDistPV0[NTrack]/F");
      Tree->Branch("TrackDistPVz", TrackDistPVz, "TrackDistPVz[NTrack]/F");
      Tree->Branch("TrackChi2NDOF", TrackChi2NDOF, "TrackChi2NDOF[NTrack]/F");
      // muons
      Tree->Branch("NMu",  &NMu, "NMu/I"); // number of Muons 
      Tree->Branch("MuPt", MuPt, "MuPt[NMu]/F"); // Muon pT
      Tree->Branch("MuEta", MuEta, "MuEta[NMu]/F"); // Muon pseudorapidity
      Tree->Branch("MuPhi", MuPhi, "MuPhi[NMu]/F"); // Muon phi
      Tree->Branch("MuC", MuC, "MuC[NMu]/F"); // Muon phi
      Tree->Branch("MuIso03", MuIso03, "MuIso03[NMu]/F"); // Muon isolation, delta_R = 0.3
      Tree->Branch("MuIso04", MuIso04, "MuIso04[NMu]/F"); // Muon isolation, delta_R = 0.4
      Tree->Branch("MuHitsValid", MuHitsValid, "MuHitsValid[NMu]/I"); // Muon valid hits number
      Tree->Branch("MuHitsPixel", MuHitsPixel, "MuHitsPixel[NMu]/I"); // Muon pixel hits number
      Tree->Branch("MuDistPV0", MuDistPV0, "MuDistPV0[NMu]/F"); // Muon distance to the primary vertex (projection on transverse plane)
      Tree->Branch("MuDistPVz", MuDistPVz, "MuDistPVz[NMu]/F"); // Muon distance to the primary vertex (z projection)
      Tree->Branch("MuTrackChi2NDOF", MuTrackChi2NDOF, "MuTrackChi2NDOF[NMu]/F"); // Muon track number of degrees of freedom
      // primary vertex
      Tree->Branch("NPV",  &NPV, "NPV/I"); // total number of primary vertices
      Tree->Branch("PVNDOF",  &PVNDOF, "PVNDOF/I"); // number of degrees of freedom of the primary vertex
      Tree->Branch("PVZ",  &PVZ, "PVZ/F"); // z component of the primary vertex
      Tree->Branch("PVRho",  &PVRho, "PVRho/F"); // rho of the primary vertex (projection on transverse plane)
   }

}

// destructor
Analyzer::~Analyzer()
{
}

// initialise event variables with needed default (zero) values; called in the beginning of each event
void Analyzer::InitBranchVars()
{
   RunNumber = 0;
   EventNumber = 0;
   NTrack = 0;
   NMu = 0;
   NPV = 0;
   PVNDOF = 0;
   PVZ = 0;
   PVRho = 0;
}

// Store event info (fill corresponding tree variables)
bool Analyzer::FillEvent(const edm::Event &iEvent)
{
   RunNumber = iEvent.id().run();
   EventNumber = iEvent.id().event();
   return true;
}

bool Analyzer::FillTrack(const edm::Handle<reco::TrackCollection> &Tracks,
   const reco::VertexCollection::const_iterator &PV)
{
   NTrack = 0;
   for(reco::TrackCollection::const_iterator iter = Tracks->begin(); iter != Tracks->end(); iter++)
   {
      TrackHitsPixel[NTrack] = 0;
      TrackHitsValid[NTrack] = 0;

      const reco::HitPattern &p = iter->hitPattern();
      for(int i = 0; i < p.numberOfHits(); i++) 
      {
         uint32_t hit = p.getHitPattern(i);
         if(p.validHitFilter(hit) && p.pixelHitFilter(hit))
            TrackHitsPixel[NTrack]++;
         if(p.validHitFilter(hit))
            TrackHitsValid[NTrack]++;
      }
      
      TrackPt[NTrack] = iter->pt();
      TrackEta[NTrack] = iter->eta();
      TrackPhi[NTrack] = iter->phi();
      TrackCharge[NTrack] = iter->charge();
      
      if(iter->ndof())
         TrackChi2NDOF[NTrack] = iter->chi2() / iter->ndof();
      
      double DX = PV->x() - iter->vx();
      double DY = PV->y() - iter->vy();
      TrackDistPV0[NTrack] = sqrt(DX * DX + DY * DY);
      TrackDistPVz[NTrack] = fabs(PV->z() - iter->vz());
      
      NTrack = NTrack + 1;
   }

   return true;
}

// muon selection
bool Analyzer::FillMuon(const edm::Handle<reco::TrackCollection> &muons,
   const reco::VertexCollection::const_iterator &pv)
{
   using namespace std;
   NMu = 0;
   NMu0 = 0;
   // loop over muons
   for(reco::TrackCollection::const_iterator it = muons->begin(); it != muons->end(); it++)
   {
      NMu0++;
      if(NMu == MaxNmu)
      {
         printf("Maximum number of muons %d reached, skipping the rest\n", MaxNmu);
         return false;
      }
      MuHitsValid[NMu] = 0;
      MuHitsPixel[NMu] = 0;
      const reco::HitPattern &p = it->hitPattern();
      for(int i = 0; i < p.numberOfHits(); i++) 
      {
         uint32_t hit = p.getHitPattern(i);
         if(p.validHitFilter(hit) && p.pixelHitFilter(hit))
            MuHitsPixel[NMu]++;
         if(p.validHitFilter(hit))
            MuHitsValid[NMu]++;
      }
      // fill three momentum (pT, eta, phi)
      MuPt[NMu] = it->pt();// * it->charge();
      MuEta[NMu] = it->eta();
      MuPhi[NMu] = it->phi();
      MuC[NMu] = it->charge();
      // fill chi2/ndof
      if(it->ndof()) MuTrackChi2NDOF[NMu] = it->chi2() / it->ndof();
      // fill distance to primary vertex
      MuDistPV0[NMu] = TMath::Sqrt(TMath::Power(pv->x() - it->vx(), 2.0) + TMath::Power(pv->y() - it->vy(), 2.0));
      MuDistPVz[NMu] = TMath::Abs(pv->z() - it->vz());
      // store muon
      NMu++;
      // determine muon sign (in the end the event will be stored only there are opposite signed leptons)
      if(it->charge() == +1)
         SignLeptonP = 1;
      if(it->charge() == -1)
         SignLeptonM = 1;
   }
   cout<<"Muons before selection: "<<NMu0<<endl;
   return true;
}

// select primary vertex
bool Analyzer::FillPrimaryVertex(const edm::Handle<reco::VertexCollection> &Vertex)
{
   // if no primary vertices in the event, return false status
   if(Vertex->size() == 0)
      return false;
   // take the first primary vertex
   reco::VertexCollection::const_iterator PV = Vertex->begin();
   // fill z and rho (projection on transverse plane)
   PVZ = PV->z();
   PVRho = TMath::Sqrt(TMath::Power(PV->x(), 2.0) + TMath::Power(PV->y(), 2.0));
   // fill number of primary veritces
   NPV = Vertex->size();
   // fill number of degrees of freedom
   PVNDOF = PV->ndof();
   // return true status
   return true;
}

// fill 4-momentum (p) with provided particle pointer
void Analyzer::FillFourMomentum(const reco::Candidate* particle, float* p)
{
   // if NULL pointer provided, initialise with default (zero)
   if(particle == NULL)
   {
      p[0] = p[1] = p[2] = p[3] = 0.0;
      return;
   }

   p[0] = particle->px();
   p[1] = particle->py();
   p[2] = particle->pz();
   p[3] = particle->mass();
}

// select MC generator level information
// (analysis specific ttbar dileptonic decay)

// ------------ method called for each event  ------------
void Analyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;

   NEvents++;
   if(NEvents % 1000 == 0)
   {
      printf("************* NEVENTS = %d K, selected = %d *************\n", NEvents / 1000, NEventsSelected);
   }

   // declare event contents
   Handle<reco::VertexCollection> hVertex;
   edm::Handle<reco::TrackCollection> hMuon;
   edm::Handle<reco::TrackCollection> hTrack;

   // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   // >>>>>>>>> event selection >>>>>>>>>
   // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   //
   // initialise event variables with default values
   InitBranchVars();
   // process generator level, if needed
   // process reco level, if needed
   if(FlagRECO)
   {
      // primary vertex
      iEvent.getByLabel(InputTagPrimaryVertex, hVertex);
      reco::VertexCollection::const_iterator Vertex = hVertex->begin();
      FillPrimaryVertex(hVertex);
   
      // tracks
      iEvent.getByLabel(InputTagTracks, hTrack);
      FillTrack(hTrack, Vertex);

      // muons
      iEvent.getByLabel(InputTagMuons, hMuon);
      FillMuon(hMuon, Vertex);
   }
   // fill event info
   FillEvent(iEvent);
   // all done: store event
   Tree->Fill();
   NEventsSelected++;
}


// ------------ method called when starting to processes a run  ------------
void Analyzer::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup)
{
}

// below is some default stuff, was not modified

// ------------ method called once each job just before starting event loop  ------------
void Analyzer::beginJob() 
{
}

// ------------ method called once each job just after ending the event loop  ------------
void Analyzer::endJob() 
{
}

// ------------ method called when ending the processing of a run  ------------
void Analyzer::endRun(const edm::Run &iRun, const edm::EventSetup &iSetup) 
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void Analyzer::beginLuminosityBlock(const edm::LuminosityBlock &iBlock, const edm::EventSetup &iSetup)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void Analyzer::endLuminosityBlock(const edm::LuminosityBlock &iBlock, const edm::EventSetup &iSetup)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Analyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
