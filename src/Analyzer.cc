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

// Centrality
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"

// Track
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

// Vertex 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Muon
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"

// Photon
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"

// PF
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

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

using namespace edm;
using namespace reco;
using namespace std;

class Analyzer : public EDAnalyzer
{
public:
   explicit Analyzer(const ParameterSet &);
   ~Analyzer();

   static void fillDescriptions(ConfigurationDescriptions &descriptions);

private:
   virtual void beginJob();
   virtual void analyze(const Event &, const EventSetup &);
   virtual void endJob();

   virtual void beginRun(const Run &, const EventSetup &);
   virtual void endRun(const Run &, const EventSetup &);
   virtual void beginLuminosityBlock(const LuminosityBlock &, const EventSetup &);
   virtual void endLuminosityBlock(const LuminosityBlock &, const EventSetup &);

   bool FillEvent(const Event &iEvent);
   bool FillTrack(const Handle<TrackCollection> &Tracks, const VertexCollection::const_iterator &PV);
   bool FillMuon(const Handle<TrackCollection> &muons, const VertexCollection::const_iterator &pv);
   bool FillPhoton(const Handle<PhotonCollection> &Photons);
   bool FillParticleFlow(const Handle<PFCandidateCollection> &PF);
   bool FillPrimaryVertex(const Handle<VertexCollection> &Vertex);
   const Candidate* GetFinalState(const Candidate* particle, const int id);
   void FillFourMomentum(const Candidate* particle, float* p);
   void InitBranchVars();

   InputTag InputTagCentrality;
   InputTag InputTagTracks;
   InputTag InputTagMuons;
   InputTag InputTagElectrons;
   InputTag InputTagBtags;
   InputTag InputTagPrimaryVertex;
   InputTag InputTagPhotons;
   InputTag InputTagPF;

   int FlagMC;
   int FlagRECO;
   int FlagGEN;
   int NEvents;
   int NEventsSelected;
   int SignLeptonP;
   int SignLeptonM;

   TFile* File;
   TTree* Tree;

   int RunNumber;
   int EventNumber;
   float Centrality;
   
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
   
   static const int MaxNPhoton = 1000;
   int NPhoton;
   float PhotonPt[MaxNPhoton];
   float PhotonEta[MaxNPhoton];
   float PhotonPhi[MaxNPhoton];
   float PhotonEcalR03[MaxNPhoton];
   float PhotonHcalR03[MaxNPhoton];
   float PhotonTrackR03[MaxNPhoton];
   float PhotonNTrackR03[MaxNPhoton];
   float PhotonPFChargedIsolation[MaxNPhoton];
   float PhotonPFNeutralIsolation[MaxNPhoton];
   float PhotonPFPhotonIsolation[MaxNPhoton];
   float PhotonE15[MaxNPhoton];
   float PhotonE25[MaxNPhoton];
   float PhotonE33[MaxNPhoton];
   float PhotonE55[MaxNPhoton];
   float PhotonSigmaEtaEta[MaxNPhoton];
   float PhotonSigmaIEtaIEta[MaxNPhoton];
   float PhotonR9[MaxNPhoton];
   float PhotonR15[MaxNPhoton];
   float PhotonR25[MaxNPhoton];
   float PhotonHOverE[MaxNPhoton];
   
   static const int MaxNPF = 10000;
   int NPF;
   float PFPt[MaxNPF];
   float PFEta[MaxNPF];
   float PFPhi[MaxNPF];
   int PFID[MaxNPF];
   int NPV;
   int PVNDOF;
   float PVZ;
   float PVRho;
};

Analyzer::Analyzer(const ParameterSet &iConfig)
{
   // for proper log files writing (immediate output)
   setbuf(stdout, NULL);

   // input tags
   InputTagCentrality    = InputTag("hiCentrality");
   InputTagMuons         = InputTag("globalMuons");
   InputTagTracks        = InputTag("generalTracks");
   InputTagPrimaryVertex = InputTag("offlinePrimaryVertices");
   InputTagPhotons       = InputTag("photons");
   InputTagPF            = InputTag("particleFlow");
   // vertex input tag used for pp collisions = "offlinePrimaryVertices"
   // PbPb might be "hiSelectedVertex"

   // read configuration parameters
   FlagMC = 0;   // iConfig.getParameter<int>("mc"); // true for MC, false for data
   FlagRECO = 1; // iConfig.getParameter<int>("reco"); // if true, RECO level processed
   FlagGEN = 0;  // iConfig.getParameter<int>("gen"); // if true, generator level processed (works only for MC)
   NEvents = 0;  // number of processed events
   NEventsSelected = 0; // number of selected events
   Service<TFileService> FileService;
   Tree = FileService->make<TTree>("EventObjects", "EventObjects"); //make output tree

   Tree->Branch("RunNumber",   &RunNumber,   "RunNumber/I");
   Tree->Branch("EventNumber", &EventNumber, "EventNumber/I");
   Tree->Branch("Centrality",  &Centrality,  "Centrality/F");

   if(FlagRECO)
   {
      Tree->Branch("NTrack",  &NTrack, "NTrack/I");
      Tree->Branch("TrackPt", TrackPt, "TrackPt[NTrack]/F");
      Tree->Branch("TrackEta", TrackEta, "TrackEta[NTrack]/F");
      Tree->Branch("TrackPhi", TrackPhi, "TrackPhi[NTrack]/F");
      Tree->Branch("TrackCharge", TrackCharge, "TrackCharge[NTrack]/I");
      Tree->Branch("TrackHitsValid", TrackHitsValid, "TrackHitsValid[NTrack]/I");
      Tree->Branch("TrackHitsPixel", TrackHitsPixel, "TrackHitsPixel[NTrack]/I");
      Tree->Branch("TrackDistPV0", TrackDistPV0, "TrackDistPV0[NTrack]/F");
      Tree->Branch("TrackDistPVz", TrackDistPVz, "TrackDistPVz[NTrack]/F");
      Tree->Branch("TrackChi2NDOF", TrackChi2NDOF, "TrackChi2NDOF[NTrack]/F");
      
      Tree->Branch("NMu",  &NMu, "NMu/I");
      Tree->Branch("MuPt", MuPt, "MuPt[NMu]/F");
      Tree->Branch("MuEta", MuEta, "MuEta[NMu]/F");
      Tree->Branch("MuPhi", MuPhi, "MuPhi[NMu]/F");
      Tree->Branch("MuC", MuC, "MuC[NMu]/F");
      Tree->Branch("MuIso03", MuIso03, "MuIso03[NMu]/F");
      Tree->Branch("MuIso04", MuIso04, "MuIso04[NMu]/F");
      Tree->Branch("MuHitsValid", MuHitsValid, "MuHitsValid[NMu]/I");
      Tree->Branch("MuHitsPixel", MuHitsPixel, "MuHitsPixel[NMu]/I");
      Tree->Branch("MuDistPV0", MuDistPV0, "MuDistPV0[NMu]/F");
      Tree->Branch("MuDistPVz", MuDistPVz, "MuDistPVz[NMu]/F");
      Tree->Branch("MuTrackChi2NDOF", MuTrackChi2NDOF, "MuTrackChi2NDOF[NMu]/F");

      Tree->Branch("NPhoton", &NPhoton, "NPhoton/I");
      Tree->Branch("PhotonPt", PhotonPt, "PhotonPt[NPhoton]/F");
      Tree->Branch("PhotonEta", PhotonEta, "PhotonEta[NPhoton]/F");
      Tree->Branch("PhotonPhi", PhotonPhi, "PhotonPhi[NPhoton]/F");
      Tree->Branch("PhotonEcalR03", PhotonEcalR03, "PhotonEcalR03[NPhoton]/F");
      Tree->Branch("PhotonHcalR03", PhotonHcalR03, "PhotonHcalR03[NPhoton]/F");
      Tree->Branch("PhotonTrackR03", PhotonTrackR03, "PhotonTrackR03[NPhoton]/F");
      Tree->Branch("PhotonNTrackR03", PhotonNTrackR03, "PhotonNTrackR03[NPhoton]/F");
      Tree->Branch("PhotonPFChargedIsolation", PhotonPFChargedIsolation, "PhotonPFChargedIsolation[NPhoton]/F");
      Tree->Branch("PhotonPFNeutralIsolation", PhotonPFNeutralIsolation, "PhotonPFNeutralIsolation[NPhoton]/F");
      Tree->Branch("PhotonPFPhotonIsolation", PhotonPFPhotonIsolation, "PhotonPFPhotonIsolation[NPhoton]/F");
      Tree->Branch("PhotonE15", PhotonE15, "PhotonE15[NPhoton]/F");
      Tree->Branch("PhotonE25", PhotonE25, "PhotonE25[NPhoton]/F");
      Tree->Branch("PhotonE33", PhotonE33, "PhotonE33[NPhoton]/F");
      Tree->Branch("PhotonE55", PhotonE55, "PhotonE55[NPhoton]/F");
      Tree->Branch("PhotonSigmaEtaEta", PhotonSigmaEtaEta, "PhotonSigmaEtaEta[NPhoton]/F");
      Tree->Branch("PhotonSigmaIEtaIEta", PhotonSigmaIEtaIEta, "PhotonSigmaIEtaIEta[NPhoton]/F");
      Tree->Branch("PhotonR9", PhotonR9, "PhotonR9[NPhoton]/F");
      Tree->Branch("PhotonR15", PhotonR15, "PhotonR15[NPhoton]/F");
      Tree->Branch("PhotonR25", PhotonR25, "PhotonR25[NPhoton]/F");
      Tree->Branch("PhotonHOverE", PhotonHOverE, "PhotonHOverE[NPhoton]/F");

      Tree->Branch("NPF", &NPF, "NPF/I");
      Tree->Branch("PFPt", PFPt, "PFPt[NPF]/F");
      Tree->Branch("PFEta", PFEta, "PFEta[NPF]/F");
      Tree->Branch("PFPhi", PFPhi, "PFPhi[NPF]/F");
      Tree->Branch("PFID", PFID, "PFID[NPF]/I");
      
      Tree->Branch("NPV",  &NPV, "NPV/I");
      Tree->Branch("PVNDOF",  &PVNDOF, "PVNDOF/I");
      Tree->Branch("PVZ",  &PVZ, "PVZ/F");
      Tree->Branch("PVRho",  &PVRho, "PVRho/F");
   }
}

Analyzer::~Analyzer()
{
}

// initialise event variables with needed default (zero) values; called in the beginning of each event
void Analyzer::InitBranchVars()
{
   RunNumber = 0;
   EventNumber = 0;
   Centrality = -1;
   NTrack = 0;
   NMu = 0;
   NPhoton = 0;
   NPF = 0;
   NPV = 0;
   PVNDOF = 0;
   PVZ = 0;
   PVRho = 0;
}

// Store event info (fill corresponding tree variables)
bool Analyzer::FillEvent(const Event &iEvent)
{
   RunNumber = iEvent.id().run();
   EventNumber = iEvent.id().event();

   edm::Handle<int> hCBin;
   iEvent.getByLabel(InputTagCentralityBin, hCBin);
   Centrality = (*hCBin) * 0.5;

   return true;
}

bool Analyzer::FillTrack(const Handle<TrackCollection> &Tracks,
   const VertexCollection::const_iterator &PV)
{
   NTrack = 0;
   for(TrackCollection::const_iterator iter = Tracks->begin(); iter != Tracks->end(); iter++)
   {
      TrackHitsPixel[NTrack] = 0;
      TrackHitsValid[NTrack] = 0;

      const HitPattern &p = iter->hitPattern();
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
bool Analyzer::FillMuon(const Handle<TrackCollection> &muons,
   const VertexCollection::const_iterator &pv)
{
   using namespace std;
   NMu = 0;
   NMu0 = 0;
   // loop over muons
   for(TrackCollection::const_iterator it = muons->begin(); it != muons->end(); it++)
   {
      NMu0++;
      if(NMu == MaxNmu)
      {
         printf("Maximum number of muons %d reached, skipping the rest\n", MaxNmu);
         return false;
      }
      MuHitsValid[NMu] = 0;
      MuHitsPixel[NMu] = 0;
      const HitPattern &p = it->hitPattern();
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

bool Analyzer::FillPhoton(const Handle<PhotonCollection> &Photons)
{
   NPhoton = Photons->size();

   for(int i = 0; i < NPhoton; i++)
   {
      const Photon &P = Photons->at(i);

      PhotonPt[i] = P.pt();
      PhotonEta[i] = P.eta();
      PhotonPhi[i] = P.phi();
      
      PhotonEcalR03[i] = P.ecalRecHitSumEtConeDR03();
      PhotonHcalR03[i] = P.hcalTowerSumEtConeDR03();
      PhotonTrackR03[i] = P.trkSumPtSolidConeDR03();
      PhotonNTrackR03[i] = P.nTrkSolidConeDR03();
      PhotonPFChargedIsolation[i] = P.chargedHadronIso();
      PhotonPFNeutralIsolation[i] = P.neutralHadronIso();
      PhotonPFPhotonIsolation[i] = P.photonIso();

      PhotonE15[i] = P.e1x5();
      PhotonE25[i] = P.e2x5();
      PhotonE33[i] = P.e3x3();
      PhotonE55[i] = P.e5x5();
      PhotonSigmaEtaEta[i] = P.sigmaEtaEta();
      PhotonSigmaIEtaIEta[i] = P.sigmaIetaIeta();
      PhotonR9[i] = P.r9();
      PhotonR15[i] = P.r1x5();
      PhotonR25[i] = P.r2x5();
      PhotonHOverE[i] = P.hadronicOverEm();
   }

   return true;
}
   
bool Analyzer::FillParticleFlow(const Handle<PFCandidateCollection> &PF)
{
   NPF = PF->size();

   for(int i = 0; i < NPF; i++)
   {
      const PFCandidate &P = PF->at(i);

      PFPt[i] = P.pt();
      PFEta[i] = P.eta();
      PFPhi[i] = P.phi();
      PFID[i] = P.particleId();
   }

   return true;
}

// select primary vertex
bool Analyzer::FillPrimaryVertex(const Handle<VertexCollection> &Vertex)
{
   // if no primary vertices in the event, return false status
   if(Vertex->size() == 0)
      return false;
   // take the first primary vertex
   VertexCollection::const_iterator PV = Vertex->begin();
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
void Analyzer::FillFourMomentum(const Candidate* particle, float* p)
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
void Analyzer::analyze(const Event &iEvent, const EventSetup &iSetup)
{
   NEvents++;
   if(NEvents % 1000 == 0)
   {
      printf("************* NEVENTS = %d K, selected = %d *************\n", NEvents / 1000, NEventsSelected);
   }

   // declare event contents
   Handle<VertexCollection> hVertex;
   Handle<TrackCollection> hMuon;
   Handle<TrackCollection> hTrack;
   Handle<PhotonCollection> hPhoton;
   Handle<PFCandidateCollection> hPF;

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
      VertexCollection::const_iterator Vertex = hVertex->begin();
      FillPrimaryVertex(hVertex);
   
      // tracks
      iEvent.getByLabel(InputTagTracks, hTrack);
      FillTrack(hTrack, Vertex);

      // muons
      iEvent.getByLabel(InputTagMuons, hMuon);
      FillMuon(hMuon, Vertex);

      // photons
      iEvent.getByLabel(InputTagPhotons, hPhoton);
      FillPhoton(hPhoton);

      // PF
      iEvent.getByLabel(InputTagPF, hPF);
      FillParticleFlow(hPF);
   }
   // fill event info
   FillEvent(iEvent);
   // all done: store event
   Tree->Fill();
   NEventsSelected++;
}


// ------------ method called when starting to processes a run  ------------
void Analyzer::beginRun(const Run &iRun, const EventSetup &iSetup)
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
void Analyzer::endRun(const Run &iRun, const EventSetup &iSetup) 
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void Analyzer::beginLuminosityBlock(const LuminosityBlock &iBlock, const EventSetup &iSetup)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void Analyzer::endLuminosityBlock(const LuminosityBlock &iBlock, const EventSetup &iSetup)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Analyzer::fillDescriptions(ConfigurationDescriptions &descriptions)
{
   ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
