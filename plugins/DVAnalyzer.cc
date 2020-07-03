// -*- C++ -*-
//
// Package:    LLPAnalyzer/DVAnalyzer
// Class:      DVAnalyzer
//
/**\class DVAnalyzer DVAnalyzer.cc LLPAnalyzer/DVAnalyzer/plugins/DVAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ang Li
//         Created:  Thu, 11 Jun 2020 14:08:05 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "TH1.h"
#include "TH2.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;
using reco::TransientTrack;

typedef std::set<reco::TrackRef> track_set;
typedef std::vector<reco::TrackRef> track_vec;

class DVAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>  {
   public:
      explicit DVAnalyzer(const edm::ParameterSet&);
      ~DVAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void beginRun(edm::Run const& , edm::EventSetup const&) override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endRun(edm::Run const& , edm::EventSetup const&) override;
      virtual void endJob() override;
      std::vector<TransientVertex> refitVtx(std::vector<reco::TransientTrack>& ttks);
      reco::TrackRefVector SharedTracks(track_set t1, track_set t2);
      std::vector<TransientTrack> combineTracks(track_set t1, track_set t2, std::map<reco::TrackRef, size_t> seed_track_ref_map, std::vector<TransientTrack> tTracks);
      track_set vertexTracks(const reco::Vertex v) const;
      track_vec vertex_track_vec(const reco::Vertex& v) const;
      bool is_track_subset(const track_set& a, const track_set& b) const;


      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::BeamSpot> bsToken_;
      edm::Handle<reco::BeamSpot> bsHandle_;
      edm::EDGetTokenT<std::vector<pat::Jet>> jetToken_;
      edm::Handle<std::vector<pat::Jet>> jetHandle_;
      edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
      edm::Handle<std::vector<pat::Muon>> muonHandle_;
      //edm::EDGetTokenT<std::vector<pat::PackedCandidate>> tracksToken_;
      //edm::Handle<std::vector<pat::PackedCandidate>> tracksHandle_;
      edm::EDGetTokenT<reco::TrackRefVector> trackRefToken_;
      edm::Handle<reco::TrackRefVector> trackRefHandle_;
      edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
      edm::Handle<std::vector<reco::Vertex>> vtxHandle_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
      edm::Handle<edm::TriggerResults> triggerResultsHandle_;

      HLTConfigProvider hltConfig_;

      std::string processName_;
      std::string triggerName_;

      std::unique_ptr<KalmanVertexFitter> fitter;

      TH1F* h_PFJet_PT_;
      TH1F* h_PFJet_eta_;
      TH1F* h_PFJet_HT_;
      TH1F* h_PFJet_N_;
      TH1F* h_muon_PT_;
      TH1F* h_muon_eta_;
      TH1F* h_track_PT_;
      TH1* h_Event_cutflow_;
      TH1* h_PFJet_cutflow_;
      TH1* h_Track_cutflow_;
      TH1* h_NTrackPerVtx_;
      TH1* h_dBV_;
      TH1* h_sigma_dBV_;
      TH2* h_vtx_XY_;
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
DVAnalyzer::DVAnalyzer(const edm::ParameterSet& iConfig)
 :
  bsToken_(consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamspotTag"))),
  jetToken_(consumes<std::vector<pat::Jet>>(iConfig.getUntrackedParameter<edm::InputTag>("jetTag"))),
  muonToken_(consumes<std::vector<pat::Muon>>(iConfig.getUntrackedParameter<edm::InputTag>("muonTag"))),
  //tracksToken_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("tracksTag"))),
  trackRefToken_(consumes<reco::TrackRefVector>(iConfig.getUntrackedParameter<edm::InputTag>("trackRefTag"))),
  vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getUntrackedParameter<edm::InputTag>("vertexTag"))),
  triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("triggerTag"))),
  processName_(iConfig.getUntrackedParameter<std::string>("processName")),
  triggerName_(iConfig.getUntrackedParameter<std::string>("triggerNameTag")),
  fitter(new KalmanVertexFitter())

{
   //now do what ever initialization is needed

}


DVAnalyzer::~DVAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DVAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // assign approprate Handles and get the object collections
  //iEvent.getByToken(genJetToken_, genJetHandle_);

  iEvent.getByToken(bsToken_, bsHandle_);
  reco::BeamSpot bs = *bsHandle_.product();
  VertexState theBeam = VertexState(bs);

  iEvent.getByToken(jetToken_, jetHandle_);
  auto jets = *jetHandle_.product();

  iEvent.getByToken(muonToken_, muonHandle_);
  auto muons = *muonHandle_.product();
  
  //iEvent.getByToken(tracksToken_, tracksHandle_);
  //auto tracks = *tracksHandle_.product();

  iEvent.getByToken(trackRefToken_, trackRefHandle_);
  auto tRef = *trackRefHandle_.product();

  iEvent.getByToken(vtxToken_, vtxHandle_);
  const auto& vertex = (*vtxHandle_)[0];

  // load Trigger paths
  iEvent.getByToken(triggerResultsToken_, triggerResultsHandle_);

  if (!triggerResultsHandle_.isValid()) {
    std::cout << "****Error in getting TriggerResults product from Event!" << std::endl;
    return;
  }

  assert(triggerResultsHandle_->size()==hltConfig_.size());
  const unsigned int ntrigs(hltConfig_.size());
  const unsigned int triggerIdx(hltConfig_.triggerIndex(triggerName_));
  assert(triggerIdx==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName_));
  if(triggerIdx>=ntrigs){
    std::cout << "Trigger path: " << triggerName_ << " not available! " << std::endl;
  }
  bool accept = triggerResultsHandle_->accept(triggerIdx);
  if(!accept){
    return;
  }

  // get transient track builder
  //ESHandle<TransientTrackBuilder> theB;
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  // loop over objects in iEvent
  
  //PFJets
  double pfJet_HT = 0;
  int pfJet_N = 0;

  h_Event_cutflow_->Fill("allEvents",1);

  for (const auto & jet:jets){
    h_PFJet_cutflow_->Fill("all jets",1);
    // PF Jet
    if(!jet.isPFJet()) continue;
    h_PFJet_cutflow_->Fill("PFJet",1);
    // pt>20GeV
    if(jet.pt()<=20) continue;
    h_PFJet_cutflow_->Fill("pT<20GeV",1);
    // |eta|<2.5
    if(std::abs(jet.eta())>2.5) continue;
    h_PFJet_cutflow_->Fill("|eta|<2.5",1);
    // # constituents > 1
    if( jet.numberOfDaughters() <= 1 ) continue;
    h_PFJet_cutflow_->Fill("#const>1",1);
    // Neutral hadron energy fraction < 0.9
    if( (jet.neutralHadronEnergyFraction()) > 0.9 ) continue;
    h_PFJet_cutflow_->Fill("f_E_NH<0.9",1);
    // Neutral EM energy fraction <0.9
    if( (jet.neutralEmEnergyFraction()) > 0.9 ) continue;
    h_PFJet_cutflow_->Fill("f_E_EM<0.9",1);
    // muon energy fraction < 0.8
    if( (jet.muonEnergyFraction()) > 0.8 ) continue;
    h_PFJet_cutflow_->Fill("f_E_MU<0.8",1);
    // cuts for |eta|<2.4
    if(std::abs(jet.eta())<2.4){
      // charged hadron energy fraction >0
      if( (jet.chargedHadronEnergyFraction()) <= 0 ) continue;
      h_PFJet_cutflow_->Fill("f_HE_c>0",1);
      // charged multiplicity>0
      if( (jet.chargedMultiplicity()) <= 0 ) continue;
      h_PFJet_cutflow_->Fill("M_c>0",1);
      // charged EM energy fraction < 0.8
      if( (jet.chargedEmEnergyFraction()) > 0.8 ) continue;
      h_PFJet_cutflow_->Fill("f_EM_c>0.8",1);
    }
    h_PFJet_PT_->Fill(jet.pt());
    h_PFJet_eta_->Fill(jet.eta());
    if(jet.pt()>=40){
      pfJet_HT += jet.pt();
    }
    ++pfJet_N;
  }
  // Event preselection:
  // HT > 1200 GeV
  if(pfJet_HT<=1200) return;
  h_Event_cutflow_->Fill("HT>1200", 1);
  h_PFJet_HT_->Fill(pfJet_HT);
  // At least four jet satisfying above criteria
  if(pfJet_N<4) return;
  h_Event_cutflow_->Fill("NJet>=4", 1);
  h_PFJet_N_->Fill(pfJet_N);

  //muons
  for (const auto & muon:muons){
    if(muon.pt()<27) continue;
    if(std::abs(muon.eta())>2.4) continue;
    if(muon::isTightMuon(muon, vertex)){
      h_muon_PT_->Fill(muon.pt());
      h_muon_eta_->Fill(muon.eta());
    }
  }

  //tracks

  std::vector<reco::TransientTrack> tTracks;
  std::map<reco::TrackRef, size_t> seed_track_ref_map;

  for (auto tk = tRef.begin(); tk!=tRef.end();++tk){
    auto track = *tk;

    // track pt > 1 GeV
    if(track->pt()<=1) continue;
    h_Track_cutflow_->Fill("pt>1GeV", 1);
    
    //rmin = 1
    if(!track->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel,1)) continue;
    h_Track_cutflow_->Fill("rmin==1", 1);

    // npxl >= 2
    if(track->hitPattern().pixelLayersWithMeasurement()<2) continue;
    h_Track_cutflow_->Fill("npxl>=2", 1);

    // nstl >= 6
    if(track->hitPattern().stripLayersWithMeasurement() < 6) continue;
    h_Track_cutflow_->Fill("nstl>=6", 1);

    // |dxy|/sigma_dxy > 4 
    if( (std::abs(track->dxy())/track->dxyError()) <= 4 ) continue;
    h_Track_cutflow_->Fill("dxy/sigma_dxy>4", 1);

    h_track_PT_->Fill(track->pt());

    reco::TransientTrack t(track, &(*bFieldHandle));
    tTracks.push_back(t);
    seed_track_ref_map[track] = tTracks.size()-1;
  }

  /*
  for (const auto & track:tracks){
    h_Track_cutflow_->Fill("all tracks", 1);

    // track pt > 1 GeV
    if(track.pt()<=1) continue;
    h_Track_cutflow_->Fill("pt>1GeV", 1);
    
    //rmin = 1
    if(track.lostInnerHits()!=-1) continue;
    h_Track_cutflow_->Fill("rmin==1", 1);

    // npxl >= 2
    if(track.pixelLayersWithMeasurement()<2) continue;
    h_Track_cutflow_->Fill("npxl>=2", 1);

    // nstl >= 6
    if(track.stripLayersWithMeasurement() < 6) continue;
    h_Track_cutflow_->Fill("nstl>=6", 1);

    // |dxy|/sigma_dxy > 4 
    if( (std::abs(track.dxy())/track.dxyError()) <= 4 ) continue;
    h_Track_cutflow_->Fill("dxy/sigma_dxy>4", 1);

    h_track_PT_->Fill(track.pt());

    if(track.hasTrackDetails()){
      reco::TransientTrack t(track.pseudoTrack(), &(*bFieldHandle));
      tTracks.push_back(t);
      seed_track_ref_map[track] = tTracks.size()-1;
      //reco::TransientTrack t = (*theB).build(track.bestTrack());
    }
  }*/

  // refit vertices

  if(tTracks.size()<2) return;
  std::vector<reco::Vertex> v_refit;
  for(auto tTrack = tTracks.cbegin(); tTrack!=tTracks.cend()-1; ++tTrack){
    const auto& tT1 = *tTrack;
    for(auto tTrack2 = tTrack+1; tTrack2!=tTracks.cend(); ++tTrack2){
      const auto& tT2 = *tTrack2;
      std::vector<reco::TransientTrack> tTrackVector({tT1, tT2});
      TransientVertex tv = fitter->vertex(tTrackVector); 
      if(tv.isValid() && (tv.normalisedChiSquared())<5){
        v_refit.push_back(reco::Vertex(tv));
      }
    }
  }
  VertexDistance3D vtxDistance;
  std::vector<reco::Vertex>::iterator v[2];
  for(v[0]=v_refit.begin(); v[0]!=v_refit.end(); ++v[0]) {
    const auto& v1 = *v[0];
    if(!v1.isValid()) {
      continue;
    }
    bool duplicate = false;
    bool merge = false;
    bool refit = false;
    track_set track[2]; //tracks in vertices
    track_set tk_rm[2]; //tracks that are going to be removed from vertices
    track[0] = vertexTracks(v1);
    if(track[0].size()<2) {
      v[0] = v_refit.erase(v[0])-1;
      continue;
    }
    //std::cout << "start with a vtx" << std::endl;
    for(v[1]=v[0]+1; v[1]!=v_refit.end(); ++v[1]){
      const auto& v2 = *v[1];
      if(!v2.isValid()) continue;

      track[1] = vertexTracks(v2);
      if(track[1].size()<2){
        v[1] = v_refit.erase(v[1])-1;
        continue;
      }

      // if one vertex have all tracks included in the other vertex
      // remove the second vertex
      if(is_track_subset(track[0],track[1])){
        duplicate = true;
        break;
      }

      auto sharedTracks = SharedTracks(track[0], track[1]);
      if(sharedTracks.size()==0) continue;
      //std::cout << "# Stracks: " << sharedTracks.size() << std::endl;
      //std::cout << "t1: " << track[0].size() << std::endl;
      //std::cout << "t2: " << track[1].size() << std::endl;
      if(vtxDistance.distance(v1, v2).significance()<4){
        // combine the two vertex tracks and refit
        merge = true;
        //std::cout << "will merge" << std::endl;
      }
      else{
        refit = true;
        //std::cout << "will refit" << std::endl;
      }
      for (auto tkIter=sharedTracks.begin(); tkIter!=sharedTracks.end(); ++tkIter) {
        auto tk = *tkIter;
        TransientTrack ttk = tTracks[seed_track_ref_map[tk]];
        auto ip3D_pair1 = IPTools::absoluteImpactParameter3D(ttk, v1);
        auto ip3D_pair2 = IPTools::absoluteImpactParameter3D(ttk, v2);
        double ip3D_1_sig = ip3D_pair1.second.significance();
        double ip3D_2_sig = ip3D_pair2.second.significance();
        ip3D_pair1.first = ip3D_pair1.first && (ip3D_pair1.second.significance()<5);
        ip3D_pair2.first = ip3D_pair2.first && (ip3D_pair2.second.significance()<5);
        bool rm_from_1 = !ip3D_pair1.first;
        bool rm_from_2 = !ip3D_pair2.first;
        if(ip3D_1_sig<1.5 && ip3D_2_sig<1.5){
          if(v1.tracksSize()>=v2.tracksSize()){
            rm_from_2 = true;
          }
          else{
            rm_from_1 = true;
          }
        }
        else if(ip3D_1_sig<=ip3D_2_sig){
          rm_from_2 = true;
        }
        else{
          rm_from_1 = true;
        }
        if(rm_from_1) tk_rm[0].insert(tk);
        if(rm_from_2) tk_rm[1].insert(tk);
        //break;
      }
      //std::cout << "tks to rm from 0: " << tk_rm[0].size() << std::endl;
      //std::cout << "tks to rm from 1: " << tk_rm[1].size() << std::endl;
      break;
    }
    if(duplicate){
      //std::cout << "duplicate" << std::endl;
      v_refit.erase(v[1]);
    }
    else if(merge){
      //std::cout << "merge"<<track[0].size() <<"  " << track[1].size() << std::endl;
      track_set combinedSet;
      for (auto tk:track[0]){
        combinedSet.insert(tk);
      }
      for (auto tk:track[1]){
        combinedSet.insert(tk);
      }
      //std::cout << "combined: " << combinedSet.size() << std::endl;
      std::vector<TransientTrack> cT;
      for (auto tk:combinedSet){
        cT.push_back(tTracks[seed_track_ref_map[tk]]);
      }
      auto v_new = refitVtx(cT);
      if(v_new.size()==1 && vertexTracks(v_new[0])==combinedSet){
        v_refit.erase(v[1]);
        *v[0] = reco::Vertex(v_new[0]);
      }
      else refit = true;
      //std::cout << "merge done" << std::endl;
    }
    if(refit){
      //std::cout << "refit" << std::endl;
      for (int i = 0; i<2; ++i){
        if(tk_rm[i].empty()) continue;
        std::vector<reco::TransientTrack> ttks;
        for (auto tk:track[i]){
          if(tk_rm[i].count(tk)==0)
            ttks.push_back(tTracks[seed_track_ref_map[tk]]);
        }
        auto v_new = refitVtx(ttks);
        if(v_new.size()==1)
          *v[i] = reco::Vertex(v_new[0]);
        else{
          v_refit.erase(v[i]);
        }
      }
      //std::cout << "refit done" << std::endl;
    }

    // start over if vertices are changed
    if(duplicate || merge || refit){
      //std::cout << "start over" << std::endl;
      v[0] = v_refit.begin()-1;
    }
  }

  for (v[0]=v_refit.begin(); v[0]!=v_refit.end(); ++ v[0]){
    const track_vec tks = vertex_track_vec(*v[0]);
    const size_t ntks = tks.size();
    if(ntks<3)
      continue;

    std::vector<TransientTrack> ttks(ntks-1);
    for(size_t i = 0; i<ntks; ++i){
      for(size_t j = 0; j<ntks; ++j){
        if(j!=i)
          ttks[j-(j>=i)] = TransientTrack(tks[j], &(*bFieldHandle));
      }
      //std::cout << "rm" << std::endl;
      reco::Vertex v_rm_track(TransientVertex(fitter->vertex(ttks)));
      const double dz = std::abs(v_rm_track.z()-v[0]->z());
      if(v_rm_track.chi2()<0 || 
          dz > 0.005){
        *v[0] = v_rm_track;
        --v[0];
        break;
      }
    }
  }

  VertexDistanceXY dBVCalc = VertexDistanceXY();
  for (auto vIter=v_refit.begin(); vIter!=v_refit.end();++vIter){
    if(vIter->tracksSize()<5) continue;
    Measurement1D dBV = dBVCalc.distance(*vIter, theBeam);
    if(dBV.value()<=0.01) continue;
    if(dBV.error()>=0.0025) continue;
    if(dBV.value()>2.09) continue;
    h_dBV_->Fill(dBV.value());
    h_sigma_dBV_->Fill(dBV.error());
    h_NTrackPerVtx_->Fill(vIter->tracksSize());
    h_vtx_XY_->Fill(vIter->x(), vIter->y());
  }


  //for(auto vIter=v_refit.begin(); vIter!=v_refit.end(); ++vIter){
  //  auto v_after = *vIter;
  //  if(! (v_after.isValid() && (v_after.totalChiSquared()/v_after.degreesOfFreedom())<5)){
  //    v_refit.erase(vIter);
  //    --vIter;
  //  }
  //  std::cout << "# tks: " << v_after.tracks().size() << std::endl;
  //}
  //std::cout << "# vtx: " << v_refit.size() << std::endl;

    //Handle<TrackCollection> tracks;
    //iEvent.getByToken(tracksToken_, tracks);
    //for(TrackCollection::const_iterator itTrack = tracks->begin();
    //    itTrack != tracks->end();
    //    ++itTrack) {
    //  // do something with track parameters, e.g, plot the charge.
    //  // int charge = itTrack->charge();
    //}

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
DVAnalyzer::beginJob()
{
  edm::Service<TFileService> fileService;
  if(!fileService) throw edm::Exception(edm::errors::Configuration, "TFileService is not registered in cfg file");
  h_PFJet_PT_ = fileService->make<TH1F>("h_PFJet_PT","h_PFJet_PT",100,0,100);
  h_PFJet_eta_ = fileService->make<TH1F>("h_PFJet_eta","h_PFJet_eta",30,-3,3);
  h_PFJet_HT_ = fileService->make<TH1F>("h_PFJet_HT","h_PFJet_HT",35,0,3500);
  h_PFJet_N_ = fileService->make<TH1F>("h_PFJet_N","h_PFJet_N",30,0,30);
  h_muon_PT_ = fileService->make<TH1F>("h_muon_PT","h_muon_PT",100,0,100);
  h_muon_eta_ = fileService->make<TH1F>("h_muon_eta","h_muon_eta",100,-5,5);
  h_track_PT_ = fileService->make<TH1F>("h_track_PT","h_track_PT",100,0,25);
  h_Event_cutflow_ = fileService->make<TH1D>("h_Event_cutflow","h_Event_cutflow",10,0,10);
  h_PFJet_cutflow_ = fileService->make<TH1D>("h_PFJet_cutflow","h_PFJet_cutflow",10,0,10);
  h_Track_cutflow_ = fileService->make<TH1D>("h_Track_cutflow","h_Track_cutflow",10,0,10);
  h_NTrackPerVtx_ = fileService->make<TH1D>("h_NTrackPerVtx","h_NTrackPerVtx",50,0,50);
  h_dBV_ = fileService->make<TH1D>("h_dBV","h_dBV",100,0,0.4);
  h_sigma_dBV_ = fileService->make<TH1D>("h_sigma_dBV","h_sigma_dBV",500,0,0.05);
  h_vtx_XY_ = fileService->make<TH2D>("h_vtx_XY","h_vtx_XY",200,-4,4,200,-4,4);
}

void DVAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if(hltConfig_.init(iRun, iSetup, processName_, changed)){
    if(changed) {
      const unsigned int n(hltConfig_.size());
      unsigned int triggerIdx(hltConfig_.triggerIndex(triggerName_));
      if (triggerIdx>=n) {
        std::cout << "TriggerName: " << triggerName_ << " not available in config!" << std::endl;
      }
    }
  }
  else {
    std::cout << "Warning, didn't find trigger process HLT,\t" << processName_ << std::endl;
  }
}

std::vector<TransientVertex> DVAnalyzer::refitVtx(std::vector<reco::TransientTrack>& ttks) {
  if (ttks.size() < 2)
    return std::vector<TransientVertex>();
  std::vector<TransientVertex> v(1, fitter->vertex(ttks));
  if (v[0].normalisedChiSquared() > 5)
    return std::vector<TransientVertex>();
  return v;
}

reco::TrackRefVector DVAnalyzer::SharedTracks(track_set t1, track_set t2){
  reco::TrackRefVector sT;
  for(auto tk:t1){
    if(t2.count(tk)>0){
      sT.push_back(tk);
    }
  }
  return sT;
}

std::vector<TransientTrack> DVAnalyzer::combineTracks(track_set t1, track_set t2, std::map<reco::TrackRef, size_t> seed_track_ref_map, std::vector<TransientTrack> tTracks) {
  track_set combinedSet;
  for (auto tk:t1){
    combinedSet.insert(tk);
  }
  for (auto tk:t2){
    combinedSet.insert(tk);
  }
  std::vector<TransientTrack> tks;
  for (auto tk:combinedSet){
    tks.push_back(tTracks[seed_track_ref_map[tk]]);
  }
  return tks;
}


track_set DVAnalyzer::vertexTracks(const reco::Vertex v) const
{
  track_set result;
  for(auto t=v.tracks_begin(); t!=v.tracks_end(); ++t){
    result.insert(t->castTo<reco::TrackRef>());
  }
  return result;
}

track_vec DVAnalyzer::vertex_track_vec(const reco::Vertex& v) const {
    track_set s = vertexTracks(v);
    return track_vec(s.begin(), s.end());
}

bool DVAnalyzer::is_track_subset(const track_set& a, const track_set& b) const {
  bool is_subset = true;
  const track_set& s = a.size() <= b.size() ? a : b;
  const track_set& l = a.size() <= b.size() ? b : a;

  for (auto t : s)
    if (l.count(t)<1){
      is_subset = false;
      break;
    }
  return is_subset;
}


void
DVAnalyzer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{}

// ------------ method called once each job just after ending the event loop  ------------
void
DVAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DVAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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

//define this as a plug-in
DEFINE_FWK_MODULE(DVAnalyzer);
