// -*- C++ -*-
//
// Package:    testCOde/trackRefProd
// Class:      trackRefProd
// 
/**\class trackRefProd trackRefProd.cc testCOde/trackRefProd/plugins/trackRefProd.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ang Li
//         Created:  Thu, 25 Jun 2020 19:02:24 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"


#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "LLPAnalyzer/Formats/interface/TrackRefMap.h"
#include "LLPAnalyzer/Formats/interface/TrackAssociation.h"
//
// class declaration
//

class trackRefProd : public edm::stream::EDProducer<> {
   public:
      explicit trackRefProd(const edm::ParameterSet&);
      ~trackRefProd();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      typedef DVAna::TrackAssociation Association;
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      bool pass_cand(const pat::PackedCandidate& cand) const {
        if (cand.charge() && cand.hasTrackDetails()) {
          if (skip_weirdos) {
            const reco::Track& tk = cand.pseudoTrack();
            union U { float f; int i; } u;
            u.f = tk.dxyError();
            const bool weirdo =
              u.i == 0x3b8c70c2 ||  // 0.0042859027
              u.i == 0x3c060959 ||  // 0.0081809396
              u.i == 0x3cd24697 ||  // 0.0256684255
              u.i == 0x3dfc7c28 ||  // 0.1232836843
              u.i == 0x3e948f67;    // 0.2901565731
            if (skip_weirdos && weirdo) return false;
          }
          return true;
        }
        return false;
      }

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      const edm::EDGetTokenT<pat::PackedCandidateCollection> packed_candidates_token;
      const edm::EDGetTokenT<reco::VertexCollection> vtx_token;

      bool skip_weirdos;
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
trackRefProd::trackRefProd(const edm::ParameterSet& iConfig):
  packed_candidates_token(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packed_candidates_src"))),
  vtx_token(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvTag"))),
  skip_weirdos(iConfig.getParameter<bool>("skip_weirdos"))
{
  //register your products
  produces<reco::TrackRefVector>().setBranchAlias("trackRefFromCand");
  produces<reco::TrackCollection>().setBranchAlias("trackRefFromCand");
  produces<DVAna::UnpackedCandidateTracksMap>();
  produces<Association>();
}


trackRefProd::~trackRefProd()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
trackRefProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  edm::Handle<pat::PackedCandidateCollection> packed_candidates;
  iEvent.getByToken(packed_candidates_token, packed_candidates);

  edm::Handle<reco::VertexCollection> primary_vertices;
  iEvent.getByToken(vtx_token, primary_vertices);

  std::vector<reco::VertexRef> vertices;
  for(size_t ivtx=0; ivtx<primary_vertices->size(); ++ivtx){
    vertices.push_back(reco::VertexRef(primary_vertices, ivtx));
  }

  auto tracks_map = std::make_unique<reco::TrackRefVector>();
  auto tracks = std::make_unique<reco::TrackCollection>();
  auto trackRef_map = std::make_unique<DVAna::UnpackedCandidateTracksMap>();
  //auto trackAsso = std::make_unique<DVAna::TrackAssociation>();
  std::unique_ptr<Association> trackAsso;
  trackAsso.reset(new Association(&iEvent.productGetter()));


  reco::TrackRefProd h_output_tracks = iEvent.getRefBeforePut<reco::TrackCollection>();

  for (size_t i = 0; i<packed_candidates->size(); ++i){
    const pat::PackedCandidate& cand = (*packed_candidates)[i];
    if(pass_cand(cand)){
      const reco::Track& tk = cand.pseudoTrack();
      tracks->push_back(tk);
      tracks_map->push_back(reco::TrackRef(h_output_tracks, tracks->size()-1));
      trackRef_map->insert(std::pair<reco::CandidatePtr, reco::TrackRef>(reco::CandidatePtr(packed_candidates, i), reco::TrackRef(h_output_tracks, tracks_map->size())));
      //also fill the track vertex association
      for(size_t iv = 0; iv<primary_vertices->size(); ++iv){
        reco::VertexRef vtxref = vertices.at(iv);
        if(cand.fromPV(iv)>1){
          trackAsso->insert(vtxref, reco::TrackRef(h_output_tracks, tracks->size()-1));
          //std::cout << "trackAsso: vtx: " << vtxref.key() << 
          //  " associated with track: " << reco::TrackRef(h_output_tracks, tracks->size()-1).key() << " with: " << cand.fromPV(iv) << std::endl;
        }
      }


    }
  }

  iEvent.put(std::move(tracks));
  iEvent.put(std::move(tracks_map));
  iEvent.put(std::move(trackRef_map));
  iEvent.put(std::move(trackAsso));

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
trackRefProd::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
trackRefProd::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
trackRefProd::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
trackRefProd::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
trackRefProd::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
trackRefProd::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
trackRefProd::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(trackRefProd);
