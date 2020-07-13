#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "LLPAnalyzer/Formats/interface/TrackRefMap.h"
#include "LLPAnalyzer/DVAnalyzer/interface/TrackRescaler.h"


class RescaleTracks : public edm::EDProducer {
  public:
    explicit RescaleTracks(const edm::ParameterSet&);

  private:
    virtual void produce(edm::Event&, const edm::EventSetup&) override;

    const edm::EDGetTokenT<reco::TrackCollection> trackToken_;

    bool doRescale_;
    DVAna::TrackRescaler rescaler_;


};

RescaleTracks::RescaleTracks(const edm::ParameterSet& iConfig):
  trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackTag"))),
  doRescale_(iConfig.getParameter<bool>("doRescale"))
{
  produces<reco::TrackCollection>();
  produces<DVAna::RescaledTrackMap>();
}

void
RescaleTracks::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackToken_, tracks);

  auto output_tracks = std::make_unique<reco::TrackCollection>();
  auto output_tracks_map = std::make_unique<DVAna::RescaledTrackMap>();

  reco::TrackRefProd h_output_tracks = iEvent.getRefBeforePut<reco::TrackCollection>();

  rescaler_.setup(doRescale_,0,0);

  for (size_t i=0; i<tracks->size(); ++i){
    reco::TrackRef tk(tracks,i);
    auto tk_n = rescaler_.scale(*tk).rescaled_tk;
    //std::cout << "before: " << std::abs(tk->dxy()/tk->dxyError()) << " after: " << std::abs(tk_n.dxy()/tk_n.dxyError());
    output_tracks->push_back(rescaler_.scale(*tk).rescaled_tk);
    output_tracks_map->insert({tk, reco::TrackRef(h_output_tracks, output_tracks->size()-1)});
  }

  iEvent.put(std::move(output_tracks));
  iEvent.put(std::move(output_tracks_map));

}


DEFINE_FWK_MODULE(RescaleTracks);
