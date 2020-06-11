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

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

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

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::PFJetCollection> pfJetToken_;
      edm::Handle<reco::PFJetCollection> pfJetHandle_;
      edm::EDGetTokenT<std::vector<reco::Muon>> muonToken_;
      edm::Handle<std::vector<reco::Muon>> muonHandle_;
      edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
      edm::Handle<std::vector<reco::Track>> tracksHandle_;
      edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
      edm::Handle<std::vector<reco::Vertex>> vtxHandle_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
      edm::Handle<edm::TriggerResults> triggerResultsHandle_;

      HLTConfigProvider hltConfig_;

      std::string processName_;
      std::string triggerName_;

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
  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

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

    Handle<TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);
    for(TrackCollection::const_iterator itTrack = tracks->begin();
        itTrack != tracks->end();
        ++itTrack) {
      // do something with track parameters, e.g, plot the charge.
      // int charge = itTrack->charge();
    }

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
}

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
