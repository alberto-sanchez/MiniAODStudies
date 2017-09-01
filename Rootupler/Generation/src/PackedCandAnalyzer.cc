// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

class PackedCandAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PackedCandAnalyzer(const edm::ParameterSet&);
      ~PackedCandAnalyzer() {}

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
};

PackedCandAnalyzer::PackedCandAnalyzer(const edm::ParameterSet& iConfig):
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
}

void PackedCandAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);
    edm::Handle<pat::PackedCandidateCollection> pfs;
    iEvent.getByToken(pfToken_, pfs);
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);

    std::vector<const reco::Candidate *> leptons;
    for (const pat::Muon &mu : *muons) leptons.push_back(&mu);
//    for (const pat::Electron &el : *electrons) leptons.push_back(&el);

    for (const reco::Candidate *lep : leptons) {
        if (lep->pt() < 5) continue;
        // initialize sums
        double charged = 0, neutral = 0, pileup  = 0;
        // now get a list of the PF candidates used to build this lepton, so to exclude them
        std::vector<reco::CandidatePtr> footprint;
        for (unsigned int i = 0, n = lep->numberOfSourceCandidatePtrs(); i < n; ++i) {
            footprint.push_back(lep->sourceCandidatePtr(i));
        }
        // now loop on pf candidates
        for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
            const pat::PackedCandidate &pf = (*pfs)[i];
            if (deltaR(pf,*lep) < 0.2) {
                // pfcandidate-based footprint removal
                if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs,i)) != footprint.end()) {
                    continue;
                }
                if (pf.charge() == 0) {
                    if (pf.pt() > 0.5) neutral += pf.pt();
                } else if (pf.fromPV() >= 2) {
                    charged += pf.pt();
                } else {
                    if (pf.pt() > 0.5) pileup += pf.pt();
                }
            }
        }
        // do deltaBeta
        double iso = charged + std::max(0.0, neutral-0.5*pileup);
	int matched = 0;
	const pat::TriggerObjectStandAloneCollection muHLTMatches1 = (dynamic_cast<const pat::Muon*>(lep))->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon10JpsiBarrel");
	const pat::TriggerObjectStandAloneCollection muHLTMatches2 = (dynamic_cast<const pat::Muon*>(lep))->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon16Jpsi");
	const pat::TriggerObjectStandAloneCollection muHLTMatches3 = (dynamic_cast<const pat::Muon*>(lep))->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon20JpsiBarrelnoCow");
	const pat::TriggerObjectStandAloneCollection muHLTMatches4 = (dynamic_cast<const pat::Muon*>(lep))->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon25Jpsis");
	const pat::TriggerObjectStandAloneCollection muHLTMatches5 = (dynamic_cast<const pat::Muon*>(lep))->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon10UpsilonBarrelnoCow");
	const pat::TriggerObjectStandAloneCollection muHLTMatches6 = (dynamic_cast<const pat::Muon*>(lep))->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon12Upsilons");
	const pat::TriggerObjectStandAloneCollection muHLTMatches7 = (dynamic_cast<const pat::Muon*>(lep))->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon14PhiBarrelnoCow");

	if (muHLTMatches1.size() > 0) matched += 1;
	if (muHLTMatches2.size() > 0) matched += 2;
	if (muHLTMatches3.size() > 0) matched += 4;
	if (muHLTMatches4.size() > 0) matched += 8;
	if (muHLTMatches5.size() > 0) matched += 16;
	if (muHLTMatches6.size() > 0) matched += 32;
	if (muHLTMatches7.size() > 0) matched += 64;

        printf("*** %-8s of pt %6.1f, eta %+4.2f, phi %+4.2f : relIso = %5.2f => matched = %d\n",
                    abs(lep->pdgId())==13 ? "muon" : "electron",
                    lep->pt(), lep->eta(), lep->phi(), iso/lep->pt(),matched);

    }

    // Let's compute the fraction of charged pt from particles with dz < 0.1 cm
    for (const pat::Jet &j :  *jets) {
        if (j.pt() < 40 || fabs(j.eta()) > 2.4) continue;
        double in = 0, out = 0; 
        for (unsigned int id = 0, nd = j.numberOfDaughters(); id < nd; ++id) {
            const pat::PackedCandidate &dau = dynamic_cast<const pat::PackedCandidate &>(*j.daughter(id));
            if (dau.charge() == 0) continue;
            (fabs(dau.dz())<0.1 ? in : out) += dau.pt();
        }
        double sum = in + out;
        printf("Jet with pt %6.1f, eta %+4.2f, beta(0.1) = %+5.3f, pileup mva disc %+.2f\n",
                j.pt(),j.eta(), sum ? in/sum : 0, j.userFloat("pileupJetId:fullDiscriminant"));
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(PackedCandAnalyzer);
