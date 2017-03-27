// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//
// class declaration
//

class MiniAODGenPartDumper2 : public edm::EDAnalyzer {
   public:
      explicit MiniAODGenPartDumper2 (const edm::ParameterSet&);
      ~MiniAODGenPartDumper2() {};
      bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      edm::EDGetTokenT<reco::GenParticleCollection > prunedGenToken_;
      edm::EDGetTokenT<pat::PackedGenParticleCollection > packedGenToken_;
};

MiniAODGenPartDumper2::MiniAODGenPartDumper2(const edm::ParameterSet& pset) {
  prunedGenToken_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
  packedGenToken_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
}


//Check recursively if any ancestor of particle is the given one
bool MiniAODGenPartDumper2::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
        if(ancestor == particle ) return true;
        for(size_t i=0;i< particle->numberOfMothers();i++)
        {
                if(isAncestor(ancestor,particle->mother(i))) return true;
        }
        return false;
}

void MiniAODGenPartDumper2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   // Pruned particles are the one containing "important" stuff
   edm::Handle<reco::GenParticleCollection > pruned;
   iEvent.getByToken(prunedGenToken_,pruned);

   // Packed particles are all the status 1, so usable to remake jets
   // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
   edm::Handle<pat::PackedGenParticleCollection > packed;
   iEvent.getByToken(packedGenToken_,packed);

   // Dump everything
   if (packed.isValid() && pruned.isValid()) {
      std::cout << "size of prunned " << pruned->size() << " size of packed " << packed->size() << std::endl;
      std::cout << "pruned particles: \n";
      for (size_t i=0; i<pruned->size(); i++) std::cout << (*pruned)[i].pdgId() << "(" << (*pruned)[i].status() 
          << "," << (*pruned)[i].numberOfMothers() << "," << (*pruned)[i].numberOfDaughters() << ") ";
      std::cout << "\n";
      std::cout << "packed particles (all status 1): \n";
      for (size_t j=0; j<packed->size(); j++) {
          const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
          if (motherInPrunedCollection == nullptr) std::cout << (*packed)[j].pdgId() << "(null,"     << (*packed)[j].numberOfMothers() << "," << (*packed)[j].numberOfDaughters() << ") ";
          else std::cout << (*packed)[j].pdgId() << "("  << motherInPrunedCollection->pdgId() << "," << (*packed)[j].numberOfMothers() << "," << (*packed)[j].numberOfDaughters() << ") "; 
      }
      std::cout << "\n";
   }
}

DEFINE_FWK_MODULE(MiniAODGenPartDumper2);
