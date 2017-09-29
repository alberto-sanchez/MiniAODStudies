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

class MiniAODDumpPrunned : public edm::EDAnalyzer {
   public:
      explicit MiniAODDumpPrunned (const edm::ParameterSet&);
      ~MiniAODDumpPrunned() override {};
      bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

   private:
      void analyze(const edm::Event&, const edm::EventSetup&) override;

      edm::EDGetTokenT<reco::GenParticleCollection > prunedGenToken_;
      edm::EDGetTokenT<pat::PackedGenParticleCollection > packedGenToken_;
};

MiniAODDumpPrunned::MiniAODDumpPrunned(const edm::ParameterSet& pset) {
  prunedGenToken_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
  packedGenToken_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
}


//Check recursively if any ancestor of particle is the given one
bool MiniAODDumpPrunned::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
        if(ancestor == particle ) return true;
        for(size_t i=0;i< particle->numberOfMothers();i++)
        {
                if(isAncestor(ancestor,particle->mother(i))) return true;
        }
        return false;
}

void MiniAODDumpPrunned::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

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
      for (size_t i=0; i<pruned->size(); i++) {
          const reco::Candidate * TheParticle = &(*pruned)[i];
          //const reco::Candidate *mother = TheParticle->mother(0);
          if (TheParticle->status() == 2) {
             //if (mother == nullptr || mother->status()!=2) {
                std::cout << i <<":";
                if (abs(TheParticle->pdgId()) == 511 || abs(TheParticle->pdgId()) == 531) std::cout << "**" << TheParticle->pdgId();
                else std::cout << TheParticle->pdgId();
                std::cout << "(" << TheParticle->numberOfMothers() 
                          << "," << TheParticle->numberOfDaughters()
                          << " --> ";
                for (size_t k=0; k< TheParticle->numberOfDaughters(); k++)
                    std::cout << TheParticle->daughter(k)->pdgId() << " ";  
             //} 
             std::cout << ") ";
          } else std::cout << i <<":"<< TheParticle->pdgId() << "(" << TheParticle->status() 
          << "," << TheParticle->numberOfMothers() << "," << TheParticle->numberOfDaughters() << ") ";
      }
      std::cout << "\n";
   }
}

DEFINE_FWK_MODULE(MiniAODDumpPrunned);
