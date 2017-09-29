// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

//
// class declaration
//

class aodBhadrons : public edm::EDAnalyzer {
public:
  explicit aodBhadrons(const edm::ParameterSet&);
  ~aodBhadrons() override {};
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  void           analyze(const edm::Event&, const edm::EventSetup&) override;
  bool                   isAncestor(const reco::Candidate*, const reco::Candidate*);
  bool                   isAncestor(int, const reco::Candidate*);
  const reco::Candidate *GetAncestor(const reco::Candidate *, int);
  bool                   IsQHadron(int,int);
  bool                   HasQDaughters(const reco::Candidate *, int);
  double                 GetLifetime(TLorentzVector, TVector3, TVector3);
  
  // Variables
  TLorentzVector b_p4;
  TVector3       b_pvtx;
  TVector3       b_dvtx;
  Double_t       b_ct;
  Int_t          b_id;
  Bool_t         b_mixed;
  TLorentzVector q_p4;
  Int_t          q_id;
  UInt_t         q_status;
  UInt_t         nb, ndau;
  UInt_t         nprun, npack, nc, ns, npr, nka, npi, nmu, nel, nga, nother;
  TTree          *gen_tree, *counter_tree;

  edm::EDGetTokenT<reco::GenParticleCollection > genGenToken_;
};

// constructors and destructor
aodBhadrons::aodBhadrons(const edm::ParameterSet& p) {
   edm::Service<TFileService> fs;
   gen_tree = fs->make<TTree>("GenTree","Tree of genParticles");
   gen_tree->Branch("b_p4",     "TLorentzVector", &b_p4);
   gen_tree->Branch("b_pvtx",   "TVector3",       &b_pvtx);
   gen_tree->Branch("b_dvtx",   "TVector3",       &b_dvtx);
   gen_tree->Branch("b_ct",     &b_ct,            "b_ct/D");
   gen_tree->Branch("b_id",     &b_id,            "b_id/I");
   gen_tree->Branch("b_mixed",  &b_mixed,         "b_mixed/O");
   gen_tree->Branch("q_p4",     "TLorentzVector", &q_p4);
   gen_tree->Branch("q_id",     &q_id,            "q_id/I");
   gen_tree->Branch("q_status", &q_status,        "q_status/i");
   gen_tree->Branch("nb",       &nb,              "nb/i");
   gen_tree->Branch("ndau",       &ndau,              "ndau/i");

   counter_tree = fs->make<TTree>("CounterTree","Tree of Counters");
   counter_tree->Branch("nb",       &nb,              "nb/i");
   counter_tree->Branch("nprun",    &nprun,           "nprun/i");
   counter_tree->Branch("npack",    &npack,           "npack/i");
   counter_tree->Branch("nc",       &nc,              "nc/i");
   counter_tree->Branch("ns",       &ns,              "ns/i");
   counter_tree->Branch("npr",      &npr,             "npr/i");
   counter_tree->Branch("nks",      &nka,             "nka/i");
   counter_tree->Branch("npi",      &npi,             "npi/i");
   counter_tree->Branch("nmu",      &nmu,             "nmu/i");
   counter_tree->Branch("nel",      &nel,             "nel/i");
   counter_tree->Branch("nga",      &nga,             "nga/i");
   counter_tree->Branch("nother",   &nother,          "nother/i");

   genGenToken_ = consumes<reco::GenParticleCollection>((edm::InputTag)"genParticles");
}

// member functions

//recursively check is a given particle is ancestor 
bool aodBhadrons::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor,particle->mother(i))) return true;
   }
   return false;
}

//recursively check is a given particle has ancestor with given pdg_id
bool aodBhadrons::isAncestor(int a_pdgId, const reco::Candidate * particle) {
   if (a_pdgId == particle->pdgId() ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(a_pdgId,particle->mother(i))) return true;
   }
   return false;
}

// recursively llok for an ancestor with a given id
const reco::Candidate *aodBhadrons::GetAncestor(const reco::Candidate * particle, int id) {
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      const reco::Candidate *ancestor =  particle->mother(i);
      if (ancestor != nullptr ) {
         if (ancestor->pdgId() == id ) return ancestor;
         return GetAncestor(ancestor,id);
      }
   }
   return nullptr;
}

// check if a given id belongs to a b hadron
bool aodBhadrons::IsQHadron(int id,int q) {
   int quark = abs(q);
   int i = abs(id) % 10000;
   int b = i / 1000;           // baryons
   int m = (i % 1000) / 100;   // mesons
   return (b==quark || m ==quark);
}

bool aodBhadrons::HasQDaughters(const reco::Candidate * particle, int q) {
   for (size_t i=0; i< particle->numberOfDaughters(); i++) {
      if (IsQHadron(particle->daughter(i)->pdgId(),q)) return true;
   }
   return false;
}

double aodBhadrons::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
   TVector3 pv_dv = decay_vtx - production_vtx;
   TVector3 b_p3  = b_p4.Vect();
   pv_dv.SetZ(0.);
   b_p3.SetZ(0.);
   Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
   return lxy*b_p4.M()/b_p3.Mag();
}

void aodBhadrons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   // gen particles are the one containing "important" stuff
   edm::Handle<reco::GenParticleCollection > gen_;
   iEvent.getByToken(genGenToken_,gen_);

   nb = nprun = npack = ndau = 0;
   nc = ns = npr = nka = npi = nmu = nel = nga = nother = 0;
   if (gen_.isValid()) {
      nprun = gen_->size();
      //npack = gen_->size();
      for (size_t k=0; k<gen_->size(); k++) if (std::abs((*gen_)[k].y())<6. && (*gen_)[k].status() ==1) npack++;
      for (size_t i=0; i<nprun; i++) {
         const reco::Candidate * bHadron = &(*gen_)[i];
         const reco::Candidate * bMother = bHadron->mother(0);
         int pdg_mother =  (bMother != nullptr) ? bMother->pdgId() : -1; 
         b_id = bHadron->pdgId();
         ndau = bHadron->numberOfDaughters();
         std::cout << b_id << "," << IsQHadron(b_id,5) << "," << pdg_mother << "," 
                   << IsQHadron(pdg_mother,5) << "," << HasQDaughters(bHadron,5) << "," << bHadron->status() << "," << ndau << std::endl;
         if (IsQHadron(b_id,5) && !HasQDaughters(bHadron,5) && bHadron->status() == 2) {
            if (ndau < 1) std::cout <<  b_id << " has stable daughters," << std::endl;
            b_p4.SetPtEtaPhiM(bHadron->pt(),bHadron->eta(),bHadron->phi(),bHadron->mass());
            nb++;
            q_id = bHadron->pdgId()>0 ? 5 : -5;
            b_mixed = (pdg_mother == -b_id);
            if (b_mixed) {
               b_pvtx.SetXYZ(bMother->vx(),bMother->vy(),bMother->vz());
               q_id *= -1;
            }  else b_pvtx.SetXYZ(bHadron->vx(),bHadron->vy(),bHadron->vz());
            const reco::Candidate * bDaugh = (ndau>0) ? bHadron->daughter(0): bHadron;
            b_dvtx.SetXYZ(bDaugh->vx(),bDaugh->vy(),bDaugh->vz());
            b_ct = GetLifetime(b_p4,b_pvtx,b_dvtx);
            const reco::Candidate *q = GetAncestor(bHadron,q_id);
            q_status = 0;
            q_p4.SetPtEtaPhiM(0.,0.,0.,0.);
            if (q != nullptr) {
               q_p4.SetPtEtaPhiM(q->pt(),q->eta(),q->phi(),q->mass());
               q_status = q->status();
            }  
            for (size_t j=0; j<ndau; j++) {
               const reco::Candidate *dau = bHadron->daughter(j);
               int d_id = dau->pdgId();
               int d_status = dau->status();
               if (IsQHadron(d_id,4) && d_status !=1) {
                  nc++;
                  for (size_t l=0; l<dau->numberOfDaughters(); l++) {
                     int gd_id = dau->daughter(l)->pdgId();
                     int gd_status = dau->daughter(l)->status();
                     if (IsQHadron(gd_id,4) && gd_status !=1) nc++;
                     else if (IsQHadron(gd_id,3) && gd_status !=1) ns++;
                  } 
               }
               else if (IsQHadron(d_id,3) && d_status !=1) {
                  ns++;
                  for (size_t l=0; l<dau->numberOfDaughters(); l++) {
                     int gd_id = dau->daughter(l)->pdgId();
                     int gd_status = dau->daughter(l)->status();
                     if (IsQHadron(gd_id,3) && gd_status !=1) ns++;
                  }
               }
            }
            for (size_t k=0; k<gen_->size(); k++) {
               const reco::Candidate *it = &(*gen_)[k];
               const reco::Candidate * motherIngenCollection = it->mother(0);
               int stable_id = abs(it->pdgId());
               if (std::abs(it->y()) < 6. && it->status() == 1 &&  motherIngenCollection != nullptr && isAncestor(bHadron,motherIngenCollection)) {
                  switch (stable_id) {
                     case 11:   nel++;    break;
                     case 13:   nmu++;    break;
                     case 22:   nga++;    break;
                     case 321:  nka++;    break;
                     case 211:  npi++;    break;
                     case 2212: npr++;    break;
                     default:   nother++; break;
                  }
               }
            }
            std::cout << "***" << bHadron->pdgId() << "," << b_mixed << "," << q_id << "," << nb << " " << nmu << "," << nc << "," << ns << std::endl;
            gen_tree->Fill();
         }
      } 
      counter_tree->Fill();
   } // if (packed.isValid() && gen.isValid())
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void aodBhadrons::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(aodBhadrons);
