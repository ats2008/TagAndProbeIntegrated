#ifndef MUONNUMBERFILTER_H
#define MUONNUMBERFILTER_H

#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/MuonSelectors.h>
#include <DataFormats/TauReco/interface/PFTau.h>
#include <DataFormats/TauReco/interface/PFTauDiscriminator.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/Math/interface/deltaR.h>
#include <DataFormats/Math/interface/LorentzVector.h>
#include <iostream>
#include <cmath>
#include <TNtuple.h>
#include <TString.h>

using namespace edm;
using namespace std;
// using namespace reco;

 
class muonNumberFilter : public edm::one::EDFilter<> {

    public:
        muonNumberFilter(const edm::ParameterSet &);
        ~muonNumberFilter();

    private:
        float ComputeMT(math::XYZTLorentzVector visP4, const reco::PFMET& met);
        bool filter(edm::Event &, edm::EventSetup const&);
        EDGetTokenT<reco::MuonCollection>  _muonTag;
        EDGetTokenT<reco::PFTauCollection>  _tausTag;
        EDGetTokenT<reco::PFMETCollection>  _metTag;
        EDGetTokenT<reco::PFTauDiscriminator>  _pfDescriminator;

        bool _useZMassCuts;
        bool _useWMassCuts;

};

muonNumberFilter::muonNumberFilter(const edm::ParameterSet & iConfig) :
_muonTag   (consumes<reco::MuonCollection> (iConfig.getParameter<InputTag>("src"))),
_tausTag  (consumes<reco::PFTauCollection>  (iConfig.getParameter<InputTag>("taus"))),
_metTag   (consumes<reco::PFMETCollection> (iConfig.getParameter<InputTag>("met"))),
//_pfDescriminator ( consumes<reco::PFTauDiscriminator> ( InputTag("hpsPFTauDiscriminationByDecayModeFinding")  ) )
_pfDescriminator ( consumes<reco::PFTauDiscriminator> ( InputTag("byVVVLooseDeepTau2018v2p5VSjet")  ) )
{

    _useZMassCuts = true;
    _useWMassCuts = true;

}

muonNumberFilter::~muonNumberFilter()
{}

bool muonNumberFilter::filter(edm::Event & iEvent, edm::EventSetup const& iSetup)
{
    Handle<reco::MuonCollection> muonHandle;
    iEvent.getByToken (_muonTag, muonHandle);

    // very strict - veto all events with > 1 muon
    // if (muonHandle->size() != 1) return false;

    int nmu = 0;
    int imu_selected=0;
    for (unsigned int imu = 0; imu < muonHandle->size(); imu++)
    {
        const reco::Muon& mu = muonHandle->at(imu);
        float pt = mu.pt();
        float iso = (mu.pfIsolationR04().sumChargedHadronPt + max(mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5 * mu.pfIsolationR04().sumPUPt, 0.0)) / pt;
//	std::cout<<"muon pt = "<<pt<<", eta = "<<mu.eta()<<", iso = "<<iso<<", isLoose = "<< muon::isLooseMuon(mu)<<endl;
        if (muon::isMediumMuon(mu) && ( pt > 24 ) && ( fabs(mu.eta()) < 2.4 ) and ( iso < 0.2 ) )  nmu+=1;
       // if (pt > 10 && fabs(mu.eta()) < 2.4 and iso < 0.3)
                
                // TODO ADD
                //'pt > 24 && abs(eta) < 2.1 ' # kinematics
                //'&& ( (pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5 * pfIsolationR04().sumPUPt, 0.0)) / pt() ) < 0.1 ' # isolation
                //'&& isMediumMuon()' # quality -- medium muon
                //# 'pt>0 && abs(eta) < 2.1'
                // ),
      imu_selected=imu;
    }

    if (nmu > 1)
      {
	cout<<"does not pass muon veto"<<endl;
	return false;
      }
    cout<<"does pass muon veto"<<endl;


    const reco::Muon& mu = muonHandle->at(imu_selected);
    //---------------------   get the met for mt computation etc. -----------------
    Handle<reco::PFMETCollection> metHandle;
    iEvent.getByToken (_metTag, metHandle);
    const reco::PFMET& met = metHandle->front();

    float mt = ComputeMT( mu.p4(), met );
    if (mt >= 30 && _useWMassCuts) return false; // reject W+jets

    // ------------------- get Taus -------------------------------
    Handle<reco::PFTauCollection> tauHandle;
    iEvent.getByToken (_tausTag, tauHandle);
    if (tauHandle->size() < 1) return false;

    Handle<reco::PFTauDiscriminator> tausDecayMode;
    //iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByDecayModeFinding"),tausDecayMode);
    iEvent.getByToken(_pfDescriminator,tausDecayMode);

    vector<pair<float, int>> tausIdxPtVec;
    int nTau=0;
    for (uint itau = 0; itau < tauHandle->size(); ++itau)
    {
        const reco::PFTau tau = (*tauHandle)[itau] ;
        math::XYZTLorentzVector pSum = mu.p4() + tau.p4();
        if (_useZMassCuts && (pSum.mass() <= 40 || pSum.mass() >= 80)) continue; // visible mass in (40, 80)
        if (deltaR(tau, mu) < 0.5) continue;

        // store iso against jet and idx 
        float isoMVA = tausDecayMode->operator[](itau).second;
        std::cout<<" for tau : "<<itau<<"  |  isoMVA : "<<isoMVA<<"\n";
        tausIdxPtVec.push_back(make_pair(isoMVA, itau));
        nTau++;
    }

    std::cout<<" GOT A PASS !! "<<nTau<<"\n";
    if (nTau < 1)  return false;


    return true;

}
float muonNumberFilter::ComputeMT(math::XYZTLorentzVector visP4, const reco::PFMET& met)
{
  math::XYZTLorentzVector METP4 (met.pt()*cos(met.phi()), met.pt()*sin(met.phi()), 0, met.pt());
  float scalSum = met.pt() + visP4.pt();

  math::XYZTLorentzVector vecSum (visP4);
  vecSum += METP4;
  float vecSumPt = vecSum.pt();
  return sqrt (scalSum*scalSum - vecSumPt*vecSumPt);
}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(muonNumberFilter);

#endif
