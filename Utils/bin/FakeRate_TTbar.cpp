// C++ includes
#include <iostream>
#include <cmath>
#include <vector>
#include <dirent.h>
#include <time.h>
#include <iomanip>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TChain.h"
#include "TLorentzVector.h"

// Bacon includes
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"

using namespace std;
using namespace baconhep;

bool passElectronID(TElectron* elec, const float (&cuts)[10], double rho );
bool passTightMuonID(TMuon* muon);
bool passLooseMuonID(TMuon* muon);

// Loadbar: x = iEvent,n = totEvents, w = width
static inline void loadbar(unsigned int x, unsigned int n, unsigned int r = 20, unsigned int w = 20)
{
    if ( x % (n/r+1) != 0 ) return;
 
    float ratio  =  x/(float)n;
    unsigned int   c      =  ratio * w;
 
    cout << setw(3) << (int)(ratio*100) << "% [";
    for (unsigned int x=0; x<c; x++) cout << "=";
    for (unsigned int x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
}

int main(int argc, char* argv[]) {    
    // ----- Variables to be set --------------
    // Selection
    bool eJet = true; // eJet = true: jet -> e
    bool tightSel = true; // tight or loose lepton selection
    float leptonPt = 20;
    
    // Binning
    const int nEtaBinsEl = 6, nPtBinsEl = 6;
    const Double_t etaBinsEl[nEtaBinsEl+1] = {0,0.5,1,1.5,1.75,2,2.5};
    const Double_t ptBinsEl[nPtBinsEl+1] = {20,30,50,80,120,180,300};
    
    // ID cuts
    const float looseElCuts[2][10] = {{0.004,0.06,0.01,0.12,0.02,0.1,0.05,0.15,1e-6,1},{0.007,0.03,0.03,0.1,0.02,0.1,0.05,0.15,1e-6,1}};    
    // veto cuts
//     const float tightElCuts[2][10] = {{0.006574,0.022868,0.010181,0.037553,0.009924,0.015310,0.131191,0.074355,1,1},{0.005681,0.032046,0.028766,0.081902,0.027261,0.147154,0.106055,0.090185,1,1}};
    const float tightElCuts[2][10] = {{0.004,0.03,0.01,0.12,0.02,0.1,0.05,0.1,1e-6,0},{0.005,0.02,0.03,0.1,0.02,0.1,0.05,0.1,1e-6,0}};
    
    // Directories
    int maxInFiles=500;
//     TString outDirPNG = "/afs/cern.ch/user/j/jlauwers/www/protected/VBS/TP/FakeRate_rebin/";
    TString outDirROOT = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/Results/";
    TString eosDir = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/";
    TString inDir = "eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/TP_ANALYSIS/BACON_TREES/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV-TP2023HGCALDR/";
//     TString inDir = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/BaconTrees/";
    
    // Constants
    const Float_t pi = 3.1416;
    
    // Verbose output
    int verbose = 2; // 0: no messages - 1: basic output and load bar - 2: calculate lepton efficiency - 3: all debug information 
    // ----------------------------------------
    
    // Parse command line parameters
    if (argc < 2) {
        cout << cout << "Studying jet -> " << (eJet?"e":"mu") << " process" << endl;
    }
    else if( string(argv[1]) == "e" ){
        eJet = true;
        if( verbose > 0 ) cout << "Studying jet -> e process" << endl;
    }
    else if( string(argv[1]) == "m" ){
        eJet = false;
        if( verbose > 0 ) cout << "Studying jet -> mu process" << endl;
    }
        
    // Set timer
    clock_t t1,t2;
    t1=clock();
    
    // Add all BaconTrees to a chain
    TChain* tree = new TChain("Events");
    DIR *dpdf;
    struct dirent *epdf;
    int nFiles = 0;
    if( verbose > 2 ) maxInFiles=1;

    dpdf = opendir(eosDir+inDir);
    if (dpdf != NULL){
        while ((epdf = readdir(dpdf))){
            string fname = epdf->d_name;
            if (fname != "." && fname != "..") {
                tree->Add(eosDir+inDir+fname);
                nFiles++;
                
                if( verbose > 2 ) cout << "Adding file: " << epdf->d_name << endl;
                if( nFiles == maxInFiles ) break;
            }
        }
    }
    else cout << "Nothing in: " << eosDir << inDir << endl;
    if( verbose > 0 ) cout << "Added " << nFiles << " files to chain." << endl;
    
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("Electron*",1);
    tree->SetBranchStatus("Muon*",1);
    tree->SetBranchStatus("Jet04*",1);
    tree->SetBranchStatus("GenParticle*",1);
    tree->SetBranchStatus("rhoIso",1);
    
    TClonesArray *fElectron = new TClonesArray("baconhep::TElectron");
    TClonesArray *fMuon = new TClonesArray("baconhep::TMuon");
    TClonesArray *fJet = new TClonesArray("baconhep::TJet");
    TClonesArray *fGenParticle = new TClonesArray("baconhep::TGenParticle");
    TEventInfo *eventInfo=0;
        
    tree->SetBranchAddress("Electron", &fElectron);
    tree->SetBranchAddress("Muon", &fMuon);
    tree->SetBranchAddress("Jet04", &fJet);
    tree->SetBranchAddress("GenParticle", &fGenParticle);
    tree->SetBranchAddress("Info", &eventInfo);
    
    TString strSel = "jet_to";
    if( eJet ) strSel += "_e";
    else strSel += "_mu";
    
    TH1::SetDefaultSumw2();
    TH2F *hJetNum, *hDenom, *hJetNumB, *hDenomB;
    TH2D *hBinCentres, *hBinCentresB;

    hJetNum= new TH2F("Jet_Numerator_"+strSel,"Jet_Numerator_"+strSel,nEtaBinsEl, etaBinsEl, nPtBinsEl, ptBinsEl);
    hDenom = new TH2F("Denominator_"+strSel,"Denominator_"+strSel,nEtaBinsEl, etaBinsEl, nPtBinsEl, ptBinsEl);
    hBinCentres = new TH2D("Pt_centre_"+strSel,"Pt_centre_"+strSel,nEtaBinsEl, etaBinsEl, nPtBinsEl, ptBinsEl);

    hJetNumB= new TH2F("Jet_Numerator_b_"+strSel,"Jet_Numerator_"+strSel,nEtaBinsEl, etaBinsEl, nPtBinsEl, ptBinsEl);
    hDenomB = new TH2F("Denominator_b_"+strSel,"Denominator_"+strSel,nEtaBinsEl, etaBinsEl, nPtBinsEl, ptBinsEl);
    hBinCentresB = new TH2D("Pt_centre_b_"+strSel,"Pt_centre_"+strSel,nEtaBinsEl, etaBinsEl, nPtBinsEl, ptBinsEl);

    TH2F *hPtMigrationB=0, *hPtMigration=0;
    TH1D *hMigrationBinCentres=0, *hBMigrationBinCentres;
    hPtMigrationB = new TH2F("Pt_migration_b_"+strSel,"Pt_migration_b_"+strSel,nPtBinsEl, ptBinsEl, 60, -50., 250.);
    hPtMigration = new TH2F("Pt_migration_"+strSel,"Pt_migration_"+strSel,nPtBinsEl, ptBinsEl, 60, -50., 250.);
    hMigrationBinCentres = new TH1D("Pt_migration_centre_"+strSel,"Pt_migration_centre_"+strSel, nPtBinsEl, ptBinsEl);
    hBMigrationBinCentres = new TH1D("Pt_migration_b_centre_"+strSel,"Pt_migration_b_centre_"+strSel, nPtBinsEl, ptBinsEl);
    
    int lepElec=0, lepMuon=0, jetElec=0, jetMuon=0, l1NotMatched=0, l2NotMatched=0, elNotFound=0;
    float elecEffDenom=0, elecEffNum=0, muonEffDenom=0, muonEffNum=0;
    int jetlep2GenMatch=0, wLept=0, bJetEvents=0, nonBJetEvents=0;
    int cJets=0;
    
    // Loop over events
    int nevents = tree->GetEntries(); // GetEntriesFast fails for chain
    if( verbose > 2 ) nevents=1000;
    if( verbose > 0 ) cout << "nevents: " << nevents << endl;
    for( int i = 0; i < nevents; ++i ) {
        if( verbose > 0 ) loadbar(i+1,nevents);
        tree->GetEntry(i);
        int nElectrons = fElectron->GetEntriesFast();
        int nMuons = fMuon->GetEntriesFast();
        int nGenParticles = fGenParticle->GetEntriesFast();
        int nJets = fJet->GetEntriesFast();
        if( verbose > 2 ) cout << "nElectrons: " << nElectrons << ", nMuons: " << nMuons << endl;
        
        // Calculate lepton efficiencies
        if( verbose > 0 ) {
            for( int iE = 0; iE < nElectrons; ++iE ) {
                TElectron *elec = (TElectron*)((*fElectron)[iE]);
                if( fabs(elec->eta) < 2.5 && elec->pt > leptonPt) {
                    elecEffDenom++;
                    bool inEndcap = fabs(elec->scEta) > 1.479;
                    bool passSel;
                    if( tightSel) passSel = passElectronID(elec, tightElCuts[inEndcap], eventInfo->rhoIso);
                    else passSel = passElectronID(elec, looseElCuts[inEndcap], eventInfo->rhoIso);
                    if( passSel ) elecEffNum++;
                }
            }
            for( int iM = 0; iM < nMuons; ++iM ) {
                TMuon *muon = (TMuon*)((*fMuon)[iM]);
                if( fabs(muon->eta) < 2.5 && muon->pt > leptonPt) {
                    muonEffDenom++;
                    bool passSel;
                    if( tightSel) passSel = passTightMuonID(muon);
                    else passSel =  passLooseMuonID(muon);
                    if( passSel ) muonEffNum++;
                }
            }
        }
        
        // -- Require genLeptons in tracker to be matched --
        vector<int> genElecIndex, genMuonIndex;
        vector<float> elecEta, elecPhi, muonEta, muonPhi; // Has to be removed from jet collection
        unsigned int nElecReq=0, nMuonReq=0;
        const unsigned int parentId = 24;
        
        // Count number of electrons and muons expected in tracker
        for( int iG = 0; iG < nGenParticles; ++iG ) {
            TGenParticle *genP = (TGenParticle*)((*fGenParticle)[iG]); 
            if( genP->status == 1 && genP->parent>=0 && fabs(genP->eta) < 2.8 ) {  // 2.8 to correct for dR
                if( abs(genP->pdgId) == 11 ) {
                    bool foundW = false, stuck = false/*, inTracker = fabs(genP->eta) < 2.5*/;
                    do {
                        if(genP->parent>=0) {
                            TGenParticle *genTemp = (TGenParticle*)((*fGenParticle)[genP->parent]);
                            genP = genTemp;
                        }
                        else stuck = true;
                        if( abs(genP->pdgId) == parentId) {
                            foundW = true;
                        }
                    } while( !stuck && !foundW ); 
                    if( foundW ) {
                        genElecIndex.push_back(iG);
                        /*if( inTracker )*/ nElecReq++;
                    }
                }
                if( abs(genP->pdgId) == 13 ) {
                    bool foundW = false, stuck = false/*, inTracker = fabs(genP->eta) < 2.5*/;
                    do {
                        if(genP->parent>=0) {
                            TGenParticle *genTemp = (TGenParticle*)((*fGenParticle)[genP->parent]);
                            genP = genTemp;
                        }
                        else stuck = true;
                        if( abs(genP->pdgId) == parentId) {
                            foundW = true;
                        }
                    } while( !stuck && !foundW ); 
                    if( foundW ) {
                        genMuonIndex.push_back(iG);
                        /*if( inTracker )*/ nMuonReq++;
                    }
                }
            }
        }
        if( verbose > 2 ) cout << "Required electrons in tracker: " << nElecReq << endl;
        if( verbose > 2 ) cout << "Required muons in tracker: " << nMuonReq << endl;
        
        vector<int> genElecIndex2 = genElecIndex, genMuonIndex2 = genMuonIndex;
    
        if( nElecReq > 0 ) {
            for( int iE = 0; iE < nElectrons; ++iE ) {
                TElectron *elec = (TElectron*)((*fElectron)[iE]);
//                 bool inEndcap = fabs(elec->scEta) > 1.479;
                bool passSel;
                if( tightSel) passSel = true/*elec->pt > leptonPt && passElectronID(elec, tightElCuts[inEndcap], eventInfo->rhoIso )*/;
                else passSel = true/*elec->pt > leptonPt && passElectronID(elec, looseElCuts[inEndcap], eventInfo->rhoIso )*/;
                if( passSel ) {
                    l1NotMatched++;
                    if( verbose > 2 ) cout << "Electron " << iE <<" passed lepton id" << endl;
                    
                    // Match with genParticle
                    for( vector<int>::iterator iG = genElecIndex.begin(); iG != genElecIndex.end(); ++iG ) {
                        TGenParticle *genP = (TGenParticle*)((*fGenParticle)[*iG]); 
                        Double_t dR = TMath::Sqrt( TMath::Power(genP->eta - elec->eta, 2) + TMath::Power(fabs(fabs(genP->phi - elec->phi)-pi)-pi, 2) );
                        bool sameCharge = genP->pdgId == (elec->q)*-11.;
                        if( verbose > 2 ) cout << "dR: " << (dR < 0.3) << ", sameCharge: " << sameCharge << endl;
                        
                        if( sameCharge && dR < 0.3 ) {
                            lepElec++;
                            if( verbose > 2 ) cout << "Matched with genElectron " << *iG << " from W" << endl;
                            elecEta.push_back(elec->eta);
                            elecPhi.push_back(elec->phi);
    //                         if( fabs(genP->eta) >= 2.5 ) nLeptonsReq++;  // corection for genLepton outside tracker en reco lepton inside
                            genElecIndex.erase(iG); // genLepton can only be matched to one reco lepton
                            wLept++;
                            l1NotMatched--;
                            break;                                            
                        }
                    }
                }
            }
            if( elecEta.size() != nElecReq ) {
                if( verbose > 2 ) cout << "Electron not found" << endl;
                elNotFound++;
                continue;
            }
        }
        
        if( nMuonReq > 0 ) { 
            for( int iM = 0; iM < nMuons; ++iM ) {
                TMuon *muon = (TMuon*)((*fMuon)[iM]);
                bool passSel/*= fabs(muon->eta) < 2.5 && ((muon->typeBits)/32)%2*/; // otherwise too many muons
                if( tightSel) passSel = true /*muon->pt > leptonPt && passTightMuonID(muon)*/;
                else passSel =  muon->pt > leptonPt && passLooseMuonID(muon);
                if( passSel ) {
                    l1NotMatched++;
                    if( verbose > 2 ) cout << "Muon " << iM << " passed lepton id" << endl;
                    
                    // Match with genParticle
                    for( vector<int>::iterator iG = genMuonIndex.begin(); iG != genMuonIndex.end(); ++iG ) {
                        if( verbose > 2 ) cout << "Genleptons left:  " << genMuonIndex.size() << endl;
                        TGenParticle *genP = (TGenParticle*)((*fGenParticle)[*iG]); 
                        Double_t dR = TMath::Sqrt( TMath::Power(genP->eta - muon->eta, 2) + TMath::Power(fabs(fabs(genP->phi - muon->phi)-pi)-pi, 2) );
                        bool sameCharge = genP->pdgId == (muon->q)*-13.;
                        if( verbose > 2 ) cout << "dR: " << (dR < 0.3) << ", sameCharge: " << sameCharge << "gen charge: " << genP->pdgId << ", reco charge: " << (muon->q)*-13 << endl;
                        
                            if( sameCharge && dR < 0.3 ) {
                            lepMuon++;
                            if( verbose > 2 ) cout << "Matched with genMuon " << *iG << " from W" << endl;
                            muonEta.push_back(muon->eta);
                            muonPhi.push_back(muon->phi);
    //                             if( fabs(genP->eta) >= 2.5 ) nLeptonsReq++; // corection for genLepton outside tracker en reco lepton inside
                            genMuonIndex.erase(iG); // genLepton can only be matched to one reco lepton
                            wLept++;
                            l1NotMatched--;
                            break;                                            
                        }
                    }
                }
            }
            if( muonEta.size() != nMuonReq ) {
                if( verbose > 2 ) cout << "Muon not found" << endl;
                continue;
            }
        }
        
        // -- Check if there is a third lepton --
        int secLeptIndex = -1; // more than 1 leptons from jets is negligible
        TMuon *muon=0; 
        TElectron *elec=0;
        
        // Jet -> e
        if( eJet ) {
            for( int iE = 0; iE < nElectrons; ++iE ) {
                TElectron *elec = (TElectron*)((*fElectron)[iE]);
                bool inEndcap = fabs(elec->scEta) > 1.479;
                bool passSel;
                if( tightSel) passSel = elec->pt > leptonPt && fabs(elec->eta) < 2.5 && passElectronID(elec, tightElCuts[inEndcap], eventInfo->rhoIso);
                else passSel = elec->pt > leptonPt && fabs(elec->eta) < 2.5 && passElectronID(elec, looseElCuts[inEndcap], eventInfo->rhoIso);
                if( passSel ) {
                    // skip previously matched electron
                    bool isWLepton = false;
                    for( unsigned int iLep=0; iLep  < elecEta.size(); ++iLep ) {
                        if( elecEta[iLep] == elec->eta || elecPhi[iLep] == elec->phi) isWLepton = true; 
                    }
                    if( isWLepton ) continue;
                    
                    // Match with genParticle to be sure it doesn't come from W    
                    for( vector<int>::iterator iG = genElecIndex2.begin(); iG != genElecIndex2.end(); ++iG ) {
                        TGenParticle *genP = (TGenParticle*)((*fGenParticle)[*iG]); 
                        Double_t dR = TMath::Sqrt( TMath::Power(genP->eta - elec->eta, 2) + TMath::Power(fabs(fabs(genP->phi - elec->phi)-pi)-pi, 2) );
//                         bool sameCharge = genP->pdgId == (elec->q)*-11.;
                        if( dR < 0.3 /*&& sameCharge*/ ) {
                            isWLepton=true;
                            break;
                        }
                    }
                    if( isWLepton ) continue;
                        
                    secLeptIndex = iE;
                    if( verbose > 2 ) cout << "Electron " << iE <<" passed lepton id" << endl;
                    jetElec++;
                        
                    break;
                }
            } 
            if( secLeptIndex >= 0 ) 
                elec = (TElectron*)((*fElectron)[secLeptIndex]);
        }
        
        // Jet -> mu
        else {
            for( int iM = 0; iM < nMuons; ++iM ) {
                TMuon *muon = (TMuon*)((*fMuon)[iM]);
                bool passSel;
                if( tightSel) passSel = muon->pt > leptonPt && fabs(muon->eta) < 2.5 && passTightMuonID(muon);
                else passSel = muon->pt > leptonPt && fabs(muon->eta) < 2.5 && passLooseMuonID(muon);
                if( passSel ) {
                    // skip previously matched muon
                    bool isWLepton = false;
                    for( unsigned int iLep=0; iLep  < muonEta.size(); ++iLep ) {
                        if( muonEta[iLep] == muon->eta && muonPhi[iLep] == muon->phi) isWLepton = true; 
                    }
                    if( isWLepton ) continue;
                    
                    // Match with genParticle to be sure it doesn't come from W                    
                    for( vector<int>::iterator iG = genMuonIndex2.begin(); iG != genMuonIndex2.end(); ++iG ) {
                        TGenParticle *genP = (TGenParticle*)((*fGenParticle)[*iG]); 
                        Double_t dR = TMath::Sqrt( TMath::Power(genP->eta - muon->eta, 2) + TMath::Power(fabs(fabs(genP->phi - muon->phi)-pi)-pi, 2) );
//                         bool sameCharge = genP->pdgId == (muon->q)*-13.;
                        if( dR < 0.3 /*&& sameCharge*/ ) {
                            isWLepton=true;
                            break;
                        }
                    }
                    if( isWLepton ) continue;
                    
                    secLeptIndex = iM;
                    if( verbose > 2 ) cout << "Muon " << iM <<" passed lepton id" << endl;
                    jetMuon++;
                        
                    break;
                }
            } 
            if( secLeptIndex >= 0 ) 
                muon = (TMuon*)((*fMuon)[secLeptIndex]);          
        }
        
        // -- Loop over jets --
        for( int iJ = 0; iJ < nJets; ++iJ ) {
            TJet *jet = (TJet*)((*fJet)[iJ]);
            
            // Skip W leptons 
            bool skipJet = false;
            for( unsigned int iLep=0; iLep  < elecEta.size(); ++iLep ) {
                 Double_t dR = TMath::Sqrt( TMath::Power(elecEta[iLep] - jet->eta, 2) + TMath::Power(fabs(fabs(elecPhi[iLep] - jet->phi)-pi)-pi, 2) );
                 if( dR < 0.3 ) skipJet = true;
            }
            for( unsigned int iLep=0; iLep  < muonEta.size(); ++iLep ) {
                 Double_t dR = TMath::Sqrt( TMath::Power(muonEta[iLep] - jet->eta, 2) + TMath::Power(fabs(fabs(muonPhi[iLep] - jet->phi)-pi)-pi, 2) );
                 if( dR < 0.3 ) skipJet = true;
            } 
            if( skipJet ) {
                if( verbose > 2 ) cout << "Skipping jet, matched to lepton." << endl;
                continue;
            }
            
            // Skip PU jets ( jets not matched to gen jet )
            if( jet->genpt == 0.0 ) {
                if( verbose > 2 ) cout << "Skipping pu jet." << endl;
                continue;
            }
            
            // Skip b-tagged jets
            if( jet->csv > 0.679 ) {
                if( verbose > 2 ) cout << "Skipping b-tagged jet." << endl;
                bJetEvents++;
                if( eJet ) continue; // muon fake rate independent of b-tag
            }
            else nonBJetEvents++;
            
            // Fill histograms
            bool bJet = false;
            if( abs(jet->mcFlavor) == 5 ) bJet = true;
            
            if(bJet) {
                hDenomB->Fill(fabs(jet->eta), jet->pt);
                int binNumber = hBinCentresB->FindBin(fabs(jet->eta),jet->pt);
                hBinCentresB->SetBinContent(binNumber, hBinCentresB->GetBinContent(binNumber) + jet->pt );
            }
            else {
                hDenom->Fill(fabs(jet->eta), jet->pt);
                int binNumber = hBinCentres->FindBin(fabs(jet->eta),jet->pt);
                hBinCentres->SetBinContent(binNumber, hBinCentres->GetBinContent(binNumber) + jet->pt );
            }
            if( secLeptIndex >= 0 ) {
                Double_t dR;
                if(eJet) dR = TMath::Sqrt( TMath::Power(elec->eta - jet->eta, 2) + TMath::Power(fabs(fabs(elec->phi - jet->phi)-pi)-pi, 2) );
                else dR = TMath::Sqrt( TMath::Power(muon->eta - jet->eta, 2) + TMath::Power(fabs(fabs(muon->phi - jet->phi)-pi)-pi, 2) );
                if( dR < 0.3 ) {
                    if(eJet) {
                        if( bJet ) {
                            hPtMigrationB->Fill(jet->pt, jet->pt - elec->pt);
                        }
                        else {
                            hPtMigration->Fill(jet->pt, jet->pt - elec->pt);
                        }
                    }
                    else {
                        if( bJet ) {
                            hPtMigrationB->Fill(jet->pt, jet->pt - muon->pt);
                        }
                        else {
                            hPtMigration->Fill(jet->pt, jet->pt - muon->pt);
                        }
                    }
                    if( bJet) {
                        int binNumber = hBMigrationBinCentres->FindBin(jet->pt);
                        hBMigrationBinCentres->SetBinContent(binNumber, hBMigrationBinCentres->GetBinContent(binNumber) + jet->pt );
                        hJetNumB->Fill(fabs(jet->eta), jet->pt);
                    }
                    else {
                        int binNumber = hMigrationBinCentres->FindBin(jet->pt);
                        hMigrationBinCentres->SetBinContent(binNumber, hMigrationBinCentres->GetBinContent(binNumber) + jet->pt );
                        hJetNum->Fill(fabs(jet->eta), jet->pt);                        
                    }
                    if( abs(jet->mcFlavor) == 4 && fabs(jet->eta) < 2.5 && jet->csv < 0.679) cJets++;
                        
                    secLeptIndex = -1; // match with only 1 jet
                    if( verbose > 2 ) cout << "Matched with jet" << endl;
                }                
            }
        }
        // Lepton should be matched
        if( secLeptIndex >= 0 ) {
            l2NotMatched++;
        }
    }    
    
    // Output
    cout << "Electron efficiency: " << elecEffNum/elecEffDenom << ", Muon efficiency: " << muonEffNum/muonEffDenom << endl;
    cout << "# elec1 passing sel and matched to gen electron: " << lepElec << ", # muons1 passing sel and matched to gen muon: " << lepMuon << endl;
    cout << "# elec1 not found: " << elNotFound << endl;
    cout << "# elec2 passing sel: " << jetElec << ", # muons2 passing sel: " << jetMuon << endl;
    cout << "lep1 not matched with W: " << l1NotMatched << ", lep2 not matched with jet: " << l2NotMatched << endl;
    cout << "lep1 matched to W: " << wLept << endl;
    cout << "lept2 matched with a gen lepton: " << jetlep2GenMatch << endl;
    cout << "b-tagged jets: " << bJetEvents << ", non b-taged jets: " << nonBJetEvents << endl;
    cout << "c-jets in numerator: " << cJets << endl;
    cout << "Fake rate: " << hJetNum->Integral()/hDenom->Integral() << endl;
    cout << "Fake rate  b-jets: " << hJetNumB->Integral()/hDenomB->Integral() << endl;
    cout << "Average fake rate: " << (hJetNumB->Integral()+hJetNum->Integral())/(hDenomB->Integral()+hDenom->Integral()) << endl;
    
    // -- Calculate and draw fake rate --
    TFile* outFile = new TFile(outDirROOT+"FakeRate_TTbar_14TeV_HGCal.root","UPDATE");
    hJetNum->Write();
    hJetNumB->Write();
    hDenom->Write();
    hDenomB->Write();
    hPtMigrationB->Write();
    hPtMigration->Write();
    
    // Bin centres 
    hBinCentres->Divide(hDenom);
    hBinCentresB->Divide(hDenomB);
    TH1D *hPtMigBins = hPtMigration->ProjectionX("_ptbinsce");
    TH1D *hPtBMigBins = hPtMigrationB->ProjectionX("B_ptbinsce");
    hMigrationBinCentres->Divide(hPtMigBins);
    hBMigrationBinCentres->Divide(hPtBMigBins);
    hBinCentres->Write();
    hMigrationBinCentres->Write();
    hBinCentresB->Write();
    hBMigrationBinCentres->Write();
    
    outFile->Delete();
    
    t2=clock();
    float diff ((float)t2-(float)t1);
    cout<< " Total runtime: " << diff/CLOCKS_PER_SEC <<endl;
}

bool passElectronID(TElectron* elec, const float (&cuts)[10], double rho ) { 
    // electron iso rho correction
    float eff_area;
    if( abs(elec->eta)<1.0 ) eff_area = 0.13;
    else if(abs(elec->eta)>1.0 && abs(elec->eta)<1.479 ) eff_area = 0.14;
    else if(abs(elec->eta)>1.479 && abs(elec->eta)<2.0 ) eff_area = 0.07;
    else if(abs(elec->eta)>2.0 && abs(elec->eta)<2.2 ) eff_area = 0.09;
    else if(abs(elec->eta)>2.2 && abs(elec->eta)<2.3 ) eff_area = 0.11;
    else if(abs(elec->eta)>2.3 && abs(elec->eta)<2.4 ) eff_area = 0.11;
    else eff_area = 0.14;
    
    return( fabs(elec->dEtaIn) < cuts[0] &&
            fabs(elec->dPhiIn) < cuts[1] &&
            elec->sieie  < cuts[2] &&
            elec->hovere < cuts[3] &&
            fabs(elec->d0) < cuts[4] &&
            fabs(elec->dz) < cuts[5] &&
            fabs((1 - elec->eoverp)/elec->ecalEnergy) < cuts[6] && 
//             ((elec->chHadIso03 + max(elec->gammaIso03+elec->neuHadIso03-0.5* elec->puIso03,0.0))/elec->pt) < cuts[7] && // electron iso dBeta
            (elec->chHadIso03 + max(elec->gammaIso03 + elec->neuHadIso03 - max(rho, 0.0)*eff_area, 0.0))/elec->pt < cuts[7] && // electron iso rho correction
            (!elec->isConv) &&
            elec->nMissingHits <= cuts[9] 
          );    
}

bool passTightMuonID(TMuon* muon) {
    return( ((muon->typeBits)/2)%2 &&
            ((muon->typeBits)/32)%2 &&
            muon->tkNchi2 < 10 && //
            muon->nValidHits > 0 &&
            muon->nMatchStn > 1 &&
            fabs(muon->d0) < 0.2 && // d0 = -dxy
            fabs(muon->dz) < 0.5 &&
            muon->nPixHits > 0 &&
            muon->nTkLayers > 5 &&
           ((muon->chHadIso04 + max(muon->gammaIso04+muon->neuHadIso04-0.5* muon->puIso04,0.0))/muon->pt) < 0.12 // muon iso  
    );  
}

// bool passLooseMuonID(TMuon* muon) {
//     return( ((muon->typeBits)/32)%2 &&
//             ( ((muon->typeBits)/2)%2 || ((muon->typeBits)/4)%2 ) &&
//             ((muon->chHadIso04 + max(muon->gammaIso04+muon->neuHadIso04-0.5* muon->puIso04,0.0))/muon->pt) < 0.2 // muon iso  
//           );
// }

bool passLooseMuonID(TMuon* muon) {
    return( ((muon->typeBits)/2)%2 &&
            ((muon->typeBits)/32)%2 &&
            muon->tkNchi2 < 10 && //
            muon->nValidHits > 0 &&
            muon->nMatchStn > 1 &&
            fabs(muon->d0) < 0.2 && // d0 = -dxy
            fabs(muon->dz) < 0.5 &&
            muon->nPixHits > 0 &&
            muon->nTkLayers > 5 &&
           ((muon->chHadIso04 + max(muon->gammaIso04+muon->neuHadIso04-0.5* muon->puIso04,0.0))/muon->pt) < 0.16 // muon iso  
    ); 
}
        
        
