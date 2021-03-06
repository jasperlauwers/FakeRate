// C++ includes
#include <iostream>
#include <cmath>
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

using namespace std;
using namespace baconhep;

bool passElectronID(TElectron* elec, const float (&cuts)[10]);
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
    bool eLepton = true; // eLepton = true: W -> e
    bool eJet = true; // eJet = true: jet -> e
    bool tightSel = true; // tight or loose lepton selection
    float leptonPt = 20;
    
    // Binning
//     const int nEtaBins = 10, nPtBins = 16;
//     const float etaBins[nEtaBins+1] = {0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5};
//     const float ptBins[nPtBins+1] = {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250};
    const int nEtaBinsEl = 6, nPtBinsEl = 10;
    const Double_t etaBinsEl[nEtaBinsEl+1] = {0,0.5,1,1.5,1.75,2,2.5};
    const Double_t ptBinsEl[nPtBinsEl+1] = {20,30,40,50,60,80,100,130,170,220,300};
    const int nEtaBinsMu = 4, nPtBinsMu = 3;
    const Double_t etaBinsMu[nEtaBinsMu+1] = {0,0.75,1.5,2,2.5};
    const Double_t ptBinsMu[nPtBinsMu+1] = {20,30,50,120};
    const Double_t ptBinsMujet[nPtBinsMu+1] = {20,40,60,160};
    
    // ID cuts
    // tightElCuts[barrel/endcap][IDcuts], CSA14 selection, conditions: 50ns, poor detector alignment
//     const float tightElCuts[2][10] = {{0.012,0.024,0.01,0.074,0.0091,0.017,0.026,0.10,1e-6,1},{0.019,0.043,0.029,0.08,0.037,0.065,0.076,0.14,1e-6,1}};
    const float looseElCuts[2][10] = {{0.0181,0.0936,0.0123,0.141,0.0166,0.54342,0.1353,0.24,1e-6,1},{0.0124,0.0642,0.035,0.1115,0.098,0.9187,0.1443,0.3529,1e-6,1}};    
    // veto cuts
    const float tightElCuts[2][10] = {{0.021,0.25,0.012,0.24,0.031,0.5,0.32,0.24,1e-6,2},{0.028,0.23,0.035,0.19,0.22,0.91,0.13,0.24,1e-6,3}}; 
    // CSA14 selection, conditions: 25ns, better detector alignment
//     const float tightElCuts[2][10] = {{0.0091,0.031,0.0106,0.0532,0.0126,0.0116,0.0609,0.1649,1e-6,1},{0.0106,0.0359,0.0305,0.0835,0.0163,0.5999,0.1126,0.2075,1e-6,1}};
//     const float looseElCuts[2][10] = {{0.0181,0.0936,0.0123,0.141,0.0166,0.54342,0.1353,0.24,1e-6,1},{0.0124,0.0642,0.035,0.1115,0.098,0.9187,0.1443,0.3529,1e-6,1}};
    
    // Directories
    int maxInFiles=5000;
    TString outDirPNG = "/afs/cern.ch/user/j/jlauwers/www/protected/VBS/TP/FakeRate_rebin/";
    TString outDirROOT = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/Results/";
    TString inDir = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/CSA14/WJetsToLNu_13TeV-madgraph-pythia8-tauola_v3/";
//     TString inDir = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/BaconTrees/";
    
    // Constants
    const Float_t pi = 3.1416;
    
    // Verbose output
    int verbose = 1; // 0: no messages - 1: basic output and load bar - 2: calculate lepton efficiency - 3: all debug information 
    // ----------------------------------------
    
    // Parse command line parameters
    if (argc < 2) {
        cout << cout << "Studying W -> " << (eLepton?"e":"mu") << ", jet -> " << (eJet?"e":"mu") << " process" << endl;
    }
    else if( string(argv[1]) == "ee" ){
        eLepton = true;
        eJet = true;
        if( verbose > 0 ) cout << "Studying W -> e, jet -> e process" << endl;
    }
    else if( string(argv[1]) == "em" ){
        eLepton = true;
        eJet = false;
        if( verbose > 0 ) cout << "Studying W -> e, jet -> mu process" << endl;
    }
    else if( string(argv[1]) == "me" ){
        eLepton = false;
        eJet = true;
        if( verbose > 0 ) cout << "Studying W -> mu, jet -> e process" << endl;
    }
    else if( string(argv[1]) == "mm" ){
        eLepton = false;
        eJet = false;
        if( verbose > 0 ) cout << "Studying W -> mu, jet -> mu process" << endl;
    }
        
    // Set timer
    clock_t t1,t2;
    t1=clock();
    
    // Add all BaconTrees to a chain
    TChain* tree = new TChain("ntupler/Events");
    DIR *dpdf;
    struct dirent *epdf;
    int nFiles = 0;
    if( verbose > 2 ) maxInFiles=1;

    dpdf = opendir(inDir);
    if (dpdf != NULL){
        while ((epdf = readdir(dpdf))){
            string fname = epdf->d_name;
            if (fname != "." && fname != "..") {
                tree->Add(inDir+fname);
                nFiles++;
                
                if( verbose > 2 ) cout << "Adding file: " << epdf->d_name << endl;
                if( nFiles == maxInFiles ) break;
            }
        }
    }
    if( verbose > 0 ) cout << "Added " << nFiles << " files to chain." << endl;
    
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("Electron*",1);
    tree->SetBranchStatus("Muon*",1);
    tree->SetBranchStatus("Jet05*",1);
    tree->SetBranchStatus("GenParticle*",1);
    
    TClonesArray *fElectron = new TClonesArray("baconhep::TElectron");
    TClonesArray *fMuon = new TClonesArray("baconhep::TMuon");
    TClonesArray *fJet = new TClonesArray("baconhep::TJet");
    TClonesArray *fGenParticle = new TClonesArray("baconhep::TGenParticle");
    tree->SetBranchAddress("Electron", &fElectron);
    tree->SetBranchAddress("Muon", &fMuon);
    tree->SetBranchAddress("Jet05", &fJet);
    tree->SetBranchAddress("GenParticle", &fGenParticle);
    
    TString strSel = "W_to";
    if( eLepton ) strSel += "_e";
    else strSel += "_mu";
    if( eJet ) strSel += "_jet_to_e";
    else strSel += "_jet_to_mu";
    
    TH1::SetDefaultSumw2();
    TH2F *hJetNum, *hLeptNum, *hDenom, *hDenomJet=0;
    TH2D *hBinCentres;
    if( eJet ) {
        hJetNum= new TH2F("Jet_Numerator_"+strSel,"Jet_Numerator_"+strSel,nEtaBinsEl, etaBinsEl, nPtBinsEl, ptBinsEl);
        hLeptNum = new TH2F("Lepton_Numerator_"+strSel,"Lepton_Numerator_"+strSel,nEtaBinsEl, etaBinsEl, nPtBinsEl, ptBinsEl);
        hDenom = new TH2F("Denominator_"+strSel,"Denominator_"+strSel,nEtaBinsEl, etaBinsEl, nPtBinsEl, ptBinsEl);
        hBinCentres = new TH2D("Pt_centre_"+strSel,"Pt_centre_"+strSel,nEtaBinsEl, etaBinsEl, nPtBinsEl, ptBinsEl);
    }
    else {
        hJetNum= new TH2F("Jet_Numerator_"+strSel,"Jet_Numerator_"+strSel,nEtaBinsMu, etaBinsMu, nPtBinsMu, ptBinsMujet);
        hLeptNum = new TH2F("Lepton_Numerator_"+strSel,"Lepton_Numerator_"+strSel,nEtaBinsMu, etaBinsMu, nPtBinsMu, ptBinsMu);
        hDenom = new TH2F("Denominator_"+strSel,"Denominator_"+strSel,nEtaBinsMu, etaBinsMu, nPtBinsMu, ptBinsMu);
        hDenomJet = new TH2F("Jet_Denominator_"+strSel,"Jet_Denominator_"+strSel,nEtaBinsMu, etaBinsMu, nPtBinsMu, ptBinsMujet);
        hBinCentres = new TH2D("Pt_centre_"+strSel,"Pt_centre_"+strSel,nEtaBinsMu, etaBinsMu, nPtBinsMu, ptBinsMujet);
    }
    TH2F *hPtMigrationCentral=0, *hPtMigrationForward=0, *hEMigrationCentral=0, *hEMigrationForward=0;
    TH1D *hMigrationBinCentres;
    if( eJet ) {
        hPtMigrationCentral = new TH2F("Pt_migration_barrel_"+strSel,"Pt_migration_barrel_"+strSel,nPtBinsEl, ptBinsEl, 60, -50., 250.);
        hPtMigrationForward = new TH2F("Pt_migration_endcap_"+strSel,"Pt_migration_endcap_"+strSel,nPtBinsEl, ptBinsEl, 60, -50., 250.);
        hEMigrationCentral = new TH2F("E_migration_barrel_"+strSel,"E_migration_barrel_"+strSel,nPtBinsEl, ptBinsEl, 60, -50., 250.);
        hEMigrationForward = new TH2F("E_migration_endcap_"+strSel,"E_migration_endcap_"+strSel,nPtBinsEl, ptBinsEl, 60, -50., 250.);
        hMigrationBinCentres = new TH1D("Pt_migration_centre_"+strSel,"Pt_migration_centre_"+strSel, nPtBinsEl, ptBinsEl);
    }
    else {
        hPtMigrationCentral = new TH2F("Pt_migration_barrel_"+strSel,"Pt_migration_"+strSel,nPtBinsMu, ptBinsMujet, 60, -50., 250.);
        hPtMigrationForward = new TH2F("Pt_migration_endcap_"+strSel,"Pt_migration_endcap_"+strSel,nPtBinsMu, ptBinsMujet, 60, -50., 250.);
        hEMigrationCentral = new TH2F("E_migration_barrel_"+strSel,"E_migration_barrel_"+strSel,nPtBinsMu, ptBinsMujet, 60, -50., 250.);
        hEMigrationForward = new TH2F("E_migration_endcap_"+strSel,"E_migration_endcap_"+strSel,nPtBinsMu, ptBinsMujet, 60, -50., 250.);
        hMigrationBinCentres = new TH1D("Pt_migration_centre_"+strSel,"Pt_migration_centre_"+strSel, nPtBinsMu, ptBinsMujet);
    }
    
    int lepElec=0, lepMuon=0, jetElec=0, jetMuon=0, l1NotMatched=0, l2NotMatched=0;
    float elecEffDenom=0, elecEffNum=0, muonEffDenom=0, muonEffNum=0;
    int jetlep2GenMatch=0, wLept=0;
    
    // Loop over events
    int nevents = tree->GetEntries(); // GetEntriesFast fails for chain
    if( verbose > 0 ) cout << "nevents: " << nevents << endl;
    for( int i = 0; i < nevents; ++i ) {
        if( verbose > 0 ) loadbar(i+1,nevents);
        tree->GetEntry(i);
        int nElectrons = fElectron->GetEntriesFast();
        int nMuons = fMuon->GetEntriesFast();
        if( verbose > 2 ) cout << "nElectrons: " << nElectrons << ", nMuons: " << nMuons << endl;
        
        // Calculate lepton efficiencies
        if( verbose > 1 ) {
            for( int iE = 0; iE < nElectrons; ++iE ) {
                TElectron *elec = (TElectron*)((*fElectron)[iE]);
                if( fabs(elec->eta) < 2.5 && elec->pt > leptonPt) {
                    elecEffDenom++;
                    bool inEndcap = fabs(elec->scEta) > 1.479;
                    bool passSel;
                    if( tightSel) passSel = passElectronID(elec, tightElCuts[inEndcap] );
                    else passSel = passElectronID(elec, looseElCuts[inEndcap] );
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
        
        // -- Select at one tight lepton with pt > leptonPt GeV that matches a genLepton --
        float leptonEta=-1, leptonPhi=-1; // Has to be removed from jet collection
        // W -> e
        if( eLepton ) {
            if( nElectrons < 1 ) continue;
            bool passFullSel = false;
            for( int iE = 0; iE < nElectrons; ++iE ) {
                TElectron *elec = (TElectron*)((*fElectron)[iE]);
                bool inEndcap = fabs(elec->scEta) > 1.479;
                bool passSel;
                if( tightSel) passSel = elec->pt > leptonPt && passElectronID(elec, tightElCuts[inEndcap] );
                else passSel = elec->pt > leptonPt && passElectronID(elec, looseElCuts[inEndcap] );
                if( passSel ) {
                    if( verbose > 2 ) cout << "Electron " << iE <<" passed lepton id" << endl;
                    l1NotMatched++;
                    
                    // Match with genParticle
                    int nGenParticles = fGenParticle->GetEntriesFast();
                    for( int iG = 0; iG < nGenParticles; ++iG ) {
                        TGenParticle *genP = (TGenParticle*)((*fGenParticle)[iG]); 
                        if( abs(genP->pdgId) == 11) {
                            Double_t dR = TMath::Sqrt( TMath::Power(genP->eta - elec->eta, 2) + TMath::Power(abs(abs(genP->phi - elec->phi)-pi)-pi, 2) );
                            if( dR < 0.3 ) {
                                lepElec++;
                                if( verbose > 2 ) cout << "Matched with genElectron" << endl;
                                
                                // Look for a corresponding neutrino to get the W mass
                                int nPdgId = (genP->pdgId==11)?-12:12;
                                for( int iG2 = 0; iG2 < nGenParticles; ++iG2 ) {
                                    TGenParticle *genP2 = (TGenParticle*)((*fGenParticle)[iG2]); 
                                    if( genP2->pdgId == nPdgId) {
                                        Double_t invMass = sqrt(2*genP->pt*genP2->pt*(cosh(genP->eta-genP2->eta)-cos(genP->phi-genP2->phi)));
                                        if( invMass>65 && invMass<95 ) {
                                            passFullSel = true;
                                            leptonEta = elec->eta;
                                            leptonPhi = elec->phi;
                                            if( verbose > 2 ) cout << "Matched with W boson" << endl;
                                            wLept++;
                                            l1NotMatched--;
                                            break;                                            
                                        }
                                    }
                                }
                                if( passFullSel ) break; // first pt genlepton
                            }
                        }
                    }
                    if( passFullSel ) break; // first pt lepton
                }
            }
            if( !passFullSel ) continue;
        }
        
        // W -> mu
        else {
            if( nMuons < 1 ) continue;
            bool passFullSel = false;
            for( int iM = 0; iM < nMuons; ++iM ) {
                TMuon *muon = (TMuon*)((*fMuon)[iM]);
                bool passSel;
                if( tightSel) passSel = muon->pt > leptonPt && passTightMuonID(muon);
                else passSel =  muon->pt > leptonPt && passLooseMuonID(muon);
                if( passSel ) {
                    if( verbose > 2 ) cout << "Muon " << iM << " passed lepton id" << endl;
                    l1NotMatched++;
                    
                    // Match with genParticle
                    int nGenParticles = fGenParticle->GetEntriesFast();
                    for( int iG = 0; iG < nGenParticles; ++iG ) {
                        TGenParticle *genP = (TGenParticle*)((*fGenParticle)[iG]); 
                        if( abs(genP->pdgId) == 13) {
                            Double_t dR = TMath::Sqrt( TMath::Power(genP->eta - muon->eta, 2) + TMath::Power(abs(abs(genP->phi - muon->phi)-pi)-pi, 2) );
                            if( dR < 0.3 ) {
                                lepMuon++;
                                if( verbose > 2 ) cout << "Matched with genMuon" << endl;
                                
                                // Look for a corresponding neutrino to get the W mass
                                int nPdgId = (genP->pdgId==13)?-14:14;
                                for( int iG2 = 0; iG2 < nGenParticles; ++iG2 ) {
                                    TGenParticle *genP2 = (TGenParticle*)((*fGenParticle)[iG2]); 
                                    if( genP2->pdgId == nPdgId) {
                                        Double_t invMass = sqrt(2*genP->pt*genP2->pt*(cosh(genP->eta-genP2->eta)-cos(genP->phi-genP2->phi)));
                                        if( invMass>65 && invMass<95 ) {
                                            passFullSel = true;
                                            leptonEta = muon->eta;
                                            leptonPhi = muon->phi;
                                            if( verbose > 2 ) cout << "Matched with W boson" << endl;
                                            wLept++;
                                            l1NotMatched--;
                                            break;                                            
                                        }
                                    }
                                }
                                if( passFullSel ) break; // first pt genlepton
                            }
                        }
                    }
                    if( passFullSel ) break; // first pt lepton
                }
            }
            if( !passFullSel ) continue;
        }
        
        // -- Check if there is a second lepton --
        int secLeptIndex = -1; // more than 1 leptons from jets is negligible
        TMuon *muon=0; 
        TElectron *elec=0;
        
        // Jet -> e
        if( eJet ) {
            for( int iE = 0; iE < nElectrons; ++iE ) {
                TElectron *elec = (TElectron*)((*fElectron)[iE]);
                bool inEndcap = fabs(elec->scEta) > 1.479;
                bool passSel;
                if( tightSel) passSel = elec->pt > leptonPt && passElectronID(elec, tightElCuts[inEndcap] );
                else passSel = elec->pt > leptonPt && passElectronID(elec, looseElCuts[inEndcap] );
                if( passSel ) {
                    if( !eLepton || (leptonEta != elec->eta && leptonPhi != elec->phi) ) { // skip previously matched electron
                        secLeptIndex = iE;
                        if( verbose > 2 ) cout << "Electron " << iE <<" passed lepton id" << endl;
                        jetElec++;
                        
                        // Match with genparticle
                        int nGenParticles = fGenParticle->GetEntriesFast();
                        for( int iG = 0; iG < nGenParticles; ++iG ) {
                            TGenParticle *genP = (TGenParticle*)((*fGenParticle)[iG]); 
                            if( abs(genP->pdgId) == 11) {
                                Double_t dR = TMath::Sqrt( TMath::Power(genP->eta - elec->eta, 2) + TMath::Power(abs(abs(genP->phi - elec->phi)-pi)-pi, 2) );
                                if( dR < 0.3 ) {
                                    if( verbose > 2 ) cout << "Jet electron matched with genElecton" << endl;
                                    jetlep2GenMatch++;
                                    break;
                                }
                            }
                        }
                        
                        break;
                    }
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
                if( tightSel) passSel = muon->pt > leptonPt && passTightMuonID(muon);
                else passSel = muon->pt > leptonPt && passLooseMuonID(muon);
                if( passSel ) {
                    if( eLepton || (leptonEta != muon->eta && leptonPhi != muon->phi) ) { // skip previously matched muon
                        secLeptIndex = iM;
                        if( verbose > 2 ) cout << "Muon " << iM <<" passed lepton id" << endl;
                        jetMuon++;
                        
                        // Match with genparticle
                        int nGenParticles = fGenParticle->GetEntriesFast();
                        for( int iG = 0; iG < nGenParticles; ++iG ) {
                            TGenParticle *genP = (TGenParticle*)((*fGenParticle)[iG]); 
                            if( abs(genP->pdgId) == 13) {
                                Double_t dR = TMath::Sqrt( TMath::Power(genP->eta - muon->eta, 2) + TMath::Power(abs(abs(genP->phi - muon->phi)-pi)-pi, 2) );
                                if( dR < 0.3 ) {
                                    if( verbose > 2 ) cout << "Jet muon matched with gen muon" << endl;
                                    jetlep2GenMatch++;
                                    break;
                                }
                            }
                        }
                        
                        break;
                    }
                }
            } 
            if( secLeptIndex >= 0 ) 
                muon = (TMuon*)((*fMuon)[secLeptIndex]);          
        }
        
        // -- Loop ever jets --
        int nJets = fJet->GetEntriesFast();
        for( int iJ = 0; iJ < nJets; ++iJ ) {
            TJet *jet = (TJet*)((*fJet)[iJ]);
            
            // Skip first lepton 
            Double_t dR = TMath::Sqrt( TMath::Power(leptonEta - jet->eta, 2) + TMath::Power(abs(abs(leptonPhi - jet->phi)-pi)-pi, 2) );
            if( dR < 0.3 ) continue;
            
            // Fill histograms
            hDenom->Fill(abs(jet->eta), jet->pt);
            if( !eJet ) hDenomJet->Fill(abs(jet->eta), jet->pt);
            int binNumber = hBinCentres->FindBin(abs(jet->eta),jet->pt);
            hBinCentres->SetBinContent(binNumber, hBinCentres->GetBinContent(binNumber) + jet->pt );
            if( secLeptIndex >= 0 ) {
                Double_t dR;
                if(eJet) dR = TMath::Sqrt( TMath::Power(elec->eta - jet->eta, 2) + TMath::Power(abs(abs(elec->phi - jet->phi)-pi)-pi, 2) );
                else dR = TMath::Sqrt( TMath::Power(muon->eta - jet->eta, 2) + TMath::Power(abs(abs(muon->phi - jet->phi)-pi)-pi, 2) );
                if( dR < 0.3 ) {
                    if(eJet) {
                        hLeptNum->Fill(abs(elec->eta), elec->pt);
                        TLorentzVector jet4V;
                        jet4V.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
                        
                        if( abs(jet->eta) < 1.5 ) {
                            hPtMigrationCentral->Fill(jet->pt, jet->pt - elec->pt);
                            hEMigrationCentral->Fill(jet4V.E(), jet->pt - elec->pt);
                        }
                        else {
                            hPtMigrationForward->Fill(jet->pt, jet->pt - elec->pt);
                            hEMigrationForward->Fill(jet4V.E(), jet->pt - elec->pt);
                        }
                    }
                    else {
                        hLeptNum->Fill(abs(muon->eta), muon->pt);
                        TLorentzVector jet4V;
                        jet4V.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
                        
                        if( abs(jet->eta) < 1.5 ) {
                            hPtMigrationCentral->Fill(jet->pt, jet->pt - muon->pt);
                            hEMigrationCentral->Fill(jet4V.E(), jet->pt - muon->pt);
                        }
                        else {
                            hPtMigrationForward->Fill(jet->pt, jet->pt - muon->pt);
                            hEMigrationForward->Fill(jet4V.E(), jet->pt - muon->pt);
                        }
                    }
                    int binNumber = hMigrationBinCentres->FindBin(jet->pt);
                    hMigrationBinCentres->SetBinContent(binNumber, hMigrationBinCentres->GetBinContent(binNumber) + jet->pt );
                        
                    hJetNum->Fill(abs(jet->eta), jet->pt);
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
    cout << "# elec2 passing sel: " << jetElec << ", # muons2 passing sel: " << jetMuon << endl;
    cout << "lep1 not matched with W: " << l1NotMatched << ", lep2 not matched with jet: " << l2NotMatched << endl;
    cout << "lep1 matched to W: " << wLept << endl;
    cout << "lept2 matched with a gen lepton: " << jetlep2GenMatch << endl;
    cout << "Total fake rate: " << hJetNum->Integral()/hDenom->Integral() << endl;
    
    // -- Calculate and draw fake rate --
    TFile* outFile = new TFile(outDirROOT+"FakeRate_centre.root","UPDATE");
    hLeptNum->Write();
    hJetNum->Write();
    hDenom->Write();
    if( !eJet )hDenomJet->Write();
    hPtMigrationCentral->Write();
    hEMigrationCentral->Write();
    hPtMigrationForward->Write();
    hEMigrationForward->Write();
    
    // Bin centres 
    if( eJet ) hBinCentres->Divide(hDenom);
    else hBinCentres->Divide(hDenomJet);
    TH1D *hPtMigBins = hPtMigrationCentral->ProjectionX("_ptbinsce");
    hPtMigBins->Add( hPtMigrationForward->ProjectionX("_ptbinsfw") );
    hMigrationBinCentres->Divide(hPtMigBins);
    hBinCentres->Write();
    hMigrationBinCentres->Write();

    TString strTitle = "Fake_rate_";
    strTitle += strSel;
    
    TH2F *hLeptFakeRate = (TH2F*) hLeptNum->Clone(strTitle+"_leptNum");
    hLeptFakeRate->Divide(hDenom);
    hLeptFakeRate->SetTitle(strTitle+"_leptNum");
    
    TH2F *hJetFakeRate = (TH2F*) hJetNum->Clone(strTitle+"_jetNum");
    hJetFakeRate->Divide(hDenom);
    hJetFakeRate->SetTitle(strTitle+"_jetNum");
    
    TCanvas *c1 = new TCanvas("c1","c1");
    TH1D *hJetFakeRate_eta =  hJetNum->ProjectionX(strTitle+"_jetNum_eta", 0, -1, "e");
    TH1D *hDenom_eta =  hDenom->ProjectionX("Denominator_jetNum_eta", 0, -1, "e");
    hJetFakeRate_eta->Divide(hDenom_eta);
    hJetFakeRate_eta->GetXaxis()->SetTitle("eta");
    hJetFakeRate_eta->Draw();
    c1->Print(outDirPNG+strTitle+"_jetNum_eta.png","png");
    hJetFakeRate_eta->Write();
    
    TCanvas *c2 = new TCanvas("c2","c2");
    TH1D *hJetFakeRate_pt =  hJetNum->ProjectionY(strTitle+"_jetNum_pt", 0, -1, "e");
    TH1D *hDenom_pt =  hDenom->ProjectionY("Denominator_jetNum_pt", 0, -1, "e");
    hJetFakeRate_pt->Divide(hDenom_pt);
    hJetFakeRate_pt->GetXaxis()->SetTitle("pt");
    hJetFakeRate_pt->Draw();
    c2->Print(outDirPNG+strTitle+"_jetNum_pt.png","png");
    hJetFakeRate_pt->Write();
    
    TCanvas *c3 = new TCanvas("c3","c3");
    hLeptFakeRate->GetXaxis()->SetTitle("eta");
    hLeptFakeRate->GetYaxis()->SetTitle("pt");
    hLeptFakeRate->Draw("TEXT");
    c3->Print(outDirPNG+strTitle+"_leptNum"+".png","png");
    hLeptFakeRate->Write();
    
    TCanvas *c4 = new TCanvas("c4","c4");
    hJetFakeRate->GetXaxis()->SetTitle("eta");
    hJetFakeRate->GetYaxis()->SetTitle("pt");
    hJetFakeRate->Draw("TEXT");
    c4->Print(outDirPNG+strTitle+"_jetNum"+".png","png");
    hJetFakeRate->Write();
    outFile->Delete();
    
    t2=clock();
    float diff ((float)t2-(float)t1);
    cout<< " Total runtime: " << diff/CLOCKS_PER_SEC <<endl;
}

bool passElectronID(TElectron* elec, const float (&cuts)[10] ) {  
    return( fabs(elec->dEtaIn) < cuts[0] &&
            fabs(elec->dPhiIn) < cuts[1] &&
            elec->sieie  < cuts[2] &&
            elec->hovere < cuts[3] &&
            fabs(elec->d0) < cuts[4] &&
            fabs(elec->dz) < cuts[5] &&
//             fabs(elec->eoverp)  < cuts[6] && // eoverp = E/p  -> missing
            ((elec->chHadIso03 + max(elec->gammaIso03+elec->neuHadIso03-0.5* elec->puIso03,0.0))/elec->pt) < cuts[7] && // electron iso
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
           ((muon->chHadIso04 + max(muon->gammaIso04+muon->neuHadIso04-0.5* muon->puIso04,0.0))/muon->pt) < 0.2 // muon iso  
    );  
}

bool passLooseMuonID(TMuon* muon) {
    return( ((muon->typeBits)/32)%2 &&
            ( ((muon->typeBits)/2)%2 || ((muon->typeBits)/4)%2 ) &&
            ((muon->chHadIso04 + max(muon->gammaIso04+muon->neuHadIso04-0.5* muon->puIso04,0.0))/muon->pt) < 0.2 // muon iso  
          );
}
        
        
