// C++ includes
#include <iostream>
#include <cmath>
#include <dirent.h>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TChain.h"

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

int main() {
    // ----- Variables to be set --------------
    // Selection
    bool eLepton = true; // eLepton = true: W -> e
    bool eJet = false; // eJet = true: jet -> e
    bool tightSel = true; // tight or loose lepton selection
    float leptonPt = 20;
//     float jetPt = 0;
    
    // Binning
    const int nEtaBins = 8, nPtBins = 7;
    const float etaBins[nEtaBins+1] = {-2.5,-2,-1.5,-1,0,1,1.5,2,2.5};
    const float ptBins[nPtBins+1] = {0,10,20,30,40,60,100,140};
    
    // ID cuts
    // tightElCuts[barrel/endcap][IDcuts], CSA14 selection, conditions: 25ns, better detector alignment
    const float tightElCuts[2][10] = {{0.0091,0.031,0.0106,0.0532,0.0126,0.0116,0.0609,0.1649,1e-6,1},{0.0106,0.0359,0.0305,0.0835,0.0163,0.5999,0.1126,0.2075,1e-6,1}};
    const float looseElCuts[2][10] = {{0.0181,0.0936,0.0123,0.141,0.0166,0.54342,0.1353,0.24,1e-6,1},{0.0124,0.0642,0.035,0.1115,0.098,0.9187,0.1443,0.3529,1e-6,1}};
    
    // Directories
    TString outDirPNG = "/afs/cern.ch/user/j/jlauwers/www/protected/VBS/TP/FakeRate/";
    TString outDirROOT = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/Results/";
//     TString inDir = "eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/CSA14/WJetsToLNu_13TeV-madgraph-pythia8-tauola_v3/";
    TString inDir = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/BaconTrees/";
    
    // Constants
    const Float_t pi = 3.1416;
    
    // Verbose output
    int verbose = 1; // 0: no messages - 3: all debug information 
    // ----------------------------------------
    
    // Add all BaconTrees to a chain
    TChain* tree = new TChain("ntupler/Events");
    DIR *dpdf;
    struct dirent *epdf;

    dpdf = opendir(inDir);
    if (dpdf != NULL){
        while ((epdf = readdir(dpdf))){
            string fname = epdf->d_name;
            if (fname != "." && fname != "..") {
                tree->Add(inDir+fname);
                if( verbose > 2 ) cout << "Adding file: " << epdf->d_name << endl;
            }
        }
    }

//     tree->Add("../../BaconTrees/baconTree_1009_1_RDq.root");
//     TFile* inFile = TFile::Open("/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/BaconTrees/baconTree_1009_1_RDq.root");
//     TTree* tree = (TTree*) inFile->Get("ntupler/Events");
    
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
    
    TH2F *hNum = new TH2F("Numerator_"+strSel,"Numerator_"+strSel,nEtaBins, etaBins, nPtBins, ptBins);
    TH2F *hDenom = new TH2F("Denominator_"+strSel,"Denominator_"+strSel,nEtaBins, etaBins, nPtBins, ptBins);
    TH1::SetDefaultSumw2();
    
    int lepElec=0, lepMuon=0, jetElec=0, jetMuon=0, l1NotMatched=0, l2NotMatched=0;
    float elecEffDenom=0, elecEffNum=0, muonEffDenom=0, muonEffNum=0;
    
    // Loop over events
    int nevents = tree->GetEntries(); // GetEntriesFast fails for chain
    if( verbose > 0 ) cout << "nevents: " << nevents << endl;
    for( int i = 0; i < nevents; ++i ) {
        tree->GetEntry(i);
        int nElectrons = fElectron->GetEntriesFast();
        int nMuons = fMuon->GetEntriesFast();
//         int nLeptons = nElectrons + nMuons;
        if( verbose > 2 ) cout << "nElectrons: " << nElectrons << ", nMuons: " << nMuons << endl;
        
        // Calculate lepton efficiencies
        if( verbose > 0 ) {
            for( int iE = 0; iE < nElectrons; ++iE ) {
                TElectron *elec = (TElectron*)((*fElectron)[iE]);
                if( fabs(elec->eta) < 2.5 && elec->pt > leptonPt) {
                    elecEffDenom++;
                    bool inEndcap = elec->scEta > 1.479;
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
        
        // -- Select at least one tight lepton with pt > leptonPt GeV that matches a genLepton --
        float leptonEta=-1, leptonPhi=-1; // Has to be removed from jet collection
        // W -> e
        if( eLepton ) {
            if( nElectrons < 1 ) continue;
            bool passFullSel = false;
            for( int iE = 0; iE < nElectrons; ++iE ) {
                TElectron *elec = (TElectron*)((*fElectron)[iE]);
                bool inEndcap = elec->scEta > 1.479;
                bool passSel;
                if( tightSel) passSel = passElectronID(elec, tightElCuts[inEndcap] ) && elec->pt > leptonPt;
                else passSel = passElectronID(elec, looseElCuts[inEndcap] ) && elec->pt > leptonPt;
                if( passSel ) {
                    if( verbose > 1 ) cout << "Electron " << iE <<" passed lepton id" << endl;
                    l1NotMatched++;
                    
                    // Match with genParticle
                    int nGenParticles = fGenParticle->GetEntriesFast();
                    for( int iG = 0; iG < nGenParticles; ++iG ) {
                        TGenParticle *genP = (TGenParticle*)((*fGenParticle)[iG]); 
                        if( abs(genP->pdgId) == 11) {
                            Double_t dR = TMath::Sqrt( TMath::Power(genP->eta - elec->eta, 2) + TMath::Power(abs(abs(genP->phi - elec->phi)-pi)-pi, 2) );
                            if( dR < 0.3 ) {
                                passFullSel = true;
                                leptonEta = elec->eta;
                                leptonPhi = elec->phi;
                                if( verbose > 1 ) cout << "Matched with genElecton" << endl;
                                lepElec++;
                                l1NotMatched--;
                                break;
                            }
                        }
                    }
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
                if( tightSel) passSel = passTightMuonID(muon) && muon->pt > leptonPt;
                else passSel =  passLooseMuonID(muon) && muon->pt > leptonPt;
                if( passSel ) {
                    if( verbose > 1 ) cout << "Muon " << iM << " passed lepton id" << endl;
                    l1NotMatched++;
                    
                    // Match with genParticle
                    int nGenParticles = fGenParticle->GetEntriesFast();
                    for( int iG = 0; iG < nGenParticles; ++iG ) {
                        TGenParticle *genP = (TGenParticle*)((*fGenParticle)[iG]); 
                        if( abs(genP->pdgId) == 13) {
                            Double_t dR = TMath::Sqrt( TMath::Power(genP->eta - muon->eta, 2) + TMath::Power(abs(abs(genP->phi - muon->phi)-pi)-pi, 2) );
                            if( dR < 0.3 ) {
                                passFullSel = true;
                                leptonEta = muon->eta;
                                leptonPhi = muon->phi;
                                if( verbose > 1 ) cout << "Matched with genMuon" << endl;
                                lepMuon++;
                                l1NotMatched--;
                                break;
                            }
                        }
                    }
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
                bool inEndcap = elec->scEta > 1.479;
                bool passSel;
                if( tightSel) passSel = passElectronID(elec, tightElCuts[inEndcap] ) && elec->pt > leptonPt;
                else passSel = passElectronID(elec, looseElCuts[inEndcap] ) && elec->pt > leptonPt;
                if( passSel ) {
                    if( eLepton && secLeptIndex == -1) secLeptIndex = -2; // skip first electron
                    else {
                        secLeptIndex = iE;
                        if( verbose > 1 ) cout << "Electron " << iE <<" passed lepton id" << endl;
                        jetElec++;
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
                if( tightSel) passSel = passTightMuonID(muon) && muon->pt > leptonPt;
                else passSel =  passLooseMuonID(muon) && muon->pt > leptonPt;
                if( passSel ) {
                    if( (!eLepton) && secLeptIndex == -1) secLeptIndex = -2; // skip first muon
                    else {
                        secLeptIndex = iM;
                        if( verbose > 1 ) cout << "Muon " << iM <<" passed lepton id" << endl;
                        jetMuon++;    
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
            hDenom->Fill(jet->eta, jet->pt);
            if( secLeptIndex >= 0 ) {
                Double_t dR;
                if(eJet) dR = TMath::Sqrt( TMath::Power(elec->eta - jet->eta, 2) + TMath::Power(abs(abs(elec->phi - jet->phi)-pi)-pi, 2) );
                else dR = TMath::Sqrt( TMath::Power(muon->eta - jet->eta, 2) + TMath::Power(abs(abs(muon->phi - jet->phi)-pi)-pi, 2) );
                if( dR < 0.3 ) {
                    hNum->Fill(jet->eta, jet->pt);
                    secLeptIndex = -1; // match with only 1 jet
                    if( verbose > 1 ) cout << "Matched with jet" << endl;
                }                
            }
        }
        // Lepton should be matched
        if( secLeptIndex >= 0 ) {
            l2NotMatched++;
//             if(eJet) hNum->Fill(elec->eta, elec->pt);
//             else hNum->Fill(muon->eta, muon->pt);
        }
    }    
    
    // Verbose output
    if( verbose > 0 ) {
        cout << "Electron efficiency: " << elecEffNum/elecEffDenom << ", Muon efficiency: " << muonEffNum/muonEffDenom << endl;
        cout << "# elec1 passing sel: " << lepElec << ", # muons1 passing sel: " << lepMuon << endl;
        cout << "# elec2 passing sel: " << jetElec << ", # muons2 passing sel: " << jetMuon << endl;
        cout << "lep1 not matched with gen: " << l1NotMatched << ", lep2 not matched with jet: " << l2NotMatched << endl;
        cout << "Total fake rate: " << hNum->Integral()/hDenom->Integral() << endl;
    }
    
    // -- Calculate and draw fake rate --
    TFile* outFile = new TFile(outDirROOT+"FakeRate.root","UPDATE");
    hNum->Write();
    hDenom->Write();
    
    TString strTitle = "Fake_rate_";
    strTitle += strSel;
    
    TH2F *hFakeRate = (TH2F*) hNum->Clone(strTitle);
    hFakeRate->Divide(hDenom);
    hFakeRate->SetTitle(strTitle);
    hFakeRate->Write();
    
    TCanvas *c1 = new TCanvas("c1","c1");
    TH1D *hFakeRate_eta =  hFakeRate->ProjectionX("_eta", 0, -1, "e");
    hFakeRate_eta->GetXaxis()->SetTitle("eta");
    hFakeRate_eta->Draw();
    c1->Print(outDirPNG+strTitle+"_eta.png","png");
    
    TCanvas *c2 = new TCanvas("c2","c2");
    TH1D *hFakeRate_pt =  hFakeRate->ProjectionY("_pt", 0, -1, "e");
    hFakeRate_pt->GetXaxis()->SetTitle("pt");
    hFakeRate_pt->Draw();
    c2->Print(outDirPNG+strTitle+"_pt.png","png");
    
    TCanvas *c3 = new TCanvas("c3","c3");
    hFakeRate->GetXaxis()->SetTitle("eta");
    hFakeRate->GetYaxis()->SetTitle("pt");
    hFakeRate->Draw("TEXT E");
    c3->Print(outDirPNG+strTitle+".png","png");
    outFile->Delete();
}

bool passElectronID(TElectron* elec, const float (&cuts)[10] ) {
    return( fabs(elec->dEtaIn) < cuts[0] &&
            fabs(elec->dPhiIn) < cuts[1] &&
            elec->sieie  < cuts[2] &&
            elec->hovere < cuts[3] &&
            fabs(elec->d0) < cuts[4] &&
            fabs(elec->dz) < cuts[5] &&
//             fabs(elec->eoverp)  < cuts[6] && // eoverp = E/p  -> missing
            ((elec->chHadIso03 + max(elec->gammaIso03+elec->neuHadIso03-0.5* elec->puIso03,0.0))/elec->pt) < cuts[7] &&
            (!elec->isConv) &&
            elec->nMissingHits <= cuts[9] );    
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
            muon->nTkLayers > 5 );  
}

bool passLooseMuonID(TMuon* muon) {
    return( ((muon->typeBits)/32)%2 &&
            ( ((muon->typeBits)/2)%2 || ((muon->typeBits)/4)%2 ) );
}
        
        