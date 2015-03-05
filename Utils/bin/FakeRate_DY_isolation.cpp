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

// Bacon includes
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"

using namespace std;
using namespace baconhep;

bool passElectronID(TElectron* elec, const float (&cuts)[10]);
bool passTightMuonID(TMuon* muon);
bool passLooseMuonID(TMuon* muon);
void passElectronID(TElectron* elec, const float (&cuts)[10], vector<bool> &passCuts );
void passTightMuonID(TMuon* muon, vector<bool> &passCuts );

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

int main() {
    // ----- Variables to be set --------------
    // Selection
    bool tightSel = false; // tight or loose lepton selection
    float leptonPt = 20;
    bool electron = true;
    
    // ID cuts
    // tightElCuts[barrel/endcap][IDcuts], CSA14 selection, conditions: 50ns, poor detector alignment
//     const float tightElCuts[2][10] = {{0.012,0.024,0.01,0.074,0.0091,0.017,0.026,0.10,1e-6,1},{0.019,0.043,0.029,0.08,0.037,0.065,0.076,0.14,1e-6,1}};
//     const float looseElCuts[2][10] = {{0.016,0.08,0.012,0.15,0.019,0.036,0.11,0.18,1e-6,1},{0.025,0.097,0.032,0.12,0.099,0.88,0.11,0.21,1e-6,1}}; 
    // veto cuts
    const float looseElCuts[2][10] = {{0.021,0.25,0.012,0.24,0.031,0.5,0.32,0.4,1e-6,2},{0.028,0.23,0.035,0.19,0.22,0.91,0.13,0.4,1e-6,3}}; 
    // CSA14 selection, conditions: 25ns, better detector alignment
//     const float tightElCuts[2][10] = {{0.0091,0.031,0.0106,0.0532,0.0126,0.0116,0.0609,0.1649,1e-6,1},{0.0106,0.0359,0.0305,0.0835,0.0163,0.5999,0.1126,0.2075,1e-6,1}};
//     const float looseElCuts[2][10] = {{0.0181,0.0936,0.0123,0.141,0.0166,0.54342,0.1353,0.24,1e-6,1},{0.0124,0.0642,0.035,0.1115,0.098,0.9187,0.1443,0.3529,1e-6,1}};
    // modified tight cuts: loose iso/ sieie, dz
    const float tightElCuts[2][10] = {{0.012,0.024,0.012,0.074,0.0091,0.5 ,0.026,0.18,1e-6,1},{0.019,0.043,0.032,0.08,0.037,0.91,0.076,0.21,1e-6,1}};
    
    // Directories
    int maxInFiles=10;
    TString outDirPNG = "/afs/cern.ch/user/j/jlauwers/www/protected/VBS/TP/FakeRate/";
    TString outDirROOT = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/Results/";
    TString inDir = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/eos/cms/store/group/upgrade/delphes/VBS_FakeRate/";
    
    // Verbose output
    int verbose = 1; // 0: no messages - 3: all debug information 
    // ----------------------------------------
    
    // use correct inDir
    if( electron ) inDir += "DYToEE_Shashlik/";
    else inDir += "DYToMM_Shashlik/";
    
    // Set timer
    clock_t t1,t2;
    t1=clock();
    
    // Add all BaconTrees to a chain
    TChain* tree = new TChain("Events");
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
    tree->SetBranchStatus("rhoIso",1);
    
    TClonesArray *fElectron = new TClonesArray("baconhep::TElectron");
    TClonesArray *fMuon = new TClonesArray("baconhep::TMuon");
    TEventInfo *eventInfo=0;
    tree->SetBranchAddress("Electron", &fElectron);
    tree->SetBranchAddress("Muon", &fMuon);
    tree->SetBranchAddress("Info", &eventInfo);
    
    TH1::SetDefaultSumw2();
    
    // Iso
    TH1F *hIsoElec = new TH1F("isoRhoCorr_elec","isoRhoCorr_elec",200,0,2);
    TH1F *hIsoMuon = new TH1F("isoRhoCorr_muon","isoRhoCorr_muon",200,0,2);
    TH1F *hCumulElec = new TH1F("cumulative_isoRhoCorr_elec","cumulative_isoRhoCorr_elec",200,0,2);
    TH1F *hCumulMuon = new TH1F("cumulative_isoRhoCorr_muon","cumulative_isoRhoCorr_muon",200,0,2);
    
    // ID cut study
    vector<bool> elecPassCuts, muonPassCuts;
    TH1F *hElecPassCuts = new TH1F("elec_pass_cut","elec_pass_cut",10,0,10);
    TH1F *hMuonPassCuts = new TH1F("muon_pass_cut","muon_pass_cut",11,0,11);
    TH1F *hElecPassOtherCuts = new TH1F("elec_pass_other_cuts","elec_pass_other_cuts",10,0,10);
    TH1F *hMuonPassOtherCuts = new TH1F("muon_pass_other_cuts","muon_pass_other_cuts",11,0,11);
    
    
    float elecEffDenom=0, elecEffNum=0, muonEffDenom=0, muonEffNum=0;
    
    // Loop over events
    int nevents = tree->GetEntries(); // GetEntriesFast fails for chain
    if( verbose > 0 ) cout << "nevents: " << nevents << endl;
    for( int i = 0; i < nevents; ++i ) {
        if( verbose > 0 ) loadbar(i+1,nevents);
        tree->GetEntry(i);
        int nElectrons = fElectron->GetEntriesFast();
        int nMuons = fMuon->GetEntriesFast();
        if( verbose > 2 ) cout << "nElectrons: " << nElectrons << ", nMuons: " << nMuons << endl;
        
        // Get rhoIso
//         cout << "rhoIso: " << eventInfo->rhoIso << endl;
        
        // Calculate lepton efficiencies
        if( verbose > 0 ) {
            for( int iE = 0; iE < nElectrons; ++iE ) {
                TElectron *elec = (TElectron*)((*fElectron)[iE]);
                if( fabs(elec->eta) < 2.5 && elec->pt > leptonPt) {
                    hIsoElec->Fill(eventInfo->rhoIso / elec->pt);
                    elecEffDenom++;
                    bool inEndcap = fabs(elec->scEta) > 1.479;
                    bool passSel;
                    if( tightSel) passSel = passElectronID(elec, tightElCuts[inEndcap] );
                    else passSel = passElectronID(elec, looseElCuts[inEndcap] );
                    if( passSel ) elecEffNum++;
                    
                    vector<bool> elecPassOtherCuts(10, true);
                    if( tightSel) passElectronID(elec, tightElCuts[inEndcap], elecPassCuts);
                    else passElectronID(elec, looseElCuts[inEndcap], elecPassCuts);
                    int i=1;
                    for( vector<bool>::const_iterator it = elecPassCuts.begin(); it != elecPassCuts.end(); ++it, ++i ) {
                        if( *it ) hElecPassCuts->SetBinContent(i, hElecPassCuts->GetBinContent(i) + 1 );
                        else {
                            int j=1;
                            for( vector<bool>::iterator it2 = elecPassOtherCuts.begin(); it2 != elecPassOtherCuts.end(); ++it2, ++j ) {
                                if ( i != j) *it2 = false;
                            }
                        }
                    }
                    
                    int j=1;
                    for( vector<bool>::iterator it2 = elecPassOtherCuts.begin(); it2 != elecPassOtherCuts.end(); ++it2, ++j ) {
                        if ( *it2 ) hElecPassOtherCuts->SetBinContent(j, hElecPassOtherCuts->GetBinContent(j) + 1 );
                    }
                    
                    elecPassOtherCuts.clear();
                    elecPassCuts.clear();
                    
                }
            }
            for( int iM = 0; iM < nMuons; ++iM ) {
                TMuon *muon = (TMuon*)((*fMuon)[iM]);
                if( fabs(muon->eta) < 2.5 && muon->pt > leptonPt) {
                    hIsoMuon->Fill(eventInfo->rhoIso / muon->pt);
                    muonEffDenom++;
                    bool passSel;
                    if( tightSel) passSel = passTightMuonID(muon);
                    else passSel =  passLooseMuonID(muon);
                    if( passSel ) muonEffNum++;
                    
                    vector<bool> muonPassOtherCuts(11, true);
                    passTightMuonID(muon, muonPassCuts);
                    int i=1;
                    for( vector<bool>::const_iterator it = muonPassCuts.begin(); it != muonPassCuts.end(); ++it, ++i ) {
                        if( *it ) hMuonPassCuts->SetBinContent(i, hMuonPassCuts->GetBinContent(i) + 1 );
                        else {
                            int j=1;
                            for( vector<bool>::iterator it2 = muonPassOtherCuts.begin(); it2 != muonPassOtherCuts.end(); ++it2, ++j ) {
                                if ( i != j) *it2 = false;
                            }
                        }
                    }
                    
                    int j=1;
                    for( vector<bool>::iterator it2 = muonPassOtherCuts.begin(); it2 != muonPassOtherCuts.end(); ++it2, ++j ) {
                        if ( *it2 ) hMuonPassOtherCuts->SetBinContent(j, hMuonPassOtherCuts->GetBinContent(j) + 1 );
                    }
                    
                    muonPassOtherCuts.clear();
                    muonPassCuts.clear();
                }
            }
        }
    }
    
    // Calculate Isolation
    Double_t integralElec, integralMuon, totElec, totMuon;
    totElec = hIsoElec->Integral(0,200+1);
    totMuon = hIsoMuon->Integral(0,200+1);
    
    for(int ibin=0; ibin<200; ibin++){
        integralElec = hIsoElec->Integral(0,ibin);
        hCumulElec->SetBinContent(ibin, integralElec/totElec);
        integralMuon = hIsoMuon->Integral(0,ibin);
        hCumulMuon->SetBinContent(ibin, integralMuon/totMuon);
    }
           
    // Verbose output
    if( verbose > 0 ) {
        cout << "Electron efficiency: " << elecEffNum/elecEffDenom << ", Muon efficiency: " << muonEffNum/muonEffDenom << endl;
        cout << "80% efficiency electrons: " << hCumulElec->GetBinCenter(hCumulElec->FindFirstBinAbove(0.8)) << endl;
        cout << "80% efficiency muons: " << hCumulMuon->GetBinCenter(hCumulMuon->FindFirstBinAbove(0.8)) << endl;
    }
    
    // -- Calculate and draw fake rate --
    TString outName = "FakeRate_DY";
    outName += (electron?"_ee":"_mm");
    outName += "_iso.root";
    
    TFile* outFile = new TFile(outDirROOT+outName,"UPDATE");
    
    // ID cut study
    hElecPassCuts->GetXaxis()->SetBinLabel(1,"dEtaIn");
    hElecPassCuts->GetXaxis()->SetBinLabel(2,"dPhiIn");
    hElecPassCuts->GetXaxis()->SetBinLabel(3,"sieie");
    hElecPassCuts->GetXaxis()->SetBinLabel(4,"H/E");
    hElecPassCuts->GetXaxis()->SetBinLabel(5,"d0");
    hElecPassCuts->GetXaxis()->SetBinLabel(6,"dZ");
    hElecPassCuts->GetXaxis()->SetBinLabel(7,"iso");
    hElecPassCuts->GetXaxis()->SetBinLabel(8,"conv");
    hElecPassCuts->GetXaxis()->SetBinLabel(9,"missingHits");
    hElecPassCuts->GetXaxis()->SetBinLabel(10,"ALL");
    hElecPassCuts->SetBinContent(10, elecEffNum/elecEffDenom);
    
    hElecPassOtherCuts->GetXaxis()->SetBinLabel(1,"dEtaIn");
    hElecPassOtherCuts->GetXaxis()->SetBinLabel(2,"dPhiIn");
    hElecPassOtherCuts->GetXaxis()->SetBinLabel(3,"sieie");
    hElecPassOtherCuts->GetXaxis()->SetBinLabel(4,"H/E");
    hElecPassOtherCuts->GetXaxis()->SetBinLabel(5,"d0");
    hElecPassOtherCuts->GetXaxis()->SetBinLabel(6,"dZ");
    hElecPassOtherCuts->GetXaxis()->SetBinLabel(7,"iso");
    hElecPassOtherCuts->GetXaxis()->SetBinLabel(8,"conv");
    hElecPassOtherCuts->GetXaxis()->SetBinLabel(9,"missingHits");
    hElecPassOtherCuts->GetXaxis()->SetBinLabel(10,"ALL");
    hElecPassOtherCuts->SetBinContent(10, elecEffNum/elecEffDenom);
    
    hMuonPassCuts->GetXaxis()->SetBinLabel(1,"isGlobal");
    hMuonPassCuts->GetXaxis()->SetBinLabel(2,"isPF");
    hMuonPassCuts->GetXaxis()->SetBinLabel(3,"chi2");
    hMuonPassCuts->GetXaxis()->SetBinLabel(4,"validHits");
    hMuonPassCuts->GetXaxis()->SetBinLabel(5,"matchedSegm");
    hMuonPassCuts->GetXaxis()->SetBinLabel(6,"do");
    hMuonPassCuts->GetXaxis()->SetBinLabel(7,"dz");
    hMuonPassCuts->GetXaxis()->SetBinLabel(8,"nPixel");
    hMuonPassCuts->GetXaxis()->SetBinLabel(9,"trackLayers");
    hMuonPassCuts->GetXaxis()->SetBinLabel(10,"iso");
    hMuonPassCuts->GetXaxis()->SetBinLabel(11,"ALL");
    hMuonPassCuts->SetBinContent(11, muonEffNum/muonEffDenom);
    
    hMuonPassOtherCuts->GetXaxis()->SetBinLabel(1,"isGlobal");
    hMuonPassOtherCuts->GetXaxis()->SetBinLabel(2,"isPF");
    hMuonPassOtherCuts->GetXaxis()->SetBinLabel(3,"chi2");
    hMuonPassOtherCuts->GetXaxis()->SetBinLabel(4,"validHits");
    hMuonPassOtherCuts->GetXaxis()->SetBinLabel(5,"matchedSegm");
    hMuonPassOtherCuts->GetXaxis()->SetBinLabel(6,"do");
    hMuonPassOtherCuts->GetXaxis()->SetBinLabel(7,"dz");
    hMuonPassOtherCuts->GetXaxis()->SetBinLabel(8,"nPixel");
    hMuonPassOtherCuts->GetXaxis()->SetBinLabel(9,"trackLayers");
    hMuonPassOtherCuts->GetXaxis()->SetBinLabel(10,"iso");
    hMuonPassOtherCuts->GetXaxis()->SetBinLabel(11,"ALL");
    hMuonPassOtherCuts->SetBinContent(11, muonEffNum/muonEffDenom);
    
    for( int i=1; i<11; ++i ) {
        if( i < 10) {
            hElecPassCuts->SetBinContent(i, hElecPassCuts->GetBinContent(i)/elecEffDenom);
            hElecPassOtherCuts->SetBinContent(i, hElecPassOtherCuts->GetBinContent(i)/elecEffDenom);
        }
        hMuonPassCuts->SetBinContent(i, hMuonPassCuts->GetBinContent(i)/muonEffDenom); 
        hMuonPassOtherCuts->SetBinContent(i, hMuonPassOtherCuts->GetBinContent(i)/muonEffDenom); 
    }
    hElecPassCuts->Write();
    hMuonPassCuts->Write();
    hElecPassOtherCuts->Write();
    hMuonPassOtherCuts->Write();

    hIsoElec->Write();
    hIsoMuon->Write();
    hCumulElec->Write();
    hCumulMuon->Write();
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

void passElectronID(TElectron* elec, const float (&cuts)[10], vector<bool> &passCuts ) {
    if( fabs(elec->dEtaIn) < cuts[0] ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( fabs(elec->dPhiIn) < cuts[1] ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( elec->sieie  < cuts[2] ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( elec->hovere < cuts[3] ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( fabs(elec->d0) < cuts[4] ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( fabs(elec->dz) < cuts[5] ) passCuts.push_back(true);
    else passCuts.push_back(false);
//             fabs(elec->eoverp)  < cuts[6] && // eoverp = E/p  -> missing
    if( ((elec->chHadIso03 + max(elec->gammaIso03+elec->neuHadIso03-0.5* elec->puIso03,0.0))/elec->pt) < cuts[7] ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( (!elec->isConv) ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( elec->nMissingHits <= cuts[9] ) passCuts.push_back(true); 
    else passCuts.push_back(false);
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
           ((muon->chHadIso04 + max(muon->gammaIso04+muon->neuHadIso04-0.5* muon->puIso04,0.0))/muon->pt) < 0.3 // muon iso  
    );  
}

void passTightMuonID(TMuon* muon, vector<bool> &passCuts ) {
    if( ((muon->typeBits)/2)%2 ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( ((muon->typeBits)/32)%2 ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( muon->tkNchi2 < 10 ) passCuts.push_back(true);
    else passCuts.push_back(false); //
    if( muon->nValidHits > 0 ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( muon->nMatchStn > 1 ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( fabs(muon->d0) < 0.2 ) passCuts.push_back(true);
    else passCuts.push_back(false); // d0 = -dxy
    if( fabs(muon->dz) < 0.5 ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( muon->nPixHits > 0 ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( muon->nTkLayers > 5 ) passCuts.push_back(true);
    else passCuts.push_back(false);
    if( ((muon->chHadIso04 + max(muon->gammaIso04+muon->neuHadIso04-0.5* muon->puIso04,0.0))/muon->pt) < 0.3 ) passCuts.push_back(true);
    else passCuts.push_back(false);
}

bool passLooseMuonID(TMuon* muon) {
    return( ((muon->typeBits)/32)%2 &&
            ( ((muon->typeBits)/2)%2 || ((muon->typeBits)/4)%2 ) &&
            ((muon->chHadIso04 + max(muon->gammaIso04+muon->neuHadIso04-0.5* muon->puIso04,0.0))/muon->pt) < 0.2 // muon iso  
          );
}
        
        