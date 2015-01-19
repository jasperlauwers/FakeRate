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
    bool eLepton = true; // eLepton = true: W -> e
    bool eJet = false; // eJet = true: jet -> e
    bool tightSel = true; // tight or loose lepton selection
    float leptonPt = 80;
    
    // ID cuts
    // tightElCuts[barrel/endcap][IDcuts], CSA14 selection, conditions: 50ns, poor detector alignment
    const float tightElCuts[2][10] = {{0.012,0.024,0.01,0.074,0.0091,0.017,0.026,0.10,1e-6,1},{0.019,0.043,0.029,0.08,0.037,0.065,0.076,0.14,1e-6,1}};
    const float looseElCuts[2][10] = {{0.0181,0.0936,0.0123,0.141,0.0166,0.54342,0.1353,0.24,1e-6,1},{0.0124,0.0642,0.035,0.1115,0.098,0.9187,0.1443,0.3529,1e-6,1}};    
    // CSA14 selection, conditions: 25ns, better detector alignment
//     const float tightElCuts[2][10] = {{0.0091,0.031,0.0106,0.0532,0.0126,0.0116,0.0609,0.1649,1e-6,1},{0.0106,0.0359,0.0305,0.0835,0.0163,0.5999,0.1126,0.2075,1e-6,1}};
//     const float looseElCuts[2][10] = {{0.0181,0.0936,0.0123,0.141,0.0166,0.54342,0.1353,0.24,1e-6,1},{0.0124,0.0642,0.035,0.1115,0.098,0.9187,0.1443,0.3529,1e-6,1}};
    
    // Directories
    int maxInFiles=10;
    TString outDirPNG = "/afs/cern.ch/user/j/jlauwers/www/protected/VBS/TP/FakeRate/";
    TString outDirROOT = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/Results/";
    TString inDir = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/CSA14/WJetsToLNu_13TeV-madgraph-pythia8-tauola_v3/";
//     TString inDir = "/afs/cern.ch/work/j/jlauwers/VBS/TP/FakeRate/BaconTrees/";
    
    // Verbose output
    int verbose = 1; // 0: no messages - 3: all debug information 
    // ----------------------------------------
    
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
    tree->SetBranchStatus("rhoIso",1);
    
    TClonesArray *fElectron = new TClonesArray("baconhep::TElectron");
    TClonesArray *fMuon = new TClonesArray("baconhep::TMuon");
    TEventInfo *eventInfo=0;
    tree->SetBranchAddress("Electron", &fElectron);
    tree->SetBranchAddress("Muon", &fMuon);
    tree->SetBranchAddress("Info", &eventInfo);
    
    TString strSel = "W_to";
    if( eLepton ) strSel += "_e";
    else strSel += "_mu";
    if( eJet ) strSel += "_jet_to_e";
    else strSel += "_jet_to_mu";
    
    TH1::SetDefaultSumw2();
    
    // Iso
    TH1F *hIsoElec = new TH1F("isoRhoCorr_elec","isoRhoCorr_elec",200,0,2);
    TH1F *hIsoMuon = new TH1F("isoRhoCorr_muon","isoRhoCorr_muon",200,0,2);
    TH1F *hCumulElec = new TH1F("cumulative_isoRhoCorr_elec","cumulative_isoRhoCorr_elec",200,0,2);
    TH1F *hCumulMuon = new TH1F("cumulative_isoRhoCorr_muon","cumulative_isoRhoCorr_muon",200,0,2);
    
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
    TFile* outFile = new TFile(outDirROOT+"FakeRate_iso.root","UPDATE");
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
        
        