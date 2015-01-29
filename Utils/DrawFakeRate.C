// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH2F.h"
#include "TCanvas.h"

void DrawFakeRate() {
    // Jet -> Electron
    TFile *inFile = new TFile("FakeRate_rebin.root","READ");
    TH2F* eeDen = (TH2F*) inFile->Get("Denominator_W_to_e_jet_to_e");
    TH2F* meDen = (TH2F*) inFile->Get("Denominator_W_to_mu_jet_to_e");
    TH2F* eeNum = (TH2F*) inFile->Get("Lepton_Numerator_W_to_e_jet_to_e");
    TH2F* meNum = (TH2F*) inFile->Get("Lepton_Numerator_W_to_mu_jet_to_e");
    TH2F* eeJetNum = (TH2F*) inFile->Get("Jet_Numerator_W_to_e_jet_to_e");
    TH2F* meJetNum = (TH2F*) inFile->Get("Jet_Numerator_W_to_mu_jet_to_e");
    
    eeNum->Add(meNum);
    eeDen->Add(meDen);
    eeJetNum->Add(meJetNum);

    eeNum->GetXaxis()->SetTitle("|eta|");
    eeNum->GetYaxis()->SetTitle("pt");
    eeNum->SetTitle("Fake rate jet->e");
    eeNum->SetMarkerSize(1.2);
    eeJetNum->GetXaxis()->SetTitle("|eta|");
    eeJetNum->GetYaxis()->SetTitle("pt");
    eeJetNum->SetTitle("Fake rate jet->e-jet");
    eeJetNum->SetMarkerSize(1.2);

    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    eeNum->Draw("text");
    c1->Print("Numerator_jet_to_e.png","png");
    eeNum->Draw("colz");
    c1->Print("Numerator_jet_to_e_colz.png","png");
    
    TCanvas *c2 = new TCanvas("c2","c2",1200,600);
    eeJetNum->Draw("text");
    c2->Print("Numerator_jet_to_ejet.png","png");
    eeJetNum->Draw("colz");
    c2->Print("Numerator_jet_to_ejet_colz.png","png");
    
    eeNum->Divide(eeDen);
    eeJetNum->Divide(eeDen);
    
    TCanvas *c3 = new TCanvas("c3","c3",1200,600);
    eeNum->Draw("text");
    c3->Print("fake_rate_jet_to_e.png","png");
    eeNum->Draw("colz");
    c3->Print("fake_rate_jet_to_e_colz.png","png");
    
    TCanvas *c4 = new TCanvas("c4","c4",1200,600);
    eeJetNum->Draw("text");
    c4->Print("fake_rate_jet_to_ejet.png","png");
    eeJetNum->Draw("colz");
    c4->Print("fake_rate_jet_to_ejet_colz.png","png");
    
    TCanvas *c5 = new TCanvas("c5","c5",1200,600);
    TH2F* eeMig = (TH2F*) inFile->Get("Pt_migration_W_to_e_jet_to_e");
    TH2F* meMig = (TH2F*) inFile->Get("Pt_migration_W_to_mu_jet_to_e");
    eeMig->Add(meMig);
    eeMig->GetXaxis()->SetTitle("jet pt");
    eeMig->GetYaxis()->SetTitle("p_{t,elec} - p_{t,jet} ");
    eeMig->SetTitle("Pt migration jet -> e");
    eeMig->Draw("text");
    c5->Print("Pt_migration_jet_to_e.png","png");
    
    
    // Jet -> Muon
    TH2F* emDen = (TH2F*) inFile->Get("Denominator_W_to_e_jet_to_mu");
    TH2F* mmDen = (TH2F*) inFile->Get("Denominator_W_to_mu_jet_to_mu");
    TH2F* emNum = (TH2F*) inFile->Get("Lepton_Numerator_W_to_e_jet_to_mu");
    TH2F* mmNum = (TH2F*) inFile->Get("Lepton_Numerator_W_to_mu_jet_to_mu");
    TH2F* emJetDen = (TH2F*) inFile->Get("Jet_Numerator_W_to_e_jet_to_mu");
    TH2F* mmJetDen = (TH2F*) inFile->Get("Jet_Numerator_W_to_mu_jet_to_mu");
    TH2F* emJetNum = (TH2F*) inFile->Get("Jet_Denominator_W_to_e_jet_to_mu");
    TH2F* mmJetNum = (TH2F*) inFile->Get("Jet_Denominator_W_to_mu_jet_to_mu");
    
    emNum->Add(mmNum);
    emDen->Add(mmDen);
    emJetNum->Add(mmJetNum);
    emJetDen->Add(mmJetDen);

    emNum->GetXaxis()->SetTitle("|eta|");
    emNum->GetYaxis()->SetTitle("pt");
    emNum->SetTitle("Fake rate jet->mu");
    emNum->SetMarkerSize(1.2);
    emJetNum->GetXaxis()->SetTitle("|eta|");
    emJetNum->GetYaxis()->SetTitle("pt");
    emJetNum->SetTitle("Fake rate jet->mu-jet");
    emJetNum->SetMarkerSize(1.2);

    TCanvas *c6 = new TCanvas("c6","c6",1200,600);
    emNum->Draw("text");
    c6->Print("Numerator_jet_to_mu.png","png");
    emNum->Draw("colz");
    c6->Print("Numerator_jet_to_mu_colz.png","png");
    
    TCanvas *c7 = new TCanvas("c7","c7",1200,600);
    emJetNum->Draw("text");
    c7->Print("Numerator_jet_to_mujet.png","png");
    emJetNum->Draw("colz");
    c7->Print("Numerator_jet_to_mujet_colz.png","png");
    
    emNum->Divide(emDen);
    emJetNum->Divide(emJetDen);
    
    TCanvas *c8 = new TCanvas("c8","c8",1200,600);
    emNum->Draw("text");
    c8->Print("fake_rate_jet_to_mu.png","png");
    emNum->Draw("colz");
    c8->Print("fake_rate_jet_to_mu_colz.png","png");
    
    TCanvas *c9 = new TCanvas("c9","c9",1200,600);
    emJetNum->Draw("text");
    c9->Print("fake_rate_jet_to_mujet.png","png");
    emJetNum->Draw("colz");
    c9->Print("fake_rate_jet_to_mujet_colz.png","png");
    
    TCanvas *c10 = new TCanvas("c10","c10",1200,600);
    TH2F* emMig = (TH2F*) inFile->Get("Pt_migration_W_to_e_jet_to_mu");
    TH2F* mmMig = (TH2F*) inFile->Get("Pt_migration_W_to_mu_jet_to_mu");
    emMig->Add(mmMig);
    emMig->GetXaxis()->SetTitle("jet pt");
    emMig->GetYaxis()->SetTitle("p_{t,elec} - p_{t,jet} ");
    emMig->SetTitle("Pt migration jet -> mu");
    emMig->Draw("text");
    c10->Print("Pt_migration_jet_to_mu.png","png");
}