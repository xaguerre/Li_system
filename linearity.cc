#include <math.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TTree.h"
#include <fstream>
#include <string>
#include <cstring>
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TH2F.h"
#include "TApplication.h"
#include "TMultiGraph.h"
#include "TFeldmanCousins.h"
#include "TGaxis.h"
#include "TLeaf.h"
using namespace std;


TH1D* spectre_charge(int om_number, int file_number ){
  TFile *file = new TFile(Form("histo_brut/histo_charge_amplitude_energie_%i.root", file_number), "READ");
  gROOT->cd();
  TH2F* charge = (TH2F*)file->Get("histo_pm_charge");

  TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  file->Close();
  return spectre_charge;
}

void txt_to_root() {

  std::ifstream parametres("/home/aguerre/Bureau/Thèse/Li_system/Resultats_txt/Resultats_charge.txt");
  TFile file("test.root","RECREATE");

  int compteur = 0;
  int i_om;
  double n_evt;
  double mean_charge;
  double sigma;
  double nbg;
  double C = 0;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("i_om", &i_om);
  Result_tree.Branch("n_evt", &n_evt);
  Result_tree.Branch("mean_charge", &mean_charge);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("nbg", &nbg);
  Result_tree.Branch("C", &C);

  while(parametres >> i_om)
  {
    parametres >> n_evt >> mean_charge >> sigma >> nbg >> C;
    Result_tree.Fill();
    compteur++;
  }

  file.cd();
  Result_tree.Write();
  file.Close();

}

void Get_Mean() {

  TFile file("Resultats_root/getMean.root","RECREATE");

  int i_om =0;
  double mean =0;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("i_om", &i_om);
  Result_tree.Branch("mean", &mean);

  TFile tree_file("histo_brut/histo_Li_system_568.root","READ");
  int om_number;
  double time;
  double charge_tree;
  double amplitude_tree;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);
  tree->SetBranchStatus("charge_tree",1);
  tree->SetBranchAddress("charge_tree", &charge_tree);
  tree->SetBranchStatus("amplitude_tree",1);
  tree->SetBranchAddress("amplitude_tree", &amplitude_tree);

  for (int om = 0; om < 712; om++) {
    double temps = 138;
    if ((om > 259 && om < 520) || (om > 583 && om < 647) || (om > 679)) {temps=386;}
    TH1D *spectre = new TH1D ("spectre_amplitude", "", 700, 0, 2300 );
    tree->Project("spectre_amplitude", "amplitude_tree", Form("om_number == %d && time > %f && time < %f", om, temps, temps+34));

    mean = spectre->GetMean();
    i_om = om;
    Result_tree.Fill();
    delete spectre;
  }

  file.cd();
  Result_tree.Write();
  file.Close();
}

int pic_number(double temps){
  int pic_number = 0;
  if (floor(temps/40) == 0 || floor(temps/40) == 6) {
    pic_number = 1;
  }
  if (floor(temps/40) == 1 || floor(temps/40) == 7) {
    pic_number = 2;
  }
  if (floor(temps/40) == 2 || floor(temps/40) == 8) {
    pic_number = 3;
  }
  if (floor(temps/40) == 3 || floor(temps/40) == 9) {
    pic_number = 4;
  }
  if (floor(temps/40) == 4 || floor(temps/40) == 10) {
    pic_number = 5;
  }
  if (floor(temps/40) == 5 || floor(temps/40) == 11) {
    pic_number = 6;
  }
  return pic_number;
}

int intensity_chooser(int pic){
  int intensity;
  if (pic == 1) {
    intensity = 70;
  }
  if (pic == 2) {
    intensity = 80;
  }
  if (pic == 3) {
    intensity = 90;
  }
  if (pic == 4) {
    intensity = 100;
  }
  if (pic == 5) {
    intensity = 110;
  }
  if (pic == 6) {
    intensity = 120;
  }
  return intensity;
}

void fit_LI_amplitude(int run_number){
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TFile file("Resultats_root/Resultat_amplitude_Li.root","RECREATE");
  std::ofstream outFile("Resultats_txt/failed_fit.txt");

  int i_om;
  double constante;
  double mean;
  double mean_error;
  double sigma;
  double error;
  int pic =0;
  double Khi2 = 0;
  int intensity = 0;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("i_om", &i_om);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("constante", &constante);
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_error", &mean_error);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("intensity", &intensity);

  double debut = 17;
  TCanvas* canvas = new TCanvas;
  TH1F comp_map("comp_map", "comp_map", 100, 0, 100);
  TFile tree_file(Form("histo_brut/histo_Li_system_%d.root", run_number), "READ");
  gROOT->cd();
  int om_number;
  double time;
  double charge_tree;
  double amplitude_tree;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);
  tree->SetBranchStatus("charge_tree",1);
  tree->SetBranchAddress("charge_tree", &charge_tree);
  tree->SetBranchStatus("amplitude_tree",1);
  tree->SetBranchAddress("amplitude_tree", &amplitude_tree);

  for(int om = 78; om < 176; om+=1)
  {
    for (double j = debut; j < debut+6*40.5; j = j+40.5)
    {
      double temps = j;
      if ((om > 259 && om < 520) || (om > 583 && om < 647) || (om > 679)) {temps=(265+j-debut);}
      TH1D *spectre = new TH1D ("spectre_amplitude", "", 700, 0, 2300 );
      tree->Project("spectre_amplitude", "amplitude_tree", Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, temps, ceil(temps+32.5)));
      std::cout << "temps = " << temps << " - " << ceil(temps+32.5) << '\n';
      if (spectre->GetEntries() < 300) {
        std::cout << "" << '\n';
        std::cout << "trop peu d'entries pour l'OM " << om << '\n';
        std::cout << "" << '\n';
        delete spectre;
      }
      else if (spectre->GetMean() > 1900) {
        std::cout << "" << '\n';
        std::cout << "the amplitude sature" << '\n';
        std::cout << "" << '\n';
        // outFile << "om = " << om << " and pic = " << floor(temps/40) <<  endl;
        delete spectre;
      }
      else if (spectre->GetMean() < 10){
        std::cout << "" << '\n';
        std::cout << "too few charge" << '\n';
        std::cout << "" << '\n';
        // outFile << "om = " << om << " and pic = " << floor(temps/40) <<  endl;
        delete spectre;
      }
      else{
        TF1 *f_Gaus = new TF1("f_Gaus", "gaus(0)", 0, 1200);
        f_Gaus->SetParNames("N_evt","mean_charge","Sigma");
        // f_Gaus->SetParameters(25, spectre->GetMean(), 100);
        f_Gaus->SetParameters(25, spectre->GetMean(), spectre->GetRMS());
        f_Gaus->SetRange(spectre->GetMean()-400, spectre->GetMean()+400);
        f_Gaus->Draw("same");
        spectre->Fit(f_Gaus, "RQ0");
        f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
        spectre->Fit(f_Gaus, "RQ0");
        f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
        spectre->Fit(f_Gaus, "RQ0");

        pic = pic_number(temps);
        Khi2 = f_Gaus->GetChisquare()/f_Gaus->GetNDF();
        i_om = om;
        constante = (f_Gaus->GetParameter(0));
        mean = (f_Gaus->GetParameter(1));
        sigma = (f_Gaus->GetParameter(2));
        mean_error = f_Gaus->GetParError(1) ;
        intensity = intensity_chooser(pic);
        Result_tree.Fill();
        spectre->Draw();
        f_Gaus->Draw("same");
        canvas->SaveAs(Form("fit/fit_Li/amplitude_Li/OM_%03d_pic_%.0d.png", om, pic));
        delete spectre;
        delete f_Gaus;
      }
      temps = 0;
    }
  }

  file.cd();
  Result_tree.Write();
  file.Close();
  outFile.close();
  return;
}

int bundle_number(int i_om){
  int bundle_number = 0;
  if ((i_om%13 > 7 && i_om/13 < 10 && i_om < 260) || (i_om > 663 && i_om < 672) || (i_om < 552 && i_om > 545) || (i_om < 536 && i_om > 529)) {
    bundle_number = 1;
  }
  else if ((i_om%13 < 8 && i_om/13 < 6 && i_om < 260) ||(i_om > 647 && i_om < 652) ||(i_om < 546 && i_om > 535) ||(i_om < 530 && i_om > 519)){
    bundle_number = 2;
  }
  if ((i_om%13 > 7 && i_om/13 > 9 && i_om < 260) || (i_om > 671 && i_om < 680) || (i_om < 568 && i_om > 561) || (i_om < 584 && i_om > 577)) {
    bundle_number = 3;
  }
  else if ((i_om%13 < 8 && i_om/13 > 5 && i_om/13 < 14 && i_om < 260) || (i_om < 660 && i_om > 651)) {
    bundle_number = 4;
  }
  else if ((i_om%13 < 8 && i_om/13 > 13 && i_om < 260) || (i_om < 664 && i_om > 659) || (i_om < 562 && i_om > 551) ||(i_om < 578 && i_om > 567)) {
    bundle_number = 5;
  }
  else if ((i_om%13 > 7 && (i_om/13-20) < 10 && i_om < 520 && i_om > 259) ||(i_om > 695 && i_om < 704) ||(i_om < 600 && i_om > 593) ||(i_om < 616 && i_om > 609)){
    bundle_number = 6;
  }
  else if ((i_om%13 > 7 && (i_om/13-20) > 9 && i_om < 520 && i_om > 259) ||(i_om > 703 && i_om < 712) ||(i_om < 648 && i_om > 641) ||(i_om < 632 && i_om > 625)){
    bundle_number = 7;
  }
  else if ((i_om%13 < 8 && (i_om/13-20) < 6 && i_om < 520 && i_om > 259) ||(i_om > 679 && i_om < 684) ||(i_om < 594 && i_om > 583) ||(i_om < 610 && i_om > 599)){
    bundle_number = 8;
  }
  else if ((i_om%13 < 8 && (i_om/13-20) > 5 && (i_om/13-20) < 14 && i_om < 520 && i_om > 259) ||(i_om > 683 && i_om < 692)){
    bundle_number = 9;
  }
  else if ((i_om%13 < 8 && (i_om/13-20) > 13 && i_om < 520 && i_om > 259) ||(i_om > 691 && i_om < 696) ||(i_om < 642 && i_om > 631) ||(i_om < 626 && i_om > 615)){
    bundle_number = 10;
  }
  return bundle_number;
}

void Li_bundle_variation_new(){
  TFile file("Resultats_root/Li_bundle_variation.root","RECREATE");
  TFile tree_file("Resultats_root/Resultat_amplitude_Li.root","READ");

  double mean;
  int i_om;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("i_om",1);
  tree->SetBranchAddress("i_om", &i_om);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);

  TH1D *Li_bundle_1 = new TH1D ("Li_bundle_1", "", 1000, 0, 2500);
  TH1D *Li_bundle_2 = new TH1D ("Li_bundle_2", "", 1000, 0, 2500);
  TH1D *Li_bundle_3 = new TH1D ("Li_bundle_3", "", 1000, 0, 2500);
  TH1D *Li_bundle_4 = new TH1D ("Li_bundle_4", "", 1000, 0, 2500);
  TH1D *Li_bundle_5 = new TH1D ("Li_bundle_5", "", 1000, 0, 2500);
  TH1D *Li_bundle_6 = new TH1D ("Li_bundle_6", "", 1000, 0, 2500);
  TH1D *Li_bundle_7 = new TH1D ("Li_bundle_7", "", 1000, 0, 2500);
  TH1D *Li_bundle_8 = new TH1D ("Li_bundle_8", "", 1000, 0, 2500);
  TH1D *Li_bundle_9 = new TH1D ("Li_bundle_9", "", 1000, 0, 2500);
  TH1D *Li_bundle_10 = new TH1D ("Li_bundle_10", "", 1000, 0, 2500);
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (bundle_number(i_om) == 1){
      Li_bundle_1->Fill(mean);
    }
    if (bundle_number(i_om) == 2){
      Li_bundle_2->Fill(mean);
    }
    if (bundle_number(i_om) == 3){
      Li_bundle_3->Fill(mean);
    }
    if (bundle_number(i_om) == 4){
      Li_bundle_4->Fill(mean);
    }
    if (bundle_number(i_om) == 5){
      Li_bundle_5->Fill(mean);
    }
    if (bundle_number(i_om) == 6){
      Li_bundle_6->Fill(mean);
    }
    if (bundle_number(i_om) == 7){
      Li_bundle_7->Fill(mean);
    }
    if (bundle_number(i_om) == 8){
      Li_bundle_8->Fill(mean);
    }
    if (bundle_number(i_om) == 9){
      Li_bundle_9->Fill(mean);
    }
    if (bundle_number(i_om) == 10){
      Li_bundle_10->Fill(mean);
    }
  }

  file.cd();
  Li_bundle_1->Write();
  Li_bundle_2->Write();
  Li_bundle_3->Write();
  Li_bundle_4->Write();
  Li_bundle_5->Write();
  Li_bundle_6->Write();
  Li_bundle_7->Write();
  Li_bundle_8->Write();
  Li_bundle_9->Write();
  Li_bundle_10->Write();
  file.Close();
}

void om_linearity_prep(string file) {

  int i_om;
  double mean;
  double mean_error;
  int intensity;
  double Khi2;
  double a;
  double b;

  TFile tree_file(Form("Resultats_root/Resultat_%s.root", file.c_str()),"READ");
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("i_om",1);
  tree->SetBranchAddress("i_om", &i_om);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("mean_error",1);
  tree->SetBranchAddress("mean_error", &mean_error);
  tree->SetBranchStatus("intensity",1);
  tree->SetBranchAddress("intensity", &intensity);

  TH3D* linearity = new TH3D("OM_linearity", "OM_linearity",
  712, 0, 712,
  200, 0, 2000,
  100, 70, 120);

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    linearity->Fill(i_om, mean, intensity);
  }

  TFile newfile("Resultats_root/linearity_charge_prep.root","RECREATE");

  newfile.cd();
  linearity->Write();
  newfile.Close();
}

void om_linearity(){
  TFile file("Resultats_root/linearity_charge_prep.root","READ");
  TH3D* om_linearity_prep = (TH3D*)file.Get("OM_linearity");
  om_linearity_prep->Draw();
  int i_om;
  double Khi2;
  double a;
  double b;

  TF1 *fit = new TF1("fit", "pol1", 60,130);

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("i_om", &i_om);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("a", &a);
  Result_tree.Branch("b", &b);

  for (int i = 0; i < 712; i++) {
    TCanvas *canvas = new TCanvas();
    om_linearity_prep->Draw();

    om_linearity_prep->GetXaxis()->SetRange(i+1, i+1);
    TH1 *polyfit = om_linearity_prep->Project3D("yz");
    polyfit->Fit("fit", "RQ");

    i_om = i;
    a = fit->GetParameter(0);
    b = fit->GetParameter(1);
    Khi2 = fit->GetChisquare()/fit->GetNDF();
    Result_tree.Fill();
    canvas->SaveAs(Form("fit_linearite/linearite_OM_%d", i));
    delete canvas;
  }
  TFile newfile("Resultats_root/linearity_charge.root","RECREATE");
  newfile.cd();
  Result_tree.Write();
  newfile.Close();

}

void ecart( ) {
  TFile tree_file("Resultats_root/linearity_charge.root","READ");

  int i_om;
  double a;
  double b;
  double Khi2;
  double mean;
  int intensity;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("a",1);
  tree->SetBranchAddress("a", &a);
  tree->SetBranchStatus("b",1);
  tree->SetBranchAddress("b", &b);

  double facteur[712][2];

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    facteur[i][0] = a;
    facteur[i][1] = b;
  }
  string file = "amplitude_Li";
  TFile tree_file2(Form("Resultats_root/Resultat_%s.root", file.c_str()),"READ");
  TTree* tree2 = (TTree*)tree_file2.Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("i_om",1);
  tree2->SetBranchAddress("i_om", &i_om);
  tree2->SetBranchStatus("mean",1);
  tree2->SetBranchAddress("mean", &mean);
  tree2->SetBranchStatus("intensity",1);
  tree2->SetBranchAddress("intensity", &intensity);

  int om, intensite, bundle;
  double ecart;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om", &om);
  Result_tree.Branch("ecart", &ecart);
  Result_tree.Branch("intensite", &intensite);
  Result_tree.Branch("bundle", &bundle);
  int compteur = 0;
  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    om = i_om;
    compteur++;

    ecart = 100*(mean - (facteur[i_om][1]*intensity + facteur[i_om][0]))/(facteur[i_om][1]*intensity + facteur[i_om][0]);
    intensite = intensity;
    bundle = bundle_number(i_om);
    std::cout << "om = " << om << " and bundle = " << bundle << '\n';
    Result_tree.Fill();
  }
  TFile newfile("Resultats_root/linearity_charge_ecart.root","RECREATE");
  newfile.cd();
  Result_tree.Write();
  newfile.Close();
}
