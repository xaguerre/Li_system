#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstring>
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TApplication.h"
#include "TMultiGraph.h"
#include "TFeldmanCousins.h"
#include "TGaxis.h"
#include "TLeaf.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <boost/filesystem.hpp>
using namespace std;

TH2F* amplitude_spectre = NULL;

void Load_spectre(int run_number){
  TFile *file = new TFile(Form("histo_brut/histo_charge_amplitude_energie_%d.root", run_number), "READ");
  gROOT->cd();
  amplitude_spectre = (TH2F*)file->Get("histo_pm_amplitude");
  // amplitude_spectre->RebinY(4);
  return;
}

TH1D* spectre_amplitude(int om_number){
  TH1D* spectre_amplitude = amplitude_spectre->ProjectionY(Form("amplitude%03d",om_number), om_number+1, om_number+1);
  return spectre_amplitude;
}

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

void fit_all_om_charge(){
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  std::ofstream outFile("Resultats_txt/Resultats_charge.txt");
  std::ifstream parametres("/home/aguerre/Bureau/Thèse/Li_system/Resultats_txt/Resultats_charge_run_547.txt");

  TFile file("Resultats_root/Charge_Li_run_547_557.root","RECREATE");

  double chi = 0;
  int i_om;
  double n_evt;
  double mean_charge;
  double mean_error;
  double sigma;
  double  nbg;
  double C = 0;
  int run_number;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("i_om", &i_om);
  Result_tree.Branch("n_evt", &n_evt);
  Result_tree.Branch("mean_charge", &mean_charge);
  Result_tree.Branch("mean_error", &mean_error);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("C", &C);
  Result_tree.Branch("nbg", &nbg);
  Result_tree.Branch("chi", chi);
  int tab[4] = {547, 551, 554, 557};

  int om_tab[712];
  double N_evt_tab[712];
  double mean_tab[712];
  double sigma_tab[712];
  double Nbg_tab[712];
  double C_tab[712];
  int compteur = 0;

  while(parametres >> om_tab[compteur])  // tant que l'on peut mettre la ligne dans "contenu"
  {
    parametres >> N_evt_tab[compteur] >> mean_tab[compteur] >> sigma_tab[compteur] >> Nbg_tab[compteur] >> C_tab[compteur];
    compteur++;
  }

  TCanvas* canvas = new TCanvas;
  canvas->SetLogy();

  TH2F mean_charge_map("histo_om_mean_charge_map", "mean_charge_map", 20, 0, 20, 13, 0, 13);
  for (int run = 0; run < 4; run++) {
    for(int om = 0; om < 520; om+=1){
      run_number = tab[run];
      TH1D* spectre_om = spectre_charge(om_tab[om], run_number);
      spectre_om->Draw();
      TF1* f_ComptonEdgePoly = new TF1 ("f_ComptonEdgePoly","[0]*(0.5*(1+TMath::Erf(([1]-x)/(TMath::Sqrt(2)*[2]))) + [3]*x + [4])", 40000, 90000);
      f_ComptonEdgePoly->SetParNames("N_evt","mean_charge","Sigma", "Nbg", "C" );

      std::cout << "om = " << om_tab[om] << " and mean = " << mean_tab[om] << " and C = " << C_tab[om] << '\n';

      f_ComptonEdgePoly->SetParameters(N_evt_tab[om], mean_tab[om], sigma_tab[om], Nbg_tab[om]/N_evt_tab[om], C_tab[om]/N_evt_tab[om]);
      f_ComptonEdgePoly->SetRange(mean_tab[om]-2.5*sigma_tab[om], mean_tab[om]+2.5*sigma_tab[om]);
      f_ComptonEdgePoly->SetParLimits(0, 0, 200);
      f_ComptonEdgePoly->SetParLimits(1, 0, 200000);
      f_ComptonEdgePoly->SetParLimits(2, 0, 10000);
      f_ComptonEdgePoly->SetParLimits(3, -10, 0);
      f_ComptonEdgePoly->SetParLimits(4, 0, 20);
      f_ComptonEdgePoly->Draw("same");
      spectre_om->Fit(f_ComptonEdgePoly, "RQ");
      f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
      spectre_om->Fit(f_ComptonEdgePoly, "RQ");
      f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
      spectre_om->Fit(f_ComptonEdgePoly, "RQ");

      i_om = om_tab[om];
      n_evt = (f_ComptonEdgePoly->GetParameter(0));
      mean_charge = (f_ComptonEdgePoly->GetParameter(1));
      sigma = (f_ComptonEdgePoly->GetParameter(2));
      nbg = (f_ComptonEdgePoly->GetParameter(3));
      C = (f_ComptonEdgePoly->GetParameter(4));
      mean_error = f_ComptonEdgePoly->GetParError(1);
      chi = (f_ComptonEdgePoly->GetChisquare()/f_ComptonEdgePoly->GetNDF());
      std::cout << "chi = " << chi << '\n';
      Result_tree.Fill();
      //mapping
      int om_col = (om % 13 );
      int om_row = (om / 13);
      mean_charge_map.SetBinContent( om_row+1, om_col+1, mean_charge);

      canvas->SaveAs(Form("fit/fit_Tl/test/charge_fit_om_%03d_run_%d.png", om_tab[om], run_number));

      outFile << i_om << "\t" << n_evt << "\t" << mean_charge << "\t" << sigma << "\t"<< nbg << "\t" << C << endl;

      delete spectre_om;
      delete f_ComptonEdgePoly;
    }
  }
    double comp[520];

    for (int i = 0; i < 520; i++) {
      Result_tree.GetEntry(i);
      comp[i] = mean_charge;
      // error1[i] = mean_error;
    }
    file.cd();

    double yaxis[4];
    double yaxis_error[4];
    double xaxis[4] = {0, 19, 38, 57};
    double xaxiserror[4] = {0.5, 0.5, 0.5, 0.5};
    for (int j = 0; j < 520; j++){
      for (int i = 0; i < 4; i++) {
        int nombre = 520*i+j;
        // int nombre = i;
        Result_tree.GetEntry(nombre);
        yaxis[i] = mean_charge/comp[j];
        // std::cout << "om = " << j << " and var = " << yaxis << '\n';
        if ((yaxis[i] < 1.1) && (yaxis[i] > 1.05)) {

          std::cout << " run == " << tab[i] << " om = " << j << " and var = " << yaxis[i] << '\n';
        }
        // yaxis_error[i] = mean_error/comp*1.0 + (mean_charge/(comp*comp))*1.0*error1;
      }
      TGraphErrors comp_map (4, xaxis, yaxis, xaxiserror, yaxis_error);
      comp_map.SetName(Form("fit_Tl_om_%d", j));
      comp_map.SetNameTitle(Form("fit_Tl_om_%d", j), Form("evolution du gain de l'OM %d", j));
      comp_map.GetXaxis()->SetTitle("Temps (h)");
      comp_map.GetYaxis()->SetTitle("Gain(t)/Gain(0)");
      comp_map.GetXaxis()->SetTimeDisplay(1);
      comp_map.SetMarkerColor(2);
      comp_map.SetMarkerStyle(34);
      comp_map.SetMarkerSize(2);


      // TCanvas* canvas2 = new TCanvas;
      // comp_map.Draw();
      // canvas2->SaveAs(Form("fit/fit_Tl/variation/charge_fit_om_%03d.png", j));
      comp_map.Write();

    }




    mean_charge_map.Write();
    Result_tree.Write();

    file.Close();
    outFile.close();

    return;
  }

double* om_gain_fit(int om, int run_number){
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  TFile file(Form("Resultats_root/Ampl_Tl_run_%d.root", run_number),"RECREATE");

  double chi = 0;
  int i_om;
  double n_evt;
  double mean_error;
  double sigma;
  double  nbg;
  double C = 0;
  double ndf = 0;
  double chin = 0;
  double mean = 0;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("i_om", &i_om);
  Result_tree.Branch("n_evt", &n_evt);
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_error", &mean_error);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("C", &C);
  Result_tree.Branch("nbg", &nbg);
  Result_tree.Branch("chi", chi);

  double* tab = new double[3];
  TCanvas* canvas = new TCanvas;
  canvas->SetLogy();
  TH1D* spectre_om = NULL;
  spectre_om = spectre_amplitude(om);
  spectre_om->Draw();
  TF1* f_ComptonEdgePoly = new TF1 ("f_ComptonEdgePoly","[0]/2.0*(1+TMath::Erf(([1]-x)/(TMath::Sqrt(2)*[2])))+[3]*x", 400, 900);
  f_ComptonEdgePoly->SetParNames("N_evt","Mean","Sigma","Nbg" );

  if ((om % 13) == 12 )        //om multiple de (13)-1
  {
    f_ComptonEdgePoly->SetParameters(120, 723, 68, 3.91e-5);
    f_ComptonEdgePoly->SetRange(600,1000);
    f_ComptonEdgePoly->Draw("same");
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
  }
  else if ((om % 13) == 0)       //om multiple de 13
  {
    f_ComptonEdgePoly->SetParameters(112, 681, 56, 1.2e-05);
    f_ComptonEdgePoly->SetRange(500,1000);
    f_ComptonEdgePoly->Draw("same");
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-1.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
  }
  else         //om normaux (8pouces)
  {
    f_ComptonEdgePoly->SetParameters(111, 609, 37, 4.19e-05);
    f_ComptonEdgePoly->SetRange(450,1000);
    f_ComptonEdgePoly->Draw("same");
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+7.5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-3.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+7.5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");

  }
canvas->SaveAs(Form("fit/fit_Tl/amplitude/amplitude_fit_om_%d_run_.png", om));

  chi = (f_ComptonEdgePoly->GetChisquare());
  std::cout << chi << '\n';
  ndf = (f_ComptonEdgePoly->GetNDF());
  chin = (f_ComptonEdgePoly->GetChisquare()/f_ComptonEdgePoly->GetNDF());

  if (chin < 1.5 && mean > 70) {
    n_evt = f_ComptonEdgePoly->GetParameter(0);
    mean = (f_ComptonEdgePoly->GetParameter(1))/2.6;
    sigma = (f_ComptonEdgePoly->GetParameter(2))/2.6;
    nbg = f_ComptonEdgePoly->GetParameter(3);
  }
  else{
    n_evt = 0;
    mean = 0;
    sigma = 0;
    nbg = 0;
  }
  delete f_ComptonEdgePoly;
  delete canvas;
  tab[0] = mean;
  tab[1] = mean_error;
  tab[2] = chin;
  return tab;
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

void fit_LI_amplitude(int run_number, int time2 = 0){
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TFile file(Form("Resultats_root/Amplitude_Li_run_%d.root", run_number),"RECREATE");
  std::ofstream outFile("Resultats_txt/failed_fit.txt");

  int i_om;
  double constante;
  double mean, time;
  double mean_error, nevent;
  double sigma;
  int pic =0;
  double Khi2 = 0;
  int intensity = 0;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("i_om", &i_om);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("constante", &constante);
  Result_tree.Branch("Amplitude", &mean);
  Result_tree.Branch("Amplitude_error", &mean_error);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("intensity", &intensity);
  Result_tree.Branch("time", &time2);
  Result_tree.Branch("nevent", &nevent);
  double debut = 0;
  TCanvas* canvas = new TCanvas;
  TH1F comp_map("comp_map", "comp_map", 100, 0, 100);
  TFile tree_file(Form("histo_brut/Li_system_%d.root", run_number), "READ");
  int om_number;
  double charge_tree;
  double amplitude_tree;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  gROOT->cd();
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);
  tree->SetBranchStatus("charge_tree",1);
  tree->SetBranchAddress("charge_tree", &charge_tree);
  tree->SetBranchStatus("amplitude_tree",1);
  tree->SetBranchAddress("amplitude_tree", &amplitude_tree);
  int n_evt = 0;
  while (n_evt < 100) {
    TH1D *spectre = new TH1D ("spectre_amplitude", "", 700, 0, 2300);
    tree->Project("spectre_amplitude", "amplitude_tree", Form("om_number == 800 && time < %f  && amplitude_tree > 10", debut));
    debut++;

    n_evt = spectre->GetEntries();
    // std::cout << "n_evt = " << n_evt << " debut = " << debut << '\n';
    delete spectre;
  }

  std::cout << "debut = " << debut << '\n';


  for(int om = 800; om < 805; om+=1)
  {
    if (om ==712) {
      om = 800;
    }
    for (double j = debut; j < debut+6*40.5; j = j+40.5)
    {
      double temps = j;
      if ((om > 259 && om < 520) || (om > 583 && om < 647) || (om > 679 && om < 712) ) {temps=(265+j-debut);}
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
        nevent = spectre->Integral();
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
        if (om < 712) {
          canvas->SaveAs(Form("fit/fit_Li/amplitude_Li/OM_%03d_pic_%d_run_%d.png", om, pic, run_number));
        }
        if (om > 712) {
          canvas->SaveAs(Form("fit/fit_Li/amplitude_Li_ref/OM_%03d_pic_%d_run_%d.png", om, pic, run_number));
        }
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

void fit_ref(int run_number, int time2 = 0) {
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TFile newfile(Form("Resultats_root/fit_ref_run_%d.root", run_number),"RECREATE");
  int om_number;
  double p0, p1, p2, p3, p4, p5, p6,time, p1_er, p4_er;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("p0", &p0);
  Result_tree.Branch("p1", &p1);
  Result_tree.Branch("p2", &p2);
  Result_tree.Branch("p3", &p3);
  Result_tree.Branch("p4", &p4);
  Result_tree.Branch("p5", &p5);
  Result_tree.Branch("p6", &p6);
  Result_tree.Branch("p4_er", &p4_er);
  Result_tree.Branch("p1_er", &p1_er);
  Result_tree.Branch("time", &time2);
  Result_tree.Branch("run_number", &run_number);

  TFile tree_file(Form("histo_brut/histo_Li_system_%d.root", run_number), "READ");
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
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);
  TCanvas* canvas = new TCanvas;

  for (int om = 800; om < 805; om++) {
    om_number = om;
    TH1D *spectre = new TH1D ("spectre_amplitude", "", 600, 0, 600) ;
    tree->Project("spectre_amplitude", "amplitude_tree", Form("om_number == %d", om));
    spectre->Draw();
    TF1* f_MultipleGaus = new TF1 ("f_MultipleGaus","[0]*(7.03*TMath::Gaus(x[0],[1],[2]*sqrt([1])) + 1.84*TMath::Gaus(x[0],([1]+([1]/976.0)*72),[2]*sqrt([1]))) + 0.54*TMath::Gaus(x[0],([1]+([1]/976.0)*84),[2]*sqrt([1])) + [3]*(1.52*TMath::Gaus(x[0],[4],[2]*sqrt([4])) + 0.44*TMath::Gaus(x[0],([4]+([1]/976.0)*73),[2]*sqrt([4])) + 0.15*TMath::Gaus(x[0],([4]+([1]/976.0)*85),[2]*sqrt([4]))) + exp(-[5]*x[0]+[6])", 0, 1200);
    f_MultipleGaus->SetParNames("P0","P1","P2","P3","P4","P5","P6");
    f_MultipleGaus->SetRange(60,180);
    f_MultipleGaus->SetParameters(312.8, 160.1, 1.675, 12798, 66.59, 0.01237, 9.224);
    f_MultipleGaus->Draw("same");
    spectre->Fit(f_MultipleGaus, "RQ0");
    p0 = f_MultipleGaus->GetParameter(0);
    p1 = f_MultipleGaus->GetParameter(1);
    p1_er = f_MultipleGaus->GetParError(1);

    std::cout << "par 1 = " << p1 << " +- " << p1_er << '\n';

    p2 = f_MultipleGaus->GetParameter(2);
    p3 = f_MultipleGaus->GetParameter(3);
    p4 = f_MultipleGaus->GetParameter(4);
    p4_er = f_MultipleGaus->GetParError(4);
    p5 = f_MultipleGaus->GetParameter(5);
    p6 = f_MultipleGaus->GetParameter(6);
    canvas->SaveAs(Form("fit/fit_ref/ref_fit_om_%d_run_%d.png", om, run_number));
    om_number = om;
    Result_tree.Fill();

    delete spectre;
  }
  newfile.cd();
  Result_tree.Write();
  newfile.Close();
}

// void fit_ref_Hea(int run_number) {
//   gStyle->SetOptFit(1);
//   gStyle->SetOptStat(1);
//   TFile newfile("Resultats_root/fit_ref.root","RECREATE");
//   int om_number;
//   double p0, p1, p2;
//   TTree Result_tree("Result_tree","");
//   Result_tree.Branch("om_number", &om_number);
//   Result_tree.Branch("p0", &p0);
//   Result_tree.Branch("p1", &p1);
//   Result_tree.Branch("p2", &p2);
//
//
//   TFile tree_file(Form("histo_brut/histo_Li_system_%d.root", run_number), "READ");
//   double time;
//   double charge_tree;
//   double amplitude_tree;
//   TTree* tree = (TTree*)tree_file.Get("Result_tree");
//   tree->SetBranchStatus("*",0);
//   tree->SetBranchStatus("om_number",1);
//   tree->SetBranchAddress("om_number", &om_number);
//   tree->SetBranchStatus("time",1);
//   tree->SetBranchAddress("time", &time);
//   tree->SetBranchStatus("charge_tree",1);
//   tree->SetBranchAddress("charge_tree", &charge_tree);
//   tree->SetBranchStatus("amplitude_tree",1);
//   tree->SetBranchAddress("amplitude_tree", &amplitude_tree);
//   TCanvas* canvas = new TCanvas;
//
//   for (int om = 800; om < 805; om++) {
//     om_number = om;
//     TH1D *spectre = new TH1D ("spectre_amplitude", "", 600, 0, 600) ;
//     tree->Project("spectre_amplitude", "amplitude_tree", Form("om_number == %d", om));
//     spectre->Draw();
//     TF1* f_Heaviside = new TF1 ("f_Heaviside","(1 + exp(x * [0])) / (1 + exp((x - [1]) / [2]))", 0, 1200);  // Heaviside for gamma
//     f_Heaviside->SetParNames("P0","P1","P2");
//     f_Heaviside->SetRange(60,180);
//     f_Heaviside->SetParameters(312.8, 160.1, 1.675);
//     f_Heaviside->Draw("same");
//     spectre->Fit(f_Heaviside, "RQ0");
//     p0 = f_Heaviside->GetParameter(0);
//     p1 = f_Heaviside->GetParameter(1);
//     p2 = f_Heaviside->GetParameter(2);
//     canvas->SaveAs(Form("fit/fit_ref/ref_fit_om_%d.png", om));
//     Result_tree.Fill();
//     delete spectre;
//   }
//   newfile.cd();
//   Result_tree.Write();
//   newfile.Close();
// }
//
// void fit_ref_tot(int run_number) {
//   gStyle->SetOptFit(1);
//   gStyle->SetOptStat(1);
//   TFile newfile("Resultats_root/fit_ref.root","RECREATE");
//   int om_number;
//   double p0, p1, p2, p3, p4, p5, p6;
//   TTree Result_tree("Result_tree","");
//   Result_tree.Branch("om_number", &om_number);
//   Result_tree.Branch("p0", &p0);
//   Result_tree.Branch("p1", &p1);
//   Result_tree.Branch("p2", &p2);
//
//
//   TFile tree_file(Form("histo_brut/histo_Li_system_%d.root", run_number), "READ");
//   double time;
//   double charge_tree;
//   double amplitude_tree;
//   TTree* tree = (TTree*)tree_file.Get("Result_tree");
//   tree->SetBranchStatus("*",0);
//   tree->SetBranchStatus("om_number",1);
//   tree->SetBranchAddress("om_number", &om_number);
//   tree->SetBranchStatus("time",1);
//   tree->SetBranchAddress("time", &time);
//   tree->SetBranchStatus("charge_tree",1);
//   tree->SetBranchAddress("charge_tree", &charge_tree);
//   tree->SetBranchStatus("amplitude_tree",1);
//   tree->SetBranchAddress("amplitude_tree", &amplitude_tree);
//   TCanvas* canvas = new TCanvas;
//
//   for (int om = 800; om < 801; om++) {
//     om_number = om;
//     TH1D *spectre = new TH1D ("spectre_amplitude", "", 600, 0, 600) ;
//     tree->Project("spectre_amplitude", "amplitude_tree", Form("om_number == %d", om));
//     spectre->Draw();
//
//     TF1* f_fit_tot = new TF1 ("f_fit_tot" , "[0]*(7.11*TMath::Gaus(x[0],[1],[2]*sqrt([1])) + 1.84*TMath::Gaus(x[0],([1]+([1]/976.0)*72),[2]*sqrt([1]))) + 0.54*TMath::Gaus(x[0],([1]+([1]/976.0)*84),[2]*sqrt([1])) + [3]*(1.52*TMath::Gaus(x[0],[4],[2]*sqrt([4])) + 0.44*TMath::Gaus(x[0],([4]+([1]/976.0)*73),[2]*sqrt([4])) + 0.15*TMath::Gaus(x[0],([4]+([1]/976.0)*85),[2]*sqrt([4]))) + exp(-[5]*x[0]+[6]) + ((1 + exp(x*[7]))/(1+exp((x-[8])/[9]))) + ((1+exp(x*[10]))/(1+exp((x-[11])/[12])))", 0, 1200);
//     f_fit_tot->SetRange(60,200);
//     double param[13] = {312.8, 160.1, 1.675, 12798, 66.59, 0.01237, 9.224, -3, -1, -1.8,-2.83, 0.2126, -213};
//     f_fit_tot->SetParameters(param);
//
//     f_fit_tot->Draw("same");
//     // spectre->Fit(f_fit_tot, "RQ0");
//     // spectre->Fit(f_fit_tot, "RQ0");
//     // spectre->Fit(f_fit_tot, "RQ0");
//     // spectre->Fit(f_fit_tot, "RQ0");
//     // spectre->Fit(f_fit_tot, "RQ0");
//     // spectre->Fit(f_fit_tot, "RQ0");
//     // spectre->Fit(f_fit_tot, "RQ0");
//     // canvas->SetLogy();
//     canvas->SaveAs(Form("fit/fit_ref/ref_fit_om_%d.png", om));
//     // Result_tree.Fill();
//     // delete spectre;
//   }
//   newfile.cd();
//   Result_tree.Write();
//   newfile.Close();
// }
//
// void Ref_variation(std::vector<int> run_number, std::vector<int> time) {
//
//   int om_number;
//   double pic;
//   double charge;
//   double amplitude;
//   int ntime = time.size();
//   TFile fit_file(Form("Resultats_root/Amplitude_Li_run_%d.root", run_number[1]),"READ");
//   TTree* Li_tree = (TTree*)fit_file.Get("Result_tree");
//   Li_tree->SetBranchStatus("*",0);
//   Li_tree->SetBranchStatus("om_number",1);
//   Li_tree->SetBranchAddress("om_number", &om_number);
//   Li_tree->SetBranchStatus("run_number",1);
//   Li_tree->SetBranchAddress("run_number", &run_number);
//   Li_tree->SetBranchStatus("pic",1);
//   Li_tree->SetBranchAddress("pic", &pic);
//   Li_tree->SetBranchStatus("charge_tree",1);
//   Li_tree->SetBranchAddress("charge_tree", &charge);
//   Li_tree->SetBranchStatus("amplitude_tree",1);
//   Li_tree->SetBranchAddress("amplitude_tree", &amplitude);
//
//   TH3D* MC_Simu = new TH3D("LI_fit"), "LI_fit"),
//                           805, 0, 805,
//                           6, 1, 6,
//                           , charge_bin_min, charge_bin_max);
//
//   for (int i = 0; i < 805; i++) {
//     if (i == 712) {
//       i = 800
//     }
//     Result_tree.GetEntry(i);
//
//     double yaxis[ntime];
//     double yaxis_error[ntime];
//     for (int j = 0; j < 717; j++){
//       for (int i = 0; i < 6; i++) {
//         int nombre = 717*i+j;
//         Result_tree.GetEntry(nombre);
//         yaxis = mean_charge/comp[j];
//         // std::cout << "om = " << j << " and var = " << yaxis << '\n';
//         if ((yaxis < 0.9) || (yaxis > 1.1)) {
//
//           std::cout << " run == " << tab[i] << " om = " << j << " and var = " << yaxis << '\n';
//         }
//         // yaxis_error[i] = mean_error/comp*1.0 + (mean_charge/(comp*comp))*1.0*error1;
//       }
//     }
//
//
//   }
// }

void file_merger(std::vector<int> run, std::vector<int> run_ref, string addfile = "") {
  double Amplitude, Amplitude_error, Khi2, p1, p4, Amplitude_1MeV, Amplitude_05MeV, Ampl_norm, Ampl1_norm, Ampl05_norm, Ampl05_norm_error, Ampl1_norm_error, p1_er, p4_er;
  int i_om, run_number, pic, time, om_number;
  double norm[717][6];
  double norm1[717];
  double norm05[717];
  double norm1_er[717];
  double norm05_er[717];
  for (size_t i = 0; i < 717; i++) {
    for (size_t j = 0; j < run.size(); j++) {
      for (size_t k = 0; k < 6; k++) {
        norm[i][j] = 0;
      }
      norm1[i] = 0;
      norm05[i] = 0;
      norm1_er[i] = 0;
      norm05_er[i] = 0;
    }
  }

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("i_om", &i_om);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("Amplitude", &Amplitude);
  Result_tree.Branch("Amplitude_error", &Amplitude_error);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("time", &time);
  Result_tree.Branch("Amplitude_1MeV", &p1);
  Result_tree.Branch("Amplitude_05MeV", &p4);
  Result_tree.Branch("Ampl_norm", &Ampl_norm);
  Result_tree.Branch("Ampl05_norm", &Ampl05_norm);
  Result_tree.Branch("Ampl05_norm_error", &Ampl05_norm_error);
  Result_tree.Branch("Ampl1_norm", &Ampl1_norm);
  Result_tree.Branch("Ampl1_norm_error", &Ampl1_norm_error);

  Ampl05_norm = 0;
  Ampl1_norm = 0;
  Amplitude_1MeV = 0;
  Amplitude_05MeV = 0;
  Ampl05_norm = 0;
  Ampl1_norm = 0;
  for (int j = 0; j < run.size(); j++) {
    TFile *file1 = new TFile(Form("Resultats_root/Amplitude_Li_run_%d.root", run[j]), "READ");
    TTree* tree = (TTree*)file1->Get("Result_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("i_om",1);
    tree->SetBranchAddress("i_om", &i_om);
    tree->SetBranchStatus("pic",1);
    tree->SetBranchAddress("pic", &pic);
    tree->SetBranchStatus("Khi2",1);
    tree->SetBranchAddress("Khi2", &Khi2);
    tree->SetBranchStatus("Amplitude",1);
    tree->SetBranchAddress("Amplitude", &Amplitude);
    tree->SetBranchStatus("Amplitude_error",1);
    tree->SetBranchAddress("Amplitude_error", &Amplitude_error);
    tree->SetBranchStatus("run_number",1);
    tree->SetBranchAddress("run_number", &run_number);
    tree->SetBranchStatus("time",1);
    tree->SetBranchAddress("time", &time);
    for (size_t i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      if (j == 0) {
        norm[i_om][pic] = Amplitude;
      }
      if (norm[i_om][pic] != 0) {
        Ampl_norm = Amplitude/norm[i_om][pic];
      }

      Result_tree.Fill();
    }
    file1->Close();
  }
  pic = 0;
  Khi2 = 0;
  Amplitude = 0;
  Amplitude_error = 0;
  Ampl_norm = 0;

  for (int j = 0; j < run_ref.size(); j++) {
    TFile *file1 = new TFile(Form("Resultats_root/fit_ref_run_%d.root", run_ref[j]), "READ");
    TTree* ref_tree = (TTree*)file1->Get("Result_tree");
    ref_tree->SetBranchStatus("*",0);
    ref_tree->SetBranchStatus("om_number",1);
    ref_tree->SetBranchAddress("om_number", &i_om);
    ref_tree->SetBranchStatus("p1",1);
    ref_tree->SetBranchAddress("p1", &p1);
    ref_tree->SetBranchStatus("p1_er",1);
    ref_tree->SetBranchAddress("p1_er", &p1_er);
    ref_tree->SetBranchStatus("p4",1);
    ref_tree->SetBranchAddress("p4", &p4);
    ref_tree->SetBranchStatus("p4_er",1);
    ref_tree->SetBranchAddress("p4_er", &p4_er);
    ref_tree->SetBranchStatus("run_number",1);
    ref_tree->SetBranchAddress("run_number", &run_number);
    ref_tree->SetBranchStatus("time",1);
    ref_tree->SetBranchAddress("time", &time);
    for (size_t i = 0; i < ref_tree->GetEntries(); i++) {
      ref_tree->GetEntry(i);
      if (j == 0) {
        norm1[i_om] = p1;
        norm1_er[i_om] = p1;
        norm05[i_om] = p4;
        norm05_er[i_om] = p1;
      }
      if (norm05[i_om] != 0) {
        Ampl05_norm = p4/norm05[i_om]*1.;
        Ampl05_norm_error = abs(p4_er/norm05[i_om]*1.) + abs(norm05_er[i_om]*p4/(norm05[i_om]*norm1[i_om]));
      }
      if (norm1[i_om] != 0) {
        Ampl1_norm = p1/norm1[i_om]*1.;
        Ampl1_norm_error = abs(p1_er/norm1[i_om]*1.) + abs(norm1_er[i_om]*p1/(norm1[i_om]*norm1[i_om])*1.);
      }
      Result_tree.Fill();
    }
    file1->Close();
  }

  if (addfile.compare("") != 0) {
    TFile *file1 = new TFile(Form("Resultats_root/%s.root", addfile.c_str()), "READ");
    TTree* tree = (TTree*)file1->Get("Result_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("i_om",1);
    tree->SetBranchAddress("i_om", &i_om);
    tree->SetBranchStatus("pic",1);
    tree->SetBranchAddress("pic", &pic);
    tree->SetBranchStatus("Khi2",1);
    tree->SetBranchAddress("Khi2", &Khi2);
    tree->SetBranchStatus("Amplitude",1);
    tree->SetBranchAddress("Amplitude", &Amplitude);
    tree->SetBranchStatus("Amplitude_error",1);
    tree->SetBranchAddress("Amplitude_error", &Amplitude_error);
    tree->SetBranchStatus("run_number",1);
    tree->SetBranchAddress("run_number", &run_number);
    tree->SetBranchStatus("time",1);
    tree->SetBranchAddress("time", &time);
    tree->SetBranchStatus("Amplitude_1MeV",1);
    tree->SetBranchAddress("Amplitude_1MeV", &Amplitude_1MeV);
    tree->SetBranchStatus("Ampl1_norm",1);
    tree->SetBranchAddress("Ampl1_norm", &Ampl1_norm);
    tree->SetBranchStatus("Ampl1_norm_error",1);
    tree->SetBranchAddress("Ampl1_norm_error", &Ampl1_norm_error);
    tree->SetBranchStatus("Amplitude_05MeV",1);
    tree->SetBranchAddress("Amplitude_05MeV", &Amplitude_05MeV);
    tree->SetBranchStatus("Ampl05_norm",1);
    tree->SetBranchAddress("Ampl05_norm", &Ampl05_norm);
    tree->SetBranchStatus("Ampl05_norm_error",1);
    tree->SetBranchAddress("Ampl05_norm_error", &Ampl05_norm_error);

    for (size_t i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      Result_tree.Fill();
    }
    file1->Close();
  }

  TFile *newfile = new TFile("Resultats_root/fit_total.root", "RECREATE");
  newfile->cd();
  Result_tree.Write();
  newfile->Close();

}

// void TGrapher() {
//       double comp[520];
//
//     for (int i = 0; i < 520; i++) {
//       Result_tree.GetEntry(i);
//       comp[i] = mean_charge;
//       // error1[i] = mean_error;
//     }
//     file.cd();
//
//     double yaxis[4];
//     double yaxis_error[4];
//     double xaxis[4] = {0, 19, 38, 57};
//     double xaxiserror[4] = {0.5, 0.5, 0.5, 0.5};
//     for (int j = 0; j < 520; j++){
//       for (int i = 0; i < 4; i++) {
//         int nombre = 520*i+j;
//         // int nombre = i;
//         Result_tree.GetEntry(nombre);
//         yaxis[i] = mean_charge/comp[j];
//         // std::cout << "om = " << j << " and var = " << yaxis << '\n';
//         if ((yaxis[i] < 1.1) && (yaxis[i] > 1.05)) {
//
//           std::cout << " run == " << tab[i] << " om = " << j << " and var = " << yaxis[i] << '\n';
//         }
//         // yaxis_error[i] = mean_error/comp*1.0 + (mean_charge/(comp*comp))*1.0*error1;
//       }
//       TGraphErrors comp_map (4, xaxis, yaxis, xaxiserror, yaxis_error);
//       comp_map.SetName(Form("fit_Tl_om_%d", j));
//       comp_map.SetNameTitle(Form("fit_Tl_om_%d", j), Form("evolution du gain de l'OM %d", j));
//       comp_map.GetXaxis()->SetTitle("Temps (h)");
//       comp_map.GetYaxis()->SetTitle("Gain(t)/Gain(0)");
//       comp_map.SetMarkerColor(2);
//       comp_map.SetMarkerStyle(34);
//       comp_map.SetMarkerSize(2);
//
//
//       // TCanvas* canvas2 = new TCanvas;
//       // comp_map.Draw();
//       // canvas2->SaveAs(Form("fit/fit_Tl/variation/charge_fit_om_%03d.png", j));
//       comp_map.Write();
//
//     }
//
//
//
//
//     mean_charge_map.Write();
//     Result_tree.Write();
//
//     file.Close();
//     outFile.close();
//
//     return;
//
// }
//

int main(int argc, char const *argv[]){
  int n_run, run, t;
  std::vector<int> run_number, time, ref_run_number, ref_time;
  int compteur = 0;
  string file;
  bool add = false;

  for(int i = 0; i<argc; i++){
    if (std::string(argv[i]) == "-add" ) {
      file = argv[i+1];
      std::cout << file << '\n';
      add = true;
    }
  }
  std::cout << "How many run do you want ?" << '\n';
  std::cin >> n_run;
  std::cout << "Write the runs you want" << '\n';
  while (compteur < n_run && cin >> run) {
    run_number.push_back(run);
    std::cout << "Write the delay between the runs (0 for the first)" << '\n';
    while (compteur < n_run && cin >> t) {
      time.push_back(t);
      break;
    }

    compteur++;
    if (compteur < n_run) {
      std::cout << "Write the runs you want" << '\n';
    }
  }
  compteur = 0;
  n_run = 0;
  // std::cout << "How many Ref run do you want ?" << '\n';
  // std::cin >> n_run;
  // std::cout << "Write the Ref OM only runs you want" << '\n';
  // while (compteur < n_run && cin >> run) {
  //   ref_run_number.push_back(run);
  //   std::cout << "Write the delay between the runs (0 for the first)" << '\n';
  //   while (compteur < n_run && cin >> t) {
  //     ref_time.push_back(t);
  //     break;
  //   }
  //   compteur++;
  //   if (compteur < n_run) {
  //     std::cout << "Write the Ref runs you want" << '\n';
  //   }
  // }
  std::cout << "Code start running" << '\n';

  for (int i = 0; i < run_number.size(); i++) {
    fit_LI_amplitude(run_number[i], time[i]);
  }
  // for (int i = 0; i < ref_run_number.size(); i++) {
  //   fit_ref(ref_run_number[i], ref_time[i]);
  // }
  if (add == false) {
    file_merger(run_number, ref_run_number);
  }
  else{
    file_merger(run_number, ref_run_number, file);
  }


  // amplitude_variation(run_number, time);

  return 0;
}
//
// void Ref_variation_tgraph(std::vector<int> run_number, std::vector<int> time) {
//
//   int om_number;
//   double pic;
//   double charge;
//   double amplitude;
//   int ntime = time.size();
//   TFile fit_file(Form("Resultats_root/fit_total.root", run_number[1]),"READ");
//   TTree* Li_tree = (TTree*)fit_file.Get("Result_tree");
//   Li_tree->SetBranchStatus("*",0);
//   Li_tree->SetBranchStatus("i_om",1);
//   Li_tree->SetBranchAddress("i_om", &i_om);
//   Li_tree->SetBranchStatus("run_number",1);
//   Li_tree->SetBranchAddress("run_number", &run_number);
//   Li_tree->SetBranchStatus("pic",1);
//   Li_tree->SetBranchAddress("pic", &pic);
//   Li_tree->SetBranchStatus("Khi2",1);
//   Li_tree->SetBranchAddress("Khi2", &Khi2);
//   Li_tree->SetBranchStatus("charge_tree",1);
//   Li_tree->SetBranchAddress("charge_tree", &charge);
//   Li_tree->SetBranchStatus("Amplitude",1);
//   Li_tree->SetBranchAddress("Amplitude", &Amplitude);
//   Li_tree->SetBranchStatus("Amplitude_error",1);
//   Li_tree->SetBranchAddress("Amplitude_error", &Amplitude_error);
//   Li_tree->SetBranchStatus("Amplitude_1MeV",1);
//   Li_tree->SetBranchAddress("Amplitude_1MeV", &Amplitude_1MeV);
//   Li_tree->SetBranchStatus("Amplitude_05MeV",1);
//   Li_tree->SetBranchAddress("Amplitude_05MeV", &Amplitude_05MeV);
//   Li_tree->SetBranchStatus("Amplitude_1MeV_error",1);
//   Li_tree->SetBranchAddress("Amplitude_1MeV_error", &Amplitude_1MeV_error);
//   Li_tree->SetBranchStatus("Amplitude_05MeV_error",1);
//   Li_tree->SetBranchAddress("Amplitude_05MeV_error", &Amplitude_05MeV_error);
//
//   double comp[717];
//
//   for (int i = 0; i < 805; i++) {
//     if (i == 712) {
//       i = 800;
//     }
//     Result_tree.GetEntry(i);
//     comp[i] = charge;
//     // error1[i] = mean_error;
//
//     for (int i = 0; i < time.size(); i++) {
//       double xaxis[i] = temps[i];
//       double xaxiserror[i] = temps[i]*0.01;
//     }
//     double yaxis[ntime];
//     double yaxis_error[ntime];
//     for (int j = 0; j < 717; j++){
//       for (int i = 0; i < 6; i++) {
//         int nombre = 717*i+j;
//         Result_tree.GetEntry(nombre);
//         yaxis = mean_charge/comp[j];
//         // std::cout << "om = " << j << " and var = " << yaxis << '\n';
//         if ((yaxis < 0.9) || (yaxis > 1.1)) {
//
//           std::cout << " run == " << tab[i] << " om = " << j << " and var = " << yaxis << '\n';
//         }
//         // yaxis_error[i] = mean_error/comp*1.0 + (mean_charge/(comp*comp))*1.0*error1;
//       }
//     }
//     TGraphErrors comp_map (4, xaxis, yaxis, xaxiserror, yaxis_error);
//
//     comp_map.SetNameTitle("fit_Tl", "Li gain evolution");
//     comp_map.GetXaxis()->SetTitle("Temps (h)");
//     comp_map.GetYaxis()->SetTitle("Gain(t)/Gain(0)");
//     comp_map.SetMarkerColor(2);
//     comp_map.SetMarkerStyle(34);
//     comp_map.SetMarkerSize(2);
//
//   }
// }
