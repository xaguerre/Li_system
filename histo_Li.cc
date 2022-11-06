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
#include <TParameter.h>
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

string namer(int om){
  string name;
  if (om == 712) {name = "Ref MW1";}
  if (om == 713) {name = "Ref MW2";}
  if (om == 714) {name = "Ref GV";}
  if (om == 715) {name = "Ref XW1";}
  if (om == 716) {name = "Ref XW2";}
  return name;
}

void txt_to_root() {

  std::ifstream parametres("/home/aguerre/Bureau/Thèse/Li_system/Resultats_txt/Resultats_charge.txt");
  TFile file("test.root","RECREATE");

  int compteur = 0;
  int om_number;
  double n_evt;
  double mean_charge;
  double sigma;
  double nbg;
  double C = 0;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("n_evt", &n_evt);
  Result_tree.Branch("mean_charge", &mean_charge);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("nbg", &nbg);
  Result_tree.Branch("C", &C);

  while(parametres >> om_number)
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

  TFile file("root/getMean.root","RECREATE");

  int i_om =0;
  double mean =0;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &i_om);
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

std::vector<double> time_measurer(int run_number){
  TFile tree_file(Form("histo_brut/Li_system_%d.root", run_number), "READ");
  int om_number;
  double time;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  gROOT->cd();
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);

  TH1D *spectre = new TH1D ("time_spectre", "", 1000, 0, 350);
  tree->Project("time_spectre", "time");
  spectre->Draw();

  int plus =0;
  int moins = 0;
  std::vector<double> time_measurer;
  // return time_measurer;
  int i;
  for (i = 0; i < 1000; i++) {
    if (spectre->GetBinContent(i) > 0) {
      plus = 1;
    }
    if (spectre->GetBinContent(i) == 0) {
      moins = 1;
    }
    if (plus - moins == 0) {
      time_measurer.push_back(spectre->GetBinCenter(i));
      plus = 0;
      moins = 0;
    }
  }
  std::vector<double> interval;
  interval.push_back(time_measurer.at(0));
  for (int j = 1; j < time_measurer.size()-2; j+=2) {
    interval.push_back((time_measurer.at(j)+time_measurer.at(j+1))/2);
  }
  delete spectre;
  return interval;
}

double* Ref_corrector(int run, string correction, double *gain_tab) {

  TFile file(Form("calcul_gain_ref/root/Merged_Fit/Fit_Ref_%s.root", correction.c_str()), "READ");

  // std::cout << Form("calcul_gain_ref/root/Merged_Fit/Fit_Ref_%s.root", correction.c_str()) << '\n';

  int run_number;
  double gain;
  TTree* tree = (TTree*)file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &gain);
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number);

  int compteur = 0;
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (run_number == run+1) {
      gain_tab[compteur] = gain;
      compteur++;
    }
  }
  file.Close();
  return gain_tab;
}

void fit_LI_amplitude(int run_number, double *ref_gain_table){
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TFile file(Form("root/Fit_Ampl/Amplitude_Li_run_%d.root", run_number),"RECREATE");

  std::vector<double> interval;
  interval = time_measurer(run_number);

  int om_number;
  double constante;
  double mean, time, start_time;
  double mean_error, nevent;
  double sigma;
  int pic =0;
  double Khi2 = 0;
  int intensity = 0;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("constante", &constante);
  Result_tree.Branch("Amplitude", &mean);
  Result_tree.Branch("Amplitude_error", &mean_error);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("intensity", &intensity);
  Result_tree.Branch("time", &start_time);
  Result_tree.Branch("nevent", &nevent);


  int debut = 0;
  TCanvas* canvas = new TCanvas;
  TH1F comp_map("comp_map", "comp_map", 100, 0, 100);
  TFile *tree_file = new TFile (Form("histo_brut/Li_system_%d.root", run_number), "READ");

  double charge_tree;
  double amplitude_tree;
  TTree* tree = (TTree*)tree_file->Get("Result_tree");
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

  TParameter<double> *param = new TParameter<double>("start_time", time);
  param = (TParameter<double>*)(tree_file->Get("start_time"));
  start_time = param->GetVal();
  std::cout << "time = " << start_time << '\n';
  gROOT->cd();

  int number = 800;
  if (run_number > 836) {
    number = 712;
  }

  for(int om = number; om < number + 5; om+=1)
  {
    if (om == 712) {
      om = number;
    }
    if ((om > 259 && om < 520) || (om > 583 && om < 647) || (om > 679 && om < 712) ) {debut = debut + (interval.size()/2);}
    for (double j = debut; j < debut + (interval.size()/2); j++)
    {

      TH1D *spectre = new TH1D ("spectre_amplitude", "", 700, 0, 2300 );
      tree->Project("spectre_amplitude", Form("amplitude_tree*%f", ref_gain_table[om - number]), Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, interval.at(j), interval.at(j+1)));

      if (om >711) {
        tree->Project("spectre_amplitude", Form("amplitude_tree*%f", ref_gain_table[om - number]), Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, interval.at(j), interval.at(j+1)));
        std::cout << "ref_gain = " << ref_gain_table[om-712] << '\n';
      }
      else tree->Project("spectre_amplitude", "amplitude_tree", Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, interval.at(j), interval.at(j+1)));

      std::cout << "temps = " << interval.at(j) << " - " << interval.at(j+1) << '\n';
      if (spectre->GetEntries() < 300) {
        std::cout << "" << '\n';
        std::cout << "trop peu d'entries pour l'OM " << om << '\n';
        std::cout << "" << '\n';
        pic = j;
        Khi2 = 0;
        om_number = om;
        mean = 0;
        mean_error = 0;
        if (om > 711 && run_number < 837) {
          om_number = om-88;
        }
        Result_tree.Fill();
        delete spectre;
      }
      else if (spectre->GetMean() > 1900) {
        std::cout << "" << '\n';
        std::cout << "the amplitude sature" << '\n';
        std::cout << "" << '\n';
        pic = j;
        Khi2 = 0;
        om_number = om;
        mean = 0;
        mean_error = 0;
        if (om > 711 && run_number < 837) {
          om_number = om-88;
        }
        Result_tree.Fill();
        delete spectre;
      }
      else if (spectre->GetMean() < 20){
        std::cout << "" << '\n';
        std::cout << "too few charge" << '\n';
        std::cout << "" << '\n';
        pic = j;
        Khi2 = 0;
        om_number = om;
        mean = 0;
        mean_error = 0;
        if (om > 711 && run_number < 837) {
          om_number = om-88;
        }
        Result_tree.Fill();
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

        pic = j;
        Khi2 = f_Gaus->GetChisquare()/f_Gaus->GetNDF();
        om_number = om;
        constante = (f_Gaus->GetParameter(0));
        mean = (f_Gaus->GetParameter(1));
        sigma = (f_Gaus->GetParameter(2));
        mean_error = f_Gaus->GetParError(1) ;
        intensity = intensity_chooser(pic);
        spectre->Draw();
        f_Gaus->Draw("same");

        if (om < 712) {
          canvas->SaveAs(Form("fit/fit_Li/amplitude_fit/om_%d/OM_%03d_pic_%d_run_%d.png", om, pic, om_number, run_number));
        }
        if (om > 711 && run_number < 837) {
          canvas->SaveAs(Form("fit/fit_Li/amplitude_fit_ref/om_%d/OM_%03d_pic_%d_run_%d.png", om-88, om -88, pic, run_number));
          om_number = om-88;
        }
        if (om > 711 && run_number > 836) {
          canvas->SaveAs(Form("fit/fit_Li/amplitude_fit_ref/om_%d/OM_%03d_pic_%d_run_%d.png", om, om, pic, run_number));
        }
        Result_tree.Fill();

        delete spectre;
        delete f_Gaus;
      }
    }
  }
  std::cout << "time = " << start_time << '\n';

  file.cd();
  Result_tree.Write();
  file.Close();
  return;
}

void fit_LI_amplitude_Ref(int run_number, double *ref_gain_table){
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TFile file(Form("root/Fit_Ampl/Amplitude_Li_run_%d.root", run_number),"RECREATE");

  std::vector<double> interval;
  interval = time_measurer(run_number);

  int om_number;
  double constante;
  double mean, time, start_time;
  double mean_error, nevent;
  double sigma;
  int pic =0;
  double Khi2 = 0;
  int intensity = 0;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("constante", &constante);
  Result_tree.Branch("Amplitude", &mean);
  Result_tree.Branch("Amplitude_error", &mean_error);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("intensity", &intensity);
  Result_tree.Branch("time", &start_time);
  Result_tree.Branch("nevent", &nevent);


  int debut = 0;
  TCanvas* canvas = new TCanvas;
  TH1F comp_map("comp_map", "comp_map", 100, 0, 100);
  TFile *tree_file = new TFile (Form("histo_brut/Li_system_%d.root", run_number), "READ");

  double charge_tree;
  double amplitude_tree;
  TTree* tree = (TTree*)tree_file->Get("Result_tree");
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

  TParameter<double> *param = new TParameter<double>("start_time", time);
  param = (TParameter<double>*)(tree_file->Get("start_time"));
  start_time = param->GetVal();
  std::cout << "time = " << start_time << '\n';
  gROOT->cd();

  int number = 800;
  if (run_number > 836) {
    number = 712;
  }

  for(int om = number; om < number + 5; om+=1)
  {
    if (om == 712) {
      om = number;
    }
    if ((om > 259 && om < 520) || (om > 583 && om < 647) || (om > 679 && om < 712) ) {debut = debut + (interval.size()/2);}
    for (double j = debut; j < debut + (interval.size()/2); j++)
    {

      TH1D *spectre = new TH1D ("spectre_amplitude", "", 700, 0, 2300 );
      tree->Project("spectre_amplitude", Form("amplitude_tree*%f", ref_gain_table[om - number]), Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, interval.at(j), interval.at(j+1)));

      if (om >711) {
        tree->Project("spectre_amplitude", Form("amplitude_tree*%f", ref_gain_table[om - number]), Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, interval.at(j), interval.at(j+1)));
        std::cout << "ref_gain = " << ref_gain_table[om-712] << '\n';
      }
      else tree->Project("spectre_amplitude", "amplitude_tree", Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, interval.at(j), interval.at(j+1)));

      std::cout << "temps = " << interval.at(j) << " - " << interval.at(j+1) << '\n';
      if (spectre->GetEntries() < 300) {
        std::cout << "" << '\n';
        std::cout << "trop peu d'entries pour l'OM " << om << '\n';
        std::cout << "" << '\n';
        pic = j;
        Khi2 = 0;
        om_number = om;
        mean = 0;
        mean_error = 0;
        if (om > 711 && run_number < 837) {
          om_number = om-88;
        }
        Result_tree.Fill();
        delete spectre;
      }
      else if (spectre->GetMean() > 1900) {
        std::cout << "" << '\n';
        std::cout << "the amplitude sature" << '\n';
        std::cout << "" << '\n';
        pic = j;
        Khi2 = 0;
        om_number = om;
        mean = 0;
        mean_error = 0;
        if (om > 711 && run_number < 837) {
          om_number = om-88;
        }
        Result_tree.Fill();
        delete spectre;
      }
      else if (spectre->GetMean() < 20){
        std::cout << "" << '\n';
        std::cout << "too few charge" << '\n';
        std::cout << "" << '\n';
        pic = j;
        Khi2 = 0;
        om_number = om;
        mean = 0;
        mean_error = 0;
        if (om > 711 && run_number < 837) {
          om_number = om-88;
        }
        Result_tree.Fill();
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

        pic = j;
        Khi2 = f_Gaus->GetChisquare()/f_Gaus->GetNDF();
        om_number = om;
        constante = (f_Gaus->GetParameter(0));
        mean = (f_Gaus->GetParameter(1));
        sigma = (f_Gaus->GetParameter(2));
        mean_error = f_Gaus->GetParError(1) ;
        intensity = intensity_chooser(pic);
        spectre->Draw();
        f_Gaus->Draw("same");

        if (om < 712) {
          canvas->SaveAs(Form("fit/fit_Li/amplitude_fit/om_%d/OM_%03d_pic_%d_run_%d.png", om, pic, om_number, run_number));
        }
        if (om > 711 && run_number < 837) {
          canvas->SaveAs(Form("fit/fit_Li/amplitude_fit_ref/om_%d/OM_%03d_pic_%d_run_%d.png", om-88, om -88, pic, run_number));
          om_number = om-88;
        }
        if (om > 711 && run_number > 836) {
          canvas->SaveAs(Form("fit/fit_Li/amplitude_fit_ref/om_%d/OM_%03d_pic_%d_run_%d.png", om, om, pic, run_number));
        }
        Result_tree.Fill();

        delete spectre;
        delete f_Gaus;
      }
    }
  }
  std::cout << "time = " << start_time << '\n';

  file.cd();
  Result_tree.Write();
  file.Close();
  return;
}


int bundle_number(int om_number){
  int bundle_number = 0;
  if ((om_number%13 > 7 && om_number/13 < 10 && om_number < 260) || (om_number > 663 && om_number < 672) || (om_number < 552 && om_number > 545) || (om_number < 536 && om_number > 529)) {
    bundle_number = 1;
  }
  else if ((om_number%13 < 8 && om_number/13 < 6 && om_number < 260) ||(om_number > 647 && om_number < 652) ||(om_number < 546 && om_number > 535) ||(om_number < 530 && om_number > 519)){
    bundle_number = 2;
  }
  if ((om_number%13 > 7 && om_number/13 > 9 && om_number < 260) || (om_number > 671 && om_number < 680) || (om_number < 568 && om_number > 561) || (om_number < 584 && om_number > 577)) {
    bundle_number = 3;
  }
  else if ((om_number%13 < 8 && om_number/13 > 5 && om_number/13 < 14 && om_number < 260) || (om_number < 660 && om_number > 651)) {
    bundle_number = 4;
  }
  else if ((om_number%13 < 8 && om_number/13 > 13 && om_number < 260) || (om_number < 664 && om_number > 659) || (om_number < 562 && om_number > 551) ||(om_number < 578 && om_number > 567)) {
    bundle_number = 5;
  }
  else if ((om_number%13 > 7 && (om_number/13-20) < 10 && om_number < 520 && om_number > 259) ||(om_number > 695 && om_number < 704) ||(om_number < 600 && om_number > 593) ||(om_number < 616 && om_number > 609)){
    bundle_number = 6;
  }
  else if ((om_number%13 > 7 && (om_number/13-20) > 9 && om_number < 520 && om_number > 259) ||(om_number > 703 && om_number < 712) ||(om_number < 648 && om_number > 641) ||(om_number < 632 && om_number > 625)){
    bundle_number = 7;
  }
  else if ((om_number%13 < 8 && (om_number/13-20) < 6 && om_number < 520 && om_number > 259) ||(om_number > 679 && om_number < 684) ||(om_number < 594 && om_number > 583) ||(om_number < 610 && om_number > 599)){
    bundle_number = 8;
  }
  else if ((om_number%13 < 8 && (om_number/13-20) > 5 && (om_number/13-20) < 14 && om_number < 520 && om_number > 259) ||(om_number > 683 && om_number < 692)){
    bundle_number = 9;
  }
  else if ((om_number%13 < 8 && (om_number/13-20) > 13 && om_number < 520 && om_number > 259) ||(om_number > 691 && om_number < 696) ||(om_number < 642 && om_number > 631) ||(om_number < 626 && om_number > 615)){
    bundle_number = 10;
  }
  return bundle_number;
}

void file_merger(std::vector<int> run, string addfile = "") {
  double Amplitude, Amplitude_error, Khi2, time;
  int om_number, run_number, pic;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("Amplitude", &Amplitude);
  Result_tree.Branch("Amplitude_error", &Amplitude_error);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("time", &time);

  for (int j = 0; j < run.size(); j++) {
    TFile *file1 = new TFile(Form("root/Fit_Ampl/Amplitude_Li_run_%d.root", run[j]), "READ");
    TTree* tree = (TTree*)file1->Get("Result_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("om_number",1);
    tree->SetBranchAddress("om_number", &om_number);
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
      Result_tree.Fill();
    }

    file1->Close();
  }
  pic = 0;
  Khi2 = 0;
  Amplitude = 0;
  Amplitude_error = 0;

  if (addfile.compare("") != 0) {
    TFile *file1 = new TFile(Form("root/%s.root", addfile.c_str()), "READ");
    TTree* tree = (TTree*)file1->Get("Result_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("om_number",1);
    tree->SetBranchAddress("om_number", &om_number);
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
      Result_tree.Fill();
    }
    file1->Close();
  }

  TFile *newfile = new TFile("root/Merged/fit_total.root", "RECREATE");
  newfile->cd();
  Result_tree.Write();
  newfile->Close();

}

void TGrapher(std::string file_name, int n_run) {

  int om_number, pic;
  double Amplitude, Amplitude_error, time;
  TFile file(Form("root/TGraph/TGraph_%s.root", file_name.c_str()), "RECREATE");

  TFile *tree_file = new TFile(Form("root/Merged/%s.root", file_name.c_str()), "READ");
  TTree* tree = (TTree*)tree_file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("pic",1);
  tree->SetBranchAddress("pic", &pic);
  tree->SetBranchStatus("Amplitude",1);
  tree->SetBranchAddress("Amplitude", &Amplitude);
  tree->SetBranchStatus("Amplitude_error",1);
  tree->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);

  std::cout << "nrun = " << n_run << '\n';

  TGraphErrors *gain_graph[5][6];
  double yaxis[n_run];
  double yaxis_error[n_run];
  double xaxis[n_run];
  double xaxis_error[n_run];

  int jmin = 712;
  int jmax =717;
  double Ampl_ref, Ampl_ref_error;
  file.cd();
  for (int j = jmin; j < jmax; j++){
    auto canvas = new TCanvas(Form("Allfit_om_%d", j),"",1600,800);
    canvas->Divide(3,2);
    for (int i = 0; i < 6; i++) {
      for (int k = 0; k < n_run; k++) {
        tree->GetEntry((jmax-jmin)*6*k + (j-jmin)*6 + i);// + (j-jmin)*i + i);
        if (k == 0) {
          Ampl_ref = Amplitude;
          Ampl_ref_error = Amplitude_error;
        }
        std::cout << "entry = " << (jmax-jmin)*6*k + i + (j-jmin)*6 << '\n';
        std::cout << "k = " << k << "; i = " << i << "; j = " << j << " and Amplitude = " <<  Amplitude << '\n';
        yaxis[k] = Amplitude/Ampl_ref;
        yaxis_error[k] = Amplitude_error/Ampl_ref + (Ampl_ref_error*Amplitude)/(Ampl_ref*Ampl_ref);
        if (k == 0) {
          yaxis_error[k] = Amplitude_error/Ampl_ref;
        }
        xaxis[k] = time;
        xaxis_error[k] = 0.01;
      }
      std::cout << "" << '\n';
      gain_graph[jmin-jmax][i] = new TGraphErrors(n_run, xaxis, yaxis, xaxis_error, yaxis_error);
      string name = to_string (j);
      if (j > 711) {
        name = namer(j);
      }
      gain_graph[jmin-jmax][i]->SetName(Form("fit_Tl_om_%s_pic_%d", name.c_str(), i));
      gain_graph[jmin-jmax][i]->SetNameTitle(Form("fit_Tl_om_%s_pic_%d", name.c_str(), i), Form("evolution du gain de l'OM %s, pic %d", name.c_str(), i));
      gain_graph[jmin-jmax][i]->GetXaxis()->SetTitle("Temps (h)");
      gain_graph[jmin-jmax][i]->GetYaxis()->SetTitle("Gain(t)/Gain(0)");
      gain_graph[jmin-jmax][i]->GetXaxis()->SetTimeDisplay(1);
      gain_graph[jmin-jmax][i]->SetMarkerColor(2);
      gain_graph[jmin-jmax][i]->SetMarkerStyle(3);
      gain_graph[jmin-jmax][i]->SetMarkerSize(2);
      gain_graph[jmin-jmax][i]->GetYaxis()->SetRangeUser(0.9, 1.1);

      canvas->cd(i+1);
      gain_graph[jmin-jmax][i]->Draw("AP");

      gain_graph[jmin-jmax][i]->Write();
      canvas->Update();

      gain_graph[jmin-jmax][i]->Write();

    }
    canvas->Update();
    canvas->Write();
    delete canvas;
  }



  file.Close();
  return;
}

void test() {
  // fit_LI_amplitude(788);
  // fit_LI_amplitude(790);
  // fit_LI_amplitude(795);
  // fit_LI_amplitude(797);
  // fit_LI_amplitude(799);
  // fit_LI_amplitude(801);
  // fit_LI_amplitude(803);
  // fit_LI_amplitude(805);
  // fit_LI_amplitude(808);
  // fit_LI_amplitude(810);
  // fit_LI_amplitude(816);
  // fit_LI_amplitude(818);
  // fit_LI_amplitude(820);
  // fit_LI_amplitude(822);
  // fit_LI_amplitude(833);
  // fit_LI_amplitude(835);
  // std::cout << "fitted" << '\n';
  std::vector<int> a;
  a.push_back(788);
  a.push_back(790);
  a.push_back(795);
  a.push_back(797);
  a.push_back(799);
  a.push_back(801);
  a.push_back(803);
  a.push_back(805);
  a.push_back(808);
  a.push_back(810);
  a.push_back(816);
  a.push_back(818);
  a.push_back(820);
  a.push_back(822);
  a.push_back(833);
  a.push_back(835);
  file_merger(a);
  std::cout << "merged" << '\n';
  TGrapher("fit_total", 16);
}

void pic_comparator(int n_run) {
  int om_number, pic, run_number;
  double Amplitude, Amplitude_error, time;
  TFile file("root/evolution_pic.root", "RECREATE");

  TFile *tree_file = new TFile("root/Merged/fit_total.root", "READ");
  TTree* tree = (TTree*)tree_file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("pic",1);
  tree->SetBranchAddress("pic", &pic);
  tree->SetBranchStatus("Amplitude",1);
  tree->SetBranchAddress("Amplitude", &Amplitude);
  tree->SetBranchStatus("Amplitude_error",1);
  tree->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number);
  gROOT->cd();

  double yaxis[n_run];
  double yaxis_error[n_run];
  double xaxis_error[n_run];
  double xaxis[n_run];
  TGraphErrors *variation_pic[5][6];
  double add[5][6];

  string name;
  double ref[5];
  double ref_err[5];
  int compteur = 0;
  for (int om = 712; om < 717; om++) {
    name = namer(om);


    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      if (pic == 3 && om_number == om) {
        if (run_number == 788) {
          ref[om-712] = Amplitude;
          ref_err[om-712] = Amplitude_error;
        }
        yaxis[compteur] = Amplitude;
        yaxis_error[compteur] = Amplitude_error;
        xaxis[compteur] = time;
        xaxis_error[compteur] = 0.01;
        compteur++;
        // std::cout << "om = " << om <<  " entry = " << compteur -1 << " and ampl = " << Amplitude << " +- " << Amplitude_error << " and time = " << time << '\n';
      }
    }
    variation_pic[om-712][2] = new TGraphErrors(n_run, xaxis, yaxis, xaxis_error, yaxis_error);
    variation_pic[om-712][2]->SetName(Form("Gain_evolution_of_the_OM_%s_of_pic_3", name.c_str()));
    variation_pic[om-712][2]->SetNameTitle(Form("Gain_evolution_of_the_OM_%s_of_pic_3", name.c_str()), Form("Gain_evolution_of_the_OM_%s_of_pic_3", name.c_str()));
    variation_pic[om-712][2]->GetXaxis()->SetTitle("Temps (h)");
    variation_pic[om-712][2]->GetYaxis()->SetTitle("Amplitude (ua)");
    // variation_pic[om-712][2]->GetYaxis()->SetRange(0,2000);
    // variation_pic[om-712][2]->GetXaxis()->SetRange(0,2);
    variation_pic[om-712][2]->GetXaxis()->SetTimeDisplay(1);
    variation_pic[om-712][2]->SetMarkerColor(2);
    variation_pic[om-712][2]->SetMarkerStyle(3);
    variation_pic[om-712][2]->SetMarkerSize(2);
    // variation_pic[om-712][2]->Draw();
    compteur = 0;
  }

  for (int om = 712; om < 717; om++) {
    for (int pique = 0; pique < 6; pique++) {
      // if (pique == 2){pique = 3;}
      name = namer(om);

      for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if (pic !=3 && om_number == om && pic == pique) {
        std::cout << "pic = " << pic << " om = " << om << " i = " << i << " amplitude = " << Amplitude << '\n';
          if (Amplitude == 0) {
            break;
          }
          if (run_number == 788) {
            yaxis[compteur] = 1;
            add[om-712][pique] = 1 - abs(Amplitude - ref[om-712])/ref[om-712];
            yaxis_error[compteur] = 0;
            xaxis[compteur] = time;
          }
          else{
            yaxis[compteur] = abs(Amplitude - ref[om-712])/ref[om-712] + add[om-712][pique];
            yaxis_error[compteur] = (Amplitude_error/ref[om - 712] + (ref_err[om-712]*Amplitude)/(ref[om - 712]*ref[om - 712]));
            std::cout << "om = " << om << " and pique = " << pique <<  " and pic = " << pic << " and run = " << run_number << '\n';
            std::cout << " and amplitude = " << abs(Amplitude - ref[om-712])/ref[om-712] << " +- " << Amplitude_error/ref[om - 712] + (ref_err[om-712]*Amplitude)/(ref[om - 712]*ref[om - 712])<< '\n';
            std::cout << "ampl = " << Amplitude << " and ref = " << ref[om-712] << '\n';
            xaxis[compteur] = time;

            xaxis_error[compteur] = 0.01;
          }
          compteur++;

        }

        variation_pic[om-712][pique] = new TGraphErrors(n_run, xaxis, yaxis, xaxis_error, yaxis_error);
        variation_pic[om-712][pique]->SetName(Form("Gain_evolution_of_the_OM_%s_of_pic_%d_in_regard_to_pic_3", name.c_str(), pique+1));
        variation_pic[om-712][pique]->SetNameTitle(Form("Gain_evolution_of_the_OM_%s_of_pic_%d_in_regard_to_pic_3", name.c_str(), pique+1), Form("Gain_evolution_of_the_OM_%s_of_pic_%d_in_regard_to_pic_3", name.c_str(), pique+1));
        variation_pic[om-712][pique]->GetXaxis()->SetTitle("Temps (h)");
        variation_pic[om-712][pique]->GetYaxis()->SetTitle("Amplitude/Amplitude pic 3");
        variation_pic[om-712][pique]->GetXaxis()->SetTimeDisplay(1);
        variation_pic[om-712][pique]->SetMarkerColor(2);
        variation_pic[om-712][pique]->SetMarkerStyle(3);
        variation_pic[om-712][pique]->SetMarkerSize(2);
        variation_pic[om-712][pique]->GetYaxis()->SetRangeUser(0.5,1.5);

      }
        compteur = 0;
    }
  }

  // for (int om = 712; om < 717; om++) {
  //   name = namer(om);
  //
  //
  //   for (int i = 0; i < tree->GetEntries(); i++) {
  //     tree->GetEntry(i);
  //     if (pic == 3 && om_number == om) {
  //       if (run_number == 788) {
  //         ref[om-712] = Amplitude;
  //         ref_err[om-712] = Amplitude_error;
  //       }
  //       yaxis[compteur] = Amplitude;
  //       yaxis_error[compteur] = Amplitude_error;
  //       xaxis[compteur] = time;
  //       xaxis_error[compteur] = 0.01;
  //       compteur++;
  //       // std::cout << "om = " << om <<  " entry = " << compteur -1 << " and ampl = " << Amplitude << " +- " << Amplitude_error << " and time = " << time << '\n';
  //     }
  //   }
  //   variation_pic[om-712][2] = new TGraphErrors(n_run, xaxis, yaxis, xaxis_error, yaxis_error);
  //   variation_pic[om-712][2]->SetName(Form("Gain_evolution_of_the_OM_%s_of_pic_3", name.c_str()));
  //   variation_pic[om-712][2]->SetNameTitle(Form("Gain_evolution_of_the_OM_%s_of_pic_3", name.c_str()), Form("Gain_evolution_of_the_OM_%s_of_pic_3", name.c_str()));
  //   variation_pic[om-712][2]->GetXaxis()->SetTitle("Temps (h)");
  //   variation_pic[om-712][2]->GetYaxis()->SetTitle("Amplitude (ua)");
  //   // variation_pic[om-712][2]->GetYaxis()->SetRange(0,2000);
  //   // variation_pic[om-712][2]->GetXaxis()->SetRange(0,2);
  //   variation_pic[om-712][2]->GetXaxis()->SetTimeDisplay(1);
  //   variation_pic[om-712][2]->SetMarkerColor(2);
  //   variation_pic[om-712][2]->SetMarkerStyle(3);
  //   variation_pic[om-712][2]->SetMarkerSize(2);
  //   // variation_pic[om-712][2]->Draw();
  //   compteur = 0;
  // }

  file.cd();
  for (int om = 712; om < 717; om++) {
    auto canvas = new TCanvas(Form("Comparison_pic_evolution_om_%d", om),"",1600,800);
    canvas->Divide(3,2);
    for (int i = 0; i < 6; i++) {
      if (i < 2) {
        canvas->cd(i+2);
        variation_pic[om-712][i]->Draw();
        canvas->Update();
      }
      if (i == 2) {
        canvas->cd(1);
        variation_pic[om-712][2]->Draw();
        canvas->Update();
      }
      if (i > 2) {
        canvas->cd(i+1);
        variation_pic[om-712][i]->Draw();
        canvas->Update();
      }
    }
    canvas->SetGridx();
    canvas->SetGridy();
    canvas->Update();
    // return;
    canvas->Write();
    delete canvas;
  }

  file.Close();
}

int main(int argc, char const *argv[]){
  int n_run, run;
  std::vector<int> run_number, ref_run_number, ref_time;
  int compteur = 0;
  string file, correction;
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
    compteur++;
    if (compteur < n_run) {
      std::cout << "Write the runs you want" << '\n';
    }
  }

  std::cout << "Correction file name ?" << '\n';
  std::cin >> correction;

  double* ref_gain_tab_base = new double[5];
  double* ref_gain_tab = new double[5];
  Ref_corrector(run_number[0], correction, ref_gain_tab_base);
  // for (size_t i = 0; i < 5; i++) {
  //   std::cout << "ref = " << ref_gain_tab_base[i] << '\n';
  // }
  // return 0;
  compteur = 0;
  n_run = 0;

  std::cout << "Code start running" << '\n';

  for (int i = 0; i < run_number.size(); i++) {
    Ref_corrector(run_number[i], correction, ref_gain_tab);
    if (i > 0) {
      for (int j = 0; j < 5; j++) {
        ref_gain_tab[j] = ref_gain_tab_base[j]/ref_gain_tab[j];
        std::cout << j << " : % = " << ref_gain_tab_base[j]/ref_gain_tab[j] << " et ref gain = " << ref_gain_tab[j] << " et base = " << ref_gain_tab_base[j] << '\n';
      }
    }
    std::cout << "ref_gain = " << ref_gain_tab[3] << '\n';

    fit_LI_amplitude(run_number[i], ref_gain_tab);
  }
  // return 0;
  // for (int i = 0; i < ref_run_number.size(); i++) {
  //   fit_ref(ref_run_number[i], ref_time[i]);
  // }
  if (add == false) {
    file_merger(run_number);
  }
  else{
    file_merger(run_number, file);
  }


  // amplitude_variation(run_number, time);

  return 0;
}
