#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooUniform.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooRealSumPdf.h"
#include "RooParamHistFunc.h"
#include "RooHistConstraint.h"
#include "RooBinSamplingPdf.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooUniform.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooRealSumPdf.h"
#include "RooParamHistFunc.h"
#include "RooHistConstraint.h"
#include "RooBinSamplingPdf.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooAddPdf.h"
#include "RooMinimizer.h"
#include "RooMinuit.h"
#include "RooFormulaVar.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TVectorD.h"
#include <iostream>
#include <memory>
using namespace RooFit;
#define _USE_MATH_DEFINES
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
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

RooRealVar charge_tree("charge_tree","charge_tree",0,2e5);

TH2F* charge_spectre = NULL;
TH2F* charge_spectre_template = NULL;

const int gain_n_bin = 10001; // Precision au 10 000 ème
const double gain_bin_min = 0.5;
const double gain_bin_max = 1.5;
const double gain_bin_width = (gain_bin_max-gain_bin_min)/(gain_n_bin-1);

void Load_spectre(int run_number){
  TFile *file = new TFile(Form("histo_brut/histo_ref_%d.root", run_number), "READ");
  gROOT->cd();
  charge_spectre = (TH2F*)file->Get("histo_pm_charge");
  return;
}

TH1D* spectre_charge(int om_number){
  TH1D* spectre_charge = charge_spectre->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  // spectre_charge->Rebin(4);
  return spectre_charge;

}

TH2F* spectre_charge_full_template(int om_number){
  TFile *file = new TFile(Form("histo_brut/modele/Modele_OM_%d.root", om_number), "READ");
  gROOT->cd();
  charge_spectre_template = (TH2F*)file->Get("Modele_Ref_OM");
  return charge_spectre_template;
}

double roofitter(TH1D* modele, TH1D* spectre_om, int om_number, double *rootab, int run_number, double gain, double compteur)
{
  using namespace RooFit;

  RooRealVar x("x", "x", 0, 200000);
  x.setBins(1024);
  RooDataHist Tl("Tl", "Tl", x, Import(*modele));

  RooHistPdf modele_pdf ("modele_pdf", "", x, Tl);
  RooDataHist spectre_data("spectre_data", "spectre_data", x, Import(*spectre_om));

  RooChi2Var RooChi2("Chi2", "Chi2", modele_pdf, spectre_data, Range(compteur, 80000), DataError(RooAbsData::Poisson));   // create the variance
  RooMinimizer *miniChi = new RooMinimizer(RooChi2);

  double Chi2 = RooChi2.getVal();
  // double Chi2 = RooChi2.getVal()/(1024. - 1);
  std::cout << "Chi2 = " << RooChi2.getVal() << '\n';
  TCanvas* can = new TCanvas;
  can->cd();
  auto frame = x.frame(Title("Fit gain simu"));

  spectre_data.plotOn(frame,DataError(RooAbsData::SumW2), DrawOption("P"));

  modele_pdf.plotOn(frame, FillColor(0));
  spectre_data.plotOn(frame, DataError(RooAbsData::SumW2), DrawOption("P"));

  // Plot model components

  modele_pdf.plotOn(frame, LineColor(kRed), Name("sum_curve"), Range(3500, 200000));

  double Tl_int = x.getVal();
  double Tl_int_error = x.getError();
  frame->GetYaxis()->SetTitle("n events");
  frame->GetXaxis()->SetTitle("charge (u.a)");
  frame ->Draw();
  can->SetLogy();
  frame->GetYaxis()->SetRangeUser(0.01, 5e6);
  // can->SaveAs("test.root");

  rootab[1] = Tl_int;
  rootab[2] = Tl_int_error;
  // rootab[0] = RooChi2.getVal()/(1024. - 1);
  rootab[0] =Chi2;
  TLatex l;
  l.SetTextFont(40);
  l.DrawLatex(90000, 80, Form("Khi2 = %.2f", Chi2));
  // return *rootab;
  if (Chi2 < 1000) {
    can->SaveAs(Form("fit/split_10min/om_%d/run_%d/fit_run_%d_om_%d_gain_%f.png", om_number, run_number, run_number, om_number, gain));
  }


  delete miniChi;
  delete can;
  delete frame;
  // delete miniLog;
  return *rootab;

}

void Fit_Ref(int run_number) {
  // Load_spectre(run_number);
  TH1::SetDefaultSumw2();

  double Chi2, gain, gain_error;
  int om_number, compteur, ndf;
  double* rootab = new double[3];

  TH1D* modele = NULL;
  TH2F* TH2modele = NULL;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("gain_error", &gain_error);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("nsrun", &compteur);
  Result_tree.Branch("ndf", &ndf);

  for (om_number = 712; om_number < 717; om_number++) {

    TH2modele = spectre_charge_full_template(om_number);
    // for (compteur = 0; compteur < 30; compteur+=1) {
    TH1D* spectre_om = NULL;
    TFile *file = new TFile(Form("histo_brut/789_10min/om_%d_file_%d.root", om_number + 88, run_number), "READ");
    std::cout << Form("histo_brut/789_10min/om_%d_file_%d.root", om_number + 88, compteur) << '\n';
    spectre_om = (TH1D*)file->Get("spectre_om");

    gROOT->cd();
    // spectre_om = spectre_charge(om_number);
    string t = Form("fit/split_10min/om_%d/run_%d", om_number, run_number);
    if (mkdir(t.c_str(), 0777) == -1)
        cout << "directory already exist" << endl;
    else
        cout << "Directory created";
    double gainmin = 0;
    double Chi2min = 1000000;

    for (int gain_count = 1; gain_count <10001; gain_count+=100) {
      double compteur = 0;
      gain = (gain_bin_min + gain_bin_width*(gain_count-1));
      // if (gain < 1.01 && gain > 0.99) {
      modele = TH2modele->ProjectionY("modele", gain_count, gain_count);
      // modele->Draw();
      // spectre_om->Draw();

      for (int i = 0; i < 40; i++) {
        modele->SetBinContent(i,0);
        modele->SetBinError(i, 0);
        spectre_om->SetBinContent(i,0);
        spectre_om->SetBinError(i, 0);
      }

      // if (om_number == 716) {
        int k = 10;
        std::cout << spectre_om->Integral(k, k) << '\n';
        while (spectre_om->Integral(k, k) > 11000) {
        // for (int i = 0; i < 17; i++) {
          modele->SetBinContent(k,0);
          modele->SetBinError(k, 0);
          spectre_om->SetBinContent(k,0);
          spectre_om->SetBinError(k, 0);
          k++;
          compteur = k;
        }
      // }
      spectre_om->Draw();
      modele->Draw("same");
      // spectre_om->Scale(1./spectre_om->Integral());
      // modele->Scale(1./modele->Integral());
      // return;

      compteur = compteur*200000/1024;
      roofitter(modele, spectre_om, om_number, rootab, run_number, gain, compteur);
      Chi2 = rootab[0];
      if (Chi2min > Chi2){
        gainmin = gain_count;
        Chi2min = Chi2;
      }
      modele->Reset();

    }
    // gainmin = 1000;
    for (int gain_count = gainmin - 200; gain_count < gainmin + 200; gain_count++) {
      gain = (gain_bin_min + gain_bin_width*(gain_count-1));
      // if (gain < 1.01 && gain > 0.99) {
      modele = TH2modele->ProjectionY("modele", gain_count, gain_count);
      // modele->Draw();
      // spectre_om->Draw();
      for (int i = 0; i < 40; i++) {
        modele->SetBinContent(i,0);
        modele->SetBinError(i, 0);
        spectre_om->SetBinContent(i,0);
        spectre_om->SetBinError(i, 0);
      }
      // if (om_number == 716) {
      int k = 10;
      while (spectre_om->Integral(k, k) > 11000) {
        // for (int i = 0; i < 17; i++) {
        modele->SetBinContent(k,0);
        modele->SetBinError(k, 0);
        spectre_om->SetBinContent(k,0);
        spectre_om->SetBinError(k, 0);
        k++;
      }

      int compte = 0;
      for (int i = 0; i < 1000; i++) {
          if (spectre_om->GetBinContent(i) > 0) {
            compte++;
          }
      }
      ndf = compte - 1;

      // }
      roofitter(modele, spectre_om, om_number, rootab, run_number, gain, k);

      Chi2 = rootab[0];
      modele->Reset();
      Result_tree.Fill();
      // }
      delete modele;
    }

    delete spectre_om;
    // }

  }
  TFile new_file(Form("root/split_10min/Complete_Fit/Fit_Ref_%d_789.root", run_number), "RECREATE");
  // TFile new_file(Form("root/Complete_Fit/Fit_Ref_%d.root", run_number), "RECREATE");
  new_file.cd();
  Result_tree.Write();
  new_file.Close();
  std::cout << "file " << Form("root/split_10min/Complete_Fit/Fit_Ref_%d_789.root", run_number) << " saved" << '\n';
  return;

}

// void minerror_calculator(string file_name) {
//   TFile file(Form("root/split_10min/Complete_Fit/%s_789.root", file_name.c_str()), "READ");
//   gROOT->cd();
//   double PolChi2, Chi2, gain, gainmin, chi2min, error_plus, error_moins;
//   int om_number, nsrun;
//
//   TTree Result_tree("Result_tree","");
//   Result_tree.Branch("om_number", &om_number);
//   Result_tree.Branch("PolChi2", &PolChi2);
//   Result_tree.Branch("gainmin", &gainmin);
//   Result_tree.Branch("chi2min", &chi2min);
//   Result_tree.Branch("nsrun", &nsrun);
//   Result_tree.Branch("error_moins", &error_moins);
//   Result_tree.Branch("error_plus", &error_plus);
//
//   TTree* tree = (TTree*)file.Get("Result_tree");
//   tree->SetBranchStatus("*",0);
//   tree->SetBranchStatus("om_number",1);
//   tree->SetBranchAddress("om_number", &om_number);
//   tree->SetBranchStatus("Chi2",1);
//   tree->SetBranchAddress("Chi2", &Chi2);
//   tree->SetBranchStatus("gain",1);
//   tree->SetBranchAddress("gain", &gain);
//   tree->SetBranchStatus("nsrun",1);
//   tree->SetBranchAddress("nsrun", &nsrun);
//
//
//   float entry = 0;
//   // for (int i = 0; i < 30; i+=1) {
//   for (int i = 0; i < 5; i++) {
//     chi2min = 10000;
//     for (int j = 0; j < tree->GetEntries(); j++) {
//       tree->GetEntry(j);
//       if (chi2min > Chi2 && om_number == i + 712) {
//         entry = j;
//         chi2min = Chi2;
//         gainmin = gain;
//       }
//     }
//     double compteur = 0;
//     double gain2 = 1e6;
//     int compareur = 0;
//     // std::cout << "gainmin = " << gainmin << '\n';
//     while (compteur < (chi2min + 9)) {
//       // std::cout << "entry = " << compteur << '\n';
//       tree->GetEntry(entry);
//       compteur = Chi2;
//
//       if (gain2 == gain) {
//         compareur++;
//       }
//       else{
//         compareur = 0;
//       }
//       if (compareur > 5) {
//         break;
//       }
//       gain2 = gain;
//       // std::cout << "gainmin = " << gain2 << '\n';
//       error_plus = abs(gainmin-gain);
//       entry++;
//     }
//
//     compareur = 0;
//     compteur = 0;
//     while (compteur < chi2min + 9) {
//       // std::cout << "entry = " << compteur << '\n';
//       tree->GetEntry(entry);
//       compteur = Chi2;
//       gain2 = gain;
//       if (gain2 == gain) {
//         compareur++;
//       }
//       else{
//         compareur = 0;
//       }
//       if (compareur >5) {
//         break;
//       }
//       error_moins = abs(gainmin-gain);
//       // std::cout << "gainmin = " << gain2 << '\n';
//       entry--;
//
//     }
//     om_number = i + 712;
//     Result_tree.Fill();
//
//   }
//
//   TFile newfile(Form("root/split_10min/Best_Fit/%s_mini_789.root", file_name.c_str()), "RECREATE");
//   newfile.cd();
//   Result_tree.Write();
//   newfile.Close();
//
// }

int ndf_counter(int om_number, int run_number) {
  TH1D* spectre_om = NULL;
  TFile *file = new TFile(Form("histo_brut/789_10min/om_%d_file_%d.root", om_number + 88, run_number), "READ");
  spectre_om = (TH1D*)file->Get("spectre_om");

  int compteur=0;
  int k = 10;
  while (spectre_om->Integral(k, k) > 11000) {
    spectre_om->SetBinContent(k,0);
    k++;
  }


  for (int i = 0; i < 409; i++) {
    spectre_om->SetBinContent(k,0);
  }



  for (int i = 0; i < 1000; i++) {
      if (spectre_om->GetBinContent(i) > 0) {
        compteur++;
      }
  }
  return compteur;
}

void minerror_calculator(string file_name, int run_number) {
  TFile file(Form("root/split_10min/Complete_Fit/%s_789.root", file_name.c_str()), "READ");
  gROOT->cd();
  double PolChi2, Chi2, gain, gainmin, chi2min, error_plus, error_moins, minChi2, mingain_tree;
  int om_number, nsrun, ndf;
  int i;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &i);
  Result_tree.Branch("PolChi2", &PolChi2);
  Result_tree.Branch("gainmin", &gainmin);
  Result_tree.Branch("chi2min", &chi2min);
  Result_tree.Branch("nsrun", &nsrun);
  Result_tree.Branch("error_moins", &error_moins);
  Result_tree.Branch("error_plus", &error_plus);
  Result_tree.Branch("mingain_tree", &mingain_tree);
  Result_tree.Branch("ndf", &ndf);

  TTree* tree = (TTree*)file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("Chi2",1);
  tree->SetBranchAddress("Chi2", &Chi2);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &gain);
  tree->SetBranchStatus("nsrun",1);
  tree->SetBranchAddress("nsrun", &nsrun);
  tree->SetBranchStatus("ndf",1);
  tree->SetBranchAddress("ndf", &ndf);

  TFile newfile(Form("root/split_10min/Best_Fit/%s_mini_789.root", file_name.c_str()), "RECREATE");
  gROOT->cd();
  float entry = 0;

  for (i = 712; i < 717; i++) {

    minChi2 = 10000;
    for (int j = 0; j < tree->GetEntries(); j++) {
      tree->GetEntry(j);
      if (om_number == i) {
        if (minChi2 > Chi2) {
          minChi2 = Chi2;
          gainmin = gain;
        }
      }
    }

    // ndf = ndf_counter(i, run_number)-1;

    mingain_tree = gainmin;
    TF1* poly = new TF1 ("poly", "pol2", gainmin - 0.01, gainmin + 0.01);
    std::cout << "gainmin = " << gainmin << '\n';
    TH2D *Chi2gain = new TH2D ("Chi2gain", "", 1000, gainmin - 0.1, gainmin + 0.1, 1000, 0, 1000);
    tree->Project("Chi2gain", "Chi2:gain", Form("om_number == %d", i));

    TCanvas* can = new TCanvas;
    can->cd();
    Chi2gain->Draw();

    poly->Draw("same");
    Chi2gain->Fit(poly, "RQ");
    PolChi2 = poly->GetChisquare();
    gainmin = poly->GetMinimumX();
    chi2min = poly->GetMinimum();
    error_plus = poly->GetX(poly->GetMinimum()+1, poly->GetMinimumX(), poly->GetMinimumX() + 0.05) - poly->GetMinimumX();
    error_moins = poly->GetMinimumX() - poly->GetX(poly->GetMinimum()+1, poly->GetMinimumX() - 0.05, poly->GetMinimumX());
    std::cout << "min = " << poly->GetMinimum() << " ; " << poly->GetMinimumX() << " + " << error_plus << " - " << error_moins << '\n';
    std::cout << "Chi2 = " << PolChi2 << '\n';
    can->SaveAs(Form("fit/split_10min/om_%d/run_%d/fit_Chi2_gain.png", i, run_number));

    Result_tree.Fill();
    delete Chi2gain;
    delete poly;
  }

    newfile.cd();
    Result_tree.Write();
    newfile.Close();

  }

  // std::vector<int> vtime = {1655219940, 1656339281, 1656427367, 1656513370, 1656596183, 1656689629, 1656772061, 1656856728, 1658395987, 1658995707, 1659445339, 1660034367, 1660633945, 1661265097, 1661783641, 1661786628, 1661931859, 1661936023};
  std::vector<int> vtime = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109};
  std::vector<int> vrun_number = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59};

  void file_merger(std::vector<int> run_number, std::vector<int> time, string previous_file_s = "") {

    TFile file(Form("root/split_10min/Merged_Fit/Fit_Ref_%d-%d_789.root", run_number.at(0), run_number.at(run_number.size()-1)), "RECREATE");
    double Chi2, gain, gain_error_plus, gain_error_moins;
    int om_number, int_run, ndf;
    double int_time;

    TTree Result_tree("Result_tree","");
    Result_tree.Branch("om_number", &om_number);
    Result_tree.Branch("Chi2", &Chi2);
    Result_tree.Branch("gain", &gain);
    Result_tree.Branch("gain_error_plus", &gain_error_plus);
    Result_tree.Branch("gain_error_moins", &gain_error_moins);
    Result_tree.Branch("run_number", &int_run);
    Result_tree.Branch("time", &int_time);
    Result_tree.Branch("ndf", &ndf);

    // if (previous_file_s.compare("") != 0){
    //   TFile previous_file(Form("root/%s.root", previous_file_s.c_str()), "READ");
    //   TTree* previous_tree = (TTree*)previous_file.Get("Result_tree");
    //   previous_tree->SetBranchStatus("*",0);
    //   previous_tree->SetBranchStatus("om_number",1);
    //   previous_tree->SetBranchAddress("om_number", &om_number);
    //   previous_tree->SetBranchStatus("Chi2",1);
    //   previous_tree->SetBranchAddress("Chi2", &Chi2);
    //   previous_tree->SetBranchStatus("gain",1);
    //   previous_tree->SetBranchAddress("gain", &gain);
    //   previous_tree->SetBranchStatus("error_moins",1);
    //   previous_tree->SetBranchAddress("error_moins", &gain_error_moins);
    //   previous_tree->SetBranchStatus("error_plus",1);
    //   previous_tree->SetBranchAddress("error_plus", &gain_error_plus);
    //   previous_tree->SetBranchStatus("run_number",1);
    //   previous_tree->SetBranchAddress("run_number", &int_run);
    //   previous_tree->SetBranchStatus("time",1);
    //   previous_tree->SetBranchAddress("time", &int_time);
    //   for (double i = 0; i < previous_tree->GetEntries(); i++) {
    //     previous_tree->GetEntry(i);
    //     Result_tree.Fill();
    //   }
    // }
    // else {
    //   for (int i = 712; i < 717; i++) {
    //     om_number = i;
    //     int_time = time.at(0);
    //     Chi2 = 0;
    //     gain = 1;
    //     gain_error_moins = 0;
    //     gain_error_plus = 0;
    //     int_run = 716;
    //     Result_tree.Fill();
    //   }
    // }

    for (int i = 0; i < run_number.size(); i++) {
      // if (i == 1) {
      //   i = 10;
      // }
      // if (i == 12) {
      //   i = 20;
      // }
      // if (i == 26) {
      //   i = 30;
      // }
      // if (i == 42) {
      //   i = 50;
      // }
      std::cout << i << '\n';
      TFile tree_file(Form("root/split_10min/Best_Fit/Fit_Ref_%d_mini_789.root", run_number.at(i)), "READ");
      // std::cout << run_number.at(i) << '\n';
      TTree* tree = (TTree*)tree_file.Get("Result_tree");
      tree->SetBranchStatus("*",0);
      tree->SetBranchStatus("om_number",1);
      tree->SetBranchAddress("om_number", &om_number);
      tree->SetBranchStatus("chi2min",1);
      tree->SetBranchAddress("chi2min", &Chi2);
      tree->SetBranchStatus("gainmin",1);
      tree->SetBranchAddress("gainmin", &gain);
      tree->SetBranchStatus("error_moins",1);
      tree->SetBranchAddress("error_moins", &gain_error_moins);
      tree->SetBranchStatus("error_plus",1);
      tree->SetBranchAddress("error_plus", &gain_error_plus);
      tree->SetBranchStatus("ndf",1);
      tree->SetBranchAddress("ndf", &ndf);

      int_run = run_number.at(i);
      int_time = time.at(i);
      std::cout << "ok" << i+1 << '\n';
      for (int j = 0; j < 5; j++) {
        tree->GetEntry(j);
        Result_tree.Fill();
      }
    }
    file.cd();
    Result_tree.Write();
    file.Close();
  }

  void TGrapher(std::string file_name, int n_run) {
    std::cout << "ok" << '\n';
    TFile file(Form("root/split_10min/TGraph/TGraph_%s_789.root", file_name.c_str()), "RECREATE");

    TFile tree_file(Form("root/split_10min/Merged_Fit/%s_789.root", file_name.c_str()), "READ");
    double int_time;
    int om_number, run_number;
    double gain;
    double gain_error_moins, gain_error_plus;
    TTree* tree = (TTree*)tree_file.Get("Result_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("om_number",1);
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchStatus("gain",1);
    tree->SetBranchAddress("gain", &gain);
    tree->SetBranchStatus("gain_error_moins",1);
    tree->SetBranchAddress("gain_error_moins", &gain_error_moins);
    tree->SetBranchStatus("gain_error_plus",1);
    tree->SetBranchAddress("gain_error_plus", &gain_error_plus);
    tree->SetBranchStatus("run_number",1);
    tree->SetBranchAddress("run_number", &run_number);
    tree->SetBranchStatus("time",1);
    tree->SetBranchAddress("time", &int_time);

    double yaxis[n_run];
    double yaxis_error_moins[n_run];
    double yaxis_error_plus[n_run];
    double xaxis[n_run];
    double xaxis_error_moins[n_run];
    double xaxis_error_plus[n_run];
    auto canvas = new TCanvas("Allfit","",1600,800);
    canvas->Divide(3,2);
    TGraphAsymmErrors *gain_graph[5]; //(n_run, xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);

    file.cd();
    int compteur = 0;
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < n_run; j++){
        tree->GetEntry(i+j*5);
        compteur++;
        // std::cout << "compteur " << i+j*5 << " and om = " << om_number << " and run = " << run_number << '\n';
        // std::cout << gain << '\n';
        yaxis[j] = gain;
        yaxis_error_moins[j] = gain_error_moins;
        yaxis_error_plus[j] = gain_error_plus;
        xaxis[j] = int_time;
        xaxis_error_plus[j] = 0.5;
        xaxis_error_moins[j] = 0.5;
      }
      // for (size_t k = 0; k < n_run; k++) {
      //   std::cout << yaxis[k] << '\n';
      // }
      gain_graph[i] = new TGraphAsymmErrors(n_run, xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);

      gain_graph[i]->SetName(Form("fit_OM_ref_%d", om_number+i));
      gain_graph[i]->SetNameTitle(Form("fit_OM_ref_%d", om_number+i), Form("Gain evolution of the OM %d with regard to the background", om_number));
      // gain_graph->GetXaxis()->SetTimeDisplay(1);
      gain_graph[i]->GetXaxis()->SetTitle("Time (min)");
      gain_graph[i]->GetYaxis()->SetTitle("Gain evolution");
      gain_graph[i]->GetYaxis()->SetRangeUser(0.9, 1.1);
      gain_graph[i]->SetMarkerColor(2);
      gain_graph[i]->SetMarkerStyle(5);
      gain_graph[i]->SetMarkerSize(2);
      canvas->cd(i+1);
      gain_graph[i]->Draw("AP");

      gain_graph[i]->Write();
      canvas->Update();
      // if (i == 4) {
      //   canvas->Write();
      // }
    }

    canvas->Update();
    canvas->Write();
    file.Close();

  }

//
// void TGrapherpart(std::string file_name, int n_run) {
//   TFile file(Form("root/split_10min/TGraph/multiTGraph_%s_789_om_%d.root", file_name.c_str(), om + 712), "RECREATE");
//
//   TFile tree_file(Form("root/split_10min/Merged_Fit/%s_789.root", file_name.c_str()), "READ");
//   double int_time;
//   int om_number, run_number;
//   double gain;
//   double gain_error_moins, gain_error_plus;
//   TTree* tree = (TTree*)tree_file.Get("Result_tree");
//   tree->SetBranchStatus("*",0);
//   tree->SetBranchStatus("om_number",1);
//   tree->SetBranchAddress("om_number", &om_number);
//   tree->SetBranchStatus("gain",1);
//   tree->SetBranchAddress("gain", &gain);
//   tree->SetBranchStatus("gain_error_moins",1);
//   tree->SetBranchAddress("gain_error_moins", &gain_error_moins);
//   tree->SetBranchStatus("gain_error_plus",1);
//   tree->SetBranchAddress("gain_error_plus", &gain_error_plus);
//   tree->SetBranchStatus("run_number",1);
//   tree->SetBranchAddress("run_number", &run_number);
//   tree->SetBranchStatus("time",1);
//   tree->SetBranchAddress("time", &int_time);
//
//   int start, end, start_tree, i;
//   double error_time;
//   // if (j == 0) {
//   //   start = 0;
//   //   end = 1;
//   //   error_time = 30;
//   //   start_tree = 0;
//   // }
//   // if (j == 1) {
//   //   start = 10;
//   //   end = 12;
//   //   error_time = 15;
//   //   start_tree = 5;
//   // }
//   // if (j == 2) {
//   //   start = 20;
//   //   end = 26;
//   //   error_time = 5;
//   //   start_tree = 15;
//   // }
//   // if (j == 3) {
//   //   start = 30;
//   //   end = 42;
//   //   error_time = 2.5;
//   //   start_tree = 45;
//   // }
//   // if (j == 4) {
//   //   start = 50;
//   //   end = 110;
//   //   error_time = 0.5;
//   //   start_tree = 105;
//   // }
//   error_time = 5;
//   start = 0;
//   end = 73;
//   std::vector<double> vyaxis, vyaxis_error_moins, vyaxis_error_plus, vxaxis, vxaxis_error_moins, vxaxis_error_plus;
//
//
//   file.cd();
//   // std::cout << "part = " << j << '\n';
//   std::cout << "end - start = " << end - start << '\n';
//
//   for (int i = 0; i < 5; i++) {
//
//     for (int k = 0; k < end - start; k++) {
//       std::cout << "k = " << k << '\n';
//       tree->GetEntry(i + start_tree + k*5);
//       std::cout << "entry = " << i + start_tree + k*5 << '\n';
//       std::cout << "gain = " << gain << '\n';
//
//       vyaxis.push_back(gain);
//       vyaxis_error_moins.push_back(gain_error_moins);
//       vyaxis_error_plus.push_back(gain_error_plus);
//       vxaxis.push_back(error_time + error_time*k*2);
//       std::cout << "time = " << vxaxis[k] << '\n';
//       vxaxis_error_plus.push_back(error_time);
//       vxaxis_error_moins.push_back(error_time);
//       std::cout << "error+ = " << vyaxis_error_plus[k] << " and error- = " << vyaxis_error_moins[k] << '\n';
//     }
//     int taille = vyaxis.size();
//     double yaxis[taille];
//     double yaxis_error_moins[taille];
//     double yaxis_error_plus[taille];
//     double xaxis[taille];
//     double xaxis_error_moins[taille];
//     double xaxis_error_plus[taille];
//
//     for (int v = 0; v < taille; v++) {
//
//       xaxis[v] = vxaxis.at(v);
//       yaxis[v] = vyaxis.at(v);
//       yaxis_error_moins[v] = vyaxis_error_moins.at(v);
//       yaxis_error_plus[v] = vyaxis_error_plus.at(v);
//       xaxis_error_moins[v] = vxaxis_error_moins.at(v);
//       xaxis_error_plus[v] = vxaxis_error_plus.at(v);
//     }
//
//
//     TGraphAsymmErrors *gain_graph = new TGraphAsymmErrors(taille, xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);
//
//
//
//     // std::cout << om_number << '\n';
//     // std::cout << "yaxis " << yaxis[0] << '\n';
//
//     gain_graph->SetName(Form("fit_OM_ref_%d", om_number));
//     gain_graph->SetNameTitle(Form("fit_OM_ref_%d", om_number), Form("Gain evolution of the OM %d with regard to the background", om_number));
//     gain_graph->GetXaxis()->SetTitle("Time");
//     gain_graph->GetYaxis()->SetTitle("Gain evolution");
//     gain_graph->GetXaxis()->SetRangeUser(-5, 1000);
//     gain_graph->SetMarkerColor(1);
//     gain_graph->SetLineColor(1);
//     gain_graph->SetMarkerStyle(5);
//     gain_graph->SetMarkerSize(2);
//
//
//     gain_graph->Draw("AP");
//     gain_graph->Write();
//     delete gain_graph
//   }
//   file.Close();
//
//
// }

void multigrapher(std::string file_name) {
  TMultiGraph *mg = new TMultiGraph();
  int om_number = 712;
  for (int i = 0; i < 5; i++) {
    TFile file(Form("root/split_10min/TGraph/multiTGraph_%s_789_part_%d.root", file_name.c_str(), i), "READ");
    TGraphAsymmErrors* graph = NULL;
    graph = (TGraphAsymmErrors*)file.Get(Form("fit_OM_ref_%d_part_%d", om_number, i));
    mg->Add(graph, "AP");
    delete graph;

  }
  // mg->Draw();
  TFile new_file("root/split_10min/TGraph/test.root", "RECREATE");
  new_file.cd();
  mg->Write();
    new_file.Close();
  }


  int main(int argc, char const *argv[]) {                            //splitted file
    int n_run, run, t;
    std::vector<int> run_number, time, ref_run_number, ref_time;
    int compteur = 0;
    std::string file;
    bool add = false;
    n_run = 74;
    ///// Create new file
    int j = 0;
    while (j < n_run ) {
      run_number.push_back(j);
      time.push_back(j);
      j++;
    }
    std::cout << "Code start running" << '\n';
    // std::cout << run_number.size() << '\n';
    // for (size_t i = 0; i < run_number.size(); i++) {
    //   std::cout << run_number[i] << '\n';
    // }
    // return 0;
    for (int i = 73; i < n_run-1; i++) {
    //   if (i == 1) {
    //     i = 10;
    //   }
    //   // std::cout <<  << '\n';
    //   if (i == 12) {
    //     i = 20;
    //   }
    //   if (i == 26) {
    //     i = 30;
    //   }
    //   if (i == 42) {
    //     i = 50;
    //   }
      Fit_Ref(run_number.at(i));
      std::cout << "Fit Ref " << i << " ok" << '\n';
    }
    // for (int i = 0; i < n_run-1; i++) {
    //   // // std::cout << i << '\n';
    //   // if (i == 1) {
    //   //   i = 10;
    //   // }
    //   // if (i == 12) {
    //   //   i = 20;
    //   // }
    //   // if (i == 26) {
    //   //   i = 30;
    //   // }
    //   // if (i == 42) {
    //   //   i = 50;
    //   // }
    //   minerror_calculator(Form("Fit_Ref_%d", run_number.at(i)), run_number.at(i));
    //   std::cout << "minerror " << i << " ok" << '\n';
    // }
    // return 0;

    std::cout << "Fit_Ref and minerror ok" << '\n';

    file_merger(run_number, vtime);
    // std::cout << "file_merger ok" << '\n';
    // for (int j = 0; j < 5; j++) {
      // for (int om_num = 0; om_num < 5; om_num++) {
        TGrapher(Form("Fit_Ref_%d-%d", run_number.at(0), run_number.at(run_number.size()-1)), n_run-1);
        std::cout << "part " << j << " done" << '\n';
      // }
    // }

    // multigrapher(Form("Fit_Ref_%d-%d", run_number.at(0), run_number.at(run_number.size()-1)));

    std::cout << "Finish !!!" << '\n';
    return 0;
  }
