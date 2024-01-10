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
#include <sys/stat.h>
#include <sys/types.h>
#include <TParameter.h>
using namespace std;

RooRealVar charge_tree("charge_tree","charge_tree",0,2e5);

TH2F* charge_spectre = NULL;
TH2D* charge_spectre_template = NULL;

const int gain_n_bin = 10001; // Precision au 10 000 ème
const double gain_bin_min = 0.5;
const double gain_bin_max = 1.5;
const double gain_bin_width = (gain_bin_max-gain_bin_min)/(gain_n_bin-1);

void Load_spectre(int run_number){
  TFile *file = new TFile(Form("../histo_brut/histo_ref_%d.root", run_number), "READ");
  gROOT->cd();
  charge_spectre = (TH2F*)file->Get("histo_pm_charge");
  return;
}

TH1D* spectre_charge(int om_number){
  TH1D* spectre_charge = charge_spectre->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  // spectre_charge->Rebin(4);
  return spectre_charge;
}

TH2D* spectre_charge_full_template(int om_number){
  TFile *file = new TFile(Form("../histo_brut/modele/Modele_OM_%d_936.root", om_number), "READ");
  gROOT->cd();
  charge_spectre_template = (TH2D*)file->Get("Modele_Ref_OM");
  return charge_spectre_template;
}

double roofitter(TH1D* modele, TH1D* spectre_om, int om_number, double *rootab, double gain, double compteur, int run_number)
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
  rootab[0] = Chi2;
  TLatex l;
  l.SetTextFont(40);
  l.DrawLatex(90000, 80, Form("Khi2 = %.2f", Chi2));
  // return *rootab;
  if (Chi2 < 3000) {
    can->SaveAs(Form("png/fit_templ/run_%d/fit_run_%d_gain_%f.png", run_number, run_number, gain));
  }


  delete miniChi;
  delete can;
  delete frame;
  // delete miniLog;
  return *rootab;

}

void Fit_Ref(int run_number) {
  Load_spectre(run_number);
  TH1::SetDefaultSumw2();

  double Chi2, gain, gain_error, ndf, time;
  int om_number, compteur;
  double* rootab = new double[3];

  TH1D* modele = NULL;
  TH2D* TH2modele = NULL;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("gain_error", &gain_error);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("nsrun", &compteur);
  Result_tree.Branch("ndf", &ndf);
  Result_tree.Branch("time", &time);

  for (om_number = 713; om_number < 714; om_number++) {


    // for (compteur = 0; compteur < 30; compteur+=1) {
    TH1D* spectre_om = NULL;
    spectre_om = spectre_charge(713);
    TFile *file = new TFile(Form("../histo_brut/histo_ref_%d.root", run_number), "READ");
    // std::cout << Form("histo_brut/histo_ref_%d.root", run_number) << '\n';
    // spectre_om = (TH1D*)file->Get("histo_pm_charge");

    time = 0;
    TParameter<double> *param = new TParameter<double>("start_time", time);
    param = (TParameter<double>*)(file->Get("start_time"));
    time = param->GetVal();
    TH2modele = spectre_charge_full_template(713);


    gROOT->cd();

    string t = Form("png/fit_templ/run_%d", run_number);
    if (mkdir(t.c_str(), 0777) == -1)
        cout << "directory already exist" << endl;
    else
        cout << "Directory created" << endl;
    double gainmin = 0;
    double Chi2min = 1000000;
    for (int gain_count = 0; gain_count <10001; gain_count+=100) {
      double compteur = 0;
      gain = (gain_bin_min + gain_bin_width*(gain_count-1));
      // if (gain < 1.01 && gain > 0.99) {
      modele = TH2modele->ProjectionY("modele", gain_count, gain_count);
      modele->Draw();
      spectre_om->Draw("same");
      for (int i = 0; i < 61; i++) {    // 25 full, 61, part
        modele->SetBinContent(i,0);
        modele->SetBinError(i, 0);
        spectre_om->SetBinContent(i,0);
        spectre_om->SetBinError(i, 0);
      }
      // if (om_number == 716) {
        int k = 10;
        std::cout << spectre_om->Integral(k, k) << '\n';
        // while (spectre_om->Integral(k, k) > 11000) {
        // // for (int i = 0; i < 17; i++) {
        //   modele->SetBinContent(k,0);
        //   modele->SetBinError(k, 0);
        //   spectre_om->SetBinContent(k,0);
        //   spectre_om->SetBinError(k, 0);
        //   k++;
        //   compteur = k;
        // }
            std::cout << "/* message */" << '\n';
      // }
      // spectre_om->Draw();
      // modele->Draw("same");
      // spectre_om->Scale(1./spectre_om->Integral());
      // modele->Scale(1./modele->Integral());
      // return;

      roofitter(modele, spectre_om, om_number, rootab, gain, compteur, run_number);
      std::cout << "Chi2 = " << rootab[0] << '\n';
      // return;
      Chi2 = rootab[0];
      if (Chi2min > Chi2){
        gainmin = gain_count;
        Chi2min = Chi2;
      }
      modele->Reset();

    }
    // gainmin = 1000;
    for (int gain_count = gainmin - 400; gain_count < gainmin + 400; gain_count++) {
      gain = (gain_bin_min + gain_bin_width*(gain_count-1));
      modele = TH2modele->ProjectionY("modele", gain_count, gain_count);
      // if (gain < 1.01 && gain > 0.99) {
      // modele = TH2modele->ProjectionY("modele", gain_count, gain_count);
      // modele->Draw();
      // spectre_om->Draw();
      for (int i = 0; i < 61; i++) {
        modele->SetBinContent(i,0);
        modele->SetBinError(i, 0);
        spectre_om->SetBinContent(i,0);
        spectre_om->SetBinError(i, 0);
      }

      // if (om_number == 716) {
      int k = 10;
      // while (spectre_om->Integral(k, k) > 11000) {
      //   // for (int i = 0; i < 17; i++) {
      //   modele->SetBinContent(k,0);
      //   modele->SetBinError(k, 0);
      //   spectre_om->SetBinContent(k,0);
      //   spectre_om->SetBinError(k, 0);
      //   k++;
      // }
      // }

      int compte = 0;
      for (int i = 0; i < 1000; i++) {
          if (spectre_om->GetBinContent(i) > 0) {
            compte++;
          }
      }
      ndf = compte - 1;

      roofitter(modele, spectre_om, om_number, rootab, gain, k, run_number);

      Chi2 = rootab[0];
      modele->Reset();
      Result_tree.Fill();
      // }
      delete modele;
    }

    delete spectre_om;
    // }

  }
  TFile new_file(Form("fit/fit_templ/Complete_Fit/Fit_Ref_%d.root", run_number), "RECREATE");
  // TFile new_file(Form("root/Complete_Fit/Fit_Ref_%d.root", run_number), "RECREATE");
  new_file.cd();
  Result_tree.Write();
  new_file.Close();
  std::cout << "file " << Form("fit/fit_templ/Complete_Fit/Fit_Ref_%d.root", run_number) << " saved" << '\n';
  return;

}

void minerror_calculator(string file_name, int run_number) {
  TFile file(Form("fit/fit_templ/Complete_Fit/%s.root", file_name.c_str()), "READ");
  gROOT->cd();
  double PolChi2, Chi2, gain, gainmin, chi2min, error_plus, error_moins, minChi2, mingain_tree, time, ndf;
  int om_number, nsrun;
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
  Result_tree.Branch("time", &time);

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
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);


  TFile newfile(Form("fit/fit_templ/Best_Fit/%s_mini.root", file_name.c_str()), "RECREATE");
  gROOT->cd();

  for (i = 713; i < 714; i++) {

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
    if (run_number == 789) {
      Chi2gain->GetYaxis()->SetRangeUser(800,3000);
    }
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
    can->SaveAs(Form("png/fit_templ/run_%d/fit_Chi2_gain.png", run_number));

    Result_tree.Fill();
    delete Chi2gain;
    delete poly;
  }

    newfile.cd();
    Result_tree.Write();
    newfile.Close();

  }

std::vector<int> vrun_number = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59};

void file_merger(std::vector<int> run_number) {


  TFile file(Form("fit/fit_templ/Merged_Fit/Fit_Ref_%d-%d.root", run_number.at(0), run_number.at(run_number.size()-1)), "RECREATE");
// TFile file(Form("root/Merged_Fit/Fit_Ref_%d-%d.root", 736, 836), "RECREATE");
  double Chi2, gain, gain_error_plus, gain_error_moins, ndf, correction;
  int om_number, int_run;
  double time;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("gain_error_plus", &gain_error_plus);
  Result_tree.Branch("gain_error_moins", &gain_error_moins);
  Result_tree.Branch("run_number", &int_run);
  Result_tree.Branch("time", &time);
  Result_tree.Branch("ndf", &ndf);

  for (int i = 0; i < run_number.size(); i++) {
    TFile tree_file(Form("fit/fit_templ/Best_Fit/Fit_Ref_%d_mini.root", run_number.at(i)), "READ");
    std::cout << run_number.at(i) << '\n';
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
    tree->SetBranchStatus("time",1);
    tree->SetBranchAddress("time", &time);

    int_run = run_number.at(i);
    std::cout << "ok " << i+1 << '\n';
    tree->GetEntry(0);
    if (i == 0)correction = gain-1;
    std::cout << "correction = " << correction << " and gain = " << gain << '\n';
    gain = gain - correction;
    Result_tree.Fill();
  }
  file.cd();
  Result_tree.Write();
  file.Close();
}

int main(int argc, char const *argv[]) {
  int n_run, run;
  std::vector<int> run_number, time, ref_run_number, ref_time;
  int compteur = 0;
  std::string file;

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

  std::cout << "Code start running" << '\n';
  // Fit_Ref(run_number.at(n_run - 1));
  for (int i = 0; i < n_run; i++) {
    Fit_Ref(run_number.at(i));
    minerror_calculator(Form("Fit_Ref_%d", run_number.at(i)), run_number.at(i));
  }

  std::cout << "Fit_Ref and minerror ok" << '\n';

  file_merger(run_number);
  std::cout << "file_merger ok" << '\n';

  std::cout << "Finish !!!" << '\n';
  return 0;
}
