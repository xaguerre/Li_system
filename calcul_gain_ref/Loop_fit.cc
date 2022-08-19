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
using namespace std;

RooRealVar charge_tree("charge_tree","charge_tree",0,2e5);

TH2F* charge_spectre = NULL;
TH2F* charge_spectre_template = NULL;

const int gain_n_bin = 2001; // Precision au 10 000 ème
const float gain_bin_min = 0.9;
const float gain_bin_max = 1.1;
const float gain_bin_width = (gain_bin_max-gain_bin_min)/(gain_n_bin-1);

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

double roofitter(TH1D* modele, TH1D* spectre_om, int om_number, double *rootab, int run_number, double gain)
{
  using namespace RooFit;
  std::cout << "/* message */" << '\n';
  RooRealVar x("x", "x", 0, 200000);
  x.setBins(1024);
  RooDataHist Tl("Tl", "Tl", x, Import(*modele));

  RooRealVar modele_rrv("modele_rrv", "Tl", 0.2, 0.1, 0.4);

  RooHistPdf modele_pdf ("modele_pdf", "", x, Tl);
  RooDataHist spectre_data("spectre_data", "spectre_data", x, Import(*spectre_om));


  RooChi2Var RooChi2("Chi2", "Chi2", modele_pdf, spectre_data, Range(3500, 60000), DataError(RooAbsData::Poisson));   // create the variance
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

  int startx = 0;
  modele_pdf.plotOn(frame, LineColor(kRed), Name("sum_curve"), Range(3500, 200000));

  double Tl_int = modele_rrv.getVal();
  double Tl_int_error = modele_rrv.getError();
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
  can->SaveAs(Form("fit/test_full/fit_om_%d_gain_%f.png"/*, om_number*/, om_number, gain));

  delete miniChi;
  delete can;
  delete frame;
  // delete miniLog;
  return *rootab;

}

void Fit_Ref(int run_number) {
  Load_spectre(run_number);
  TH1::SetDefaultSumw2();

  double time, charge_tree, amplitude_tree, Chi2, gain, gain_error;
  int om_number, compteur;
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
  for (om_number = 712; om_number < 717; om_number++) {
    TH2modele = spectre_charge_full_template(om_number);
    for (compteur = 1; compteur < 60; compteur++) {
    TH1D* spectre_om = NULL;
    TFile *file = new TFile(Form("histo_brut/716/file_%d.root", compteur+165), "READ");
    // TFile *file = new TFile(Form("histo_brut/histo_ref_%d.root", run_number), "READ");
    gROOT->cd();
    spectre_om = spectre_charge(om_number);
    for (int gain_count = 1; gain_count <2001; gain_count++) {
      gain = (gain_bin_min + gain_bin_width*(gain_count-1));
      // if (gain < 1.01 && gain > 0.99) {
        modele = TH2modele->ProjectionY("modele", gain_count, gain_count);
        modele->Draw();
        spectre_om->Draw();

        for (int i = 0; i < 17; i++) {
          modele->SetBinContent(i,0);
          // modele->SetBinError(i, 0);
          spectre_om->SetBinContent(i,0);
          // spectre_om->SetBinError(i, 0);
        }
        std::cout << "/* message */" << '\n';
        roofitter(modele, spectre_om, om_number, rootab, run_number, gain);
        Chi2 = rootab[0];
        modele->Reset();
        Result_tree.Fill();
      // }
    }

    delete spectre_om;
    }

  }
  // TFile new_file(Form("root/Fit_Ref_%d.root", run_number), "RECREATE");
  TFile new_file(Form("root/Complete_Fit/Fit_Ref_%d.root", run_number), "RECREATE");
  new_file.cd();
  Result_tree.Write();
  new_file.Close();

  return;

}

void minimum_calculator(string file_name) {
  TFile file(Form("root/Complete_Fit/%s", file_name.c_str()), "READ");
  gROOT->cd();
  double PolChi2, Chi2, gain, mingain, minChi2, gainmin, chi2min;
  int om_number, nsrun;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("PolChi2", &PolChi2);
  Result_tree.Branch("gainmin", &gainmin);
  Result_tree.Branch("chi2min", &chi2min);
  Result_tree.Branch("nsrun", &nsrun);

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

  for (int i = 1; i < 30; i++) {
    minChi2 = 10000;
    for (int j = 0; j < tree->GetEntries(); j++) {
      tree->GetEntry(j);
      if (nsrun == i) {
        if (minChi2 > Chi2) {
          minChi2 = Chi2;
          mingain = gain;
        }
      }
    }

    TF1* poly = new TF1 ("poly", "pol2", mingain - 0.005, mingain + 0.005);
    std::cout << "mingain = " << mingain << '\n';
    TH2D *Chi2gain = new TH2D ("Chi2gain", "", 1000, 0.99, 1.01, 1000, 0, 1000);
    tree->Project("Chi2gain", "Chi2:gain", Form("nsrun == %d", i));
    TCanvas* can = new TCanvas;
    can->cd();
    Chi2gain->Draw();

    poly->Draw("same");
    Chi2gain->Fit(poly, "RQ");

    PolChi2 = poly->GetChisquare();
    std::cout << "min = " << poly->GetMinimum() << " ; " << poly->GetMinimumX() << '\n';
    gainmin = poly->GetMinimumX();
    chi2min = poly->GetMinimum();
    nsrun = i;
    Result_tree.Fill();

    can->SaveAs(Form("fit/fit_Chi2/fit_om_%d_part_%d.png", om_number, i));

    delete can;
    delete Chi2gain;
    delete poly;
  }

  TFile newfile(Form("root/Best_Fit/%s_mini.root", file_name.c_str(), "RECREATE");
  newfile.cd();
  Result_tree.Write();
  newfile.Close();

}

std::vector<int> vtime = {1655219940, 1656339281, 1656427367, 1656513370, 1656596183, 1656689629, 1656772061, 1656856728, 1658395987, 1658995707, 1659445339, 1660034367};
std::vector<int> vrun_number = {736, 737, 738, 739, 740, 741, 742};

void file_merger(std::vector<int> run_number, std::vector<int> time, string previous_file_s = "") {
  TFile file(Form("root/Merged_Fit/Fit_Ref_%d-%d.root", run_number.at(0), run_number.at(run_number.size()-1)), "RECREATE");

  double Chi2, gain, gain_error;
  int om_number, int_run;
  double int_time;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("gain_error", &gain_error);
  Result_tree.Branch("run_number", &int_run);
  Result_tree.Branch("time", &int_time);

  if (previous_file_s.compare("") != 0){
    TFile previous_file(Form("root/%s.root", previous_file_s.c_str()), "READ");
    TTree* previous_tree = (TTree*)previous_file.Get("Result_tree");
    previous_tree->SetBranchStatus("*",0);
    previous_tree->SetBranchStatus("om_number",1);
    previous_tree->SetBranchAddress("om_number", &om_number);
    previous_tree->SetBranchStatus("Chi2",1);
    previous_tree->SetBranchAddress("Chi2", &Chi2);
    previous_tree->SetBranchStatus("gain",1);
    previous_tree->SetBranchAddress("gain", &gain);
    previous_tree->SetBranchStatus("gain_error",1);
    previous_tree->SetBranchAddress("gain_error", &gain_error);
    previous_tree->SetBranchStatus("run_number",1);
    previous_tree->SetBranchAddress("run_number", &int_run);
    previous_tree->SetBranchStatus("time",1);
    previous_tree->SetBranchAddress("time", &int_time);
    for (double i = 0; i < previous_tree->GetEntries(); i++) {
      previous_tree->GetEntry(i);
      Result_tree.Fill();
    }
  }
  else {
    for (int i = 712; i < 717; i++) {
      om_number = i;
      int_time = time.at(0);
      Chi2 = 0;
      gain = 1;
      gain_error = 0;
      int_run = 716;
      Result_tree.Fill();
    }
  }

  for (int i = 0; i < run_number.size(); i++) {
    TFile tree_file(Form("root/Best_Fit/Fit_Ref_mini_%d.root", run_number.at(i)), "READ");

    TTree* tree = (TTree*)tree_file.Get("Result_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("om_number",1);
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchStatus("Chi2",1);
    tree->SetBranchAddress("Chi2", &Chi2);
    tree->SetBranchStatus("gain",1);
    tree->SetBranchAddress("gain", &gain);
    tree->SetBranchStatus("gain_error",1);
    tree->SetBranchAddress("gain_error", &gain_error);

    int_run = run_number.at(i);
    int_time = time.at(i+1);
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

  TFile file(Form("root/TGraph/TGraph_%s.root", file_name.c_str()), "RECREATE");

  TFile tree_file(Form("root/Merged_Fit/%s.root", file_name.c_str()), "READ");
  double int_time;
  int om_number, run_number;
  double gain;
  double gain_error;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &gain);
  tree->SetBranchStatus("gain_error",1);
  tree->SetBranchAddress("gain_error", &gain_error);
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &int_time);

  double yaxis[n_run];
  double yaxis_error[n_run];
  double xaxis[n_run];
  double xaxiserror[n_run];

  file.cd();
  int compteur = 0;
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < n_run; j++){
      tree->GetEntry(i+j*5);
      compteur++;
      std::cout << "compteur " << i+j*5 << " and om = " << om_number << " and run = " << run_number << '\n';
      std::cout << gain << '\n';
      yaxis[j] = gain;
      yaxis_error[j] = 0.001;
      xaxis[j] = int_time;
      xaxiserror[j] = 1;
    }
    for (size_t k = 0; k < n_run; k++) {
        std::cout << yaxis[k] << '\n';
    }

    TGraphErrors gain_graph(n_run, xaxis, yaxis, xaxiserror, yaxis_error);
    gain_graph.SetName(Form("fit_OM_ref_%d", om_number));
    gain_graph.SetNameTitle(Form("fit_OM_ref_%d", om_number), Form("Gain evolution of the OM %d", om_number));
    gain_graph.GetXaxis()->SetTimeDisplay(1);
    gain_graph.GetXaxis()->SetTitle("Time");
    gain_graph.GetYaxis()->SetTitle("Gain evolution");
    gain_graph.SetMarkerColor(2);
    gain_graph.SetMarkerStyle(34);
    gain_graph.SetMarkerSize(2);


    // TCanvas* canvas2 = new TCanvas;
    // gain_graph.Draw();
    // canvas2->SaveAs(Form("fit/fit_Tl/variation/charge_fit_om_%03d.png", j));
    gain_graph.Write();

  }
  file.Close();
}

void comparator(){

  TFile file("variation_gain/variation_gain_716.root", "READ");
  gROOT->cd();
  double Chi2, gain, mingain, minChi2, gainmin, measured_gain, fitted_gain;
  double measured_min[29], fitted_min[29];
  int om_number, nsrun;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("measuredgain", &measured_gain);
  Result_tree.Branch("fittedgain", &fitted_gain);
  Result_tree.Branch("nsrun", &nsrun);

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


  TH2D *Measuredgain = new TH2D ("Measuredgain", "", 1000, 0, 29, 1000, 0.99, 1.01);
  for (int i = 1; i < 30; i++) {
    minChi2 = 10000;
    for (int j = 0; j < tree->GetEntries(); j++) {
      tree->GetEntry(j);
      if (nsrun == i) {
        if (minChi2 > Chi2) {
          minChi2 = Chi2;
          mingain = gain;
        }
      }
      measured_min[i] = mingain;
    }
    Measuredgain->Fill(i, measured_min[i]);
    std::cout << i << '\n';
  }

  TFile file2("root/Fittedgain.root", "READ");
  gROOT->cd();
  TTree* tree2 = (TTree*)file2.Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);
  tree2->SetBranchStatus("gainmin",1);
  tree2->SetBranchAddress("gainmin", &gainmin);
  tree2->SetBranchStatus("nsrun",1);
  tree2->SetBranchAddress("nsrun", &nsrun);
  TH2D *Tcomparator = new TH2D ("Chi2gain", "", 1000, 0, 29, 1000, -0.01, 0.01);
  TH2D *Fittedgain = new TH2D ("Fittedgain", "", 1000, 0, 29, 1000,  0.99, 1.01);
  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    Tcomparator->Fill(nsrun, measured_min[i]-gainmin);
    Fittedgain->Fill(nsrun, gainmin);

    measured_gain = measured_min[i];
    fitted_gain = gainmin;
    Result_tree.Fill();

  }
  Tcomparator->Draw();
  Tcomparator->SetMarkerStyle(2);
  TCanvas* c = new TCanvas;
  Measuredgain->Draw();
  Measuredgain->SetMarkerStyle(2);
  Measuredgain->SetMarkerColor(4);
  Fittedgain->Draw("same");
  Fittedgain->SetMarkerStyle(4);
  Fittedgain->SetMarkerColor(2);

  TFile newfile("root/comparator.root", "RECREATE");
  newfile.cd();
  Tcomparator->Write();
  Measuredgain->Write();
  Fittedgain->Write();
  Result_tree.Write();
  newfile.Close();
  return;
}



int main(int argc, char const *argv[]) {
  int n_run, run, t;
  std::vector<int> run_number, time, ref_run_number, ref_time;
  int compteur = 0;
  std::string file;
  bool add = false;

  for(int i = 0; i<argc; i++){
    if (std::string(argv[i]) == "-add" ) {
      file = argv[i+1];
      std::cout << file << '\n';
      add = true;
    }
  }

  if (add == true) {                  ///// Add run to the other runs
    string old_run;
    std::cout << "To what file do you want to add the new run(s)?" << '\n';
    std::cin >> old_run;

    std::cout << "How many run do you want to add?" << '\n';
    std::cin >> n_run;
    std::cout << "Write the run(s) you want" << '\n';
    while (compteur < n_run && cin >> run) {
      run_number.push_back(run);
      std::cout << "Write the delay between the run(s) (run time beginning)" << '\n';
      while (compteur < n_run && cin >> t) {
        time.push_back(t);
        break;
      }
      compteur++;
      if (compteur < n_run) {
        std::cout << "Write the runs you want" << '\n';
      }
    }
    std::cout << "Code start running" << '\n';

    for (int i = 0; i < n_run; i++) {
      Fit_Ref(run_number.at(i));
      minimum_calculator(Form("Fit_Ref_%d", run_number.at(i)));
    }
    std::cout << "Fit_Ref ok" << '\n';

    file_merger(run_number, time, old_run);
    std::cout << "file_merger ok" << '\n';

    TGrapher(Form("Fit_ref_716-%d.root", run_number.at(run_number.size())), n_run);

  }


  else {                            ///// Create new file
    std::cout << "How many run do you want ?" << '\n';
    std::cin >> n_run;
    std::cout << "Write the runs you want" << '\n';
    while (compteur < n_run && cin >> run) {
      run_number.push_back(run);
      // std::cout << "Write the delay between the runs (0 for the first)" << '\n';
      // while (compteur < n_run && cin >> t) {
      //   time.push_back(t);
      //   break;
      // }
      compteur++;
      if (compteur < n_run) {
        std::cout << "Write the runs you want" << '\n';
      }
    }
    std::cout << "Code start running" << '\n';

    for (int i = 0; i < n_run; i++) {
      Fit_Ref(run_number.at(i));
      minimum_calculator(Form("Fit_Ref_%d", run_number.at(i)));
    }
    std::cout << "Fit_Ref ok" << '\n';

    file_merger(run_number, vtime);
    std::cout << "file_merger ok" << '\n';

    TGrapher(Form("Fit_ref_%d-%d.root", run_number.at(0), run_number.at(run_number.size())), n_run);
  }

  return 0;
}
