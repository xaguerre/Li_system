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
#include "RooParametricBinningStepFunction.h"
using namespace std;

RooRealVar charge_tree("charge_tree","charge_tree",0,2e5);

RooDataSet *SetTable[200][5];

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

TH1D* spectre_charge_full(int om_number){
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

double roofitter(TH1D* modele, TH1D* spectre_om, int om_number, double *rootab, int run_number, double gain, int compteur)
{
  using namespace RooFit;

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
  l.DrawLatex(90000, 80, Form("Khi2/NDF = %.2f", Chi2));
  // return *rootab;
  can->SaveAs(Form("fit/variation_run_time_716/fit_om_%d_gain_%f_%d.png"/*, om_number*/, om_number, gain, compteur));

  delete miniChi;
  delete can;
  delete frame;
  // delete miniLog;
  return *rootab;

}

void Fit_Ref(int run_number) {
  Load_spectre(run_number);
  // TH1::SetDefaultSumw2();

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

  for (om_number = 712; om_number < 713; om_number++) {
    TH2modele = spectre_charge_full_template(om_number);
    for (compteur = 1; compteur < 60; compteur++) {
      TH1D* spectre_om = NULL;
      TFile *file = new TFile(Form("histo_brut/716/time_file_%d.root", compteur), "READ");
      // TFile *file = new TFile("histo_brut/histo_ref_716.root", "READ");
      gROOT->cd();
      spectre_om = (TH1D*)file->Get("spectre_om");
      std::cout << "bin content = " << spectre_om->GetBinContent(317) << " and error = " << spectre_om->GetBinError(317) << '\n';

      for (int gain_count = 1; gain_count <2001; gain_count++) {
        gain = (gain_bin_min + gain_bin_width*(gain_count-1));
        // gain_count = 1001;
        if (gain < 1.01 && gain > 0.99) {
          modele = TH2modele->ProjectionY("modele", gain_count, gain_count);

          for (int i = 0; i < 17; i++) {
            modele->SetBinContent(i,0);
            // modele->SetBinError(i, 0);
            spectre_om->SetBinContent(i,0);
            // spectre_om->SetBinError(i, 0);
          }

          roofitter(modele, spectre_om, om_number, rootab, run_number, gain, compteur);
          // return;
          Chi2 = rootab[0];
          // gain = rootab[1];
          // gain_error = rootab[2];
          modele->Reset();
          Result_tree.Fill();
        }
      }
      // return;

      delete spectre_om;
    }

  }
  TFile new_file(Form("root/time_variation_%d.root", run_number), "RECREATE");
  new_file.cd();
  Result_tree.Write();
  new_file.Close();

  return;

}


void test(/* arguments */) {
  TFile *file = new TFile("test/a/test.root", "READ");
  gROOT->cd();
  RooDataSet* test = (RooDataSet*)file->Get("spectre_set4;1");
  RooDataSet* test2 = (RooDataSet*)file->Get("spectre_set4;2");
  auto frame = charge_tree.frame(Title("Fit gain simu"));
  test->plotOn(frame, DataError(RooAbsData::SumW2), DrawOption("P"),  LineColor(kGreen));
  test2->plotOn(frame, DataError(RooAbsData::SumW2), DrawOption("P"),LineColor(kBlue));

  TFile *file2 = new TFile("test/a/roodataset_run_716_gain_0.900-0.910.root", "READ");
  RooDataSet* test3 = (RooDataSet*)file->Get("spectre_set4");
  test3->plotOn(frame, DataError(RooAbsData::SumW2), DrawOption("P"));
  frame->Draw();
}

void RooDataSetter(TTree *tree, double gain_min, double gain_max, double pas, int om_number, int run_number){

  TStopwatch twach;
  twach.Start();
  charge_tree.setBins(1024);
  double charge, time;
  int om;
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("charge_tree",1);
  tree->SetBranchAddress("charge_tree", &charge);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);

  for (int j = 0; j < (gain_max-gain_min)/pas; j++) {
   SetTable[j][om_number-712] = new RooDataSet(TString("spectre_set")+=j, TString("spectre_set")+=j, RooArgSet(charge_tree)); //, Import("charge_tree"), Cut("charge_tree < 200000 && charge_tree > 3500 && om == 800 && time < 3600"));
  }

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (i%1000000 == 0) {
      std::cout << "/* entry = " << i << '\n';

    }
    for (int j = 0; j < (gain_max-gain_min)/pas; j++) {
      if (charge < 200000 && charge > 3500 && om == 800) {
        charge_tree = charge*(gain_min+j*pas);
        SetTable[j][om_number-712]->add(RooArgSet(charge_tree));
      }
    }
  }
  TFile new_file(Form("roodataset/roodataset_run_%d_gain_%.3f-%.3f.root", run_number, gain_min, gain_max), "RECREATE");
  new_file.cd();
  for (int j = 0; j < (gain_max-gain_min)/pas; j++) {
    SetTable[j][om_number - 712]->Write();
  }
  new_file.Close();

  twach.Stop();
  twach.Print();
  return;
}

void Looper(string file, double pas, double gain_min, double gain_max, int om_number, int run_number){

  TFile tree_file(Form("histo_brut/%s.root", file.c_str()), "READ");
  gROOT->cd();
  TTree* tree = (TTree*)tree_file.Get("Result_tree");


  for (double i = gain_min; i < gain_max; i+=pas*10) {
    RooDataSetter(tree, i, i+pas*10, pas, om_number, run_number);
    std::cout << "gain =" << i << '\n';
  }
}

double roofitter_dataset(TH1D* modele, int om_number, double *rootab, int run_number, double gain, RooDataSet *spectre_set)
{
  using namespace RooFit;
  charge_tree.setBins(1024);
  RooRealVar x("x", "x", 0, 200000);
  x.setBins(1024);
  charge_tree.setBins(1024);
  RooDataHist data("data", "data", charge_tree, Import(*modele));
  RooRealVar modele_rrv("modele_rrv", "data", 1, 0, 10000);
  RooHistPdf modele_pdf ("modele_pdf", "", charge_tree, data);
  // std::cout << "/* message */" << '\n';
  // auto aframe = charge_tree.frame(Title("Fit gain simu"));
  //
  // spectre_set->plotOn(aframe, DataError(RooAbsData::SumW2), DrawOption("P"));
  //
  // modele_pdf.plotOn(aframe, FillColor(0));
  // aframe->Draw();
  // return *rootab;



  RooNLLVar RooLogL("LogL", "LogL", modele_pdf, *spectre_set, Range(3500, 50000));   // create the variance
  RooMinimizer *miniChi = new RooMinimizer(RooLogL);

  // double ConvChi = miniChi->minimize("Minuit", "Migrad");  //   Create the LogL/LogL
  // for (double i = 3515.; i < 3516; i+=0.0001) {
  //   charge_tree.setVal(i);
  //   std::cout << "i = " << i << " and " <<  modele_pdf.getVal() << '\n';
  // }

  double LogL = RooLogL.getVal();
  // double LogL = RooLogL.getVal()/(1024. - 1);
  std::cout << "LogL = " << RooLogL.getVal() << '\n';
  TCanvas* can = new TCanvas;
  can->cd();
  auto frame = charge_tree.frame(Title("Fit gain simu"));

  spectre_set->plotOn(frame, DataError(RooAbsData::SumW2), DrawOption("P"));

  // modele_pdf.plotOn(frame, FillColor(0));
  spectre_set->plotOn(frame, DataError(RooAbsData::SumW2), DrawOption("P"));

  // Plot model components

  modele_pdf.plotOn(frame, LineColor(kRed), Name("sum_curve"), Range(3500, 200000), Precision(1e-5));

  double data_int = charge_tree.getVal();
  double data_int_error = charge_tree.getError();
  frame->GetYaxis()->SetTitle("n events");
  frame->GetXaxis()->SetTitle("charge (u.a)");
  frame ->Draw();
  can->SetLogy();
  frame->GetYaxis()->SetRangeUser(0.01, 5e6);
  // can->SaveAs("test.root");

  rootab[1] = data_int;
  rootab[2] = data_int_error;
  rootab[0] = LogL;
  TLatex l;
  l.SetTextFont(40);
  l.DrawLatex(90000, 80, Form("Khi2/NDF = %.2f", LogL));
  // return *rootab;
  can->SaveAs(Form("fit/om_%d_loop/fit_om_%d_gain_%f.png", om_number, om_number, gain));

  delete miniChi;
  delete can;
  delete frame;
  // delete miniLog;
  return *rootab;

}

void Dataset_Fit_Ref(int run_number) {
  Load_spectre(run_number);
  TH1::SetDefaultSumw2();

  double time, charge_tree, amplitude_tree, Chi2, gain, gain_error;
  int om_number, compteur;
  double* rootab = new double[3];

  TH1D* modele = NULL;
  TH2F* TH2modele = NULL;
  TFile tree_file(Form("histo_brut/Li_system_%d.root", run_number), "READ");
  gROOT->cd();
  TTree* tree = (TTree*)tree_file.Get("Result_tree");

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("gain_error", &gain_error);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("nsrun", &compteur);

  double pas = (1.1-0.9)/10000;
  // RooDataSet *SetTable = NULL;

  TFile dataset_file("roodataset/full.root", "READ");
  // TFile dataset_file("roodataset/roodataset_run_716_gain_1.000-1.010.root");
  gROOT->cd();

  for (om_number = 712; om_number < 713; om_number++) {
    TH2modele = spectre_charge_full_template(om_number);
    // RooDataSetter(tree, 0.9, 1.1, 0.001, om_number);
    // for (compteur = 1; compteur < 60; compteur++) {

    int count = 0;
    int count2 = 1;
      for (int gain_count = 1; gain_count <201; gain_count++) {
        gain = (gain_bin_min + gain_bin_width*(gain_count-1));
        // if (gain < 1.1 && gain > 0.9) {

        if ((int)(gain*100)%10 == 0) {
          count++;
        }
        // count = 5;
        // count2 = 1;

        RooDataSet* spectre_set = (RooDataSet*)dataset_file.Get(Form("spectre_set%d;%d", count, count2));
        if (count2 == 20) {
          count2=0;
        }


          modele = TH2modele->ProjectionY("modele", gain_count, gain_count);


          for (int i = 0; i < 16; i++) {
            modele->SetBinContent(i,0);
            modele->SetBinError(i, 0);
          }

          roofitter_dataset(modele, om_number, rootab, run_number, gain, spectre_set);
          // return;
          Chi2 = rootab[0];
          // gain = rootab[1];
          // gain_error = rootab[2];
          modele->Reset();
          Result_tree.Fill();
        // }
      }

      // file->Close();

    // }

  }
  TFile new_file(Form("root/test_%d.root", run_number), "RECREATE");
  new_file.cd();
  Result_tree.Write();
  new_file.Close();

  return;

}





std::vector<int> vtime = {1655219940, 1656339281, 1656427367, 1656513370, 1656596183, 1656689629, 1656772061, 1656856728};
std::vector<int> vrun_number = {736, 737, 738, 739, 740, 741, 742};

void file_merger(std::vector<int> run_number, std::vector<int> time, string previous_file_s = "") {
  TFile file(Form("root/Fit_Ref_%d-%d.root", run_number.at(0), run_number.at(run_number.size()-1)), "RECREATE");

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
    TFile tree_file(Form("root/Fit_Ref_%d.root", run_number.at(i)), "READ");

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

  TFile file(Form("root/TGraph_%s.root", file_name.c_str()), "RECREATE");

  TFile tree_file(Form("root/%s.root", file_name.c_str()), "READ");
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
      yaxis_error[j] = gain_error;
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
    std::cout << "Code start running" << '\n';

    for (int i = 0; i < n_run; i++) {
      Fit_Ref(run_number.at(i));
    }
    std::cout << "Fit_Ref ok" << '\n';

    file_merger(run_number, time);
    std::cout << "file_merger ok" << '\n';

    TGrapher(Form("Fit_ref_%d-%d.root", run_number.at(0), run_number.at(run_number.size())), n_run);
  }

  return 0;
}
