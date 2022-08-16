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

TH2F* charge_spectre = NULL;
TH2F* charge_spectre_template = NULL;

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

void Load_spectre_template(){
  TFile *file = new TFile("histo_brut/histo_ref_716.root", "READ");
  gROOT->cd();
  charge_spectre_template = (TH2F*)file->Get("histo_pm_charge");
  return;
}

TH1D* spectre_charge_full_template(int om_number){
  TH1D* spectre_charge = charge_spectre_template->ProjectionY(Form("charge_modele%03d",om_number), om_number+1, om_number+1);
  // spectre_charge->Rebin(4);
  return spectre_charge;

}

double roofitter(TH1D* modele, TH1D* spectre_om, int om_number, double *rootab, int run_number)
{
  using namespace RooFit;
  RooRealVar* x = new RooRealVar("x","x",0,2e5);
  RooDataHist spectre_data("spectre_data", "spectre_data", *x, Import(*spectre_om));
  x->setBins(1024);

  Int_t nbins(1024);
  TArrayD limits(nbins+1);

  RooRealVar *gain=new RooRealVar("gain","gain",1,0.5,1.5);
  gain->setBins(10000);
  RooRealVar *binValue[nbins];
  RooRealVar *binLimValue[nbins+1];
  RooFormulaVar* binLim[nbins+1];//0 = new RooRealVar("binHeight0","bin 0 Value",0.1,0.0,1.0);

  RooArgList* list = new RooArgList("list");
  RooArgList* listBin = new RooArgList("listBin");
  gain->setVal(0.1);
  for(int i=0 ;i<nbins+1;i++)
  {
    binLimValue[i]= new RooRealVar(TString("binLimValue")+=i,TString("binLimValue")+=i, modele->GetBinLowEdge(i));

    binLim[i]= new RooFormulaVar(TString("binLim")+=i,TString("binLim")+=i,"x[0]*x[1]",RooArgList(*binLimValue[i],*gain));
    listBin->add(*binLim[i]);
    if (i==nbins){continue;}
    binValue[i]= new RooRealVar(TString("binValue")+=i,TString("binValue")+=i, modele->GetBinContent(i));
    list->add(*binValue[i]); // one less bin than limits
  }

  RooParametricBinningStepFunction  *aPdf = new RooParametricBinningStepFunction("aPdf", "PSF", *x, *list, *listBin, nbins);
  RooPlot* xFrame= x->frame();

  // cout<<"val="<<aPdf->getVal()<<endl;
  spectre_data.plotOn(xFrame, MarkerSize(0.1), DataError(RooAbsData::SumW2), DrawOption("P"));

  // aPdf->plotOn(xFrame,LineColor(kGreen+2));

  RooChi2Var RooChi2("Chi2", "Chi2", *aPdf, spectre_data);
  RooMinuit miniChi(RooChi2);
  miniChi.migrad();
  miniChi.hesse();
  RooFitResult* R;
  RooFitResult* R2;
  R = miniChi.save();
  R2 = miniChi.save();
  for (int i = 0; i < 20; i++) {
    gain->randomize();
    miniChi.migrad();
    miniChi.hesse();
    R2 = miniChi.save();
    std::cout << "R2 minNLL = " << R2->minNll() << '\n';
    if (R2->minNll() < R->minNll() ) {
      R = R2;
    }
  }

  miniChi.hesse();
  auto* a1= new TCanvas();
  TH1D* Chisto = (TH1D*)RooChi2.createHistogram("Chisto", *gain);
  Chisto->Sumw2(kFALSE);
  Chisto->Draw();
  a1->SaveAs("test.root");

  auto* b = new TCanvas();
  RooPlot* gainFrame= gain->frame();
  RooChi2.plotOn(gainFrame);
  gainFrame->Draw();

  b->SaveAs("Chi2.root");

  double Chi2 = RooChi2.getVal()/(1024 - 1);
  std::cout << "Chi2 = " << Chi2 << " or Chi2 = " << Chisto->GetMinimum() << '\n';
  // std::cout << "gain = " << gain->getVal() << " +- " << gain->getError() << '\n';
  std::cout << "gain = " << ((RooRealVar*)R->floatParsFinal().find(*gain))->getVal() << " +- " <<  ((RooRealVar*)R->floatParsFinal().find(*gain))->getPropagatedError(*R) << '\n';

  gain->setVal(((RooRealVar*)R->floatParsFinal().find(*gain))->getVal());
  aPdf->plotOn(xFrame,LineColor(kRed+2));
  auto* c1= new TCanvas();
  xFrame->Draw();
  xFrame->GetYaxis()->SetRangeUser(0.001,2e6);
  c1->SetLogy();
  c1->SaveAs(Form("fit/om_%d/run_%d.root", om_number, run_number));

  rootab[0] = Chisto->GetMinimum();
  rootab[1] = ((RooRealVar*)R->floatParsFinal().find(*gain))->getVal();
  rootab[2] = 0.0013;
  return *rootab;
  delete aPdf;
  delete Chisto;
  delete xFrame;
  // delete spectre_data;
  delete c1;
  delete a1;
  delete b;
  return *rootab;

}

void Fit_Ref(int run_number) {
  Load_spectre(run_number);
  Load_spectre_template();
  TH1::SetDefaultSumw2();

  double time, charge_tree, amplitude_tree, Chi2, gain, gain_error;
  int om_number;
  double* rootab = new double[3];

  TH1D* modele = NULL;

  // TFile tree_file(Form("histo_brut/Li_system_%d.root", run_number), "READ");
  // gROOT->cd();
  // TTree* tree = (TTree*)tree_file.Get("Result_tree");
  // tree->SetBranchStatus("*",0);
  // tree->SetBranchStatus("om_number",1);
  // tree->SetBranchAddress("om_number", &om_number);
  // tree->SetBranchStatus("time",1);
  // tree->SetBranchAddress("time", &time);
  // tree->SetBranchStatus("charge_tree",1);
  // tree->SetBranchAddress("charge_tree", &charge_tree);
  // tree->SetBranchStatus("amplitude_tree",1);
  // tree->SetBranchAddress("amplitude_tree", &amplitude_tree);

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("gain_error", &gain_error);
  Result_tree.Branch("run_number", &run_number);

  for (om_number = 712; om_number < 713; om_number++) {

    modele = spectre_charge_full_template(om_number);
    // TH1D *spectre_om = new TH1D ("spectre_om", "", 1024, 0, 200000);
    // tree->Project("spectre_om", "charge_tree", Form("om_number == %d && time < 300 && Entry$ > 120e6 ", om_number+88));
    TH1D* spectre_om = NULL;
    spectre_om = spectre_charge_full(om_number);

    for (int i = 0; i < 50; i++) {
      modele->SetBinContent(i,0);
      modele->SetBinError(i, 0);
      spectre_om->SetBinContent(i,0);
      spectre_om->SetBinError(i, 0);
    }

    roofitter(modele, spectre_om, om_number, rootab, run_number);
    return;
    Chi2 = rootab[0];
    gain = rootab[1];
    gain_error = rootab[2];



    modele->Reset();
    Result_tree.Fill();
    delete spectre_om;

  }
  // TFile new_file(Form("root/Fit_Ref_%d.root", run_number), "RECREATE");
  // new_file.cd();
  // Result_tree.Write();
  // new_file.Close();

  // return;

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
