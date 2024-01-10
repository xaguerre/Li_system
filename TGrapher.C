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

void TGrapher(std::string file_name, std::string path, int n_run) {
  std::cout << "ok" << '\n';
  TFile file(Form("%s/TGraph/TGraph_%s.root", path.c_str(), file_name.c_str()), "RECREATE");

  TFile tree_file(Form("%s/Merged_Fit/%s.root", path.c_str(), file_name.c_str()), "READ");
  double time;
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
  tree->SetBranchAddress("time", &time);

  double yaxis[n_run];
  double yaxis_error_moins[n_run];
  double yaxis_error_plus[n_run];
  double xaxis[n_run];
  double xaxis_error_moins[n_run];
  double xaxis_error_plus[n_run];
  auto canvas = new TCanvas;

  file.cd();
  int compteur = 0;
  for (int j = 0; j < n_run; j++){
    tree->GetEntry(i+j*5);
    compteur++;
    yaxis[j] = gain;
    yaxis_error_moins[j] = gain_error_moins;
    yaxis_error_plus[j] = gain_error_plus;
    xaxis[j] = time;
    xaxis_error_plus[j] = 0.00001;
    xaxis_error_moins[j] = 0.00001;
  }

  gain_graph = new TGraphAsymmErrors(n_run, xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);

  gain_graph->SetName(Form("fit_OM_ref_%d", om_number));
  gain_graph->SetNameTitle(Form("fit_OM_ref_%d", om_number), Form("Gain evolution of the OM %d with regard to the background", om_number));
  gain_graph->GetXaxis()->SetTimeDisplay(1);
  gain_graph->GetXaxis()->SetTitle("Time");
  gain_graph->GetYaxis()->SetTitle("Gain evolution");
  gain_graph->GetYaxis()->SetRangeUser(0.9, 1.1);
  gain_graph->SetMarkerColor(2);
  gain_graph->SetMarkerStyle(5);
  gain_graph->SetMarkerSize(2);
  gain_graph->Draw("AP");

  gain_graph->Write();

}

canvas->Write();
  file.Close();

}
