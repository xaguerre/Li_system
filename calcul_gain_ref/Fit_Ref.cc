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
#include <RooParametricBinningStepFunction.h>

TH2F* charge_spectre = NULL;

void Load_spectre(){
  TFile *file = new TFile("../histo_brut/histo_ref_714.root", "READ");
  gROOT->cd();
  charge_spectre = (TH2F*)file->Get("histo_pm_charge");
  return;
}

TH1D* spectre_charge_full(int om_number){
  TH1D* spectre_charge = charge_spectre->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  // spectre_charge->Rebin(4);
  return spectre_charge;

}


void roofitter(TH1D* modele, TH1D* spectre_om, int om_number)
{
	using namespace RooFit;
  RooRealVar* x = new RooRealVar("x","x",0,2e5);
  // modele->Scale(1./modele->Integral());
  // spectre_om->Scale(1./spectre_om->Integral());
  RooDataHist spectre_data("spectre_data", "spectre_data", *x, Import(*spectre_om));
  x->setBins(1024);

  Int_t nbins(1024);
	TArrayD limits(nbins+1);

	RooRealVar *gain=new RooRealVar("gain","gain",1,0.5,1.5);
  gain->setBins(1000);
	RooRealVar *binValue[nbins];
	RooRealVar *binLimValue[nbins+1];
	RooFormulaVar* binLim[nbins+1];//0 = new RooRealVar("binHeight0","bin 0 Value",0.1,0.0,1.0);

	RooArgList* list = new RooArgList("list");
	RooArgList* listBin = new RooArgList("listBin");
  gain->setVal(0.1);
	for(int i=0 ;i<nbins+1;i++)
	{
		// limits[i] = 0.0+i*0.5; //etc...
		binLimValue[i]= new RooRealVar(TString("binLimValue")+=i,TString("binLimValue")+=i, modele->GetBinLowEdge(i));

		binLim[i]= new RooFormulaVar(TString("binLim")+=i,TString("binLim")+=i,"x[0]*x[1]",RooArgList(*binLimValue[i],*gain));
		listBin->add(*binLim[i]);
		if (i==nbins){continue;}
		binValue[i]= new RooRealVar(TString("binValue")+=i,TString("binValue")+=i, modele->GetBinContent(i));
		list->add(*binValue[i]); // one less bin than limits
	}

	RooParametricBinningStepFunction  *aPdf = new RooParametricBinningStepFunction("aPdf", "PSF", *x, *list, *listBin, nbins);
  RooPlot* xFrame= x->frame();

	cout<<"val="<<aPdf->getVal()<<endl;
  spectre_data.plotOn(xFrame, MarkerSize(0.1), DataError(RooAbsData::SumW2), DrawOption("P"));

  aPdf->plotOn(xFrame,LineColor(kGreen+2));

  RooChi2Var RooChi2("Chi2", "Chi2", *aPdf, spectre_data);
  RooMinuit miniChi(RooChi2);
  // RooMinimizer *miniChi = new RooMinimizer(RooChi2);
  miniChi.migrad();
  miniChi.hesse();
  RooFitResult* R;
  RooFitResult* R2;
  R = miniChi.save();
  R2 = miniChi.save();
  for (int i = 0; i < 10; i++) {
    gain->randomize();
    // miniChi->minimize("Minuit", "Migrad");  //   Create the Chi2/LogL
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



  Chisto->Draw();
  // Chisto->GetXaxis()->SetRangeUser(0.9,1.1);
  // Chisto->GetYaxis()->SetRangeUser(0,3);
  // a1->SetLogy();
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
  // auto* c1=new TCanvas();
	// RooPlot* xFrame= x->frame();
  // spectre_data.plotOn(xFrame);
  // xFrame->Draw();
  // c1->SaveAs("test.png");
  //
  //
  //
  // auto* c2=new TCanvas();
	// RooPlot* xFrame2= x->frame();
  // xFrame2->Draw();

  // c2->SaveAs("test2.png");
  gain->setVal(((RooRealVar*)R->floatParsFinal().find(*gain))->getVal());
	aPdf->plotOn(xFrame,LineColor(kRed+2));
	auto* c1= new TCanvas();
	xFrame->Draw();
  xFrame->GetYaxis()->SetRangeUser(0.001,2e6);
  c1->SetLogy();
  c1->SaveAs(Form("fit/fit_ref_%d.png",om_number));

  // delete miniChi;
  delete xFrame;
  delete c1;
  return;

}



void Fit_Ref() {
  Load_spectre();
  TH1::SetDefaultSumw2();

  double time, charge_tree, amplitude_tree;
  int om_number;

  TFile tree_file("../histo_brut/Li_system_714.root", "READ");
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
  int n_evt = 0;
  TH1D* modele = NULL;

  // auto* can = new TCanvas();
  // can->SetLogy();
  for (int om = 712; om < 713; om++) {

    // TH1D *modele = new TH1D ("modele", "", 1024, 0, 200000);
    modele = spectre_charge_full(om);


    TH1D *spectre_om = new TH1D ("spectre_om", "", 1024, 0, 200000);
        // spectre_om = spectre_charge_full(om);
    // spectre_om->Draw();
    tree->Project("spectre_om", "charge_tree", Form("om_number == %d && time < 300 && Entry$ > 120e6 ", om+88));
    //
    for (int i = 0; i < 50; i++) {
      modele->SetBinContent(i,0);
      modele->SetBinError(i, 0);
      spectre_om->SetBinContent(i,0);
      spectre_om->SetBinError(i, 0);
    }
    // for (int i = 0; i > 50; i++) {
    //   modele->SetBinContent(i,0);
    //   spectre_om->SetBinContent(i,0);
    //   spectre_om->SetBinError(i, 0);
    //   modele->SetBinError(i, 0);
    // }

    roofitter(modele, spectre_om, om);

    modele->Reset();
    delete spectre_om;

  }


  delete modele;


}
