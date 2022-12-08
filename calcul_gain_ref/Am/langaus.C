#include "TStyle.h"
#include "TMath.h"
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

TH2F* charge_spectre = NULL;
TH1D* charge_spectre_template = NULL;

void Load_spectre(int run_number){
  TFile *file = new TFile(Form("../histo_brut/histo_ref_%d.root", run_number), "READ");
  gROOT->cd();
  charge_spectre = (TH2F*)file->Get("histo_pm_charge_50");
  return;
}

TH1D* spectre_charge(int om_number){
  TH1D* spectre_charge = charge_spectre->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  // spectre_charge->Rebin(4);
  return spectre_charge;
}

double langau_func (double *x, double *par)
{
   const double LANDAU_MPV = par[1];
   const double LANDAU_WIDTH = par[0];
   const double TOT_INT = par[2];
   const double GAUS_WIDTH = par[3];
   const double EXP_CONST = par[4];
   const double EXP_LAMBDA = par[5];
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      double np = 100.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;


      // MP shift correction
      mpc = LANDAU_MPV - mpshift * LANDAU_WIDTH;

      // Range of convolution integral
      xlow = x[0] - sc * GAUS_WIDTH;
      xupp = x[0] + sc * GAUS_WIDTH;

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(mpc - xx,0,LANDAU_WIDTH) / LANDAU_WIDTH ;
         sum += fland * TMath::Gaus(x[0],xx,GAUS_WIDTH);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(mpc - xx,0,LANDAU_WIDTH) / LANDAU_WIDTH ;
         sum += fland * TMath::Gaus(x[0],xx,GAUS_WIDTH);
      }

      return (TOT_INT * step * sum * invsq2pi / GAUS_WIDTH) + EXP_CONST*TMath::Exp(-x[0]/EXP_LAMBDA) ;
}

void langaus(int run_number) {
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");

  Load_spectre(run_number);

  TFile file(Form("fit/fit_Am/Am_fit_%d.root", run_number), "RECREATE");

  double mean, mean_error;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_error", &mean_error);
  gROOT->cd();

  TH1D* spectre_om = NULL;
  spectre_om = spectre_charge(713);

  TCanvas* canvas = new TCanvas;
  TF1 *f_langaus = new TF1 ("f_langaus", langau_func, 6000, 10000, 6);
  // f_langaus->SetParameters(240, 8700, 5.4e8, 600, 3.3e5, 6300);
  f_langaus->SetParameters(240, 8700, 1.5e7, 600, 10000, 6300);

  // spectre_om->Draw();
  // f_langaus->Draw();
  // return;


  spectre_om->Fit("f_langaus", "R");
  f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
  spectre_om->Fit("f_langaus", "R");
  f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
  spectre_om->Fit("f_langaus", "R");
  f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));

  spectre_om->Draw();
  f_langaus->Draw("lsame");
  TF1 *f_exp = new TF1 ("f_exp", "[0]*TMath::Exp(-x/[1])", 0, 20000);
  f_exp->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
  f_exp->SetParameters(f_langaus->GetParameter(4), f_langaus->GetParameter(5));
  f_exp->Draw("same");
  f_exp->SetLineColor(kRed);
  f_exp->SetLineStyle(kDashed);
  mean = f_langaus->GetParameter(1);
  mean_error = f_langaus->GetParError(1);
  Result_tree.Fill();

  canvas->SaveAs(Form("png/fit_Am/fit_run_%d.png", run_number));
  // delete canvas;
  // delete f_langaus;
  // delete f_exp;
  // delete spectre_om;

  file.cd();
  Result_tree.Write();
  file.Close();

}
