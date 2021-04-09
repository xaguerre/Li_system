#include <math.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TTree.h"
#include <fstream>
#include <string>
#include <cstring>
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TH2F.h"
#include "TApplication.h"
#include "TMultiGraph.h"
#include "TFeldmanCousins.h"
#include "TGaxis.h"
#include "TLeaf.h"
using namespace std;


TH1D* spectre_charge(int om_number, int file_number )
{
    TFile *file = new TFile(Form("histo_brut/histo_charge_amplitude_energie_%i.root", file_number), "READ");
    gROOT->cd();
    TH2F* charge = (TH2F*)file->Get("histo_pm_charge");

    TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
    file->Close();
    return spectre_charge;
}


void fit_all_om_charge()
  {
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TFile file("Resultats_root/Resultat_charge.root","RECREATE");

    int i_om;
    double n_evt;
    double mean_charge;
    double mean_error;
    double sigma;
    double  nbg;
    double lambda;
    int run_number;
    double comp;
    double error;

    TTree Result_tree("Result_tree","");
    Result_tree.Branch("run_number", &run_number);
    Result_tree.Branch("i_om", &i_om);
    Result_tree.Branch("n_evt", &n_evt);
    Result_tree.Branch("mean_charge", &mean_charge);
    Result_tree.Branch("mean_error", &mean_error);
    Result_tree.Branch("sigma", &sigma);
    Result_tree.Branch("nbg", &nbg);
    int tab[4] = {547, 551, 554, 557};

    TCanvas* canvas = new TCanvas;
    canvas->SetLogy();

    TH2F mean_charge_map("histo_om_mean_charge_map", "mean_charge_map", 20, 0, 20, 13, 0, 13);
    for (int i = 0; i < 4; i++) {
      run_number = tab[i];

      for(int om = 1; om<259; om+=2577)
      {
        TH1D* spectre_om = spectre_charge(om, tab[i]);
        spectre_om->Draw();
        TF1* f_ComptonEdgePoly = new TF1 ("f_ComptonEdgePoly","[0]/2.0*(1+TMath::Erf(([1]-x)/(TMath::Sqrt(2)*[2])))+[3]*x", 40000, 90000);
        f_ComptonEdgePoly->SetParNames("N_evt","mean_charge","Sigma","Nbg" );

        if (om == 1)        //om multiple de (13)-1
          {
            f_ComptonEdgePoly->SetParameters(120, 72367, 6893, 3.91e-5);
            f_ComptonEdgePoly->SetRange(45000,80000);
            f_ComptonEdgePoly->Draw("same");
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
            f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
            f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          }
        if (om == 258)        //om multiple de (13)-1
          {
            f_ComptonEdgePoly->SetParameters(120, 72367, 6893, 3.91e-5);
            f_ComptonEdgePoly->SetRange(50000,80000);
            f_ComptonEdgePoly->Draw("same");
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
            f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
            f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          }
          if (om == 8)        //om multiple de (13)-1
            {
              f_ComptonEdgePoly->SetParameters(120, 72367, 6893, 3.91e-5);
              f_ComptonEdgePoly->SetRange(45000,80000);
              f_ComptonEdgePoly->Draw("same");
              spectre_om->Fit(f_ComptonEdgePoly, "RQ");
              f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
              spectre_om->Fit(f_ComptonEdgePoly, "RQ");
              f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
              spectre_om->Fit(f_ComptonEdgePoly, "RQ");
            }
          // if ((om % 13) == 12 )        //om multiple de (13)-1
          // {
          //   f_ComptonEdgePoly->SetParameters(120, 72367, 6893, 3.91e-5);
          //   f_ComptonEdgePoly->SetRange(6000,100000);
          //   f_ComptonEdgePoly->Draw("same");
          //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          //   f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
          //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          //   f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
          //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          // }
          // else if ((om % 13) == 0)       //om multiple de 13
          // {
          //   f_ComptonEdgePoly->SetParameters(112, 68168, 5604, 1.2e-05);
          //   f_ComptonEdgePoly->SetRange(50000,100000);
          //   f_ComptonEdgePoly->Draw("same");
          //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          //   f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
          //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          //   f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-1.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
          //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          // }
          // else         //om normaux (8pouces)
          // {
          //   f_ComptonEdgePoly->SetParameters(111, 60978, 3787, 4.19e-05);
          //   f_ComptonEdgePoly->SetRange(55000,100000);
          //   f_ComptonEdgePoly->Draw("same");
          //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          //   f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+7.5*f_ComptonEdgePoly->GetParameter(2));
          //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          //   f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-3.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+7.5*f_ComptonEdgePoly->GetParameter(2));
          //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          // }


            i_om = om;
            n_evt = (f_ComptonEdgePoly->GetParameter(0));
            mean_charge = (f_ComptonEdgePoly->GetParameter(1));
            sigma = (f_ComptonEdgePoly->GetParameter(2));
            nbg = (f_ComptonEdgePoly->GetParameter(3));
            mean_error = f_ComptonEdgePoly->GetParError(1) ;
            //mapping
            int om_col = (om % 13 );
            int om_row = (om / 13);
            mean_charge_map.SetBinContent( om_row+1, om_col+1, mean_charge);
            Result_tree.Fill();

            canvas->SaveAs(Form("fit/fit_Tl/charge_fit_om_%03d_run_%d.png", om, run_number));

            delete spectre_om;
            delete f_ComptonEdgePoly;
        }
      }

  TH1F comp_map("comp_map", "comp_map", 100, 0, 100);

  Result_tree.GetEntry(1);
  comp = mean_charge;
  double error1 = mean_error;
  std::cout << comp << '\n';
  for (int i = 0; i < Result_tree.GetEntries(); i++) {
    Result_tree.GetEntry(i);
    if (i_om == 1) {
      comp_map.SetBinContent(i*19, mean_charge/comp*1.0);
      error = mean_error/comp*1.0 + (mean_charge/(comp*comp))*1.0*error1;
      std::cout << error << '\n';
      comp_map.SetBinError(i*19,error);
      comp_map.GetYaxis()->SetRangeUser(0.9, 1.1);
    }
  }

  file.cd();
  comp_map.Write();
  mean_charge_map.Write();
  Result_tree.Write();

  file.Close();

  return;
}

void fit_LI()
  {
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TFile file("Resultats_root/Resultat_charge_Li.root","RECREATE");

    int i_om;
    double constante;
    double mean;
    double mean_error;
    double sigma;
    int run_number;

    double error;

    double count = -1;

    TTree Result_tree("Result_tree","");
    Result_tree.Branch("run_number", &run_number);
    Result_tree.Branch("i_om", &i_om);
    Result_tree.Branch("constante", &constante);
    Result_tree.Branch("mean", &mean);
    Result_tree.Branch("mean_error", &mean_error);
    Result_tree.Branch("sigma", &sigma);
    Result_tree.Branch("count", &count);

    int tab[4] = {549, 552, 555, 559};
    double debut[4] = {50e9, 50e9, 55e9, 55e9};

    TCanvas* canvas = new TCanvas;
    canvas->SetLogy();

    TH1F comp_map("comp_map", "comp_map", 100, 0, 100);
    TH2F mean_charge_map("histo_om_mean_charge_map", "mean_charge_map", 20, 0, 20, 13, 0, 13);
    for (int i = 0; i < 4; i++) {
      run_number = tab[i];
      TFile tree_file(Form("histo_brut/histo_Li_system_%d.root", run_number),"READ");
      count++;
      int om_number;
      double time;
      double charge_tree;
      TTree* tree = (TTree*)tree_file.Get("Result_tree");
      tree->SetBranchStatus("*",0);
      tree->SetBranchStatus("om_number",1);
      tree->SetBranchAddress("om_number", &om_number);
      tree->SetBranchStatus("time",1);
      tree->SetBranchAddress("time", &time);
      tree->SetBranchStatus("charge_tree",1);
      tree->SetBranchAddress("charge_tree", &charge_tree);

      for(int om = 1; om<259; om+=259)
      {
        // for (double j = debut[i]; j < 250e9; j+=38e9)
        for (double j = debut[i]; j < debut[i]+6*38e9; j+=38e9)
        {
          TF1* f_ComptonEdgePoly = new TF1 ("f_ComptonEdgePoly","[0]/2.0*(1+TMath::Erf(([1]-x)/(TMath::Sqrt(2)*[2])))+[3]*x", 40000, 90000);

          tree->Draw("charge_tree", Form("om_number == 1 && time < %f && time > %f", j, j-40e9));
          TH1D *spectre = (TH1D*)gPad->GetPrimitive("htemp");
          TF1 *f_Gaus = new TF1("f1", "[0]*TMath::Gaus(x,[1],[2])");
          f_Gaus->SetParNames("N_evt","mean_charge","Sigma");

          if (om == 1)        //om multiple de (13)-1
          {
            f_Gaus->SetParameters(25, spectre->GetMean(1), 2100);
            f_Gaus->SetRange(0,300e9);
            f_Gaus->Draw("same");
            spectre->Fit(f_Gaus, "RQ");
            f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
            spectre->Fit(f_Gaus, "RQ");
            f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
            spectre->Fit(f_Gaus, "RQ");
          }
          // if (om == 258)        //om multiple de (13)-1
          // {
          //   f_Gaus->SetParameters(355, spectre->GetMean(1), 1400);
          //   f_Gaus->SetRange(0,200000);
          //   f_Gaus->Draw("same");
          //   spectre->Fit(f_Gaus, "RQ");
          //   f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
          //   spectre->Fit(f_Gaus, "RQ");
          //   f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
          //   spectre->Fit(f_Gaus, "RQ");
          // };

          i_om = om;
          constante = (f_Gaus->GetParameter(0));
          mean = (f_Gaus->GetParameter(1));
          sigma = (f_Gaus->GetParameter(2));
          mean_error = f_Gaus->GetParError(1) ;
          //mapping
          int om_col = (om % 13 );
          int om_row = (om / 13);
          mean_charge_map.SetBinContent( om_row+1, om_col+1, mean);

          Result_tree.Fill();
          canvas->SaveAs(Form("fit/fit_Li/charge_Li_om_%03d_run_%i_time_%f.png", om, tab[i], j));

          delete spectre;
          delete f_Gaus;
        }
      }
    }
    TH1F comp_map_0("comp_map_1", "comp_map_1", 100, 0, 100);
    TH1F comp_map_1("comp_map_2", "comp_map_2", 100, 0, 100);
    TH1F comp_map_2("comp_map_3", "comp_map_3", 100, 0, 100);
    TH1F comp_map_3("comp_map_4", "comp_map_4", 100, 0, 100);
    TH1F comp_map_4("comp_map_5", "comp_map_5", 100, 0, 100);
    TH1F comp_map_5("comp_map_6", "comp_map_6", 100, 0, 100);
    double comp[6];
    double error1[6];
    for (int i = 0; i < 6; i++) {
      Result_tree.GetEntry(i);
      comp[i] = mean;
      error1[i] = mean_error;
      std::cout << comp[i] << '\n';
    }
    for (int i = 0; i < Result_tree.GetEntries(); i++) {
      Result_tree.GetEntry(i);
      if (i_om == 1) {
        if ((i%6) == 0) {
          comp_map_0.SetBinContent(i*19/6, mean/comp[0]*1.0);
          std::cout << mean/comp[0]*1.0 << '\n';
          error = mean_error/comp[0]*1.0 + (mean/(comp[0]*comp[0]))*1.0*error1[0];
          comp_map_0.SetBinError(i*19/6,error);
          comp_map_0.GetYaxis()->SetRangeUser(0.9, 1.1);
        }
        else if((i%6) == 1) {
          comp_map_1.SetBinContent((i-1)*19/6, mean/comp[1]*1.0);
          std::cout << mean/comp[1]*1.0 << '\n';
          error = mean_error/comp[1]*1.0 + (mean/(comp[1]*comp[1]))*1.0*error1[1];
          comp_map_1.SetBinError((i-1)*19/6,error);
          comp_map_1.GetYaxis()->SetRangeUser(0.9, 1.1);
        }
        else if((i%6) == 2) {
          comp_map_2.SetBinContent((i-2)*19/6, mean/comp[2]*1.0);
          std::cout << mean/comp[2]*1.0 << '\n';
          error = mean_error/comp[2]*1.0 + (mean/(comp[2]*comp[2]))*1.0*error1[2];
          comp_map_2.SetBinError((i-2)*19/6,error);
          comp_map_2.GetYaxis()->SetRangeUser(0.9, 1.1);
        }
        else if((i%6) == 3) {
          comp_map_3.SetBinContent((i-3)*19/6, mean/comp[3]*1.0);
          std::cout << mean/comp[3]*1.0 << '\n';
          error = mean_error/comp[3]*1.0 + (mean/(comp[3]*comp[3]))*1.0*error1[3];
          comp_map_3.SetBinError((i-3)*19/6,error);
          comp_map_3.GetYaxis()->SetRangeUser(0.9, 1.1);
        }
        else if((i%6) == 4) {
          comp_map_4.SetBinContent((i-4)*19/6, mean/comp[4]*1.0);
          std::cout << mean/comp[4]*1.0 << '\n';
          error = mean_error/comp[4]*1.0 + (mean/(comp[4]*comp[4]))*1.0*error1[4];
          comp_map_4.SetBinError((i-4)*19/6,error);
          comp_map_4.GetYaxis()->SetRangeUser(0.9, 1.1);
        }
        else if((i%6) == 5) {
          comp_map_5.SetBinContent((i-5)*19/6, mean/comp[5]*1.0);
          std::cout << mean/comp[5]*1.0 << '\n';
          error = mean_error/comp[5]*1.0 + (mean/(comp[5]*comp[5]))*1.0*error1[5];
          comp_map_5.SetBinError((i-5)*19/6,error);
          comp_map_5.GetYaxis()->SetRangeUser(0.9, 1.1);
        }
      }
    }

    file.cd();
    mean_charge_map.Write();
    comp_map_0.Write();
    comp_map_1.Write();
    comp_map_2.Write();
    comp_map_3.Write();
    comp_map_4.Write();
    comp_map_5.Write();
    Result_tree.Write();
    file.Close();
    return;
}
