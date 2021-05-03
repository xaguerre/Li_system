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

void txt_to_root() {

  std::ifstream parametres("/home/aguerre/Bureau/Thèse/Li_system/Resultats_txt/Resultats_charge.txt");
  TFile file("test.root","RECREATE");

  int compteur = 0;
  int i_om;
  double n_evt;
  double mean_charge;
  double sigma;
  double nbg;
  double C = 0;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("i_om", &i_om);
  Result_tree.Branch("n_evt", &n_evt);
  Result_tree.Branch("mean_charge", &mean_charge);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("nbg", &nbg);
  Result_tree.Branch("C", &C);

  while(parametres >> i_om)
  {
    parametres >> n_evt >> mean_charge >> sigma >> nbg >> C;
    Result_tree.Fill();
    compteur++;
  }

file.cd();
Result_tree.Write();
file.Close();

}


void Test(){
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
    double error;
    double C = 0;

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
    for (int i = 0; i <1; i++) {
      run_number = tab[i];

      for(int om = 91; om<648; om+=666)
      {
        TH1D* spectre_om = spectre_charge(om, tab[i]);
        spectre_om->Draw();
        TF1* f_ComptonEdgePoly = new TF1 ("f_ComptonEdgePoly","[0]/2.0*(1+TMath::Erf(([1]-x)/(TMath::Sqrt(2)*[2])))+[3]*x + [4]", 40000, 90000);
        f_ComptonEdgePoly->SetParNames("N_evt","mean_charge","Sigma", "Nbg", "C" );
        f_ComptonEdgePoly->SetParLimits(0, 0, 200);
        f_ComptonEdgePoly->SetParLimits(1, 0, 200000);
        f_ComptonEdgePoly->SetParLimits(2, 0, 10000);
        f_ComptonEdgePoly->SetParLimits(3, -10, 0);
        f_ComptonEdgePoly->SetParLimits(4, 0, 20);

          if ((om % 13) == 12 )        //om multiple de (13)-1
          {
            f_ComptonEdgePoly->SetParameters(120, 60000, 6893, -3.91e-5, 11);
            f_ComptonEdgePoly->SetRange(60000-2.5*6893,90000+2.5*6893);
            f_ComptonEdgePoly->Draw("same");
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
            f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
            f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          }
          else if ((om % 13) == 0)       //om multiple de 13
          {
            f_ComptonEdgePoly->SetParameters(27, 63000, 4104, -2.2e-04, 5);
            f_ComptonEdgePoly->SetRange(63000-2.5*2500,68000+2.5*6893);
            f_ComptonEdgePoly->Draw("same");
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
            f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
            f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          }
          else         //om normaux (8pouces)
          {
            f_ComptonEdgePoly->SetParameters(111, 57000, 3787, -3.19e-05, 5);
            f_ComptonEdgePoly->SetRange(57000-2.5*3787,100000);
            f_ComptonEdgePoly->Draw("same");
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
            f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+7.5*f_ComptonEdgePoly->GetParameter(2));
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
            f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-3.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+7.5*f_ComptonEdgePoly->GetParameter(2));
            spectre_om->Fit(f_ComptonEdgePoly, "RQ");
          }

        i_om = om;
        n_evt = (f_ComptonEdgePoly->GetParameter(0));
        mean_charge = (f_ComptonEdgePoly->GetParameter(1));
        sigma = (f_ComptonEdgePoly->GetParameter(2));
        nbg = (f_ComptonEdgePoly->GetParameter(3));
        C = (f_ComptonEdgePoly->GetParameter(4));
        mean_error = f_ComptonEdgePoly->GetParError(1) ;
        //mapping
        int om_col = (om % 13 );
        int om_row = (om / 13);
        mean_charge_map.SetBinContent( om_row+1, om_col+1, mean_charge);
        Result_tree.Fill();

        canvas->SaveAs(Form("fit/fit_Tl/charge_fit_om_%03d_run_%d.png", om, run_number));

      }
    }



}

void fit_all_om_charge()
  {
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    std::ofstream outFile("Resultats_txt/Resultats_charge.txt");
    std::ifstream parametres("/home/aguerre/Bureau/Thèse/Li_system/Resultats_txt/Resultats_charge_run_547.txt");

    TFile file("Resultats_root/Resultat_charge.root","RECREATE");

    double chi = 0;
    int i_om;
    double n_evt;
    double mean_charge;
    double mean_error;
    double sigma;
    double  nbg;
    double lambda;
    int run_number;
    double error;
    double C = 0;

    TTree Result_tree("Result_tree","");
    Result_tree.Branch("run_number", &run_number);
    Result_tree.Branch("i_om", &i_om);
    Result_tree.Branch("n_evt", &n_evt);
    Result_tree.Branch("mean_charge", &mean_charge);
    Result_tree.Branch("mean_error", &mean_error);
    Result_tree.Branch("sigma", &sigma);
    Result_tree.Branch("C", &C);
    Result_tree.Branch("nbg", &nbg);
    Result_tree.Branch("chi", chi);
    int tab[4] = {547, 551, 554, 557};

    int run_tab[712];
    int om_tab[712];
    double N_evt_tab[712];
    double mean_tab[712];
    double sigma_tab[712];
    double Nbg_tab[712];
    double C_tab[712];
    int compteur = 0;

    while(parametres >> om_tab[compteur])  // tant que l'on peut mettre la ligne dans "contenu"
    {
      parametres >> N_evt_tab[compteur] >> mean_tab[compteur] >> sigma_tab[compteur] >> Nbg_tab[compteur] >> C_tab[compteur];
      compteur++;
    }

    TCanvas* canvas = new TCanvas;
    canvas->SetLogy();

    int k = 3;
    TH2F mean_charge_map("histo_om_mean_charge_map", "mean_charge_map", 20, 0, 20, 13, 0, 13);
      for(int om = 0; om < 520; om+=1)
      {
        int run_number = tab[k];
        TH1D* spectre_om = spectre_charge(om_tab[om], run_number);
        spectre_om->Draw();
        TF1* f_ComptonEdgePoly = new TF1 ("f_ComptonEdgePoly","[0]*(0.5*(1+TMath::Erf(([1]-x)/(TMath::Sqrt(2)*[2]))) + [3]*x + [4])", 40000, 90000);
        f_ComptonEdgePoly->SetParNames("N_evt","mean_charge","Sigma", "Nbg", "C" );

        std::cout << "om = " << om_tab[om] << " and mean = " << mean_tab[om] << " and C = " << C_tab[om] << '\n';

        f_ComptonEdgePoly->SetParameters(N_evt_tab[om], mean_tab[om], sigma_tab[om], Nbg_tab[om]/N_evt_tab[om], C_tab[om]/N_evt_tab[om]);
        f_ComptonEdgePoly->SetRange(mean_tab[om]-2.5*sigma_tab[om], mean_tab[om]+2.5*sigma_tab[om]);
        f_ComptonEdgePoly->SetParLimits(0, 0, 200);
        f_ComptonEdgePoly->SetParLimits(1, 0, 200000);
        f_ComptonEdgePoly->SetParLimits(2, 0, 10000);
        f_ComptonEdgePoly->SetParLimits(3, -10, 0);
        f_ComptonEdgePoly->SetParLimits(4, 0, 20);
        f_ComptonEdgePoly->Draw("same");
        spectre_om->Fit(f_ComptonEdgePoly, "RQ");
        f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
        spectre_om->Fit(f_ComptonEdgePoly, "RQ");
        f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
        spectre_om->Fit(f_ComptonEdgePoly, "RQ");

        i_om = om_tab[om];
        n_evt = (f_ComptonEdgePoly->GetParameter(0));
        mean_charge = (f_ComptonEdgePoly->GetParameter(1));
        sigma = (f_ComptonEdgePoly->GetParameter(2));
        nbg = (f_ComptonEdgePoly->GetParameter(3));
        C = (f_ComptonEdgePoly->GetParameter(4));
        mean_error = f_ComptonEdgePoly->GetParError(1);
        chi = (f_ComptonEdgePoly->GetChisquare()/f_ComptonEdgePoly->GetNDF());
        std::cout << "chi = " << chi << '\n';
        Result_tree.Fill();
        //mapping
        int om_col = (om % 13 );
        int om_row = (om / 13);
        mean_charge_map.SetBinContent( om_row+1, om_col+1, mean_charge);

        canvas->SaveAs(Form("fit/fit_Tl/charge_fit_om_%03d_run_%d.png", om_tab[om], run_number));

        outFile << i_om << "\t" << n_evt << "\t" << mean_charge << "\t" << sigma << "\t"<< nbg << "\t" << C << endl;

        delete spectre_om;
        delete f_ComptonEdgePoly;
        }

  //   double comp[520];
  //
  //   for (int i = 0; i < 520; i++) {
  //     Result_tree.GetEntry(i);
  //     comp[i] = mean_charge;
  //     // error1[i] = mean_error;
  //   }
  //
  // double yaxis;
  // double yaxis_error[4];
  // double xaxis[4] = {0, 19, 38, 57};
  // double xaxiserror[4] = {0.5, 0.5, 0.5, 0.5};
  //     for (int j = 0; j < 520; j++){
  //       for (int i = 0; i < 4; i++) {
  //         int nombre = 520*i+j;
  //         Result_tree.GetEntry(nombre);
  //         yaxis = mean_charge/comp[j];
  //         // std::cout << "om = " << j << " and var = " << yaxis << '\n';
  //         if ((yaxis < 0.9) || (yaxis > 1.1)) {
  //
  //           std::cout << " run == " << tab[i] << " om = " << j << " and var = " << yaxis << '\n';
  //         }
  //       // yaxis_error[i] = mean_error/comp*1.0 + (mean_charge/(comp*comp))*1.0*error1;
  //     }
  //   }
  // TGraphErrors comp_map (4, xaxis, yaxis, xaxiserror, yaxis_error);
  //
  // comp_map.SetNameTitle("fit_Tl", "evolution du gain de l'OM 258");
  // comp_map.GetXaxis()->SetTitle("Temps (h)");
  // comp_map.GetYaxis()->SetTitle("Gain(t)/Gain(0)");
  // comp_map.SetMarkerColor(2);
  // comp_map.SetMarkerStyle(34);
  // comp_map.SetMarkerSize(2);
  //
  file.cd();
  // comp_map.Write();
  mean_charge_map.Write();
  Result_tree.Write();

  file.Close();
  outFile.close();

  return;
}

void fit_LI_amplitude()
  {
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TFile file("Resultats_root/Resultat_amplitude_Li.root","RECREATE");
    std::ofstream outFile("Resultats_txt/failed_fit.txt");
    std::ofstream outFile2("Resultats_txt/failed_fit.txt");

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

    double debut = 27.5;

    TCanvas* canvas = new TCanvas;

    TH1F comp_map("comp_map", "comp_map", 100, 0, 100);

    TFile tree_file("histo_brut/histo_Li_system_565.root","READ");
    count++;
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

    for(int om = 1; om < 712; om+=10)
    {
      for (double j = debut; j < debut+6*39.8; j+=39.8)
      {
        double temps = j;
        if ((om > 259 && om < 520) || (om > 583 && om < 647) || (om > 679)) {temps+=(6*39.8+1*10);}
        TH1D *spectre = new TH1D ("spectre_amplitude", "", 700, 0, 2300 );
        tree->Project("spectre_amplitude", "amplitude_tree", Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, temps, temps+39.8));

        if (spectre->GetEntries() < 100) {
          std::cout << spectre->GetEntries() << '\n';
          std::cout << "" << '\n';
          std::cout << "trop peu d'entries pour l'OM " << om << '\n';
          std::cout << "" << '\n';
          delete spectre;
        }
        else if (spectre->GetMean() > 1900) {
          std::cout << "" << '\n';
          std::cout << "the amplitude sature" << '\n';
          std::cout << "" << '\n';
          outFile2 << om << "\t" << temps << "\t" << temps+39.8 << endl;
          delete spectre;
        }
        else{
          TF1 *f_Gaus = new TF1("f1", "[0]*TMath::Gaus(x,[1],[2])");
          f_Gaus->SetParNames("N_evt","mean_charge","Sigma");
          // f_Gaus->SetParameters(25, spectre->GetMean(), 100);
          f_Gaus->SetParameters(25, spectre->GetMean(), spectre->GetRMS());
          f_Gaus->SetRange(spectre->GetMean()-300, spectre->GetMean()+300);
          f_Gaus->Draw("same");
          spectre->Fit(f_Gaus, "RQ0");
          f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
          spectre->Fit(f_Gaus, "RQ0");
          f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
          spectre->Fit(f_Gaus, "RQ0");

          if (f_Gaus->GetChisquare()/f_Gaus->GetNDF() > 2 || f_Gaus->GetChisquare()/f_Gaus->GetNDF() < 0.5) {
            outFile << om << "\t" << temps << "\t" << temps+39.8 << endl;
          }

          i_om = om;
          constante = (f_Gaus->GetParameter(0));
          mean = (f_Gaus->GetParameter(1));
          sigma = (f_Gaus->GetParameter(2));
          mean_error = f_Gaus->GetParError(1) ;
          int om_col = (om % 13 );
          int om_row = (om / 13);

          Result_tree.Fill();
          spectre->Draw();
          f_Gaus->Draw("same");
          canvas->SaveAs(Form("fit/fit_Li/amplitude_Li/OM_%03d_time_%.2f.png", om, temps));
          delete spectre;
          delete f_Gaus;
        }
      }
    }

    file.cd();
    Result_tree.Write();
    file.Close();
    outFile.close();
    return;
  }

void Li_bundle_variation()
{
  TFile file("Resultats_root/Li_bundle_variation.root","RECREATE");
  TFile tree_file("Resultats_root/Resultat_amplitude_Li.root","READ");

  double mean;
  int i_om;
  TTree* tree = (TTree*)tree_file.Get("Resultat_amplitude_Li");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("i_om",1);
  tree->SetBranchAddress("i_om", &i_om);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);

  TH1D *Li_bundle_1 = new TH1D ("Li_bundle_1", "", 2500, 0, 2500);
  tree->Project("mean", "(i_om%13 > 7 && i_om/13 > 9 && i_om < 260) || (i_om > 663 && i_om < 672) || (i_om < 552 && i_om > 545) ||(i_om < 536 && i_om > 529) )");
  std::cout << "Li bundle 1 mean = " << Li_bundle_1->GetMean(1) << '\n';
  tree->Scan("mean", "(mean/Li_bundle_1->GetMean(1) > 1.2 && mean/Li_bundle_1->GetMean(1) < 0.8 &&i_om%13 > 7 && i_om/13 > 9 && i_om < 260) || (i_om > 663 && i_om < 672) || (i_om < 552 && i_om > 545) ||(i_om < 536 && i_om > 529) )");

  TH1D *Li_bundle_2 = new TH1D ("Li_bundle_2", "", 2500, 0, 2500);
  tree->Project("mean", "(i_om%13 < 7 && i_om/13 > 13 && i_om < 260) || (i_om > 652 && i_om < 647) || (i_om < 546 && i_om > 535) ||(i_om < 562 && i_om > 551) )");
  std::cout << "Li bundle 2 mean = " << Li_bundle_2->GetMean(1) << '\n';
  tree->Scan("mean", "(mean/Li_bundle_2->GetMean(1) > 1.2 && mean/Li_bundle_2->GetMean(1) < 0.8 &&i_om%13 > 7 && i_om/13 > 9 && i_om < 260) || (i_om > 663 && i_om < 672) || (i_om < 552 && i_om > 545) ||(i_om < 536 && i_om > 529) )");

  TH1D *Li_bundle_3 = new TH1D ("Li_bundle_3", "", 2500, 0, 2500 );
  tree->Project("mean", "(i_om%13 > 7 && i_om/13 < 10 && i_om < 260) || (i_om > 659 && i_om < 664) || (i_om < 578 && i_om > 567) ||(i_om < 530 && i_om > 519) )");
  std::cout << "Li bundle 3 mean = " << Li_bundle_3->GetMean(1) << '\n';
  tree->Scan("mean", "(mean/Li_bundle_3->GetMean(1) > 1.2 && mean/Li_bundle_3->GetMean(1) < 0.8 &&i_om%13 > 7 && i_om/13 > 9 && i_om < 260) || (i_om > 663 && i_om < 672) || (i_om < 552 && i_om > 545) ||(i_om < 536 && i_om > 529) )");

  TH1D *Li_bundle_4 = new TH1D ("Li_bundle_4", "", 2500, 0, 2500 );
  tree->Project("mean", "(i_om%13 < 7 && i_om/13 > 5 && i_om < 14 && i_om < 260) || (i_om < 660 && i_om > 651) )");
  std::cout << "Li bundle 4 mean = " << Li_bundle_4->GetMean(1) << '\n';
  tree->Scan("mean", "(mean/Li_bundle_4->GetMean(1) > 1.2 && mean/Li_bundle_4->GetMean(1) < 0.8 &&i_om%13 > 7 && i_om/13 > 9 && i_om < 260) || (i_om > 663 && i_om < 672) || (i_om < 552 && i_om > 545) ||(i_om < 536 && i_om > 529) )");

  TH1D *Li_bundle_5 = new TH1D ("Li_bundle_5", "", 2500, 0, 2500 );
  tree->Project("mean", "(i_om%13 < 7 && i_om/13 < 6 && i_om < 260) || (i_om > 659 && i_om < 664) || (i_om < 578 && i_om > 567) ||(i_om < 562 && i_om > 551) )");
  std::cout << "Li bundle 5 mean = " << Li_bundle_5->GetMean(1) << '\n';
  tree->Scan("mean", "(mean/Li_bundle_5->GetMean(1) > 1.2 && mean/Li_bundle_5->GetMean(1) < 0.8 &&i_om%13 > 7 && i_om/13 > 9 && i_om < 260) || (i_om > 663 && i_om < 672) || (i_om < 552 && i_om > 545) ||(i_om < 536 && i_om > 529) )");

  TH1D *Li_bundle_6 = new TH1D ("Li_bundle_6", "", 2500, 0, 2500 );
  tree->Project("mean", "(i_om%13 > 7 && i_om/13 < 10 && i_om < 520 && i_om > 259) || (i_om > 695 && i_om < 704) || (i_om < 600 && i_om > 593) ||(i_om < 616 && i_om > 609) )");
  std::cout << "Li bundle 6 mean = " << Li_bundle_6->GetMean(1) << '\n';
  tree->Scan("mean", "(mean/Li_bundle_6->GetMean(1) > 1.2 && mean/Li_bundle_6->GetMean(1) < 0.8 &&i_om%13 > 7 && i_om/13 > 9 && i_om < 260) || (i_om > 663 && i_om < 672) || (i_om < 552 && i_om > 545) ||(i_om < 536 && i_om > 529) )");

  TH1D *Li_bundle_7 = new TH1D ("Li_bundle_7", "", 2500, 0, 2500 );
  tree->Project("mean", "(i_om%13 > 7 && i_om/13 > 9 && i_om < 520 && i_om > 259) || (i_om > 703 && i_om < 712) || (i_om < 648 && i_om > 641) ||(i_om < 632 && i_om > 625) )");
  std::cout << "Li bundle 7 mean = " << Li_bundle_7->GetMean(1) << '\n';
  tree->Scan("mean", "(mean/Li_bundle_7->GetMean(1) > 1.2 && mean/Li_bundle_7->GetMean(1) < 0.8 &&i_om%13 > 7 && i_om/13 > 9 && i_om < 260) || (i_om > 663 && i_om < 672) || (i_om < 552 && i_om > 545) ||(i_om < 536 && i_om > 529) )");

  TH1D *Li_bundle_8 = new TH1D ("Li_bundle_8", "", 2500, 0, 2500 );
  tree->Project("mean", "(i_om%13 < 7 && i_om/13 < 6 && i_om < 520 && i_om > 259) || (i_om > 679 && i_om < 684) || (i_om < 594 && i_om > 583) ||(i_om < 610 && i_om > 599) )");
  std::cout << "Li bundle 8 mean = " << Li_bundle_8->GetMean(1) << '\n';
  tree->Scan("mean", "(mean/Li_bundle_8->GetMean(1) > 1.2 && mean/Li_bundle_8->GetMean(1) < 0.8 &&i_om%13 > 7 && i_om/13 > 9 && i_om < 260) || (i_om > 663 && i_om < 672) || (i_om < 552 && i_om > 545) ||(i_om < 536 && i_om > 529) )");

  TH1D *Li_bundle_9 = new TH1D ("Li_bundle_9", "", 2500, 0, 2500 );
  tree->Project("mean", "(i_om%13 < 7 && i_om/13 > 5 && i_om < 14 && i_om < 520 && i_om > 259) || (i_om > 683 && i_om < 691) )");
  std::cout << "Li bundle 9 mean = " << Li_bundle_9->GetMean(1) << '\n';
  tree->Scan("mean", "(mean/Li_bundle_9->GetMean(1) > 1.2 && mean/Li_bundle_9->GetMean(1) < 0.8 &&i_om%13 > 7 && i_om/13 > 9 && i_om < 260) || (i_om > 663 && i_om < 672) || (i_om < 552 && i_om > 545) ||(i_om < 536 && i_om > 529) )");

  TH1D *Li_bundle_10 = new TH1D ("Li_bundle_10", "", 2500, 0, 2500 );
  tree->Project("mean", "(i_om%13 < 7 && i_om/13 >13 && i_om < 520 && i_om > 259) || (i_om > 691 && i_om < 696) || (i_om < 642 && i_om > 631) ||(i_om < 626 && i_om > 615) )");
  std::cout << "Li bundle 10 mean = " << Li_bundle_10->GetMean(1) << '\n';
  tree->Scan("mean", "(mean/Li_bundle_10->GetMean(1) > 1.2 && mean/Li_bundle_10->GetMean(1) < 0.8 &&i_om%13 > 7 && i_om/13 > 9 && i_om < 260) || (i_om > 663 && i_om < 672) || (i_om < 552 && i_om > 545) ||(i_om < 536 && i_om > 529) )");

  file.cd();
  Li_bundle_1->Write();
  Li_bundle_2->Write();
  Li_bundle_3->Write();
  Li_bundle_4->Write();
  Li_bundle_5->Write();
  Li_bundle_6->Write();
  Li_bundle_7->Write();
  Li_bundle_8->Write();
  Li_bundle_9->Write();
  Li_bundle_10->Write();
  file.Close();
}


// void variation_Li() {
//   gStyle->SetOptFit(1);
//   gStyle->SetOptStat(0);
//   TH1::SetDefaultSumw2();
//   TH2::SetDefaultSumw2();
//
//   TFile file("Resultats_root/Resultat_Li_variation.root","RECREATE");
//   TFile tree_file("Resultats_root/Resultat_charge_Li.root","RECREATE");
//
//   TTree* tree = (TTree*)tree_file.Get("Result_tree");
//   tree->SetBranchStatus("*",0);
//   tree->SetBranchStatus("om_number",1);
//   tree->SetBranchAddress("om_number", &om_number);
//   tree->SetBranchStatus("time",1);
//   tree->SetBranchAddress("time", &time);
//   tree->SetBranchStatus("charge_tree",1);
//   tree->SetBranchAddress("charge_tree", &charge_tree);
//   tree->SetBranchStatus("amplitude_tree",1);
//   tree->SetBranchAddress("amplitude_tree", &amplitude_tree);
//
//   double comp[6];
//   double error1[6];
//   for (int i = 0; i < 6; i++) {
//     Result_tree.GetEntry(i);
//     comp[i] = mean;
//     error1[i] = mean_error;
//   }
//
//   double yaxis_0[4];
//   double yaxis_1[4];
//   double yaxis_2[4];
//   double yaxis_3[4];
//   double yaxis_4[4];
//   double yaxis_5[4];
//   double yaxis_error_0[4];
//   double yaxis_error_1[4];
//   double yaxis_error_2[4];
//   double yaxis_error_3[4];
//   double yaxis_error_4[4];
//   double yaxis_error_5[4];
//   double xaxis[4] = {0, 19, 38, 57};
//   double xaxiserror[4] = {0.5, 0.5, 0.5, 0.5};
//
//   for (int i = 0; i < 4; i++) {
//     Result_tree.GetEntry(i*6);
//     yaxis_0[i] = mean/comp[0]*1.0;
//     yaxis_error_0[i] = mean_error/comp[0]*1.0 + (mean/(comp[0]*comp[0]))*1.0*error1[0];
//   }
//   for (int i = 0; i < 4; i++) {
//     Result_tree.GetEntry(i*6+1);
//     yaxis_1[i] = mean/comp[1]*1.0;
//     yaxis_error_1[i] = mean_error/comp[1]*1.0 + (mean/(comp[1]*comp[1]))*1.0*error1[1];
//   }
//   for (int i = 0; i < 4; i++) {
//     Result_tree.GetEntry(i*6+2);
//     yaxis_2[i] = mean/comp[2]*1.0;
//     yaxis_error_2[i] = mean_error/comp[2]*1.0 + (mean/(comp[2]*comp[2]))*1.0*error1[2];
//   }
//   for (int i = 0; i < 4; i++) {
//     Result_tree.GetEntry(i*6+3);
//     yaxis_3[i] = mean/comp[3]*1.0;
//     yaxis_error_3[i] = mean_error/comp[3]*1.0 + (mean/(comp[3]*comp[3]))*1.0*error1[3];
//   }
//   for (int i = 0; i < 4; i++) {
//     Result_tree.GetEntry(i*6+4);
//     yaxis_4[i] = mean/comp[4]*1.0;
//     yaxis_error_4[i] = mean_error/comp[4]*1.0 + (mean/(comp[4]*comp[4]))*1.0*error1[4];
//   }
//   for (int i = 0; i < 4; i++) {
//     Result_tree.GetEntry(i*6+5);
//     yaxis_5[i] = mean/comp[5]*1.0;
//     yaxis_error_5[i] = mean_error/comp[5]*1.0 + (mean/(comp[5]*comp[5]))*1.0*error1[5];
//   }
//
//   TGraphErrors gr (4, xaxis, yaxis_0, xaxiserror, yaxis_error_0);
//   TGraphErrors comp_map_1 (4, xaxis, yaxis_1, xaxiserror, yaxis_error_1);
//   TGraphErrors comp_map_2 (4, xaxis, yaxis_2, xaxiserror, yaxis_error_2);
//   TGraphErrors comp_map_3 (4, xaxis, yaxis_3, xaxiserror, yaxis_error_3);
//   TGraphErrors comp_map_4 (4, xaxis, yaxis_4, xaxiserror, yaxis_error_4);
//   TGraphErrors comp_map_5 (4, xaxis, yaxis_5, xaxiserror, yaxis_error_5);
//
//   gr.SetNameTitle("pic1", "evolution du gain de l'OM 258");
//   gr.GetXaxis()->SetTitle("Temps (h)");
//   gr.GetYaxis()->SetTitle("Gain(t)/Gain(0)");
//   gr.SetMarkerColor(4);
//   gr.SetMarkerStyle(21);
//   comp_map_1.SetMarkerColor(3);
//   comp_map_1.SetMarkerStyle(21);
//   comp_map_1.SetNameTitle("pic2", "evolution du gain de l'OM 258");
//   comp_map_2.SetMarkerColor(2);
//   comp_map_2.SetMarkerStyle(21);
//   comp_map_2.SetNameTitle("pic3", "evolution du gain de l'OM 258");
//   comp_map_3.SetMarkerColor(1);
//   comp_map_3.SetMarkerStyle(21);
//   comp_map_3.SetNameTitle("pic4", "evolution du gain de l'OM 258");
//   comp_map_4.SetMarkerColor(5);
//   comp_map_4.SetMarkerStyle(21);
//   comp_map_4.SetNameTitle("pic5", "evolution du gain de l'OM 258");
//   comp_map_5.SetMarkerColor(6);
//   comp_map_5.SetMarkerStyle(21);
//   comp_map_5.SetNameTitle("pic6", "evolution du gain de l'OM 258");
//
//   file.cd();
//   gr.Write();
//   comp_map_1.Write();
//   comp_map_2.Write();
//   comp_map_3.Write();
//   comp_map_4.Write();
//   comp_map_5.Write();
//   file.Close();
//   return;
// }


// void fit_LI_ampl()
//   {
//     gStyle->SetOptFit(1);
//     gStyle->SetOptStat(0);
//     TH1::SetDefaultSumw2();
//     TH2::SetDefaultSumw2();
//
//     TFile file("Resultats_root/Resultat_charge_Li.root","RECREATE");
//
//     int i_om;
//     double constante;
//     double mean;
//     double mean_error;
//     double sigma;
//     int run_number;
//     double error;
//     double count = -1;
//
//     TTree Result_tree("Result_tree","");
//     Result_tree.Branch("run_number", &run_number);
//     Result_tree.Branch("i_om", &i_om);
//     Result_tree.Branch("constante", &constante);
//     Result_tree.Branch("mean", &mean);
//     Result_tree.Branch("mean_error", &mean_error);
//     Result_tree.Branch("sigma", &sigma);
//     Result_tree.Branch("count", &count);
//
//     int tab[4] = {549, 552, 555, 559};
//     double debut[4] = {10, 50, 55, 55};
//     int lim_min[6] = {50, 50, 95, 130, 170, 210};
//     int lim_max[6] = {48, 88, 120, 168, 210, 250};
//
//     TCanvas* canvas = new TCanvas;
//     canvas->SetLogy();
//
//     TH1F comp_map("comp_map", "comp_map", 100, 0, 100);
//     TH2F mean_charge_map("histo_om_mean_charge_map", "mean_charge_map", 20, 0, 20, 13, 0, 13);
//     for (int i = 0; i < 4; i++) {
//       // i = 0
//       run_number = tab[i];
//       TFile tree_file(Form("histo_brut/histo_Li_system_%d.root", run_number),"READ");
//       count++;
//       int om_number;
//       double time;
//       double charge_tree;
//       double amplitude_tree;
//       TTree* tree = (TTree*)tree_file.Get("Result_tree");
//       tree->SetBranchStatus("*",0);
//       tree->SetBranchStatus("om_number",1);
//       tree->SetBranchAddress("om_number", &om_number);
//       tree->SetBranchStatus("time",1);
//       tree->SetBranchAddress("time", &time);
//       tree->SetBranchStatus("charge_tree",1);
//       tree->SetBranchAddress("charge_tree", &charge_tree);
//       tree->SetBranchStatus("amplitude_tree",1);
//       tree->SetBranchAddress("amplitude_tree", &amplitude_tree);
//
//
//
//       for(int om = 1; om<259; om+=259)
//       {
//         for (double j = debut[i]; j < debut[i]+6*38; j+=38)
//
//         // for (int j = 0; j < 6; j++)
//         {
//           TF1* f_ComptonEdgePoly = new TF1 ("f_ComptonEdgePoly","[0]/2.0*(1+TMath::Erf(([1]-x)/(TMath::Sqrt(2)*[2])))+[3]*x", 40000, 90000);
//           // int min = lim_min[j];
//           // int max = lim_max[j];
//           // tree->Draw("amplitude_tree ", Form("om_number == 1 && time > %d && time < %d", lim_min[j], lim_max[j]));
//           tree->Draw("amplitude_tree", Form("om_number == 1 && time < %f && time > %f", j, j-40));
//
//           TH1D *spectre = (TH1D*)gPad->GetPrimitive("htemp");
//           spectre->Draw();
//
//           TF1 *f_Gaus = new TF1("f1", "[0]*TMath::Gaus(x,[1],[2])");
//           f_Gaus->SetParNames("N_evt","mean_charge","Sigma");
//
//           if (om == 1)        //om multiple de (13)-1
//           {
//             std::cout << "mean = " << spectre->GetMean(1) << '\n';
//             f_Gaus->SetParameters(25, spectre->GetMean(1), 2100);
//             f_Gaus->SetRange(spectre->GetMean(1)-300, spectre->GetMean(1)+300);
//             f_Gaus->Draw("same");
//             spectre->Fit(f_Gaus, "RQ");
//             f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
//             spectre->Fit(f_Gaus, "RQ");
//             f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
//             spectre->Fit(f_Gaus, "RQ");
//           }
//
//           i_om = om;
//           constante = (f_Gaus->GetParameter(0));
//           mean = (f_Gaus->GetParameter(1));
//           sigma = (f_Gaus->GetParameter(2));
//           mean_error = f_Gaus->GetParError(1) ;
//           //mapping
//           int om_col = (om % 13 );
//           int om_row = (om / 13);
//           mean_charge_map.SetBinContent( om_row+1, om_col+1, mean);
//
//           Result_tree.Fill();
//           canvas->SaveAs(Form("fit/fit_Li/amplitude_Li/OM_%03d_run_%i_time_%f.png", om, tab[i], j));
//
//           delete spectre;
//           delete f_Gaus;
//         }
//       }
//     }
//
//     return;
//
//     double comp[6];
//     double error1[6];
//     for (int i = 0; i < 6; i++) {
//       Result_tree.GetEntry(i);
//       comp[i] = mean;
//       error1[i] = mean_error;
//     }
//
//     double yaxis_0[4];
//     double yaxis_1[4];
//     double yaxis_2[4];
//     double yaxis_3[4];
//     double yaxis_4[4];
//     double yaxis_5[4];
//     double yaxis_error_0[4];
//     double yaxis_error_1[4];
//     double yaxis_error_2[4];
//     double yaxis_error_3[4];
//     double yaxis_error_4[4];
//     double yaxis_error_5[4];
//     double xaxis[4] = {0, 19, 38, 57};
//     double xaxiserror[4] = {0.5, 0.5, 0.5, 0.5};
//
//     for (int i = 0; i < 4; i++) {
//       Result_tree.GetEntry(i*6);
//       yaxis_0[i] = mean/comp[0]*1.0;
//       yaxis_error_0[i] = mean_error/comp[0]*1.0 + (mean/(comp[0]*comp[0]))*1.0*error1[0];
//     }
//     for (int i = 0; i < 4; i++) {
//       Result_tree.GetEntry(i*6+1);
//       yaxis_1[i] = mean/comp[1]*1.0;
//       yaxis_error_1[i] = mean_error/comp[1]*1.0 + (mean/(comp[1]*comp[1]))*1.0*error1[1];
//     }
//     for (int i = 0; i < 4; i++) {
//       Result_tree.GetEntry(i*6+2);
//       yaxis_2[i] = mean/comp[2]*1.0;
//       yaxis_error_2[i] = mean_error/comp[2]*1.0 + (mean/(comp[2]*comp[2]))*1.0*error1[2];
//     }
//     for (int i = 0; i < 4; i++) {
//       Result_tree.GetEntry(i*6+3);
//       yaxis_3[i] = mean/comp[3]*1.0;
//       yaxis_error_3[i] = mean_error/comp[3]*1.0 + (mean/(comp[3]*comp[3]))*1.0*error1[3];
//     }
//     for (int i = 0; i < 4; i++) {
//       Result_tree.GetEntry(i*6+4);
//       yaxis_4[i] = mean/comp[4]*1.0;
//       yaxis_error_4[i] = mean_error/comp[4]*1.0 + (mean/(comp[4]*comp[4]))*1.0*error1[4];
//     }
//     for (int i = 0; i < 4; i++) {
//       Result_tree.GetEntry(i*6+5);
//       yaxis_5[i] = mean/comp[5]*1.0;
//       yaxis_error_5[i] = mean_error/comp[5]*1.0 + (mean/(comp[5]*comp[5]))*1.0*error1[5];
//     }
//
//     TGraphErrors gr (4, xaxis, yaxis_0, xaxiserror, yaxis_error_0);
//     TGraphErrors comp_map_1 (4, xaxis, yaxis_1, xaxiserror, yaxis_error_1);
//     TGraphErrors comp_map_2 (4, xaxis, yaxis_2, xaxiserror, yaxis_error_2);
//     TGraphErrors comp_map_3 (4, xaxis, yaxis_3, xaxiserror, yaxis_error_3);
//     TGraphErrors comp_map_4 (4, xaxis, yaxis_4, xaxiserror, yaxis_error_4);
//     TGraphErrors comp_map_5 (4, xaxis, yaxis_5, xaxiserror, yaxis_error_5);
//
//     gr.SetNameTitle("pic1", "evolution du gain de l'OM 258");
//     gr.GetXaxis()->SetTitle("Temps (h)");
//     gr.GetYaxis()->SetTitle("Gain(t)/Gain(0)");
//     gr.SetMarkerColor(4);
//     gr.SetMarkerStyle(21);
//     comp_map_1.SetMarkerColor(3);
//     comp_map_1.SetMarkerStyle(21);
//     comp_map_1.SetNameTitle("pic2", "evolution du gain de l'OM 258");
//     comp_map_2.SetMarkerColor(2);
//     comp_map_2.SetMarkerStyle(21);
//     comp_map_2.SetNameTitle("pic3", "evolution du gain de l'OM 258");
//     comp_map_3.SetMarkerColor(1);
//     comp_map_3.SetMarkerStyle(21);
//     comp_map_3.SetNameTitle("pic4", "evolution du gain de l'OM 258");
//     comp_map_4.SetMarkerColor(5);
//     comp_map_4.SetMarkerStyle(21);
//     comp_map_4.SetNameTitle("pic5", "evolution du gain de l'OM 258");
//     comp_map_5.SetMarkerColor(6);
//     comp_map_5.SetMarkerStyle(21);
//     comp_map_5.SetNameTitle("pic6", "evolution du gain de l'OM 258");
//
//     file.cd();
//     mean_charge_map.Write();
//     gr.Write();
//     comp_map_1.Write();
//     comp_map_2.Write();
//     comp_map_3.Write();
//     comp_map_4.Write();
//     comp_map_5.Write();
//     Result_tree.Write();
//     file.Close();
//     return;
// }

// void fit_Ref_OM()
//   {
//     gStyle->SetOptFit(1);
//     gStyle->SetOptStat(0);
//     TH1::SetDefaultSumw2();
//     TH2::SetDefaultSumw2();
//
//     TFile file("Resultats_root/Resultat_charge_Ref.root","RECREATE");
//
//     int i_om;
//     double constante;
//     double mean;
//     double mean_error;
//     double sigma;
//     int run_number;
//     double error;
//     double count = -1;
//
//     TTree Result_tree("Result_tree","");
//     Result_tree.Branch("run_number", &run_number);
//     Result_tree.Branch("i_om", &i_om);
//     Result_tree.Branch("constante", &constante);
//     Result_tree.Branch("mean", &mean);
//     Result_tree.Branch("mean_error", &mean_error);
//     Result_tree.Branch("sigma", &sigma);
//     Result_tree.Branch("count", &count);
//
//     int tab[4] = {547, 552, 555, 559};
//     double debut[4] = {50e9, 50e9, 55e9, 55e9};
//
//     TCanvas* canvas = new TCanvas;
//     canvas->SetLogy();
//
//     TH1F comp_map("comp_map", "comp_map", 100, 0, 100);
//     TH2F mean_charge_map("histo_om_mean_charge_map", "mean_charge_map", 20, 0, 20, 13, 0, 13);
//     for (int i = 0; i < 4; i++) {
//       run_number = tab[i];
//       TFile tree_file(Form("histo_brut/histo_Li_system_%d.root", run_number),"READ");
//       count++;
//       int om_number;
//       double time;
//       double charge_tree;
//       TTree* tree = (TTree*)tree_file.Get("Result_tree");
//       tree->SetBranchStatus("*",0);
//       tree->SetBranchStatus("om_number",1);
//       tree->SetBranchAddress("om_number", &om_number);
//       tree->SetBranchStatus("time",1);
//       tree->SetBranchAddress("time", &time);
//       tree->SetBranchStatus("charge_tree",1);
//       tree->SetBranchAddress("charge_tree", &charge_tree);
//
//       for(int om = 800; om<801; om+=259)
//       {
//         // for (double j = debut[i]; j < 250e9; j+=38e9)
//         // for (double j = debut[i]; j < debut[i]+6*38e9; j+=38e9)
//         // {
//
//         double convert = par[1]/976.0;
//         double deltae = par[2]*sqrt(par[1]);
//         double deltae2 = par[2]*sqrt(par[4]);
//
//         par[0]*(7.11*TMath::Gaus(x[0],par[1],deltae) + 1.84*TMath::Gaus(x[0],(par[1]+convert*72),deltae) +
//         0.441*TMath::Gaus(x[0],(par[1]+convert*84),deltae))
//         +
//         par[3]*(1.548*TMath::Gaus(x[0],par[4],deltae2) + 0.429*TMath::Gaus(x[0],(par[4]+convert*73),deltae2) +
//         0.1057*TMath::Gaus(x[0],(par[4]+convert*85),deltae2))
//         +
//         exp(-par[5]*x[0]+par[6]);
//
//
//           tree->Draw("charge_tree", Form("om_number == 1 && time < %f && time > %f", j, j-40e9));
//           TH1D *spectre = (TH1D*)gPad->GetPrimitive("htemp");
//           TF1 *f_Gaus = new TF1("f1", "[0]*TMath::Gaus(x,[1],[2])");
//           f_Gaus->SetParNames("N_evt","mean_charge","Sigma");
//
//           if (om == 33)        //om multiple de (13)-1
//           {
//             f_Gaus->SetParameters(25, spectre->GetMean(1), 2100);
//             f_Gaus->SetRange(0,300e9);
//             f_Gaus->Draw("same");
//             spectre->Fit(f_Gaus, "RQ");
//             f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
//             spectre->Fit(f_Gaus, "RQ");
//             f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
//             spectre->Fit(f_Gaus, "RQ");
//           }
//
//           i_om = om;
//           constante = (f_Gaus->GetParameter(0));
//           mean = (f_Gaus->GetParameter(1));
//           sigma = (f_Gaus->GetParameter(2));
//           mean_error = f_Gaus->GetParError(1) ;
//           //mapping
//           int om_col = (om % 13 );
//           int om_row = (om / 13);
//           mean_charge_map.SetBinContent( om_row+1, om_col+1, mean);
//
//           Result_tree.Fill();
//           canvas->SaveAs(Form("fit/fit_Li/charge_Li_om_%03d_run_%i_time_%f.png", om, tab[i], j));
//
//           delete spectre;
//           delete f_Gaus;
//         // }
//       }
//     }
//
//     return;
// }
