#define _USE_MATH_DEFINES
#include <fstream>
#include <sstream>
#include <fstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TVector.h"
#include "TVector3.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TFractionFitter.h>
#include <TLegend.h>
using namespace std;

const int gain_n_bin = 10001; // Precision au 10 000 ème
const float gain_bin_min = 0.5;
const float gain_bin_max = 1.5;
const float gain_bin_width = (gain_bin_max-gain_bin_min)/(gain_n_bin-1);

const int amplitude_n_bin = 1024;
const float amplitude_bin_min = 0e-05;
const float amplitude_bin_max = 2e5;
const float amplitude_bin_width = (amplitude_bin_max-amplitude_bin_min)/(amplitude_n_bin-1);

TH2D* MC_Simu(int om){

  TFile *file = new TFile("../Li_system_714.root", "READ");

  double charge;
  int om_num;

  TRandom3 rando;

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_num);
  tree->SetBranchStatus("charge_tree",1);
  tree->SetBranchAddress("charge_tree", &charge);

  TH2D* MC_Simu_ubc = new TH2D("Modele_Ref_OM", "Modele_Ref_OM",
  gain_n_bin, gain_bin_min - (gain_bin_width/2), gain_bin_max + (gain_bin_width/2),
  amplitude_n_bin, amplitude_bin_min, amplitude_bin_max);
  for (int i = 0; i < tree->GetEntries()*8/10; i++) {
  // for (int i = 0; i < 10000000; i++) {

    if (i%100000 == 0) {
      std::cout << "/* entry = " << i << '\n';
    }
    double E_kolmo =0;
    tree->GetEntry(i);

    // if (om_num_new->at(j) < 520 && (om_num_new->at(j)%13)!=0 && (om_num_new->at(j)%13)!=12) {           //apply to take only the wanted OM (MW8, MW5, ...)
    if (om_num ==  om && charge > 0) {
      for(float gain = gain_bin_min; gain <= gain_bin_max; gain += gain_bin_width) {
        E_kolmo = charge*gain;
        MC_Simu_ubc->Fill(gain, E_kolmo);
      }

    }
  }



  TFile *newfile = new TFile(Form("Modele_OM_%d_714.root", om), "RECREATE");
  newfile->cd();
  MC_Simu_ubc->Write();
  newfile->Close();

  return MC_Simu_ubc;
}

TH2D* MC_Simu_936(int om){

  TFile *file = new TFile("../Li_system_936_bin.root", "READ");

  double charge;
  int om_num;

  TRandom3 rando;

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_num);
  tree->SetBranchStatus("charge_tree",1);
  tree->SetBranchAddress("charge_tree", &charge);

  TH2D* MC_Simu_ubc = new TH2D("Modele_Ref_OM", "Modele_Ref_OM",
  gain_n_bin, gain_bin_min - (gain_bin_width/2), gain_bin_max + (gain_bin_width/2),
  amplitude_n_bin, amplitude_bin_min, amplitude_bin_max);
  for (int i = 0; i < tree->GetEntries(); i++) {
  // for (int i = 0; i < 10000000; i++) {

    if (i%100000 == 0) {
      std::cout << "/* entry = " << i << '\n';
    }
    double E_kolmo =0;
    tree->GetEntry(i);

    // if (om_num_new->at(j) < 520 && (om_num_new->at(j)%13)!=0 && (om_num_new->at(j)%13)!=12) {           //apply to take only the wanted OM (MW8, MW5, ...)
    if (om_num ==  om && charge > 0) {
      for(float gain = gain_bin_min; gain <= gain_bin_max; gain += gain_bin_width) {
        E_kolmo = charge*gain;
        MC_Simu_ubc->Fill(gain, E_kolmo);
      }

    }
  }


  TFile *newfile = new TFile(Form("Modele_OM_%d_936.root", om), "RECREATE");
  newfile->cd();
  MC_Simu_ubc->Write();
  newfile->Close();

  return MC_Simu_ubc;
}


int main(int argc, char const *argv[])
{
  for (int i = 800; i < 801; i++) {
    MC_Simu(i);
  }
  // for (int i = 712; i < 717; i++) {
  //   MC_Simu_936(i);
  // }
  return 0;
}
