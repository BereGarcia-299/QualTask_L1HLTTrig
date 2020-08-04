#include <assert.h>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <tgmath.h>
#include <iostream>
#include <TChain.h>
#include <TLegend.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"
#include <TAttMarker.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TF1.h"
#include "TLine.h"
#include <fstream>

using namespace std;

float pi_dos = 360;
float pi_mi = 180;
float original = 3.14159;

float delta_phi(float phi_uno,float phi_dos){                                                                                                                                   

  float final_phi = fmod(abs(phi_uno - phi_dos), 2*original);
  float phi = final_phi > original  ? 2*original - final_phi : final_phi;
  return phi;

}


//----LogBins Function
void log_bins(const float x_low, const float x_high, const int num_bins, double *bins){

  int const size = num_bins+1;
  float logbins[26] = {};
  
  *(bins+0)= x_low;
  *(bins+ num_bins) = x_high;


  logbins[0] = TMath::Log10(x_low);
  logbins[num_bins] = TMath::Log10(x_high);
  

  float interval = (logbins[25] - logbins[0])/num_bins;

  for(int index = 1; index < num_bins; index++){
    logbins[index] = logbins[0] + index*interval;
    *(bins +index) = TMath::Power(10, logbins[index]);
    cout << "This went into  your bins_new array " << *(bins+index) << endl;
  }

  return;
}



int Fitting_Functions_TurnOnPlots(){


double bins_new[26]={};

double *ptr = bins_new; 

log_bins(1,40,25,ptr);

for (int i = 0; i < 26; i++)
{
  cout << "Array index number " << i << " saved this number: " << bins_new[i] << endl;
}


//--------------Reading In Plots
TFile *file = new TFile("/Users/berenicegarcia/Desktop/ATLAS/Qual-Task/histograms.root","READ");

ifstream points_L1_10GeV;
points_L1_10GeV.open("/Users/berenicegarcia/Desktop/ATLAS/Qual-Task/points_L1_10GeV.txt");

ifstream point_L1_10GeV_edcaps;
point_L1_10GeV_edcaps.open("/Users/berenicegarcia/Desktop/ATLAS/Qual-Task/point_L1_10GeV_edcaps.txt");


float x_cor[35] = {};
float y_cor[35] = {};

for(int iCord_barrel = 0; iCord_barrel < 34; iCord_barrel++) // To get you all the lines.
    {
	    points_L1_10GeV >> x_cor[iCord_barrel] >> y_cor[iCord_barrel];  //REading In The x and y Values

    }


float x_cor_endcaps[35] = {};
float y_cor_endcaps[35] = {};



for(int iCord_endcaps = 0; iCord_endcaps < 34; iCord_endcaps++) // To get you all the lines.
    {
	    point_L1_10GeV_edcaps >> x_cor_endcaps[iCord_endcaps] >> y_cor_endcaps[iCord_endcaps];  //REading In The x and y Values

    }





//---------------------------------------//
//----------ALl Offline Photons----------//
//---------------------------------------//
TH1D *all_photons = (TH1D*) file->Get("all_photons"); 
TH1D *all_photons_barrel = (TH1D*) file->Get("all_photons_barrel");
TH1D *all_photons_endcaps = (TH1D*) file->Get("all_photons_endcaps");

//---------------------------------------//
//--------L1EM RoI pT Distributions------//
//---------------------------------------//


TH1D *photons_wL1EM_15GeV = (TH1D*) file->Get("photons_wL1EM_15GeV");
TH1D *photons_wL1EM_25GeV = (TH1D*) file->Get("photons_wL1EM_25GeV");
TH1D *photons_wL1EM_35GeV = (TH1D*) file->Get("photons_wL1EM_35GeV");

TH1D *photons_wL1EM_10GeV_barrel = (TH1D*) file->Get("photons_wL1EM_10GeV_barrel");
TH1D *photons_wL1EM_15GeV_barrel = (TH1D*) file->Get("photons_wL1EM_15GeV_barrel");
TH1D *photons_wL1EM_25GeV_barrel = (TH1D*) file->Get("photons_wL1EM_25GeV_barrel");
TH1D *photons_wL1EM_35GeV_barrel = (TH1D*) file->Get("photons_wL1EM_35GeV_barrel");

TH1D *photons_wL1EM_10GeV_endcaps = (TH1D*) file->Get("photons_wL1EM_10GeV_endcaps");
TH1D *photons_wL1EM_15GeV_endcaps = (TH1D*) file->Get("photons_wL1EM_15GeV_endcaps");
TH1D *photons_wL1EM_25GeV_endcaps = (TH1D*) file->Get("photons_wL1EM_25GeV_endcaps");
TH1D *photons_wL1EM_35GeV_endcaps = (TH1D*) file->Get("photons_wL1EM_35GeV_endcaps");

//-----SumET Distribution
TH1D *barrel_wL1EM_10GeV_SumET_dis_50 = (TH1D*) file->Get("barrel_wL1EM_10GeV_SumET_dis_50");
TH1D *barrel_wL1EM_10GeV_SumET_dis_10 = (TH1D*) file->Get("barrel_wL1EM_10GeV_SumET_dis_10");

TH1D *endcaps_wL1EM_10GeV_SumET_dis_50 = (TH1D*) file->Get("endcaps_wL1EM_10GeV_SumET_dis_50");
TH1D *endcaps_wL1EM_10GeV_SumET_dis_10 = (TH1D*) file->Get("endcaps_wL1EM_10GeV_SumET_dis_10");

//--------------------------------------//
//--------HLT pT Distribution-----------//
//--------------------------------------//
TH1D *photons_w_HLT_15 = (TH1D*) file->Get("photons_w_HLT_15");
TH1D *photons_w_HLT_25 = (TH1D*) file->Get("photons_w_HLT_25");
TH1D *photons_w_HLT_35 = (TH1D*) file->Get("photons_w_HLT_35");

TH1D *photons_w_HLT_15_barrel = (TH1D*) file->Get("photons_w_HLT_15_barrel");
TH1D *photons_w_HLT_25_barrel = (TH1D*) file->Get("photons_w_HLT_25_barrel");
TH1D *photons_w_HLT_35_barrel = (TH1D*) file->Get("photons_w_HLT_35_barrel");

TH1D *photons_w_HLT_15_endcaps = (TH1D*) file->Get("photons_w_HLT_15_endcaps");
TH1D *photons_w_HLT_25_endcaps = (TH1D*) file->Get("photons_w_HLT_25_endcaps");
TH1D *photons_w_HLT_35_endcaps = (TH1D*) file->Get("photons_w_HLT_35_endcaps");

//--------------------------------------------------//
//------Photons Failed TO Fire A HLT Trigger--------//
//--------------------------------------------------//
TH1D *photons_nofire_HLT_pt = (TH1D*) file->Get("photons_nofire_HLT_pt");
TH1D *photphotons_wL1EM_15GeV_barrelons_nofire_HLT_eta =(TH1D*) file->Get("photons_nofire_HLT_eta");
TH1D *photons_nofire_HLT_sumET = (TH1D*) file->Get("photons_nofire_HLT_sumET");

//---------------------------------------------------------//
//------L1EM Efficiencies Dependent on Centrality----------//
//---------------------------------------------------------//
//--------Barrel & Endcaps
TH1D *all_photons_010 = (TH1D*) file->Get("all_photons_010"); 
TH1D *all_photons_1030 = (TH1D*) file->Get("all_photons_1030");
TH1D *all_photons_3050 = (TH1D*) file->Get("all_photons_3050");
TH1D *all_photons_5080 = (TH1D*) file->Get("all_photons_5080");

TH1D *photon_w_L1EM_10GeV_010 = (TH1D*) file->Get("photon_w_L1EM_10GeV_010");
TH1D *photon_w_L1EM_10GeV_1030 = (TH1D*) file->Get("photon_w_L1EM_10GeV_1030");
TH1D *photon_w_L1EM_10GeV_3050 = (TH1D*) file->Get("photon_w_L1EM_10GeV_3050");
TH1D *photon_w_L1EM_10GeV_5080 = (TH1D*) file->Get("photon_w_L1EM_10GeV_5080");

TH1D *photon_w_L1EM_15GeV_010 = (TH1D*) file->Get("photon_w_L1EM_15GeV_010");
TH1D *photon_w_L1EM_15GeV_1030 = (TH1D*) file->Get("photon_w_L1EM_15GeV_1030");
TH1D *photon_w_L1EM_15GeV_3050 = (TH1D*) file->Get("photon_w_L1EM_15GeV_3050");
TH1D *photon_w_L1EM_15GeV_5080 = (TH1D*) file->Get("photon_w_L1EM_15GeV_5080");

TH1D *photon_w_L1EM_25GeV_010 = (TH1D*) file->Get("photon_w_L1EM_25GeV_010");
TH1D *photon_w_L1EM_25GeV_1030 = (TH1D*) file->Get("photon_w_L1EM_25GeV_1030");
TH1D *photon_w_L1EM_25GeV_3050 = (TH1D*) file->Get("photon_w_L1EM_25GeV_3050");
TH1D *photon_w_L1EM_25GeV_5080 = (TH1D*) file->Get("photon_w_L1EM_25GeV_5080");

TH1D *photon_w_L1EM_35GeV_010 = (TH1D*) file->Get("photon_w_L1EM_35GeV_010");
TH1D *photon_w_L1EM_35GeV_1030 = (TH1D*) file->Get("photon_w_L1EM_35GeV_1030");
TH1D *photon_w_L1EM_35GeV_3050 = (TH1D*) file->Get("photon_w_L1EM_35GeV_3050");
TH1D *photon_w_L1EM_35GeV_5080 = (TH1D*) file->Get("photon_w_L1EM_35GeV_5080");
//------Barrel
TH1D *all_photons_010_b = (TH1D*) file->Get("all_photons_010_b"); 
TH1D *all_photons_1030_b = (TH1D*) file->Get("all_photons_1030_b");
TH1D *all_photons_3050_b = (TH1D*) file->Get("all_photons_3050_b");
TH1D *all_photons_5080_b = (TH1D*) file->Get("all_photons_5080_b");


TH1D *photon_w_L1EM_10GeV_010_b = (TH1D*) file->Get("photon_w_L1EM_10GeV_010_b");
TH1D *photon_w_L1EM_10GeV_1030_b = (TH1D*) file->Get("photon_w_L1EM_10GeV_1030_b");
TH1D *photon_w_L1EM_10GeV_3050_b = (TH1D*) file->Get("photon_w_L1EM_10GeV_3050_b");
TH1D *photon_w_L1EM_10GeV_5080_b = (TH1D*) file->Get("photon_w_L1EM_10GeV_5080_b");

TH1D *photon_w_L1EM_15GeV_010_b = (TH1D*) file->Get("photon_w_L1EM_15GeV_010_b");
TH1D *photon_w_L1EM_15GeV_1030_b = (TH1D*) file->Get("photon_w_L1EM_15GeV_1030_b");
TH1D *photon_w_L1EM_15GeV_3050_b = (TH1D*) file->Get("photon_w_L1EM_15GeV_3050_b");
TH1D *photon_w_L1EM_15GeV_5080_b = (TH1D*) file->Get("photon_w_L1EM_15GeV_5080_b");

TH1D *photon_w_L1EM_25GeV_010_b = (TH1D*) file->Get("photon_w_L1EM_25GeV_010_b"); 
TH1D *photon_w_L1EM_25GeV_1030_b = (TH1D*) file->Get("photon_w_L1EM_25GeV_1030_b");
TH1D *photon_w_L1EM_25GeV_3050_b = (TH1D*) file->Get("photon_w_L1EM_25GeV_3050_b");
TH1D *photon_w_L1EM_25GeV_5080_b = (TH1D*) file->Get("photon_w_L1EM_25GeV_5080_b");

TH1D *photon_w_L1EM_35GeV_010_b = (TH1D*) file->Get("photon_w_L1EM_35GeV_010_b");
TH1D *photon_w_L1EM_35GeV_1030_b = (TH1D*) file->Get("photon_w_L1EM_35GeV_1030_b");
TH1D *photon_w_L1EM_35GeV_3050_b = (TH1D*) file->Get("photon_w_L1EM_35GeV_3050_b");
TH1D *photon_w_L1EM_35GeV_5080_b = (TH1D*) file->Get("photon_w_L1EM_35GeV_5080_b");

//----Endcaps
TH1D *all_photons_010_e = (TH1D*) file->Get("all_photons_010_e"); 
TH1D *all_photons_1030_e = (TH1D*) file->Get("all_photons_1030_e");
TH1D *all_photons_3050_e = (TH1D*) file->Get("all_photons_3050_e");
TH1D *all_photons_5080_e = (TH1D*) file->Get("all_photons_5080_e");

TH1D *photon_w_L1EM_10GeV_010_e = (TH1D*) file->Get("photon_w_L1EM_10GeV_010_e");
TH1D *photon_w_L1EM_10GeV_1030_e = (TH1D*) file->Get("photon_w_L1EM_10GeV_1030_e");
TH1D *photon_w_L1EM_10GeV_3050_e = (TH1D*) file->Get("photon_w_L1EM_10GeV_3050_e");
TH1D *photon_w_L1EM_10GeV_5080_e = (TH1D*) file->Get("photon_w_L1EM_10GeV_5080_e");

TH1D *photon_w_L1EM_15GeV_010_e = (TH1D*) file->Get("photon_w_L1EM_15GeV_010_e");
TH1D *photon_w_L1EM_15GeV_1030_e = (TH1D*) file->Get("photon_w_L1EM_15GeV_1030_e");
TH1D *photon_w_L1EM_15GeV_3050_e = (TH1D*) file->Get("photon_w_L1EM_15GeV_3050_e");
TH1D *photon_w_L1EM_15GeV_5080_e = (TH1D*) file->Get("photon_w_L1EM_15GeV_5080_e");

TH1D *photon_w_L1EM_25GeV_010_e = (TH1D*) file->Get("photon_w_L1EM_25GeV_010_e");
TH1D *photon_w_L1EM_25GeV_1030_e =(TH1D*) file->Get("photon_w_L1EM_25GeV_1030_e");
TH1D *photon_w_L1EM_25GeV_3050_e =(TH1D*) file->Get("photon_w_L1EM_25GeV_3050_e");
TH1D *photon_w_L1EM_25GeV_5080_e = (TH1D*) file->Get("photon_w_L1EM_25GeV_5080_e");
 
TH1D *photon_w_L1EM_35GeV_010_e = (TH1D*) file->Get("photon_w_L1EM_35GeV_010_e");
TH1D *photon_w_L1EM_35GeV_1030_e = (TH1D*) file->Get("photon_w_L1EM_35GeV_1030_e");
TH1D *photon_w_L1EM_35GeV_3050_e = (TH1D*) file->Get("photon_w_L1EM_35GeV_3050_e");
TH1D *photon_w_L1EM_35GeV_5080_e = (TH1D*) file->Get("photon_w_L1EM_35GeV_5080_e");

//---------------------------------------------------------------------------------//
//-------------------HLT Efficiencies Dependent on Centrality----------------------//
//---------------------------------------------------------------------------------//

//----Barrel & EndCaps
TH1D *photon_w_HLT_15GeV_010 = (TH1D*) file->Get("photon_w_HLT_15GeV_010");
TH1D *photon_w_HLT_15GeV_1030 = (TH1D*) file->Get("photon_w_HLT_15GeV_1030");
TH1D *photon_w_HLT_15GeV_3050 = (TH1D*) file->Get("photon_w_HLT_15GeV_3050");
TH1D *photon_w_HLT_15GeV_5080 = (TH1D*) file->Get("photon_w_HLT_15GeV_5080");

TH1D *photon_w_HLT_25GeV_010 = (TH1D*) file->Get("photon_w_HLT_25GeV_010");
TH1D *photon_w_HLT_25GeV_1030 = (TH1D*) file->Get("photon_w_HLT_25GeV_1030");
TH1D *photon_w_HLT_25GeV_3050 = (TH1D*) file->Get("photon_w_HLT_25GeV_3050");
TH1D *photon_w_HLT_25GeV_5080= (TH1D*) file->Get("photon_w_HLT_25GeV_5080");

TH1D *photon_w_HLT_35GeV_010 = (TH1D*) file->Get("photon_w_HLT_35GeV_010");
TH1D *photon_w_HLT_35GeV_1030 = (TH1D*) file->Get("photon_w_HLT_35GeV_1030");
TH1D *photon_w_HLT_35GeV_3050 = (TH1D*) file->Get("photon_w_HLT_35GeV_3050");
TH1D *photon_w_HLT_35GeV_5080 = (TH1D*) file->Get("photon_w_HLT_35GeV_5080");
  
//----Barrel
TH1D *photon_w_HLT_15GeV_010_b = (TH1D*) file->Get("photon_w_HLT_15GeV_010_b");
TH1D *photon_w_HLT_15GeV_1030_b = (TH1D*) file->Get("photon_w_HLT_15GeV_1030_b"); 
TH1D *photon_w_HLT_15GeV_3050_b = (TH1D*) file->Get("photon_w_HLT_15GeV_3050_b");
TH1D *photon_w_HLT_15GeV_5080_b = (TH1D*) file->Get("photon_w_HLT_15GeV_5080_b");

TH1D *photon_w_HLT_25GeV_010_b = (TH1D*) file->Get("photon_w_HLT_25GeV_010_b");
TH1D *photon_w_HLT_25GeV_1030_b = (TH1D*) file->Get("photon_w_HLT_25GeV_1030_b");
TH1D *photon_w_HLT_25GeV_3050_b = (TH1D*) file->Get("photon_w_HLT_25GeV_3050_b");
TH1D *photon_w_HLT_25GeV_5080_b = (TH1D*) file->Get("photon_w_HLT_25GeV_5080_b");

TH1D *photon_w_HLT_35GeV_010_b = (TH1D*) file->Get("photon_w_HLT_35GeV_010_b");
TH1D *photon_w_HLT_35GeV_1030_b = (TH1D*) file->Get("photon_w_HLT_35GeV_1030_b");
TH1D *photon_w_HLT_35GeV_3050_b = (TH1D*) file->Get("photon_w_HLT_35GeV_3050_b");
TH1D *photon_w_HLT_35GeV_5080_b = (TH1D*) file->Get("photon_w_HLT_35GeV_5080_b");
  

//----Endcaps
TH1D *photon_w_HLT_15GeV_010_e = (TH1D*) file->Get("photon_w_HLT_15GeV_010_e");
TH1D *photon_w_HLT_15GeV_1030_e = (TH1D*) file->Get("photon_w_HLT_15GeV_1030_e");
TH1D *photon_w_HLT_15GeV_3050_e = (TH1D*) file->Get("photon_w_HLT_15GeV_3050_e");
TH1D *photon_w_HLT_15GeV_5080_e = (TH1D*) file->Get("photon_w_HLT_15GeV_5080_e");

TH1D *photon_w_HLT_25GeV_010_e = (TH1D*) file->Get("photon_w_HLT_25GeV_010_e");
TH1D *photon_w_HLT_25GeV_1030_e = (TH1D*) file->Get("photon_w_HLT_25GeV_1030_e");
TH1D *photon_w_HLT_25GeV_3050_e = (TH1D*) file->Get("photon_w_HLT_25GeV_3050_e");
TH1D *photon_w_HLT_25GeV_5080_e = (TH1D*) file->Get("photon_w_HLT_25GeV_5080_e");

TH1D *photon_w_HLT_35GeV_010_e = (TH1D*) file->Get("photon_w_HLT_35GeV_010_e");
TH1D *photon_w_HLT_35GeV_1030_e = (TH1D*) file->Get("photon_w_HLT_35GeV_1030_e");
TH1D *photon_w_HLT_35GeV_3050_e = (TH1D*) file->Get("photon_w_HLT_35GeV_3050_e");
TH1D *photon_w_HLT_35GeV_5080_e  = (TH1D*) file->Get("photon_w_HLT_35GeV_5080_e");
//**************************************************************************************//

//-------Startinig Plotting Code HERE-------//

//------L1 Efficiencies
//------0-10% Centrality
//------Barrel & Endcaps

//------Creating Canvas & Legend
TCanvas *new_l1em_eff_centrality = new TCanvas("new_l1em_eff_centrality","new_l1em_eff_centrality",600,500);
TLegend *new_lg_l1em_cent = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_010_15GeV = new TGraphAsymmErrors();
L1EM_eff_010_15GeV->SetName("L1EM_eff_010_15GeV");                                                                                                                           
L1EM_eff_010_15GeV->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_010_15GeV->SetMarkerStyle(24);
L1EM_eff_010_15GeV->SetTitle("L1EM Efficiency");
L1EM_eff_010_15GeV->SetLineColor(kRed);
L1EM_eff_010_15GeV->SetMarkerColor(kRed);
L1EM_eff_010_15GeV->BayesDivide(photon_w_L1EM_15GeV_010,all_photons_010);
L1EM_eff_010_15GeV->SetMaximum(1.1);
L1EM_eff_010_15GeV->SetMinimum(0);

TGraphAsymmErrors *L1EM_eff_010_10GeV = new TGraphAsymmErrors();
L1EM_eff_010_10GeV->SetName("L1EM_eff_010_10GeV");
L1EM_eff_010_10GeV->SetMarkerStyle(27);
L1EM_eff_010_10GeV->SetMarkerColor(kCyan+1);
L1EM_eff_010_10GeV->SetLineColor(kCyan+1);
L1EM_eff_010_10GeV->BayesDivide(photon_w_L1EM_10GeV_010,all_photons_010);

TGraphAsymmErrors *L1EM_eff_010_25GeV = new TGraphAsymmErrors();
L1EM_eff_010_25GeV->SetName("L1EM_eff_010_25GeV");
L1EM_eff_010_25GeV->SetMarkerStyle(25);
L1EM_eff_010_25GeV->SetMarkerColor(kBlue);
L1EM_eff_010_25GeV->SetLineColor(kBlue);
L1EM_eff_010_25GeV->BayesDivide(photon_w_L1EM_25GeV_010,all_photons_010);

TGraphAsymmErrors *L1EM_eff_010_35GeV = new TGraphAsymmErrors();
L1EM_eff_010_35GeV->SetName("L1EM_eff_010_35GeV");
L1EM_eff_010_35GeV->SetMarkerStyle(26);
L1EM_eff_010_35GeV->SetMarkerColor(kMagenta);
L1EM_eff_010_35GeV->SetLineColor(kMagenta);
L1EM_eff_010_35GeV->BayesDivide(photon_w_L1EM_35GeV_010,all_photons_010);


TPaveText *l1em_photon_0 = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_0->SetTextSize(0.03);
l1em_photon_0->SetFillColor(0);
l1em_photon_0->SetBorderSize(0);
l1em_photon_0->SetShadowColor(0);
l1em_photon_0->SetTextAlign(21);
l1em_photon_0->AddText("");
l1em_photon_0->AddText("Pb+Pb");
l1em_photon_0->AddText("0-10%");
l1em_photon_0->AddText("Tight Photons");
l1em_photon_0->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_0->AddText("DeltaR < 0.15");
l1em_photon_0->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_010_10GeV = new TF1("fit_010_10GeV", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_010_15GeV = new TF1("fit_010_15GeV", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_010_25GeV = new TF1("fit_010_25GeV","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_010_35GeV = new TF1("fit_010_35GeV","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);

//----Guessing the fit parameters implemeted into our fit function
fit_010_10GeV->SetParameter(0,10.0);		
fit_010_10GeV->SetParameter(1,2.0);

fit_010_15GeV->SetParameter(0,15.0);		
fit_010_15GeV->SetParameter(1,2.0);

fit_010_25GeV->SetParameter(0,25.0);
fit_010_25GeV->SetParameter(1,2.0);

fit_010_35GeV->SetParameter(0,35.0);
fit_010_35GeV->SetParameter(1,2.0);


L1EM_eff_010_10GeV->Fit("fit_010_10GeV","M Q N","",0,photon_w_L1EM_10GeV_010->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_010->GetXaxis()->GetNbins()+1));
L1EM_eff_010_15GeV->Fit("fit_010_15GeV","M Q N","",0,photon_w_L1EM_15GeV_010->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_010->GetXaxis()->GetNbins()+1));
L1EM_eff_010_25GeV->Fit("fit_010_25GeV", "M Q N","",0,photon_w_L1EM_25GeV_010->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_010->GetXaxis()->GetNbins()+1));
L1EM_eff_010_35GeV->Fit("fit_010_35GeV","M Q N","",0,photon_w_L1EM_35GeV_010->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_010->GetXaxis()->GetNbins()+1));

fit_010_10GeV->SetLineStyle(2);
fit_010_10GeV->SetLineColor(kCyan+1);
fit_010_10GeV->SetMarkerColor(kCyan+1);

fit_010_15GeV->SetLineStyle(2);
fit_010_15GeV->SetLineColor(kRed);
fit_010_15GeV->SetMarkerColor(kRed);

fit_010_25GeV->SetLineStyle(2);
fit_010_25GeV->SetLineColor(kBlue);
fit_010_25GeV->SetMarkerColor(kBlue);

fit_010_35GeV->SetLineStyle(2);
fit_010_35GeV->SetLineColor(kMagenta);
fit_010_35GeV->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100 = new TLine(0.4,1,169,1);
line_100->SetLineColor(kBlack);
line_100->SetLineWidth(1);




L1EM_eff_010_15GeV->Draw("A P");

line_100->Draw("same");

L1EM_eff_010_10GeV->Draw("P same");
L1EM_eff_010_25GeV->Draw("P same");
L1EM_eff_010_35GeV->Draw("P same");

fit_010_10GeV->Draw("same");
fit_010_15GeV->Draw("same");
fit_010_25GeV->Draw("same");
fit_010_35GeV->Draw("same");


l1em_photon_0->Draw("same");

new_lg_l1em_cent->AddEntry("L1EM_eff_010_10GeV","L1 RoI p_{T} > 10 GeV","pe");
new_lg_l1em_cent->AddEntry("L1EM_eff_010_15GeV","L1 RoI p_{T} > 15 GeV","pe");
new_lg_l1em_cent->AddEntry("L1EM_eff_010_25GeV","L1 RoI p_{T} > 25 GeV","pe");
new_lg_l1em_cent->AddEntry("L1EM_eff_010_35GeV","L1 RoI p_{T} > 35 GeV","pe");



new_lg_l1em_cent->Draw("same");



//------10-30% Centrality
//------Barrel & Endcaps

//------Creating Canvas & Legend
TCanvas *new_l1em_eff_centrality_1030 = new TCanvas("new_l1em_eff_centrality_1030","new_l1em_eff_centrality_1030",600,500);
TLegend *new_lg_l1em_cent_1030 = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_1030_15GeV = new TGraphAsymmErrors();
L1EM_eff_1030_15GeV->SetName("L1EM_eff_1030_15GeV");                                                                                                                           
L1EM_eff_1030_15GeV->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_1030_15GeV->SetMarkerStyle(24);
L1EM_eff_1030_15GeV->SetTitle("L1EM Efficiency");
L1EM_eff_1030_15GeV->SetLineColor(kRed);
L1EM_eff_1030_15GeV->SetMarkerColor(kRed);
L1EM_eff_1030_15GeV->BayesDivide(photon_w_L1EM_15GeV_1030,all_photons_1030);
L1EM_eff_1030_15GeV->SetMaximum(1.1);
L1EM_eff_1030_15GeV->SetMinimum(0);

TGraphAsymmErrors *L1EM_eff_1030_10GeV = new TGraphAsymmErrors();
L1EM_eff_1030_10GeV->SetName("L1EM_eff_1030_10GeV");
L1EM_eff_1030_10GeV->SetMarkerStyle(27);
L1EM_eff_1030_10GeV->SetMarkerColor(kCyan+1);
L1EM_eff_1030_10GeV->SetLineColor(kCyan+1);
L1EM_eff_1030_10GeV->BayesDivide(photon_w_L1EM_10GeV_1030,all_photons_1030);

TGraphAsymmErrors *L1EM_eff_1030_25GeV = new TGraphAsymmErrors();
L1EM_eff_1030_25GeV->SetName("L1EM_eff_1030_25GeV");
L1EM_eff_1030_25GeV->SetMarkerStyle(25);
L1EM_eff_1030_25GeV->SetMarkerColor(kBlue);
L1EM_eff_1030_25GeV->SetLineColor(kBlue);
L1EM_eff_1030_25GeV->BayesDivide(photon_w_L1EM_25GeV_1030,all_photons_1030);

TGraphAsymmErrors *L1EM_eff_1030_35GeV = new TGraphAsymmErrors();
L1EM_eff_1030_35GeV->SetName("L1EM_eff_1030_35GeV");
L1EM_eff_1030_35GeV->SetMarkerStyle(26);
L1EM_eff_1030_35GeV->SetMarkerColor(kMagenta);
L1EM_eff_1030_35GeV->SetLineColor(kMagenta);
L1EM_eff_1030_35GeV->BayesDivide(photon_w_L1EM_35GeV_1030,all_photons_1030);


TPaveText *l1em_photon_1030 = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_1030->SetTextSize(0.03);
l1em_photon_1030->SetFillColor(0);
l1em_photon_1030->SetBorderSize(0);
l1em_photon_1030->SetShadowColor(0);
l1em_photon_1030->SetTextAlign(21);
l1em_photon_1030->AddText("");
l1em_photon_1030->AddText("Pb+Pb");
l1em_photon_1030->AddText("10-30%");
l1em_photon_1030->AddText("Tight Photons");
l1em_photon_1030->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_1030->AddText("DeltaR < 0.15");
l1em_photon_1030->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_1030_10GeV = new TF1("fit_1030_10GeV", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_10GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_1030->GetXaxis()->GetNbins()+1));
TF1* fit_1030_15GeV = new TF1("fit_1030_15GeV", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_15GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_1030->GetXaxis()->GetNbins()+1));
TF1* fit_1030_25GeV = new TF1("fit_1030_25GeV","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_25GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_1030->GetXaxis()->GetNbins()+1));
TF1* fit_1030_35GeV = new TF1("fit_1030_35GeV","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_35GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_1030->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_1030_10GeV->SetParameter(0,10.0);		
fit_1030_10GeV->SetParameter(1,2.0);

fit_1030_15GeV->SetParameter(0,15.0);		
fit_1030_15GeV->SetParameter(1,2.0);

fit_1030_25GeV->SetParameter(0,25.0);
fit_1030_25GeV->SetParameter(1,2.0);

fit_1030_35GeV->SetParameter(0,35.0);
fit_1030_35GeV->SetParameter(1,2.0);


L1EM_eff_1030_10GeV->Fit("fit_1030_10GeV","M Q N","",0,photon_w_L1EM_10GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_1030->GetXaxis()->GetNbins()+1));
L1EM_eff_1030_15GeV->Fit("fit_1030_15GeV","M Q N","",0,photon_w_L1EM_15GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_1030->GetXaxis()->GetNbins()+1));
L1EM_eff_1030_25GeV->Fit("fit_1030_25GeV", "M Q N","",0,photon_w_L1EM_25GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_1030->GetXaxis()->GetNbins()+1));
L1EM_eff_1030_35GeV->Fit("fit_1030_35GeV","M Q N","",0,photon_w_L1EM_35GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_1030->GetXaxis()->GetNbins()+1));


fit_1030_10GeV->SetLineStyle(2);
fit_1030_10GeV->SetLineColor(kCyan+1);
fit_1030_10GeV->SetMarkerColor(kCyan+1);

fit_1030_15GeV->SetLineStyle(2);
fit_1030_15GeV->SetLineColor(kRed);
fit_1030_15GeV->SetMarkerColor(kRed);

fit_1030_25GeV->SetLineStyle(2);
fit_1030_25GeV->SetLineColor(kBlue);
fit_1030_25GeV->SetMarkerColor(kBlue);

fit_1030_35GeV->SetLineStyle(2);
fit_1030_35GeV->SetLineColor(kMagenta);
fit_1030_35GeV->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_1 = new TLine(0.5,1,165,1);
line_100_1->SetLineColor(kBlack);
line_100_1->SetLineWidth(1);




L1EM_eff_1030_15GeV->Draw("A P");

line_100_1->Draw("same");

L1EM_eff_1030_10GeV->Draw("P same");
L1EM_eff_1030_25GeV->Draw("P same");
L1EM_eff_1030_35GeV->Draw("P same");

fit_1030_10GeV->Draw("same");
fit_1030_15GeV->Draw("same");
fit_1030_25GeV->Draw("same");
fit_1030_35GeV->Draw("same");


l1em_photon_1030->Draw("same");

new_lg_l1em_cent_1030->AddEntry("L1EM_eff_1030_10GeV","L1 RoI p_{T} > 10 GeV","pe");
new_lg_l1em_cent_1030->AddEntry("L1EM_eff_1030_15GeV","L1 RoI p_{T} > 15 GeV","pe");
new_lg_l1em_cent_1030->AddEntry("L1EM_eff_1030_25GeV","L1 RoI p_{T} > 25 GeV","pe");
new_lg_l1em_cent_1030->AddEntry("L1EM_eff_1030_35GeV","L1 RoI p_{T} > 35 GeV","pe");



new_lg_l1em_cent_1030->Draw("same");

//------30-50% Centrality
//------Barrel & Endcaps

//------Creating Canvas & Legend
TCanvas *new_l1em_eff_centrality_3050 = new TCanvas("new_l1em_eff_centrality_3050","new_l1em_eff_centrality_3050",600,500);
TLegend *new_lg_l1em_cent_3050 = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_3050_15GeV = new TGraphAsymmErrors();
L1EM_eff_3050_15GeV->SetName("L1EM_eff_3050_15GeV");                                                                                                                           
L1EM_eff_3050_15GeV->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_3050_15GeV->SetMarkerStyle(24);
L1EM_eff_3050_15GeV->SetTitle("L1EM Efficiency");
L1EM_eff_3050_15GeV->SetLineColor(kRed);
L1EM_eff_3050_15GeV->SetMarkerColor(kRed);
L1EM_eff_3050_15GeV->BayesDivide(photon_w_L1EM_15GeV_3050,all_photons_3050);
L1EM_eff_3050_15GeV->SetMaximum(1.1);
L1EM_eff_3050_15GeV->SetMinimum(0);

TGraphAsymmErrors *L1EM_eff_3050_10GeV = new TGraphAsymmErrors();
L1EM_eff_3050_10GeV->SetName("L1EM_eff_3050_10GeV");
L1EM_eff_3050_10GeV->SetMarkerStyle(27);
L1EM_eff_3050_10GeV->SetMarkerColor(kCyan+1);
L1EM_eff_3050_10GeV->SetLineColor(kCyan+1);
L1EM_eff_3050_10GeV->BayesDivide(photon_w_L1EM_10GeV_3050,all_photons_3050);

TGraphAsymmErrors *L1EM_eff_3050_25GeV = new TGraphAsymmErrors();
L1EM_eff_3050_25GeV->SetName("L1EM_eff_3050_25GeV");
L1EM_eff_3050_25GeV->SetMarkerStyle(25);
L1EM_eff_3050_25GeV->SetMarkerColor(kBlue);
L1EM_eff_3050_25GeV->SetLineColor(kBlue);
L1EM_eff_3050_25GeV->BayesDivide(photon_w_L1EM_25GeV_3050,all_photons_3050);

TGraphAsymmErrors *L1EM_eff_3050_35GeV = new TGraphAsymmErrors();
L1EM_eff_3050_35GeV->SetName("L1EM_eff_3050_35GeV");
L1EM_eff_3050_35GeV->SetMarkerStyle(26);
L1EM_eff_3050_35GeV->SetMarkerColor(kMagenta);
L1EM_eff_3050_35GeV->SetLineColor(kMagenta);
L1EM_eff_3050_35GeV->BayesDivide(photon_w_L1EM_35GeV_3050,all_photons_3050);


TPaveText *l1em_photon_3050 = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_3050->SetTextSize(0.03);
l1em_photon_3050->SetFillColor(0);
l1em_photon_3050->SetBorderSize(0);
l1em_photon_3050->SetShadowColor(0);
l1em_photon_3050->SetTextAlign(21);
l1em_photon_3050->AddText("");
l1em_photon_3050->AddText("Pb+Pb");
l1em_photon_3050->AddText("30-50%");
l1em_photon_3050->AddText("Tight Photons");
l1em_photon_3050->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_3050->AddText("DeltaR < 0.15");
l1em_photon_3050->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_3050_10GeV = new TF1("fit_3050_10GeV", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_10GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_10GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_3050->GetXaxis()->GetNbins()+1));
TF1* fit_3050_15GeV = new TF1("fit_3050_15GeV", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_15GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_15GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_3050->GetXaxis()->GetNbins()+1));
TF1* fit_3050_25GeV = new TF1("fit_3050_25GeV","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_25GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_25GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_3050->GetXaxis()->GetNbins()+1));
TF1* fit_3050_35GeV = new TF1("fit_3050_35GeV","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_35GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_35GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_3050->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_3050_10GeV->SetParameter(0,10.0);		
fit_3050_10GeV->SetParameter(1,2.0);

fit_3050_15GeV->SetParameter(0,15.0);		
fit_3050_15GeV->SetParameter(1,2.0);

fit_3050_25GeV->SetParameter(0,25.0);
fit_3050_25GeV->SetParameter(1,2.0);

fit_3050_35GeV->SetParameter(0,35.0);
fit_3050_35GeV->SetParameter(1,2.0);


L1EM_eff_3050_10GeV->Fit("fit_3050_10GeV","M Q N","",photon_w_L1EM_10GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_10GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_3050->GetXaxis()->GetNbins()+1));
L1EM_eff_3050_15GeV->Fit("fit_3050_15GeV","M Q N","",photon_w_L1EM_15GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_15GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_3050->GetXaxis()->GetNbins()+1));
L1EM_eff_3050_25GeV->Fit("fit_3050_25GeV", "M Q N","",photon_w_L1EM_25GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_25GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_3050->GetXaxis()->GetNbins()+1));
L1EM_eff_3050_35GeV->Fit("fit_3050_35GeV","M Q N","",photon_w_L1EM_35GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_35GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_3050->GetXaxis()->GetNbins()+1));

fit_3050_10GeV->SetLineStyle(2);
fit_3050_10GeV->SetLineColor(kCyan+1);
fit_3050_10GeV->SetMarkerColor(kCyan+1);

fit_3050_15GeV->SetLineStyle(2);
fit_3050_15GeV->SetLineColor(kRed);
fit_3050_15GeV->SetMarkerColor(kRed);

fit_3050_25GeV->SetLineStyle(2);
fit_3050_25GeV->SetLineColor(kBlue);
fit_3050_25GeV->SetMarkerColor(kBlue);

fit_3050_35GeV->SetLineStyle(2);
fit_3050_35GeV->SetLineColor(kMagenta);
fit_3050_35GeV->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_2 = new TLine(2.0,1,153,1);
line_100_2->SetLineColor(kBlack);
line_100_2->SetLineWidth(1);




L1EM_eff_3050_15GeV->Draw("A P");
line_100_2->Draw("same");

L1EM_eff_3050_10GeV->Draw("P same");
L1EM_eff_3050_25GeV->Draw("P same");
L1EM_eff_3050_35GeV->Draw("P same");

fit_3050_10GeV->Draw("same");
fit_3050_15GeV->Draw("same");
fit_3050_25GeV->Draw("same");
fit_3050_35GeV->Draw("same");


l1em_photon_3050->Draw("same");

new_lg_l1em_cent_3050->AddEntry("L1EM_eff_3050_10GeV","L1 RoI p_{T} > 10 GeV","pe");
new_lg_l1em_cent_3050->AddEntry("L1EM_eff_3050_15GeV","L1 RoI p_{T} > 15 GeV","pe");
new_lg_l1em_cent_3050->AddEntry("L1EM_eff_3050_25GeV","L1 RoI p_{T} > 25 GeV","pe");
new_lg_l1em_cent_3050->AddEntry("L1EM_eff_3050_35GeV","L1 RoI p_{T} > 35 GeV","pe");



new_lg_l1em_cent_3050->Draw("same");


//------50-80% Centrality
//------Barrel & Endcaps

//------Creating Canvas & Legend
TCanvas *new_l1em_eff_centrality_5080 = new TCanvas("new_l1em_eff_centrality_5080","new_l1em_eff_centrality_5080",600,500);
TLegend *new_lg_l1em_cent_5080 = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_5080_15GeV = new TGraphAsymmErrors();
L1EM_eff_5080_15GeV->SetName("L1EM_eff_5080_15GeV");                                                                                                                           
L1EM_eff_5080_15GeV->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_5080_15GeV->SetMarkerStyle(24);
L1EM_eff_5080_15GeV->SetTitle("L1EM Efficiency");
L1EM_eff_5080_15GeV->SetLineColor(kRed);
L1EM_eff_5080_15GeV->SetMarkerColor(kRed);
L1EM_eff_5080_15GeV->BayesDivide(photon_w_L1EM_15GeV_5080,all_photons_5080);
L1EM_eff_5080_15GeV->SetMaximum(1.1);
L1EM_eff_5080_15GeV->SetMinimum(0);

TGraphAsymmErrors *L1EM_eff_5080_10GeV = new TGraphAsymmErrors();
L1EM_eff_5080_10GeV->SetName("L1EM_eff_5080_10GeV");
L1EM_eff_5080_10GeV->SetMarkerStyle(27);
L1EM_eff_5080_10GeV->SetMarkerColor(kCyan+1);
L1EM_eff_5080_10GeV->SetLineColor(kCyan+1);
L1EM_eff_5080_10GeV->BayesDivide(photon_w_L1EM_10GeV_5080,all_photons_5080);

TGraphAsymmErrors *L1EM_eff_5080_25GeV = new TGraphAsymmErrors();
L1EM_eff_5080_25GeV->SetName("L1EM_eff_5080_25GeV");
L1EM_eff_5080_25GeV->SetMarkerStyle(25);
L1EM_eff_5080_25GeV->SetMarkerColor(kBlue);
L1EM_eff_5080_25GeV->SetLineColor(kBlue);
L1EM_eff_5080_25GeV->BayesDivide(photon_w_L1EM_25GeV_5080,all_photons_5080);

TGraphAsymmErrors *L1EM_eff_5080_35GeV = new TGraphAsymmErrors();
L1EM_eff_5080_35GeV->SetName("L1EM_eff_5080_35GeV");
L1EM_eff_5080_35GeV->SetMarkerStyle(26);
L1EM_eff_5080_35GeV->SetMarkerColor(kMagenta);
L1EM_eff_5080_35GeV->SetLineColor(kMagenta);
L1EM_eff_5080_35GeV->BayesDivide(photon_w_L1EM_35GeV_5080,all_photons_5080);


TPaveText *l1em_photon_5080 = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_5080->SetTextSize(0.03);
l1em_photon_5080->SetFillColor(0);
l1em_photon_5080->SetBorderSize(0);
l1em_photon_5080->SetShadowColor(0);
l1em_photon_5080->SetTextAlign(21);
l1em_photon_5080->AddText("");
l1em_photon_5080->AddText("Pb+Pb");
l1em_photon_5080->AddText("50-80%");
l1em_photon_5080->AddText("Tight Photons");
l1em_photon_5080->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_5080->AddText("DeltaR < 0.15");
l1em_photon_5080->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_5080_10GeV = new TF1("fit_5080_10GeV", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_10GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_10GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_5080->GetXaxis()->GetNbins()+1));
TF1* fit_5080_15GeV = new TF1("fit_5080_15GeV", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_15GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_15GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_5080->GetXaxis()->GetNbins()+1));
TF1* fit_5080_25GeV = new TF1("fit_5080_25GeV","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_25GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_25GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_5080->GetXaxis()->GetNbins()+1));
TF1* fit_5080_35GeV = new TF1("fit_5080_35GeV","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_35GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_35GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_5080->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_5080_10GeV->SetParameter(0,10.0);		
fit_5080_10GeV->SetParameter(1,2.0);

fit_5080_15GeV->SetParameter(0,15.0);		
fit_5080_15GeV->SetParameter(1,2.0);

fit_5080_25GeV->SetParameter(0,25.0);
fit_5080_25GeV->SetParameter(1,2.0);

fit_5080_35GeV->SetParameter(0,35.0);
fit_5080_35GeV->SetParameter(1,2.0);


L1EM_eff_5080_10GeV->Fit("fit_5080_10GeV","M Q N","",photon_w_L1EM_10GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_10GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_5080->GetXaxis()->GetNbins()+1));
L1EM_eff_5080_15GeV->Fit("fit_5080_15GeV","M Q N","",photon_w_L1EM_15GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_15GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_5080->GetXaxis()->GetNbins()+1));
L1EM_eff_5080_25GeV->Fit("fit_5080_25GeV", "M Q N","",photon_w_L1EM_25GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_25GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_5080->GetXaxis()->GetNbins()+1));
L1EM_eff_5080_35GeV->Fit("fit_5080_35GeV","M Q N","",photon_w_L1EM_35GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_35GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_5080->GetXaxis()->GetNbins()+1));

fit_5080_10GeV->SetLineStyle(2);
fit_5080_10GeV->SetLineColor(kCyan+1);
fit_5080_10GeV->SetMarkerColor(kCyan+1);

fit_5080_15GeV->SetLineStyle(2);
fit_5080_15GeV->SetLineColor(kRed);
fit_5080_15GeV->SetMarkerColor(kRed);

fit_5080_25GeV->SetLineStyle(2);
fit_5080_25GeV->SetLineColor(kBlue);
fit_5080_25GeV->SetMarkerColor(kBlue);

fit_5080_35GeV->SetLineStyle(2);
fit_5080_35GeV->SetLineColor(kMagenta);
fit_5080_35GeV->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_3 = new TLine(2.0,1,200,1);
line_100_3->SetLineColor(kBlack);
line_100_3->SetLineWidth(1);




L1EM_eff_5080_15GeV->Draw("A P");
line_100_3->Draw("same");

L1EM_eff_5080_10GeV->Draw("P same");
L1EM_eff_5080_25GeV->Draw("P same");
L1EM_eff_5080_35GeV->Draw("P same");

fit_5080_10GeV->Draw("same");
fit_5080_15GeV->Draw("same");
fit_5080_25GeV->Draw("same");
fit_5080_35GeV->Draw("same");


l1em_photon_5080->Draw("same");

new_lg_l1em_cent_5080->AddEntry("L1EM_eff_5080_10GeV","L1 RoI p_{T} > 10 GeV","pe");
new_lg_l1em_cent_5080->AddEntry("L1EM_eff_5080_15GeV","L1 RoI p_{T} > 15 GeV","pe");
new_lg_l1em_cent_5080->AddEntry("L1EM_eff_5080_25GeV","L1 RoI p_{T} > 25 GeV","pe");
new_lg_l1em_cent_5080->AddEntry("L1EM_eff_5080_35GeV","L1 RoI p_{T} > 35 GeV","pe");



new_lg_l1em_cent_5080->Draw("same");

//-------Summary Plots
//-------Barrel & Endcaps
//-------This is for the parameter that tells us what is our pT at 50% efficiency
TCanvas *summary_new = new TCanvas("summary_new","summary_new",600,500);
TLegend *new_legend = new TLegend(0.6,0.5,0.85,0.6);

//-----0-10%
TH1D *summary_010_0 = new TH1D("summary_010_0","",4,0,40);
summary_010_0->SetMarkerStyle(22);
summary_010_0->SetMarkerColor(kBlack);
summary_010_0->SetLineColor(kBlack);
summary_010_0->SetName("summary_010_0");
summary_010_0->GetXaxis()->SetTitle("L1 p_{T} [GeV]");
summary_010_0->GetYaxis()->SetTitle("Value at 50 Percent Efficiency [GeV]");

summary_010_0->SetBinContent(1,fit_010_10GeV->GetParameter(0));
summary_010_0->SetBinError(1,fit_010_10GeV->GetParError(0));
summary_010_0->SetBinContent(2,fit_010_15GeV->GetParameter(0));
summary_010_0->SetBinError(2,fit_010_15GeV->GetParError(0));
summary_010_0->SetBinContent(3,fit_010_25GeV->GetParameter(0));
summary_010_0->SetBinError(3,fit_010_25GeV->GetParError(0));
summary_010_0->SetBinContent(4,fit_010_35GeV->GetParameter(0));
summary_010_0->SetBinError(4,fit_010_35GeV->GetParError(0));

//------10-30%
TH1D *summary_1030_0 = new TH1D("summary_1030_0","",4,0,40);
summary_1030_0->SetMarkerStyle(21);
summary_1030_0->SetMarkerColor(kRed);
summary_1030_0->SetLineColor(kRed);
summary_1030_0->SetName("summary_1030_0");

summary_1030_0->SetBinContent(1,fit_1030_10GeV->GetParameter(0));
summary_1030_0->SetBinError(1,fit_1030_10GeV->GetParError(0));
summary_1030_0->SetBinContent(2,fit_1030_15GeV->GetParameter(0));
summary_1030_0->SetBinError(2,fit_1030_15GeV->GetParError(0));
summary_1030_0->SetBinContent(3,fit_1030_25GeV->GetParameter(0));
summary_1030_0->SetBinError(3,fit_1030_25GeV->GetParError(0));
summary_1030_0->SetBinContent(4,fit_1030_35GeV->GetParameter(0));
summary_1030_0->SetBinError(4,fit_1030_35GeV->GetParError(0));

//------30-50%
TH1D *summary_3050_0 = new TH1D("summary_3050_0","",4,0,40);
summary_3050_0->SetMarkerStyle(33);
summary_3050_0->SetMarkerColor(kMagenta);
summary_3050_0->SetLineColor(kMagenta);
summary_3050_0->SetName("summary_3050_0");


summary_3050_0->SetBinContent(1,fit_3050_10GeV->GetParameter(0));
summary_3050_0->SetBinError(1,fit_3050_10GeV->GetParError(0));
summary_3050_0->SetBinContent(2,fit_3050_15GeV->GetParameter(0));
summary_3050_0->SetBinError(2,fit_3050_15GeV->GetParError(0));
summary_3050_0->SetBinContent(3,fit_3050_25GeV->GetParameter(0));
summary_3050_0->SetBinError(3,fit_3050_25GeV->GetParError(0));
summary_3050_0->SetBinContent(4,fit_3050_35GeV->GetParameter(0));
summary_3050_0->SetBinError(4,fit_3050_35GeV->GetParError(0));

//------50-80%
TH1D *summary_5080_0 = new TH1D("summary_5080_0","",4,0,40);
summary_5080_0->SetMarkerStyle(20);
summary_5080_0->SetMarkerColor(kBlue);
summary_5080_0->SetLineColor(kBlue);
summary_5080_0->SetName("summary_5080_0");

summary_5080_0->SetBinContent(1,fit_5080_10GeV->GetParameter(0));
summary_5080_0->SetBinError(1,fit_5080_10GeV->GetParError(0));
summary_5080_0->SetBinContent(2,fit_5080_15GeV->GetParameter(0));
summary_5080_0->SetBinError(2,fit_5080_15GeV->GetParError(0));
summary_5080_0->SetBinContent(3,fit_5080_25GeV->GetParameter(0));
summary_5080_0->SetBinError(3,fit_5080_25GeV->GetParError(0));
summary_5080_0->SetBinContent(4,fit_5080_35GeV->GetParameter(0));
summary_5080_0->SetBinError(4,fit_5080_35GeV->GetParError(0));


TPaveText *l1_summary_0 = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1_summary_0->SetTextSize(0.03);
l1_summary_0->SetFillColor(0);
l1_summary_0->SetBorderSize(0);
l1_summary_0->SetShadowColor(0);
l1_summary_0->SetTextAlign(21);
l1_summary_0->AddText("");
l1_summary_0->AddText("Pb+Pb");
l1_summary_0->AddText("Tight Photons");
l1_summary_0->AddText("p^{#gamma}_{T} > 15 GeV");
l1_summary_0->AddText("HLT Loose Photons");
l1_summary_0->AddText("DeltaR < 0.15");
l1_summary_0->SetTextSize(0.036);




summary_010_0->Draw();
summary_1030_0->Draw("same");
summary_3050_0->Draw("same");
summary_5080_0->Draw("same");

l1_summary_0->Draw("same");


new_legend->AddEntry("summary_010_0","0-10%","pe");
new_legend->AddEntry("summary_1030_0","10-30%","pe");
new_legend->AddEntry("summary_3050_0","30-50%","pe");
new_legend->AddEntry("summary_5080_0","50-80%","pe");
new_legend->Draw("same");

//-------Summary Plots For The Width Parameter
TCanvas *summary_new_2 = new TCanvas("summary_new_2","summary_new_2",600,500);
TLegend *new_legend_2 = new TLegend(0.6,0.5,0.85,0.6);

//-----0-10%
TH1D *summary_010_1 = new TH1D("summary_010_1","",4,0,40);
summary_010_1->SetMarkerStyle(22);
summary_010_1->SetMarkerColor(kBlack);
summary_010_1->SetLineColor(kBlack);
summary_010_1->SetName("summary_010_1");
summary_010_1->GetXaxis()->SetTitle("L1 p_{T} [GeV]");
summary_010_1->GetYaxis()->SetTitle("Width");

summary_010_1->SetBinContent(1,fit_1030_10GeV->GetParameter(1));
summary_010_1->SetBinError(1,fit_1030_10GeV->GetParError(1));
summary_010_1->SetBinContent(2,fit_010_15GeV->GetParameter(1));
summary_010_1->SetBinError(2,fit_010_15GeV->GetParError(1));
summary_010_1->SetBinContent(3,fit_010_25GeV->GetParameter(1));
summary_010_1->SetBinError(3,fit_010_25GeV->GetParError(1));
summary_010_1->SetBinContent(4,fit_010_35GeV->GetParameter(1));
summary_010_1->SetBinError(4,fit_010_35GeV->GetParError(1));

//------10-30%
TH1D *summary_1030_1 = new TH1D("summary_1030_1","",4,0,40);
summary_1030_1->SetMarkerStyle(21);
summary_1030_1->SetMarkerColor(kRed);
summary_1030_1->SetLineColor(kRed);
summary_1030_1->SetName("summary_1030_1");

summary_1030_1->SetBinContent(1,fit_1030_10GeV->GetParameter(1));
summary_1030_1->SetBinError(1,fit_1030_10GeV->GetParError(1));
summary_1030_1->SetBinContent(2,fit_1030_15GeV->GetParameter(1));
summary_1030_1->SetBinError(2,fit_1030_15GeV->GetParError(1));
summary_1030_1->SetBinContent(3,fit_1030_25GeV->GetParameter(1));
summary_1030_1->SetBinError(3,fit_1030_25GeV->GetParError(1));
summary_1030_1->SetBinContent(4,fit_1030_35GeV->GetParameter(1));
summary_1030_1->SetBinError(4,fit_1030_35GeV->GetParError(1));

//------30-50%
TH1D *summary_3050_1 = new TH1D("summary_3050_1","",4,0,40);
summary_3050_1->SetMarkerStyle(33);
summary_3050_1->SetMarkerColor(kMagenta);
summary_3050_1->SetLineColor(kMagenta);
summary_3050_1->SetName("summary_3050_1");


summary_3050_1->SetBinContent(1,fit_3050_10GeV->GetParameter(1));
summary_3050_1->SetBinError(1,fit_3050_10GeV->GetParError(1));
summary_3050_1->SetBinContent(2,fit_3050_15GeV->GetParameter(1));
summary_3050_1->SetBinError(2,fit_3050_15GeV->GetParError(1));
summary_3050_1->SetBinContent(3,fit_3050_25GeV->GetParameter(1));
summary_3050_1->SetBinError(3,fit_3050_25GeV->GetParError(1));
summary_3050_1->SetBinContent(4,fit_3050_35GeV->GetParameter(1));
summary_3050_1->SetBinError(4,fit_3050_35GeV->GetParError(1));

//------50-80%
TH1D *summary_5080_1 = new TH1D("summary_5080_1","",4,0,40);
summary_5080_1->SetMarkerStyle(20);
summary_5080_1->SetMarkerColor(kBlue);
summary_5080_1->SetLineColor(kBlue);
summary_5080_1->SetName("summary_5080_1");

summary_5080_1->SetBinContent(1,fit_5080_10GeV->GetParameter(1));
summary_5080_1->SetBinError(1,fit_5080_10GeV->GetParError(1));
summary_5080_1->SetBinContent(2,fit_5080_15GeV->GetParameter(1));
summary_5080_1->SetBinError(2,fit_5080_15GeV->GetParError(1));
summary_5080_1->SetBinContent(3,fit_5080_25GeV->GetParameter(1));
summary_5080_1->SetBinError(3,fit_5080_25GeV->GetParError(1));
summary_5080_1->SetBinContent(4,fit_5080_35GeV->GetParameter(1));
summary_5080_1->SetBinError(4,fit_5080_35GeV->GetParError(1));


TPaveText *l1_summary_1 = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1_summary_1->SetTextSize(0.03);
l1_summary_1->SetFillColor(0);
l1_summary_1->SetBorderSize(0);
l1_summary_1->SetShadowColor(0);
l1_summary_1->SetTextAlign(21);
l1_summary_1->AddText("");
l1_summary_1->AddText("Pb+Pb");
l1_summary_1->AddText("Tight Photons");
l1_summary_1->AddText("p^{#gamma}_{T} > 15 GeV");
l1_summary_1->AddText("HLT Loose Photons");
l1_summary_1->AddText("DeltaR < 0.15");
l1_summary_1->SetTextSize(0.036);




summary_010_1->Draw();
summary_1030_1->Draw("same");
summary_3050_1->Draw("same");
summary_5080_1->Draw("same");

l1_summary_1->Draw("same");


new_legend_2->AddEntry("summary_010_1","0-10%","pe");
new_legend_2->AddEntry("summary_1030_1","10-30%","pe");
new_legend_2->AddEntry("summary_3050_1","30-50%","pe");
new_legend_2->AddEntry("summary_5080_1","50-80%","pe");
new_legend_2->Draw("same");


//------This summary plot let's us know when we have reached 99% efficiency 
//----0-10% Centrality 
TCanvas *canv_summary_99eff = new TCanvas("canv_summary_99eff","canv_summary_99eff",600,500);
TLegend *summary_99eff_leg = new TLegend(0.6,0.5,0.85,0.6);

TH1D *summary_010_99eff = new TH1D("summary_010_99eff","",3,10,40);
summary_010_99eff->SetMarkerStyle(22);
summary_010_99eff->SetMarkerColor(kBlack);
summary_010_99eff->SetLineColor(kBlack);
summary_010_99eff->SetName("summary_010_99eff");

summary_010_99eff->GetXaxis()->SetTitle("L1 p_{T} [GeV]");
summary_010_99eff->GetYaxis()->SetTitle("Offline pT Value At 99 Percent Efficiency");

//summary_010_1->SetBinContent(1,fit_010_15GeV->GetParameter(0) + (fit_010_15GeV->GetParameter(1)*TMath::ErfInverse(0.99)));
cout << endl << endl;
cout << "This is the pT value of your offline photon at 99 Percent efficiency: " << fit_010_15GeV->GetParameter(0) + (fit_010_15GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << endl;
cout << endl << endl;


summary_010_99eff->SetBinContent(1,fit_010_15GeV->GetParameter(0) + (fit_010_15GeV->GetParameter(1)*TMath::ErfInverse(0.99)));
summary_010_99eff->SetBinContent(2,fit_010_25GeV->GetParameter(0) + (fit_010_25GeV->GetParameter(1)*TMath::ErfInverse(0.99)));
summary_010_99eff->SetBinContent(3,fit_010_35GeV->GetParameter(0) + (fit_010_25GeV->GetParameter(1)*TMath::ErfInverse(0.99)));

//---10-30% Centrality 
TH1D *summary_1030_99eff = new TH1D("summary_1030_99eff","",3,10,40);
summary_1030_99eff->SetMarkerStyle(21);
summary_1030_99eff->SetMarkerColor(kRed);
summary_1030_99eff->SetLineColor(kRed);
summary_1030_99eff->SetName("summary_1030_99eff");

summary_1030_99eff->SetBinContent(1,fit_1030_15GeV->GetParameter(0) + (fit_1030_15GeV->GetParameter(1)*TMath::ErfInverse(0.99)));
summary_1030_99eff->SetBinContent(2,fit_1030_25GeV->GetParameter(0) + (fit_1030_25GeV->GetParameter(1)*TMath::ErfInverse(0.99)));
summary_1030_99eff->SetBinContent(2,fit_1030_35GeV->GetParameter(0) + (fit_1030_35GeV->GetParameter(1)*TMath::ErfInverse(0.99)));

//----30-50% Centrality
TH1D *summary_3050_99eff = new TH1D("summary_3050_99eff","",3,10,40);
summary_3050_99eff->SetMarkerStyle(33);
summary_3050_99eff->SetMarkerColor(kMagenta);
summary_3050_99eff->SetLineColor(kMagenta);
summary_3050_99eff->SetName("summary_3050_99eff");

summary_3050_99eff->SetBinContent(1,fit_3050_15GeV->GetParameter(0) + (fit_3050_15GeV->GetParameter(1)*TMath::ErfInverse(0.99)));
summary_3050_99eff->SetBinContent(2,fit_3050_25GeV->GetParameter(0) + (fit_3050_25GeV->GetParameter(1)*TMath::ErfInverse(0.99)));
summary_3050_99eff->SetBinContent(2,fit_3050_35GeV->GetParameter(0) + (fit_3050_35GeV->GetParameter(1)*TMath::ErfInverse(0.99)));

//----50-80% Centrality
TH1D *summary_5080_99eff = new TH1D("summary_5080_99eff","",3,10,40);
summary_5080_99eff->SetMarkerStyle(20);
summary_5080_99eff->SetLineColor(kBlue);
summary_5080_99eff->SetMarkerColor(kBlue);
summary_5080_99eff->SetName("summary_5080_99eff");

summary_010_99eff->Draw("P");
summary_1030_99eff->Draw("P same");
summary_3050_99eff->Draw("P same");
summary_5080_99eff->Draw("P same");

summary_5080_99eff->SetBinContent(1,fit_5080_15GeV->GetParameter(0) + (fit_5080_15GeV->GetParameter(1)*TMath::ErfInverse(0.99)));
summary_5080_99eff->SetBinContent(2,fit_5080_25GeV->GetParameter(0) + (fit_5080_25GeV->GetParameter(1)*TMath::ErfInverse(0.99)));
summary_5080_99eff->SetBinContent(3,fit_5080_35GeV->GetParameter(0) + (fit_5080_35GeV->GetParameter(1)*TMath::ErfInverse(0.99)));

summary_99eff_leg->AddEntry("summary_010_99eff","0-10%","p");
summary_99eff_leg->AddEntry("summary_1030_99eff","10-30%","p");
summary_99eff_leg->AddEntry("summary_3050_99eff","30-50%","p");
summary_99eff_leg->AddEntry("summary_5080_99eff","30-80%","p");

summary_99eff_leg->Draw("same");


//------Print Out & Paste To Overleaf To Create table
cout << endl << endl;
cout << "******Table Print out for OverLeaf******" << endl << endl;
cout << "This is for 15 GeV Threshold " << endl;

cout << "\\begin{center}" << endl;
cout << " \\begin{tabular}{||c c c c||}" << endl;
cout << " \\hline" << endl;
cout << "  Centrality w/ Threshold At 15 GeV & 50\\% Efficiency & Width & 99\\% Efficiency \\\\ [0.5ex]" << endl;
cout << " \\hline\\hline " << endl;
cout << " 0-10\\% & " << fit_010_15GeV->GetParameter(0) << " $\\pm$ " << fit_010_15GeV->GetParError(0) << " & " << fit_010_15GeV->GetParameter(1) << " $\\pm$ " << fit_010_15GeV->GetParError(1) << " & " << fit_010_15GeV->GetParameter(0) + (fit_010_15GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << " \\\\" << endl;
cout << " \\hline" << endl;
cout << " 10-30\\% & " << fit_1030_15GeV->GetParameter(0) << " $\\pm$ " << fit_1030_15GeV->GetParError(0) << " & " << fit_1030_15GeV->GetParameter(1) << " $\\pm$ " << fit_1030_15GeV->GetParError(1)  << " & " << fit_1030_15GeV->GetParameter(0) + (fit_1030_15GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << " \\\\" << endl;
cout << " \\hline" << endl;
cout << " 30-50\\% & " << fit_3050_15GeV->GetParameter(0) << " $\\pm$ " << fit_3050_15GeV->GetParError(0) << " & " << fit_3050_15GeV->GetParameter(1) << " $\\pm$ " << fit_3050_15GeV->GetParError(1) <<  " & " << fit_3050_15GeV->GetParameter(0) + (fit_3050_15GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << " \\\\" << endl; 
cout <<  "\\hline" << endl;
cout << " 50-80\\% & " << fit_5080_15GeV->GetParameter(0) << " $\\pm$ " << fit_5080_15GeV->GetParError(0) << " & " << fit_5080_15GeV->GetParameter(1) << " $\\pm$ " << fit_5080_15GeV->GetParError(1) << " & " << fit_5080_15GeV->GetParameter(0) + (fit_5080_15GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << " \\\\ [1ex]" << endl;
cout << " \\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << "\\end{center}" << endl;
cout << endl << endl;


cout << "This is for the 25 GeV Threshold " << endl;
cout << "\\begin{center}" << endl;
cout << " \\begin{tabular}{||c c c c||}" << endl;
cout << " \\hline" << endl;
cout << "  Centrality w/ Threshold At 25 GeV & 50\\% Efficiency & Width & 99\\% Efficiency \\\\ [0.5ex]" << endl;
cout << " \\hline\\hline " << endl;
cout << " 0-10\\% & " << fit_010_25GeV->GetParameter(0) << " $\\pm$ " << fit_010_25GeV->GetParError(0) << " & " << fit_010_25GeV->GetParameter(1) << " $\\pm$ " << fit_010_25GeV->GetParError(1) << " & " << fit_010_25GeV->GetParameter(0) + (fit_010_25GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << " \\\\" << endl;
cout << " \\hline" << endl;
cout << " 10-30\\% & " << fit_1030_25GeV->GetParameter(0) << " $\\pm$ " << fit_1030_25GeV->GetParError(0) <<" & " << fit_1030_25GeV->GetParameter(1) << " $\\pm$ " << fit_1030_25GeV->GetParError(1) << " & " << fit_1030_25GeV->GetParameter(0) + (fit_1030_25GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << " \\\\" << endl;
cout << " \\hline" << endl;
cout << " 30-50\\% & " << fit_3050_25GeV->GetParameter(0) << " $\\pm$ " << fit_3050_25GeV->GetParError(0) << " & " << fit_3050_25GeV->GetParameter(1) << " $\\pm$ " << fit_3050_25GeV->GetParError(1) << " & " << fit_3050_25GeV->GetParameter(0) + (fit_3050_25GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << " \\\\" << endl; 
cout <<  "\\hline" << endl;
cout << " 50-80\\% & " << fit_5080_25GeV->GetParameter(0) <<  " $\\pm$ " << fit_5080_25GeV->GetParError(0) << " & " << fit_5080_25GeV->GetParameter(1) << " $\\pm$ " << fit_5080_25GeV->GetParError(1) <<" & " << fit_5080_25GeV->GetParameter(0) + (fit_5080_25GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << " \\\\ [1ex]" << endl;
cout << " \\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << "\\end{center}" << endl;
cout << endl << endl;


cout << "This is for the 35 GeV Threshold " << endl;
cout << "\\begin{center}" << endl;
cout << " \\begin{tabular}{||c c c c||}" << endl;
cout << " \\hline" << endl;
cout << "  Centrality w/ Threshold At 35 GeV & 50\\% Efficiency & Width & 99\\% Efficiency \\\\ [0.5ex]" << endl;
cout << " \\hline\\hline " << endl;
cout << " 0-10\\% & " << fit_010_35GeV->GetParameter(0) << " $\\pm$ " << fit_010_35GeV->GetParError(0) << " & " << fit_010_35GeV->GetParameter(1) << " $\\pm$ " << fit_010_35GeV->GetParError(1) << " & " << fit_010_35GeV->GetParameter(0) + (fit_010_35GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << " \\\\" << endl;
cout << " \\hline" << endl;
cout << " 10-30\\% & " << fit_1030_35GeV->GetParameter(0) << " $\\pm$ " << fit_1030_35GeV->GetParError(0)  << " & " << fit_1030_35GeV->GetParameter(1) << " $\\pm$ " << fit_1030_35GeV->GetParError(1) << " & " << fit_1030_35GeV->GetParameter(0) + (fit_1030_35GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << " \\\\" << endl;
cout << " \\hline" << endl;
cout << " 30-50\\% & " << fit_3050_35GeV->GetParameter(0) <<  " $\\pm$ " << fit_3050_35GeV->GetParError(0) << " & " << fit_3050_35GeV->GetParameter(1) << " $\\pm$ " << fit_3050_35GeV->GetParError(1) << " & " << fit_3050_35GeV->GetParameter(0) + (fit_3050_35GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << " \\\\" << endl; 
cout <<  "\\hline" << endl;
cout << " 50-80\\% & " << fit_5080_35GeV->GetParameter(0) << " $\\pm$ " << fit_5080_35GeV->GetParError(0) << " & " << fit_5080_35GeV->GetParameter(1) << " $\\pm$ " << fit_5080_35GeV->GetParError(1) <<  " & " << fit_5080_35GeV->GetParameter(0) + (fit_5080_35GeV->GetParameter(1)*TMath::ErfInverse(0.99)) << " \\\\ [1ex]" << endl;
cout << " \\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << "\\end{center}" << endl;
cout << endl << endl;




//------L1 Efficiencies
TCanvas *new_l1em_eff_barrel = new TCanvas("new_l1em_eff_barrel","new_l1em_eff_barrel",600,500);
TLegend *new_l1em_eff_barrel_lg = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_15GeV_b = new TGraphAsymmErrors();
L1EM_eff_15GeV_b->SetName("L1EM_eff_15GeV_b");                                                                                                                           
L1EM_eff_15GeV_b->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_15GeV_b->SetMarkerStyle(20);
L1EM_eff_15GeV_b->SetTitle("L1EM Efficiency");
L1EM_eff_15GeV_b->SetLineColor(kRed);
L1EM_eff_15GeV_b->SetMarkerColor(kRed);
L1EM_eff_15GeV_b->BayesDivide(photons_wL1EM_15GeV_barrel,all_photons_barrel);
L1EM_eff_15GeV_b->SetMaximum(1.1);
L1EM_eff_15GeV_b->SetMinimum(0);


TGraphAsymmErrors *L1EM_eff_10GeV_b = new TGraphAsymmErrors();
L1EM_eff_10GeV_b->SetName("L1EM_eff_10GeV_b");
L1EM_eff_10GeV_b->SetMarkerStyle(20);
L1EM_eff_10GeV_b->SetMarkerColor(kGreen+2);
L1EM_eff_10GeV_b->SetLineColor(kGreen+2);
L1EM_eff_10GeV_b->BayesDivide(photons_wL1EM_10GeV_barrel,all_photons_barrel);

TGraphAsymmErrors *L1EM_eff_25GeV_b = new TGraphAsymmErrors();
L1EM_eff_25GeV_b->SetName("L1EM_eff_25GeV_b");
L1EM_eff_25GeV_b->SetMarkerStyle(20);
L1EM_eff_25GeV_b->SetMarkerColor(kBlue);
L1EM_eff_25GeV_b->SetLineColor(kBlue);
L1EM_eff_25GeV_b->BayesDivide(photons_wL1EM_25GeV_barrel,all_photons_barrel);


TGraphAsymmErrors *L1EM_eff_35GeV_b = new TGraphAsymmErrors();
L1EM_eff_35GeV_b->SetName("L1EM_eff_35GeV_b");
L1EM_eff_35GeV_b->SetMarkerStyle(20);
L1EM_eff_35GeV_b->SetMarkerColor(kMagenta);
L1EM_eff_35GeV_b->SetLineColor(kMagenta);
L1EM_eff_35GeV_b->BayesDivide(photons_wL1EM_35GeV_barrel,all_photons_barrel);

TPaveText *l1em_photon_0_barrel = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_0_barrel->SetTextSize(0.03);
l1em_photon_0_barrel->SetFillColor(0);
l1em_photon_0_barrel->SetBorderSize(0);
l1em_photon_0_barrel->SetShadowColor(0);
l1em_photon_0_barrel->SetTextAlign(21);
l1em_photon_0_barrel->AddText("");
l1em_photon_0_barrel->AddText("Pb+Pb 2018");
l1em_photon_0_barrel->AddText("Tight Photons");
l1em_photon_0_barrel->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_0_barrel->AddText("DeltaR < 0.15");
l1em_photon_0_barrel->AddText("|#eta| < 1.37");
l1em_photon_0_barrel->SetTextSize(0.036);



//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_10GeV_KH_b = new TF1("fit_10GeV_KH_b","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_10GeV_b = new TF1("fit_10GeV_b", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_15GeV_b = new TF1("fit_15GeV_b", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_25GeV_b = new TF1("fit_25GeV_b","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_35GeV_b = new TF1("fit_35GeV_b","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);

//----Guessing the fit parameters implemeted into our fit function
fit_10GeV_KH_b->SetParameter(0,10.0);
fit_10GeV_KH_b->SetParameter(1,2.0);


fit_10GeV_b->SetParameter(0,10.0);		
fit_10GeV_b->SetParameter(1,2.0);

fit_15GeV_b->SetParameter(0,15.0);		
fit_15GeV_b->SetParameter(1,2.0);

fit_25GeV_b->SetParameter(0,25.0);
fit_25GeV_b->SetParameter(1,2.0);

fit_35GeV_b->SetParameter(0,35.0);
fit_35GeV_b->SetParameter(1,2.0);


L1EM_eff_10GeV_b->Fit("fit_10GeV_b","M Q N","",0,40);
L1EM_eff_15GeV_b->Fit("fit_15GeV_b","M Q N","",0,40);
L1EM_eff_25GeV_b->Fit("fit_25GeV_b", "M Q N","",0,40);
L1EM_eff_35GeV_b->Fit("fit_35GeV_b","M Q N","",0,40);

fit_10GeV_KH_b->SetLineStyle(2);
fit_10GeV_KH_b->SetLineWidth(5);
fit_10GeV_KH_b->SetLineColor(kBlack);
fit_10GeV_KH_b->SetMarkerColor(kBlack);

fit_10GeV_b->SetLineStyle(2);
fit_10GeV_b->SetLineWidth(5);
fit_10GeV_b->SetLineColor(kGreen+2);
fit_10GeV_b->SetMarkerColor(kGreen+2);

fit_15GeV_b->SetLineStyle(2);
fit_15GeV_b->SetLineWidth(5);
fit_15GeV_b->SetLineColor(kRed);
fit_15GeV_b->SetMarkerColor(kRed);

fit_25GeV_b->SetLineStyle(2);
fit_25GeV_b->SetLineWidth(5);
fit_25GeV_b->SetLineColor(kBlue);
fit_25GeV_b->SetMarkerColor(kBlue);

fit_35GeV_b->SetLineStyle(2);
fit_35GeV_b->SetLineWidth(5);
fit_35GeV_b->SetLineColor(kMagenta);
fit_35GeV_b->SetMarkerColor(kMagenta);

TGraph *kh_datapoints_eff = new TGraph (35, x_cor, y_cor); //----Khurt's Data Points For Pb+Pb 2018 Data
kh_datapoints_eff->SetMarkerStyle(20);
kh_datapoints_eff->SetMarkerColor(kBlack);
kh_datapoints_eff->SetLineColor(kBlack);

kh_datapoints_eff->Fit(fit_10GeV_KH_b);

//-----y=1 line
TLine *line_100_barrel = new TLine(0.4,1,43.5,1);
line_100_barrel->SetLineColor(kBlack);
line_100_barrel->SetLineWidth(1);




L1EM_eff_15GeV_b->Draw("A P");
line_100_barrel->Draw("same");

L1EM_eff_10GeV_b->Draw("P same");
L1EM_eff_25GeV_b->Draw("P same");
L1EM_eff_35GeV_b->Draw("P same");

fit_10GeV_b->Draw("same");
fit_15GeV_b->Draw("same");
fit_25GeV_b->Draw("same");
fit_35GeV_b->Draw("same");
kh_datapoints_eff->Draw("P same");

l1em_photon_0_barrel->Draw("same");

new_l1em_eff_barrel_lg->AddEntry("L1EM_eff_10GeV_b","L1 RoI p_{T} > 10 GeV","pe");
new_l1em_eff_barrel_lg->AddEntry("L1EM_eff_15GeV_b","L1 RoI p_{T} > 15 GeV","pe");
new_l1em_eff_barrel_lg->AddEntry("L1EM_eff_25GeV_b","L1 RoI p_{T} > 25 GeV","pe");
new_l1em_eff_barrel_lg->AddEntry("L1EM_eff_35GeV_b","L1 RoI p_{T} > 35 GeV","pe");

new_l1em_eff_barrel_lg->AddEntry("kh_datapoints_eff","L1 RoI > 10 GeV (KH)", "pe");



new_l1em_eff_barrel_lg->Draw("same");


//---------------------------//
//----SumET Distribution-----//
//---------------------------//

//---Barrel (This Is Using The MB and L1EM10 Triggers)
TCanvas *new_SumET_dis_b = new TCanvas("new_SumET_dis_50_b","new_SumET_dis_50_b",600,500);
TLegend *legend_b  = new TLegend(0.6,0.5,0.85,0.6);

barrel_wL1EM_10GeV_SumET_dis_10->SetMarkerStyle(20);
barrel_wL1EM_10GeV_SumET_dis_10->SetMarkerColor(kBlack);
barrel_wL1EM_10GeV_SumET_dis_10->SetLineColor(kBlack);
barrel_wL1EM_10GeV_SumET_dis_10->GetXaxis()->SetTitle("SumET [GeV]");
barrel_wL1EM_10GeV_SumET_dis_10->GetYaxis()->SetTitle("Counts");

barrel_wL1EM_10GeV_SumET_dis_50->SetMarkerStyle(20);
barrel_wL1EM_10GeV_SumET_dis_50->SetMarkerColor(kRed);
barrel_wL1EM_10GeV_SumET_dis_50->SetLineColor(kRed);

barrel_wL1EM_10GeV_SumET_dis_10->Draw();
barrel_wL1EM_10GeV_SumET_dis_50->Draw("same");

legend_b->AddEntry("barrel_wL1EM_10GeV_SumET_dis_10","HLT_noalg_L1EM10","l");
legend_b->AddEntry("barrel_wL1EM_10GeV_SumET_dis_50","HLT_noalg_eb_L1TE50","l");

legend_b->Draw("same");

TPaveText *l1em_trig_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_trig_b->SetTextSize(0.03);
l1em_trig_b->SetFillColor(0);
l1em_trig_b->SetBorderSize(0);
l1em_trig_b->SetShadowColor(0);
l1em_trig_b->SetTextAlign(21);
l1em_trig_b->AddText("");
l1em_trig_b->AddText("Pb+Pb 2018");
l1em_trig_b->AddText("Tight Photons");
l1em_trig_b->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_trig_b->AddText("DeltaR < 0.15");
l1em_trig_b->AddText("|#eta| < 1.37");
l1em_trig_b->SetTextSize(0.036);

l1em_trig_b->Draw("same");


//---Endcaps
TCanvas *new_SumET_dis_10 = new TCanvas("new_SumET_dis_10","new_SumET_dis_10",600,500);
TLegend *legend_e = new TLegend(0.6,0.5,0.85,0.6);

endcaps_wL1EM_10GeV_SumET_dis_10->SetMarkerStyle(20);
endcaps_wL1EM_10GeV_SumET_dis_10->SetMarkerColor(kBlack);
endcaps_wL1EM_10GeV_SumET_dis_10->SetLineColor(kBlack);
endcaps_wL1EM_10GeV_SumET_dis_10->GetXaxis()->SetTitle("SumET [GeV]");
endcaps_wL1EM_10GeV_SumET_dis_10->GetYaxis()->SetTitle("Counts");

endcaps_wL1EM_10GeV_SumET_dis_50->SetMarkerStyle(20);
endcaps_wL1EM_10GeV_SumET_dis_50->SetMarkerColor(kRed);
endcaps_wL1EM_10GeV_SumET_dis_50->SetLineColor(kRed);

endcaps_wL1EM_10GeV_SumET_dis_10->Draw();
endcaps_wL1EM_10GeV_SumET_dis_50->Draw("same");

legend_e->AddEntry("endcaps_wL1EM_10GeV_SumET_dis_10","HLT_noalg_L1EM10","l");
legend_e->AddEntry("endcaps_wL1EM_10GeV_SumET_dis_50","HLT_noalg_eb_L1TE50","l");

legend_e->Draw("same");

TPaveText *l1em_trig_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_trig_e->SetTextSize(0.03);
l1em_trig_e->SetFillColor(0);
l1em_trig_e->SetBorderSize(0);
l1em_trig_e->SetShadowColor(0);
l1em_trig_e->SetTextAlign(21);
l1em_trig_e->AddText("");
l1em_trig_e->AddText("Pb+Pb 2018");
l1em_trig_e->AddText("Tight Photons");
l1em_trig_e->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_trig_e->AddText("DeltaR < 0.15");
l1em_trig_e->AddText("1.56 < |#eta| < 2.37");
l1em_trig_e->SetTextSize(0.036);

l1em_trig_e->Draw("same");


//-----------------------------//
//------L1 Efficiencies--------//
TCanvas *new_l1em_eff_endcaps = new TCanvas("new_l1em_eff_endcaps","new_l1em_eff_endcaps",600,500);
TLegend *new_l1em_eff_endcaps_lg = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_15GeV_endcaps = new TGraphAsymmErrors();
L1EM_eff_15GeV_endcaps->SetName("L1EM_eff_15GeV_endcaps");                                                                                                                           
L1EM_eff_15GeV_endcaps->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_15GeV_endcaps->SetMarkerStyle(20);
L1EM_eff_15GeV_endcaps->SetTitle("L1EM Efficiency");
L1EM_eff_15GeV_endcaps->SetLineColor(kRed);
L1EM_eff_15GeV_endcaps->SetMarkerColor(kRed);
L1EM_eff_15GeV_endcaps->BayesDivide(photons_wL1EM_15GeV_endcaps,all_photons_endcaps);
L1EM_eff_15GeV_endcaps->SetMaximum(1.1);
L1EM_eff_15GeV_endcaps->SetMinimum(0);


TGraphAsymmErrors *L1EM_eff_10GeV_endcaps = new TGraphAsymmErrors();
L1EM_eff_10GeV_endcaps->SetName("L1EM_eff_10GeV_endcaps");
L1EM_eff_10GeV_endcaps->SetMarkerStyle(20);
L1EM_eff_10GeV_endcaps->SetMarkerColor(kGreen+2);
L1EM_eff_10GeV_endcaps->SetLineColor(kGreen+2);
L1EM_eff_10GeV_endcaps->BayesDivide(photons_wL1EM_10GeV_endcaps,all_photons_endcaps);

TGraphAsymmErrors *L1EM_eff_25GeV_endcaps = new TGraphAsymmErrors();
L1EM_eff_25GeV_endcaps->SetName("L1EM_eff_25GeV_endcaps");
L1EM_eff_25GeV_endcaps->SetMarkerStyle(20);
L1EM_eff_25GeV_endcaps->SetMarkerColor(kBlue);
L1EM_eff_25GeV_endcaps->SetLineColor(kBlue);
L1EM_eff_25GeV_endcaps->BayesDivide(photons_wL1EM_25GeV_endcaps,all_photons_endcaps);


TGraphAsymmErrors *L1EM_eff_35GeV_endcaps = new TGraphAsymmErrors();
L1EM_eff_35GeV_endcaps->SetName("L1EM_eff_35GeV_endcaps");
L1EM_eff_35GeV_endcaps->SetMarkerStyle(20);
L1EM_eff_35GeV_endcaps->SetMarkerColor(kMagenta);
L1EM_eff_35GeV_endcaps->SetLineColor(kMagenta);
L1EM_eff_35GeV_endcaps->BayesDivide(photons_wL1EM_35GeV_endcaps,all_photons_endcaps);

TPaveText *l1em_photon_0_endcaps = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_0_endcaps->SetTextSize(0.03);
l1em_photon_0_endcaps->SetFillColor(0);
l1em_photon_0_endcaps->SetBorderSize(0);
l1em_photon_0_endcaps->SetShadowColor(0);
l1em_photon_0_endcaps->SetTextAlign(21);
l1em_photon_0_endcaps->AddText("");
l1em_photon_0_endcaps->AddText("Pb+Pb 2018");
l1em_photon_0_endcaps->AddText("Tight Photons");
l1em_photon_0_endcaps->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_0_endcaps->AddText("DeltaR < 0.15");
l1em_photon_0_endcaps->AddText("1.56 < |#eta| < 2.37");
l1em_photon_0_endcaps->SetTextSize(0.036);



//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_10GeV_KH_e = new  TF1("fit_10GeV_KH_e","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_10GeV_e = new TF1("fit_10GeV_e", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_15GeV_e = new TF1("fit_15GeV_e", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_25GeV_e = new TF1("fit_25GeV_e","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_35GeV_e = new TF1("fit_35GeV_e","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);

//----Guessing the fit parameters implemeted into our fit function
fit_10GeV_KH_e->SetParameter(0,10.0);
fit_10GeV_KH_e->SetParameter(1,2.0);

fit_10GeV_e->SetParameter(0,10.0);		
fit_10GeV_e->SetParameter(1,2.0);

fit_15GeV_e->SetParameter(0,15.0);		
fit_15GeV_e->SetParameter(1,2.0);

fit_25GeV_e->SetParameter(0,25.0);
fit_25GeV_e->SetParameter(1,2.0);

fit_35GeV_e->SetParameter(0,35.0);
fit_35GeV_e->SetParameter(1,2.0);


L1EM_eff_10GeV_endcaps->Fit("fit_10GeV_e","M Q N","",0,40);
L1EM_eff_15GeV_endcaps->Fit("fit_15GeV_e","M Q N","",0,40);
L1EM_eff_25GeV_endcaps->Fit("fit_25GeV_e", "M Q N","",0,40);
L1EM_eff_35GeV_endcaps->Fit("fit_35GeV_e","M Q N","",0,40);

fit_10GeV_KH_e->SetLineStyle(2);
fit_10GeV_KH_e->SetLineWidth(5);
fit_10GeV_KH_e->SetLineColor(kBlack);
fit_10GeV_KH_e->SetMarkerColor(kBlack);

fit_10GeV_e->SetLineStyle(2);
fit_10GeV_e->SetLineWidth(5);
fit_10GeV_e->SetLineColor(kGreen+2);
fit_10GeV_e->SetMarkerColor(kGreen+2);

fit_15GeV_e->SetLineStyle(2);
fit_15GeV_e->SetLineWidth(5);
fit_15GeV_e->SetLineColor(kRed);
fit_15GeV_e->SetMarkerColor(kRed);

fit_25GeV_e->SetLineStyle(2);
fit_25GeV_e->SetLineWidth(5);
fit_25GeV_e->SetLineColor(kBlue);
fit_25GeV_e->SetMarkerColor(kBlue);

fit_35GeV_e->SetLineStyle(2);
fit_35GeV_e->SetLineWidth(5);
fit_35GeV_e->SetLineColor(kMagenta);
fit_35GeV_e->SetMarkerColor(kMagenta);

TGraph *kh_datapoints_endcaps = new TGraph (35, x_cor_endcaps, y_cor_endcaps); //----Khurt's Data Points For Pb+Pb 2018 Data
kh_datapoints_endcaps->SetMarkerStyle(20);
kh_datapoints_endcaps->SetMarkerColor(kBlack);
kh_datapoints_endcaps->SetLineColor(kBlack);

kh_datapoints_endcaps->Fit(fit_10GeV_KH_e);

//-----y=1 line
TLine *line_100_endcaps = new TLine(0.4,1,43.5,1);
line_100_endcaps->SetLineColor(kBlack);
line_100_endcaps->SetLineWidth(1);




L1EM_eff_15GeV_endcaps->Draw("A P");
line_100_endcaps->Draw("same");

L1EM_eff_10GeV_endcaps->Draw("P same");
L1EM_eff_25GeV_endcaps->Draw("P same");
L1EM_eff_35GeV_endcaps->Draw("P same");

fit_10GeV_e->Draw("same");
fit_15GeV_e->Draw("same");
fit_25GeV_e->Draw("same");
fit_35GeV_e->Draw("same");
kh_datapoints_endcaps->Draw("P same");

l1em_photon_0_endcaps->Draw("same");

new_l1em_eff_endcaps_lg->AddEntry("L1EM_eff_10GeV_endcaps","L1 RoI p_{T} > 10 GeV","pe");
new_l1em_eff_endcaps_lg->AddEntry("L1EM_eff_15GeV_endcaps","L1 RoI p_{T} > 15 GeV","pe");
new_l1em_eff_endcaps_lg->AddEntry("L1EM_eff_25GeV_endcaps","L1 RoI p_{T} > 25 GeV","pe");
new_l1em_eff_endcaps_lg->AddEntry("L1EM_eff_35GeV_endcaps","L1 RoI p_{T} > 35 GeV","pe");

new_l1em_eff_endcaps_lg->AddEntry("kh_datapoints_endcaps","L1 RoI > 10 GeV (KH)", "pe");



new_l1em_eff_endcaps_lg->Draw("same");



//------0-10% Centrality
//------Barrel

//------Creating Canvas & Legend
TCanvas *new_l1em_eff_centrality_b = new TCanvas("new_l1em_eff_centrality_b","new_l1em_eff_centrality_b",600,500);
TLegend *new_lg_l1em_cent_b = new TLegend(0.6,0.5,0.85,0.6);



//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_010_15GeV_b = new TGraphAsymmErrors();
L1EM_eff_010_15GeV_b->SetName("L1EM_eff_010_15GeV_b");                                                                                                                           
L1EM_eff_010_15GeV_b->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_010_15GeV_b->SetMarkerStyle(24);
L1EM_eff_010_15GeV_b->SetTitle("L1EM Efficiency");
L1EM_eff_010_15GeV_b->SetLineColor(kRed);
L1EM_eff_010_15GeV_b->SetMarkerColor(kRed);
L1EM_eff_010_15GeV_b->BayesDivide(photon_w_L1EM_15GeV_010_b,all_photons_010_b);
L1EM_eff_010_15GeV_b->SetMaximum(1.1);
L1EM_eff_010_15GeV_b->SetMinimum(0);

TGraphAsymmErrors *L1EM_eff_010_10GeV_b = new TGraphAsymmErrors();
L1EM_eff_010_10GeV_b->SetName("L1EM_eff_010_10GeV_b");
L1EM_eff_010_10GeV_b->SetMarkerStyle(27);
L1EM_eff_010_10GeV_b->SetMarkerColor(kCyan+1);
L1EM_eff_010_10GeV_b->SetLineColor(kCyan+1);
L1EM_eff_010_10GeV_b->BayesDivide(photon_w_L1EM_10GeV_010_b,all_photons_010_b);


TGraphAsymmErrors *L1EM_eff_010_25GeV_b = new TGraphAsymmErrors();
L1EM_eff_010_25GeV_b->SetName("L1EM_eff_010_25GeV_b");
L1EM_eff_010_25GeV_b->SetMarkerStyle(25);
L1EM_eff_010_25GeV_b->SetMarkerColor(kBlue);
L1EM_eff_010_25GeV_b->SetLineColor(kBlue);
L1EM_eff_010_25GeV_b->BayesDivide(photon_w_L1EM_25GeV_010_b,all_photons_010_b);

TGraphAsymmErrors *L1EM_eff_010_35GeV_b = new TGraphAsymmErrors();
L1EM_eff_010_35GeV_b->SetName("L1EM_eff_010_35GeV_b");
L1EM_eff_010_35GeV_b->SetMarkerStyle(26);
L1EM_eff_010_35GeV_b->SetMarkerColor(kMagenta);
L1EM_eff_010_35GeV_b->SetLineColor(kMagenta);
L1EM_eff_010_35GeV_b->BayesDivide(photon_w_L1EM_35GeV_010_b,all_photons_010_b);


TPaveText *l1em_photon_0_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_0_b->SetTextSize(0.03);
l1em_photon_0_b->SetFillColor(0);
l1em_photon_0_b->SetBorderSize(0);
l1em_photon_0_b->SetShadowColor(0);
l1em_photon_0_b->SetTextAlign(21);
l1em_photon_0_b->AddText("");
l1em_photon_0_b->AddText("Pb+Pb");
l1em_photon_0_b->AddText("0-10%");
l1em_photon_0_b->AddText("Tight Photons");
l1em_photon_0_b->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_0_b->AddText("DeltaR < 0.15");
l1em_photon_0_b->AddText("|#eta| < 1.37");
l1em_photon_0_b->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_010_10GeV_b = new TF1("fit_010_10GeV_b", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_010_15GeV_b = new TF1("fit_010_15GeV_b", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_010_25GeV_b = new TF1("fit_010_25GeV_b","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);
TF1* fit_010_35GeV_b = new TF1("fit_010_35GeV_b","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,40);

//----Guessing the fit parameters implemeted into our fit function
fit_010_10GeV_b->SetParameter(0,10.0);		
fit_010_10GeV_b->SetParameter(1,2.0);

fit_010_15GeV_b->SetParameter(0,15.0);		
fit_010_15GeV_b->SetParameter(1,2.0);

fit_010_25GeV_b->SetParameter(0,25.0);
fit_010_25GeV_b->SetParameter(1,2.0);

fit_010_35GeV_b->SetParameter(0,35.0);
fit_010_35GeV_b->SetParameter(1,2.0);


L1EM_eff_010_10GeV_b->Fit("fit_010_10GeV_b","M Q N","",0,photon_w_L1EM_10GeV_010_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_010_b->GetXaxis()->GetNbins()+1));
L1EM_eff_010_15GeV_b->Fit("fit_010_15GeV_b","M Q N","",0,photon_w_L1EM_15GeV_010_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_010_b->GetXaxis()->GetNbins()+1));
L1EM_eff_010_25GeV_b->Fit("fit_010_25GeV_b", "M Q N","",0,photon_w_L1EM_25GeV_010_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_010_b->GetXaxis()->GetNbins()+1));
L1EM_eff_010_35GeV_b->Fit("fit_010_35GeV_b","M Q N","",0,photon_w_L1EM_35GeV_010_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_010_b->GetXaxis()->GetNbins()+1));

fit_010_10GeV_b->SetLineStyle(2);
fit_010_10GeV_b->SetLineColor(kCyan+1);
fit_010_10GeV_b->SetMarkerColor(kCyan+1);

fit_010_15GeV_b->SetLineStyle(2);
fit_010_15GeV_b->SetLineColor(kRed);
fit_010_15GeV_b->SetMarkerColor(kRed);

fit_010_25GeV_b->SetLineStyle(2);
fit_010_25GeV_b->SetLineColor(kBlue);
fit_010_25GeV_b->SetMarkerColor(kBlue);

fit_010_35GeV_b->SetLineStyle(2);
fit_010_35GeV_b->SetLineColor(kMagenta);
fit_010_35GeV_b->SetMarkerColor(kMagenta);

TGraph *kh_datapoints = new TGraph (35, x_cor, y_cor); //----Khurt's Data Points For Pb+Pb 2018 Data
kh_datapoints->SetMarkerStyle(21);
kh_datapoints->SetMarkerColor(kBlack);
kh_datapoints->SetLineColor(kBlack);





//-----y=1 line
TLine *line_100_b = new TLine(0.4,1,169,1);
line_100_b->SetLineColor(kBlack);
line_100_b->SetLineWidth(1);




L1EM_eff_010_15GeV_b->Draw("A P");
line_100_b->Draw("same");

L1EM_eff_010_10GeV_b->Draw("P same");
L1EM_eff_010_25GeV_b->Draw("P same");
L1EM_eff_010_35GeV_b->Draw("P same");

fit_010_10GeV_b->Draw("same");
fit_010_15GeV_b->Draw("same");
fit_010_25GeV_b->Draw("same");
fit_010_35GeV_b->Draw("same");
kh_datapoints->Draw("P same");

l1em_photon_0_b->Draw("same");

new_lg_l1em_cent_b->AddEntry("L1EM_eff_010_10GeV_b","L1 RoI p_{T} > 10 GeV","pe");
new_lg_l1em_cent_b->AddEntry("L1EM_eff_010_15GeV_b","L1 RoI p_{T} > 15 GeV","pe");
new_lg_l1em_cent_b->AddEntry("L1EM_eff_010_25GeV_b","L1 RoI p_{T} > 25 GeV","pe");
new_lg_l1em_cent_b->AddEntry("L1EM_eff_010_35GeV_b","L1 RoI p_{T} > 35 GeV","pe");



new_lg_l1em_cent_b->Draw("same");



//------10-30% Centrality
//------Barrel 

//------Creating Canvas & Legend
TCanvas *new_l1em_eff_centrality_1030_b = new TCanvas("new_l1em_eff_centrality_1030_b","new_l1em_eff_centrality_1030_b",600,500);
TLegend *new_lg_l1em_cent_1030_b = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_1030_15GeV_b = new TGraphAsymmErrors();
L1EM_eff_1030_15GeV_b->SetName("L1EM_eff_1030_15GeV_b");                                                                                                                           
L1EM_eff_1030_15GeV_b->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_1030_15GeV_b->SetMarkerStyle(24);
L1EM_eff_1030_15GeV_b->SetTitle("L1EM Efficiency");
L1EM_eff_1030_15GeV_b->SetLineColor(kRed);
L1EM_eff_1030_15GeV_b->SetMarkerColor(kRed);
L1EM_eff_1030_15GeV_b->BayesDivide(photon_w_L1EM_15GeV_1030_b,all_photons_1030_b);
L1EM_eff_1030_15GeV_b->SetMaximum(1.1);
L1EM_eff_1030_15GeV_b->SetMinimum(0);

TGraphAsymmErrors *L1EM_eff_1030_10GeV_b = new TGraphAsymmErrors();
L1EM_eff_1030_10GeV_b->SetName("L1EM_eff_1030_10GeV_b");
L1EM_eff_1030_10GeV_b->SetMarkerStyle(27);
L1EM_eff_1030_10GeV_b->SetMarkerColor(kCyan+1);
L1EM_eff_1030_10GeV_b->SetLineColor(kCyan+1);
L1EM_eff_1030_10GeV_b->BayesDivide(photon_w_L1EM_10GeV_1030_b,all_photons_1030_b);


TGraphAsymmErrors *L1EM_eff_1030_25GeV_b = new TGraphAsymmErrors();
L1EM_eff_1030_25GeV_b->SetName("L1EM_eff_1030_25GeV_b");
L1EM_eff_1030_25GeV_b->SetMarkerStyle(25);
L1EM_eff_1030_25GeV_b->SetMarkerColor(kBlue);
L1EM_eff_1030_25GeV_b->SetLineColor(kBlue);
L1EM_eff_1030_25GeV_b->BayesDivide(photon_w_L1EM_25GeV_1030_b,all_photons_1030_b);

TGraphAsymmErrors *L1EM_eff_1030_35GeV_b = new TGraphAsymmErrors();
L1EM_eff_1030_35GeV_b->SetName("L1EM_eff_1030_35GeV_b");
L1EM_eff_1030_35GeV_b->SetMarkerStyle(26);
L1EM_eff_1030_35GeV_b->SetMarkerColor(kMagenta);
L1EM_eff_1030_35GeV_b->SetLineColor(kMagenta);
L1EM_eff_1030_35GeV_b->BayesDivide(photon_w_L1EM_35GeV_1030_b,all_photons_1030_b);


TPaveText *l1em_photon_1030_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_1030_b->SetTextSize(0.03);
l1em_photon_1030_b->SetFillColor(0);
l1em_photon_1030_b->SetBorderSize(0);
l1em_photon_1030_b->SetShadowColor(0);
l1em_photon_1030_b->SetTextAlign(21);
l1em_photon_1030_b->AddText("");
l1em_photon_1030_b->AddText("Pb+Pb");
l1em_photon_1030_b->AddText("10-30%");
l1em_photon_1030_b->AddText("Tight Photons");
l1em_photon_1030_b->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_1030_b->AddText("DeltaR < 0.15");
l1em_photon_1030_b->AddText("|#eta| < 1.37");
l1em_photon_1030_b->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_1030_10GeV_b = new TF1("fit_1030_10GeV_b", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_10GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_1030_b->GetXaxis()->GetNbins()+1));
TF1* fit_1030_15GeV_b = new TF1("fit_1030_15GeV_b", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_15GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_1030_b->GetXaxis()->GetNbins()+1));
TF1* fit_1030_25GeV_b = new TF1("fit_1030_25GeV_b","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_25GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_1030_b->GetXaxis()->GetNbins()+1));
TF1* fit_1030_35GeV_b = new TF1("fit_1030_35GeV_b","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_35GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_1030_b->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_1030_10GeV_b->SetParameter(0,10.0);		
fit_1030_10GeV_b->SetParameter(1,2.0);

fit_1030_15GeV_b->SetParameter(0,15.0);		
fit_1030_15GeV_b->SetParameter(1,2.0);

fit_1030_25GeV_b->SetParameter(0,25.0);
fit_1030_25GeV_b->SetParameter(1,2.0);

fit_1030_35GeV_b->SetParameter(0,35.0);
fit_1030_35GeV_b->SetParameter(1,2.0);


L1EM_eff_1030_10GeV_b->Fit("fit_1030_10GeV_b","M Q N","",0,photon_w_L1EM_10GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_1030_b->GetXaxis()->GetNbins()+1));
L1EM_eff_1030_15GeV_b->Fit("fit_1030_15GeV_b","M Q N","",0,photon_w_L1EM_15GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_1030_b->GetXaxis()->GetNbins()+1));
L1EM_eff_1030_25GeV_b->Fit("fit_1030_25GeV_b", "M Q N","",0,photon_w_L1EM_25GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_1030_b->GetXaxis()->GetNbins()+1));
L1EM_eff_1030_35GeV_b->Fit("fit_1030_35GeV_b","M Q N","",0,photon_w_L1EM_35GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_1030_b->GetXaxis()->GetNbins()+1));

fit_1030_10GeV_b->SetLineStyle(2);
fit_1030_10GeV_b->SetLineColor(kCyan+1);
fit_1030_10GeV_b->SetMarkerColor(kCyan+1);

fit_1030_15GeV_b->SetLineStyle(2);
fit_1030_15GeV_b->SetLineColor(kRed);
fit_1030_15GeV_b->SetMarkerColor(kRed);

fit_1030_25GeV_b->SetLineStyle(2);
fit_1030_25GeV_b->SetLineColor(kBlue);
fit_1030_25GeV_b->SetMarkerColor(kBlue);

fit_1030_35GeV_b->SetLineStyle(2);
fit_1030_35GeV_b->SetLineColor(kMagenta);
fit_1030_35GeV_b->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_1_b = new TLine(0.5,1,165,1);
line_100_1_b->SetLineColor(kBlack);
line_100_1_b->SetLineWidth(1);




L1EM_eff_1030_15GeV_b->Draw("A P");
line_100_1_b->Draw("same");

L1EM_eff_1030_10GeV_b->Draw("P same");
L1EM_eff_1030_25GeV_b->Draw("P same");
L1EM_eff_1030_35GeV_b->Draw("P same");

fit_1030_10GeV_b->Draw("same");
fit_1030_15GeV_b->Draw("same");
fit_1030_25GeV_b->Draw("same");
fit_1030_35GeV_b->Draw("same");


l1em_photon_1030_b->Draw("same");

new_lg_l1em_cent_1030_b->AddEntry("L1EM_eff_1030_10GeV_b","L1 RoI p_{T} > 10 GeV","pe");
new_lg_l1em_cent_1030_b->AddEntry("L1EM_eff_1030_15GeV_b","L1 RoI p_{T} > 15 GeV","pe");
new_lg_l1em_cent_1030_b->AddEntry("L1EM_eff_1030_25GeV_b","L1 RoI p_{T} > 25 GeV","pe");
new_lg_l1em_cent_1030_b->AddEntry("L1EM_eff_1030_35GeV_b","L1 RoI p_{T} > 35 GeV","pe");



new_lg_l1em_cent_1030_b->Draw("same");

//------30-50% Centrality
//------Barrel & Endcaps

//------Creating Canvas & Legend
TCanvas *new_l1em_eff_centrality_3050_b = new TCanvas("new_l1em_eff_centrality_3050_b","new_l1em_eff_centrality_3050_b",600,500);
TLegend *new_lg_l1em_cent_3050_b = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_3050_15GeV_b = new TGraphAsymmErrors();
L1EM_eff_3050_15GeV_b->SetName("L1EM_eff_3050_15GeV_b");                                                                                                                           
L1EM_eff_3050_15GeV_b->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_3050_15GeV_b->SetMarkerStyle(24);
L1EM_eff_3050_15GeV_b->SetTitle("L1EM Efficiency");
L1EM_eff_3050_15GeV_b->SetLineColor(kRed);
L1EM_eff_3050_15GeV_b->SetMarkerColor(kRed);
L1EM_eff_3050_15GeV_b->BayesDivide(photon_w_L1EM_15GeV_3050_b,all_photons_3050_b);
L1EM_eff_3050_15GeV_b->SetMaximum(1.1);
L1EM_eff_3050_15GeV_b->SetMinimum(0);

TGraphAsymmErrors *L1EM_eff_3050_10GeV_b = new TGraphAsymmErrors();
L1EM_eff_3050_10GeV_b->SetName("L1EM_eff_3050_10GeV_b");
L1EM_eff_3050_10GeV_b->SetMarkerStyle(27);
L1EM_eff_3050_10GeV_b->SetMarkerColor(kCyan+1);
L1EM_eff_3050_10GeV_b->SetLineColor(kCyan+1);
L1EM_eff_3050_10GeV_b->BayesDivide(photon_w_L1EM_10GeV_3050_b,all_photons_3050_b);

TGraphAsymmErrors *L1EM_eff_3050_25GeV_b = new TGraphAsymmErrors();
L1EM_eff_3050_25GeV_b->SetName("L1EM_eff_3050_25GeV_b");
L1EM_eff_3050_25GeV_b->SetMarkerStyle(25);
L1EM_eff_3050_25GeV_b->SetMarkerColor(kBlue);
L1EM_eff_3050_25GeV_b->SetLineColor(kBlue);
L1EM_eff_3050_25GeV_b->BayesDivide(photon_w_L1EM_25GeV_3050_b,all_photons_3050_b);

TGraphAsymmErrors *L1EM_eff_3050_35GeV_b = new TGraphAsymmErrors();
L1EM_eff_3050_35GeV_b->SetName("L1EM_eff_3050_35GeV_b");
L1EM_eff_3050_35GeV_b->SetMarkerStyle(26);
L1EM_eff_3050_35GeV_b->SetMarkerColor(kMagenta);
L1EM_eff_3050_35GeV_b->SetLineColor(kMagenta);
L1EM_eff_3050_35GeV_b->BayesDivide(photon_w_L1EM_35GeV_3050_b,all_photons_3050_b);


TPaveText *l1em_photon_3050_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_3050_b->SetTextSize(0.03);
l1em_photon_3050_b->SetFillColor(0);
l1em_photon_3050_b->SetBorderSize(0);
l1em_photon_3050_b->SetShadowColor(0);
l1em_photon_3050_b->SetTextAlign(21);
l1em_photon_3050_b->AddText("");
l1em_photon_3050_b->AddText("Pb+Pb");
l1em_photon_3050_b->AddText("30-50%");
l1em_photon_3050_b->AddText("Tight Photons");
l1em_photon_3050_b->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_3050_b->AddText("DeltaR < 0.15");
l1em_photon_3050_b->AddText("|#eta| < 1.37");
l1em_photon_3050_b->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_3050_10GeV_b = new TF1("fit_3050_10GeV_b", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_10GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_10GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_3050_b->GetXaxis()->GetNbins()+1));
TF1* fit_3050_15GeV_b = new TF1("fit_3050_15GeV_b", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_15GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_15GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_3050_b->GetXaxis()->GetNbins()+1));
TF1* fit_3050_25GeV_b = new TF1("fit_3050_25GeV_b","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_25GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_25GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_3050_b->GetXaxis()->GetNbins()+1));
TF1* fit_3050_35GeV_b = new TF1("fit_3050_35GeV_b","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_35GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_35GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_3050_b->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_3050_10GeV_b->SetParameter(0,10.0);		
fit_3050_10GeV_b->SetParameter(1,2.0);

fit_3050_15GeV_b->SetParameter(0,15.0);		
fit_3050_15GeV_b->SetParameter(1,2.0);

fit_3050_25GeV_b->SetParameter(0,25.0);
fit_3050_25GeV_b->SetParameter(1,2.0);

fit_3050_35GeV_b->SetParameter(0,35.0);
fit_3050_35GeV_b->SetParameter(1,2.0);


L1EM_eff_3050_10GeV_b->Fit("fit_3050_10GeV_b","M Q N","",photon_w_L1EM_10GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_10GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_3050_b->GetXaxis()->GetNbins()+1));
L1EM_eff_3050_15GeV_b->Fit("fit_3050_15GeV_b","M Q N","",photon_w_L1EM_15GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_15GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_3050_b->GetXaxis()->GetNbins()+1));
L1EM_eff_3050_25GeV_b->Fit("fit_3050_25GeV_b", "M Q N","",photon_w_L1EM_25GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_25GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_3050_b->GetXaxis()->GetNbins()+1));
L1EM_eff_3050_35GeV_b->Fit("fit_3050_35GeV_b","M Q N","",photon_w_L1EM_35GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_35GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_3050_b->GetXaxis()->GetNbins()+1));

fit_3050_10GeV_b->SetLineStyle(2);
fit_3050_10GeV_b->SetLineColor(kCyan+1);
fit_3050_10GeV_b->SetMarkerColor(kCyan+1);

fit_3050_15GeV_b->SetLineStyle(2);
fit_3050_15GeV_b->SetLineColor(kRed);
fit_3050_15GeV_b->SetMarkerColor(kRed);

fit_3050_25GeV_b->SetLineStyle(2);
fit_3050_25GeV_b->SetLineColor(kBlue);
fit_3050_25GeV_b->SetMarkerColor(kBlue);

fit_3050_35GeV_b->SetLineStyle(2);
fit_3050_35GeV_b->SetLineColor(kMagenta);
fit_3050_35GeV_b->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_2_b = new TLine(2.0,1,153,1);
line_100_2_b->SetLineColor(kBlack);
line_100_2_b->SetLineWidth(1);




L1EM_eff_3050_15GeV_b->Draw("A P");
line_100_2_b->Draw("same");

L1EM_eff_3050_10GeV_b->Draw("P same");
L1EM_eff_3050_25GeV_b->Draw("P same");
L1EM_eff_3050_35GeV_b->Draw("P same");

fit_3050_10GeV_b->Draw("same");
fit_3050_15GeV_b->Draw("same");
fit_3050_25GeV_b->Draw("same");
fit_3050_35GeV_b->Draw("same");


l1em_photon_3050_b->Draw("same");

new_lg_l1em_cent_3050_b->AddEntry("L1EM_eff_3050_10GeV_b","L1 RoI p_{T} > 10 GeV","pe");
new_lg_l1em_cent_3050_b->AddEntry("L1EM_eff_3050_15GeV_b","L1 RoI p_{T} > 15 GeV","pe");
new_lg_l1em_cent_3050_b->AddEntry("L1EM_eff_3050_25GeV_b","L1 RoI p_{T} > 25 GeV","pe");
new_lg_l1em_cent_3050_b->AddEntry("L1EM_eff_3050_35GeV_b","L1 RoI p_{T} > 35 GeV","pe");



new_lg_l1em_cent_3050_b->Draw("same");


//------50-80% Centrality
//------Barrel 

//------Creating Canvas & Legend
TCanvas *new_l1em_eff_centrality_5080_b = new TCanvas("new_l1em_eff_centrality_5080_b","new_l1em_eff_centrality_5080_b",600,500);
TLegend *new_lg_l1em_cent_5080_b = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_5080_15GeV_b = new TGraphAsymmErrors();
L1EM_eff_5080_15GeV_b->SetName("L1EM_eff_5080_15GeV_b");                                                                                                                           
L1EM_eff_5080_15GeV_b->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_5080_15GeV_b->SetMarkerStyle(24);
L1EM_eff_5080_15GeV_b->SetTitle("L1EM Efficiency");
L1EM_eff_5080_15GeV_b->SetLineColor(kRed);
L1EM_eff_5080_15GeV_b->SetMarkerColor(kRed);
L1EM_eff_5080_15GeV_b->BayesDivide(photon_w_L1EM_15GeV_5080_b,all_photons_5080_b);
L1EM_eff_5080_15GeV_b->SetMaximum(1.1);
L1EM_eff_5080_15GeV_b->SetMinimum(0);

TGraphAsymmErrors *L1EM_eff_5080_10GeV_b = new TGraphAsymmErrors();
L1EM_eff_5080_10GeV_b->SetName("L1EM_eff_5080_10GeV_b");
L1EM_eff_5080_10GeV_b->SetMarkerStyle(27);
L1EM_eff_5080_10GeV_b->SetMarkerColor(kCyan+1);
L1EM_eff_5080_10GeV_b->SetLineColor(kCyan+1);
L1EM_eff_5080_10GeV_b->BayesDivide(photon_w_L1EM_10GeV_5080_b,all_photons_5080_b);

TGraphAsymmErrors *L1EM_eff_5080_25GeV_b = new TGraphAsymmErrors();
L1EM_eff_5080_25GeV_b->SetName("L1EM_eff_5080_25GeV_b");
L1EM_eff_5080_25GeV_b->SetMarkerStyle(25);
L1EM_eff_5080_25GeV_b->SetMarkerColor(kBlue);
L1EM_eff_5080_25GeV_b->SetLineColor(kBlue);
L1EM_eff_5080_25GeV_b->BayesDivide(photon_w_L1EM_25GeV_5080_b,all_photons_5080_b);

TGraphAsymmErrors *L1EM_eff_5080_35GeV_b = new TGraphAsymmErrors();
L1EM_eff_5080_35GeV_b->SetName("L1EM_eff_5080_35GeV_b");
L1EM_eff_5080_35GeV_b->SetMarkerStyle(26);
L1EM_eff_5080_35GeV_b->SetMarkerColor(kMagenta);
L1EM_eff_5080_35GeV_b->SetLineColor(kMagenta);
L1EM_eff_5080_35GeV_b->BayesDivide(photon_w_L1EM_35GeV_5080_b,all_photons_5080_b);


TPaveText *l1em_photon_5080_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_5080_b->SetTextSize(0.03);
l1em_photon_5080_b->SetFillColor(0);
l1em_photon_5080_b->SetBorderSize(0);
l1em_photon_5080_b->SetShadowColor(0);
l1em_photon_5080_b->SetTextAlign(21);
l1em_photon_5080_b->AddText("");
l1em_photon_5080_b->AddText("Pb+Pb");
l1em_photon_5080_b->AddText("50-80%");
l1em_photon_5080_b->AddText("Tight Photons");
l1em_photon_5080_b->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_5080_b->AddText("DeltaR < 0.15");
l1em_photon_5080_b->AddText("|#eta| < 1.37");
l1em_photon_5080_b->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_5080_10GeV_b = new TF1("fit_5080_10GeV_b", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_10GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_10GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_5080_b->GetXaxis()->GetNbins()+1));
TF1* fit_5080_15GeV_b = new TF1("fit_5080_15GeV_b", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_15GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_15GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_5080_b->GetXaxis()->GetNbins()+1));
TF1* fit_5080_25GeV_b = new TF1("fit_5080_25GeV_b","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_25GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_25GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_5080_b->GetXaxis()->GetNbins()+1));
TF1* fit_5080_35GeV_b = new TF1("fit_5080_35GeV_b","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_35GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_35GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_5080_b->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_5080_10GeV_b->SetParameter(0,10.0);		
fit_5080_10GeV_b->SetParameter(1,2.0);


fit_5080_15GeV_b->SetParameter(0,15.0);		
fit_5080_15GeV_b->SetParameter(1,2.0);

fit_5080_25GeV_b->SetParameter(0,25.0);
fit_5080_25GeV_b->SetParameter(1,2.0);

fit_5080_35GeV_b->SetParameter(0,35.0);
fit_5080_35GeV_b->SetParameter(1,2.0);


L1EM_eff_5080_10GeV_b->Fit("fit_5080_10GeV_b","M Q N","",photon_w_L1EM_10GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_10GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_5080_b->GetXaxis()->GetNbins()+1));
L1EM_eff_5080_15GeV_b->Fit("fit_5080_15GeV_b","M Q N","",photon_w_L1EM_15GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_15GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_5080_b->GetXaxis()->GetNbins()+1));
L1EM_eff_5080_25GeV_b->Fit("fit_5080_25GeV_b", "M Q N","",photon_w_L1EM_25GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_25GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_5080_b->GetXaxis()->GetNbins()+1));
L1EM_eff_5080_35GeV_b->Fit("fit_5080_35GeV_b","M Q N","",photon_w_L1EM_35GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_35GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_5080_b->GetXaxis()->GetNbins()+1));


fit_5080_10GeV_b->SetLineStyle(2);
fit_5080_10GeV_b->SetLineColor(kCyan+1);
fit_5080_10GeV_b->SetMarkerColor(kCyan+1);

fit_5080_15GeV_b->SetLineStyle(2);
fit_5080_15GeV_b->SetLineColor(kRed);
fit_5080_15GeV_b->SetMarkerColor(kRed);

fit_5080_25GeV_b->SetLineStyle(2);
fit_5080_25GeV_b->SetLineColor(kBlue);
fit_5080_25GeV_b->SetMarkerColor(kBlue);

fit_5080_35GeV_b->SetLineStyle(2);
fit_5080_35GeV_b->SetLineColor(kMagenta);
fit_5080_35GeV_b->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_3_b = new TLine(2.0,1,200,1);
line_100_3_b->SetLineColor(kBlack);
line_100_3_b->SetLineWidth(1);




L1EM_eff_5080_15GeV_b->Draw("A P");
line_100_3_b->Draw("same");

L1EM_eff_5080_10GeV_b->Draw("P same");
L1EM_eff_5080_25GeV_b->Draw("P same");
L1EM_eff_5080_35GeV_b->Draw("P same");

fit_5080_10GeV_b->Draw("same");
fit_5080_15GeV_b->Draw("same");
fit_5080_25GeV_b->Draw("same");
fit_5080_35GeV_b->Draw("same");


l1em_photon_5080_b->Draw("same");

new_lg_l1em_cent_5080_b->AddEntry("L1EM_eff_5080_10GeV_b","L1 RoI p_{T} > 10 GeV","pe");
new_lg_l1em_cent_5080_b->AddEntry("L1EM_eff_5080_15GeV_b","L1 RoI p_{T} > 15 GeV","pe");
new_lg_l1em_cent_5080_b->AddEntry("L1EM_eff_5080_25GeV_b","L1 RoI p_{T} > 25 GeV","pe");
new_lg_l1em_cent_5080_b->AddEntry("L1EM_eff_5080_35GeV_b","L1 RoI p_{T} > 35 GeV","pe");



new_lg_l1em_cent_5080_b->Draw("same");

//-------Summary Plots
//-------Barrel 
//-------This is for the parameter 0
TCanvas *summary_new_b = new TCanvas("summary_new_b","summary_new_b",600,500);
TLegend *new_legend_b = new TLegend(0.6,0.5,0.85,0.6);

//-----0-10%
TH1D *summary_010_0_b = new TH1D("summary_010_0_b","",3,10,40);
summary_010_0_b->SetMarkerStyle(22);
summary_010_0_b->SetMarkerColor(kBlack);
summary_010_0_b->SetLineColor(kBlack);
summary_010_0_b->SetName("summary_010_0_b");
summary_010_0_b->GetXaxis()->SetTitle("L1 p_{T} [GeV]");
summary_010_0_b->GetYaxis()->SetTitle("Value at 50 Percent Efficiency [GeV]");

summary_010_0_b->SetBinContent(1,fit_010_15GeV_b->GetParameter(0));
summary_010_0_b->SetBinError(1,fit_010_15GeV_b->GetParError(0));
summary_010_0_b->SetBinContent(2,fit_010_25GeV_b->GetParameter(0));
summary_010_0_b->SetBinError(2,fit_010_25GeV_b->GetParError(0));
summary_010_0_b->SetBinContent(3,fit_010_35GeV_b->GetParameter(0));
summary_010_0_b->SetBinError(3,fit_010_35GeV_b->GetParError(0));

//------10-30%
TH1D *summary_1030_0_b = new TH1D("summary_1030_0_b","",3,10,40);
summary_1030_0_b->SetMarkerStyle(21);
summary_1030_0_b->SetMarkerColor(kRed);
summary_1030_0_b->SetLineColor(kRed);
summary_1030_0_b->SetName("summary_1030_0_b");

summary_1030_0_b->SetBinContent(1,fit_1030_15GeV_b->GetParameter(0));
summary_1030_0_b->SetBinError(1,fit_1030_15GeV_b->GetParError(0));
summary_1030_0_b->SetBinContent(2,fit_1030_25GeV_b->GetParameter(0));
summary_1030_0_b->SetBinError(2,fit_1030_25GeV_b->GetParError(0));
summary_1030_0_b->SetBinContent(3,fit_1030_35GeV_b->GetParameter(0));
summary_1030_0_b->SetBinError(3,fit_1030_35GeV_b->GetParError(0));

//------30-50%
TH1D *summary_3050_0_b = new TH1D("summary_3050_0_b","",3,10,40);
summary_3050_0_b->SetMarkerStyle(33);
summary_3050_0_b->SetMarkerColor(kMagenta);
summary_3050_0_b->SetLineColor(kMagenta);
summary_3050_0_b->SetName("summary_3050_0_b");

summary_3050_0_b->SetBinContent(1,fit_3050_15GeV_b->GetParameter(0));
summary_3050_0_b->SetBinError(1,fit_3050_15GeV_b->GetParError(0));
summary_3050_0_b->SetBinContent(2,fit_3050_25GeV_b->GetParameter(0));
summary_3050_0_b->SetBinError(2,fit_3050_25GeV_b->GetParError(0));
summary_3050_0_b->SetBinContent(3,fit_3050_35GeV_b->GetParameter(0));
summary_3050_0_b->SetBinError(3,fit_3050_35GeV_b->GetParError(0));

//------50-80%
TH1D *summary_5080_0_b = new TH1D("summary_5080_0_b","",3,10,40);
summary_5080_0_b->SetMarkerStyle(20);
summary_5080_0_b->SetMarkerColor(kBlue);
summary_5080_0_b->SetLineColor(kBlue);
summary_5080_0_b->SetName("summary_5080_0_b");

summary_5080_0_b->SetBinContent(1,fit_5080_15GeV_b->GetParameter(0));
summary_5080_0_b->SetBinError(1,fit_5080_15GeV_b->GetParError(0));
summary_5080_0_b->SetBinContent(2,fit_5080_25GeV_b->GetParameter(0));
summary_5080_0_b->SetBinError(2,fit_5080_25GeV_b->GetParError(0));
summary_5080_0_b->SetBinContent(3,fit_5080_35GeV_b->GetParameter(0));
summary_5080_0_b->SetBinError(3,fit_5080_35GeV_b->GetParError(0));


TPaveText *l1_summary_0_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1_summary_0_b->SetTextSize(0.03);
l1_summary_0_b->SetFillColor(0);
l1_summary_0_b->SetBorderSize(0);
l1_summary_0_b->SetShadowColor(0);
l1_summary_0_b->SetTextAlign(21);
l1_summary_0_b->AddText("");
l1_summary_0_b->AddText("Pb+Pb");
l1_summary_0_b->AddText("Tight Photons");
l1_summary_0_b->AddText("p^{#gamma}_{T} > 15 GeV");
l1_summary_0_b->AddText("HLT Loose Photons");
l1_summary_0_b->AddText("DeltaR < 0.15");
l1_summary_0_b->AddText("|#eta| < 1.37");
l1_summary_0_b->SetTextSize(0.036);




summary_010_0_b->Draw();
summary_1030_0_b->Draw("same");
summary_3050_0_b->Draw("same");
summary_5080_0_b->Draw("same");

l1_summary_0_b->Draw("same");


new_legend_b->AddEntry("summary_010_0_b","0-10%","pe");
new_legend_b->AddEntry("summary_1030_0_b","10-30%","pe");
new_legend_b->AddEntry("summary_3050_0_b","30-50%","pe");
new_legend_b->AddEntry("summary_5080_0_b","50-80%","pe");
new_legend_b->Draw("same");

//-------This is for the parameter 1
TCanvas *summary_new_2_b = new TCanvas("summary_new_2_b","summary_new_2_b",600,500);
TLegend *new_legend_2_b = new TLegend(0.6,0.5,0.85,0.6);

//-----0-10%
TH1D *summary_010_1_b = new TH1D("summary_010_1_b","",3,10,40);
summary_010_1_b->SetMarkerStyle(22);
summary_010_1_b->SetMarkerColor(kBlack);
summary_010_1_b->SetLineColor(kBlack);
summary_010_1_b->SetName("summary_010_1_b");
summary_010_1_b->GetXaxis()->SetTitle("L1 p_{T} [GeV]");
summary_010_1_b->GetYaxis()->SetTitle("Width");

summary_010_1_b->SetBinContent(1,fit_010_15GeV_b->GetParameter(1));
summary_010_1_b->SetBinError(1,fit_010_15GeV_b->GetParError(1));
summary_010_1_b->SetBinContent(2,fit_010_25GeV_b->GetParameter(1));
summary_010_1_b->SetBinError(2,fit_010_25GeV_b->GetParError(1));
summary_010_1_b->SetBinContent(3,fit_010_35GeV_b->GetParameter(1));
summary_010_1_b->SetBinError(3,fit_010_35GeV_b->GetParError(1));

//------10-30%
TH1D *summary_1030_1_b = new TH1D("summary_1030_1_b","",3,10,40);
summary_1030_1_b->SetMarkerStyle(21);
summary_1030_1_b->SetMarkerColor(kRed);
summary_1030_1_b->SetLineColor(kRed);
summary_1030_1_b->SetName("summary_1030_1_b");

summary_1030_1_b->SetBinContent(1,fit_1030_15GeV_b->GetParameter(1));
summary_1030_1_b->SetBinError(1,fit_1030_15GeV_b->GetParError(1));
summary_1030_1_b->SetBinContent(2,fit_1030_25GeV_b->GetParameter(1));
summary_1030_1_b->SetBinError(2,fit_1030_25GeV_b->GetParError(1));
summary_1030_1_b->SetBinContent(3,fit_1030_35GeV_b->GetParameter(1));
summary_1030_1_b->SetBinError(3,fit_1030_35GeV_b->GetParError(1));

//------30-50%
TH1D *summary_3050_1_b = new TH1D("summary_3050_1_b","",3,10,40);
summary_3050_1_b->SetMarkerStyle(33);
summary_3050_1_b->SetMarkerColor(kMagenta);
summary_3050_1_b->SetLineColor(kMagenta);
summary_3050_1_b->SetName("summary_3050_1_b");

summary_3050_1_b->SetBinContent(1,fit_3050_15GeV_b->GetParameter(1));
summary_3050_1_b->SetBinError(1,fit_3050_15GeV_b->GetParError(1));
summary_3050_1_b->SetBinContent(2,fit_3050_25GeV_b->GetParameter(1));
summary_3050_1_b->SetBinError(2,fit_3050_25GeV_b->GetParError(1));
summary_3050_1_b->SetBinContent(3,fit_3050_35GeV_b->GetParameter(1));
summary_3050_1_b->SetBinError(3,fit_3050_35GeV_b->GetParError(1));

//------50-80%
TH1D *summary_5080_1_b = new TH1D("summary_5080_1_b","",3,10,40);
summary_5080_1_b->SetMarkerStyle(20);
summary_5080_1_b->SetMarkerColor(kBlue);
summary_5080_1_b->SetLineColor(kBlue);
summary_5080_1_b->SetName("summary_5080_1_b");

summary_5080_1_b->SetBinContent(1,fit_5080_15GeV_b->GetParameter(1));
summary_5080_1_b->SetBinError(1,fit_5080_15GeV_b->GetParError(1));
summary_5080_1_b->SetBinContent(2,fit_5080_25GeV_b->GetParameter(1));
summary_5080_1_b->SetBinError(2,fit_5080_25GeV_b->GetParError(1));
summary_5080_1_b->SetBinContent(3,fit_5080_35GeV_b->GetParameter(1));
summary_5080_1_b->SetBinError(3,fit_5080_35GeV_b->GetParError(1));


TPaveText *l1_summary_1_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1_summary_1_b->SetTextSize(0.03);
l1_summary_1_b->SetFillColor(0);
l1_summary_1_b->SetBorderSize(0);
l1_summary_1_b->SetShadowColor(0);
l1_summary_1_b->SetTextAlign(21);
l1_summary_1_b->AddText("");
l1_summary_1_b->AddText("Pb+Pb");
l1_summary_1_b->AddText("Tight Photons");
l1_summary_1_b->AddText("p^{#gamma}_{T} > 15 GeV");
l1_summary_1_b->AddText("HLT Loose Photons");
l1_summary_1_b->AddText("DeltaR < 0.15");
l1_summary_1_b->AddText("|#eta| < 1.37");
l1_summary_1_b->SetTextSize(0.036);




summary_010_1_b->Draw();
summary_1030_1_b->Draw("same");
summary_3050_1_b->Draw("same");
summary_5080_1_b->Draw("same");

l1_summary_1_b->Draw("same");


new_legend_2_b->AddEntry("summary_010_1_b","0-10%","pe");
new_legend_2_b->AddEntry("summary_1030_1_b","10-30%","pe");
new_legend_2_b->AddEntry("summary_3050_1_b","30-50%","pe");
new_legend_2_b->AddEntry("summary_5080_1_b","50-80%","pe");
new_legend_2_b->Draw("same");


//------L1 Efficiencies
//------0-10% Centrality
//------Endcaps

//------Creating Canvas & Legend
TCanvas *new_l1em_eff_centrality_e = new TCanvas("new_l1em_eff_centrality_e","new_l1em_eff_centrality_e",600,500);
TLegend *new_lg_l1em_cent_e = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_010_15GeV_e = new TGraphAsymmErrors();
L1EM_eff_010_15GeV_e->SetName("L1EM_eff_010_15GeV_e");                                                                                                                           
L1EM_eff_010_15GeV_e->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_010_15GeV_e->SetMarkerStyle(24);
L1EM_eff_010_15GeV_e->SetTitle("L1EM Efficiency");
L1EM_eff_010_15GeV_e->SetLineColor(kRed);
L1EM_eff_010_15GeV_e->SetMarkerColor(kRed);
L1EM_eff_010_15GeV_e->BayesDivide(photon_w_L1EM_15GeV_010_e,all_photons_010_e);
L1EM_eff_010_15GeV_e->SetMaximum(1.1);
L1EM_eff_010_15GeV_e->SetMinimum(0);


TGraphAsymmErrors *L1EM_eff_010_10GeV_e = new TGraphAsymmErrors();
L1EM_eff_010_10GeV_e->SetName("L1EM_eff_010_10GeV_e");
L1EM_eff_010_10GeV_e->SetMarkerStyle(27);
L1EM_eff_010_10GeV_e->SetMarkerColor(kCyan+1);
L1EM_eff_010_10GeV_e->SetLineColor(kCyan+1);
L1EM_eff_010_10GeV_e->BayesDivide(photon_w_L1EM_10GeV_010_e,all_photons_010_e);

TGraphAsymmErrors *L1EM_eff_010_25GeV_e = new TGraphAsymmErrors();
L1EM_eff_010_25GeV_e->SetName("L1EM_eff_010_25GeV_e");
L1EM_eff_010_25GeV_e->SetMarkerStyle(25);
L1EM_eff_010_25GeV_e->SetMarkerColor(kBlue);
L1EM_eff_010_25GeV_e->SetLineColor(kBlue);
L1EM_eff_010_25GeV_e->BayesDivide(photon_w_L1EM_25GeV_010_e,all_photons_010_e);

TGraphAsymmErrors *L1EM_eff_010_35GeV_e = new TGraphAsymmErrors();
L1EM_eff_010_35GeV_e->SetName("L1EM_eff_010_35GeV_e");
L1EM_eff_010_35GeV_e->SetMarkerStyle(26);
L1EM_eff_010_35GeV_e->SetMarkerColor(kMagenta);
L1EM_eff_010_35GeV_e->SetLineColor(kMagenta);
L1EM_eff_010_35GeV_e->BayesDivide(photon_w_L1EM_35GeV_010_e,all_photons_010_e);


TPaveText *l1em_photon_0_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_0_e->SetTextSize(0.03);
l1em_photon_0_e->SetFillColor(0);
l1em_photon_0_e->SetBorderSize(0);
l1em_photon_0_e->SetShadowColor(0);
l1em_photon_0_e->SetTextAlign(21);
l1em_photon_0_e->AddText("");
l1em_photon_0_e->AddText("Pb+Pb");
l1em_photon_0_e->AddText("0-10%");
l1em_photon_0_e->AddText("Tight Photons");
l1em_photon_0_e->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_0_e->AddText("DeltaR < 0.15");
l1em_photon_0_e->AddText("1.56 < |#eta| < 2.37");
l1em_photon_0_e->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_010_10GeV_e = new TF1("fit_010_10GeV_e", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_10GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_010_e->GetXaxis()->GetNbins()+1));
TF1* fit_010_15GeV_e = new TF1("fit_010_15GeV_e", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_15GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_010_e->GetXaxis()->GetNbins()+1));
TF1* fit_010_25GeV_e = new TF1("fit_010_25GeV_e","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_25GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_010_e->GetXaxis()->GetNbins()+1));
TF1* fit_010_35GeV_e = new TF1("fit_010_35GeV_e","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_35GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_010_e->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_010_10GeV_e->SetParameter(0,10.0);		
fit_010_10GeV_e->SetParameter(1,2.0);

fit_010_15GeV_e->SetParameter(0,15.0);		
fit_010_15GeV_e->SetParameter(1,2.0);

fit_010_25GeV_e->SetParameter(0,25.0);
fit_010_25GeV_e->SetParameter(1,2.0);

fit_010_35GeV_e->SetParameter(0,35.0);
fit_010_35GeV_e->SetParameter(1,2.0);


L1EM_eff_010_10GeV_e->Fit("fit_010_10GeV_e","M Q N","",0,photon_w_L1EM_10GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_010_e->GetXaxis()->GetNbins()+1));
L1EM_eff_010_15GeV_e->Fit("fit_010_15GeV_e","M Q N","",0,photon_w_L1EM_15GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_010_e->GetXaxis()->GetNbins()+1));
L1EM_eff_010_25GeV_e->Fit("fit_010_25GeV_e", "M Q N","",0,photon_w_L1EM_25GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_010_e->GetXaxis()->GetNbins()+1));
L1EM_eff_010_35GeV_e->Fit("fit_010_35GeV_e","M Q N","",0,photon_w_L1EM_35GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_010_e->GetXaxis()->GetNbins()+1));



fit_010_10GeV_e->SetLineStyle(2);
fit_010_10GeV_e->SetLineColor(kCyan+1);
fit_010_10GeV_e->SetMarkerColor(kCyan+1);

fit_010_15GeV_e->SetLineStyle(2);
fit_010_15GeV_e->SetLineColor(kRed);
fit_010_15GeV_e->SetMarkerColor(kRed);

fit_010_25GeV_e->SetLineStyle(2);
fit_010_25GeV_e->SetLineColor(kBlue);
fit_010_25GeV_e->SetMarkerColor(kBlue);

fit_010_35GeV_e->SetLineStyle(2);
fit_010_35GeV_e->SetLineColor(kMagenta);
fit_010_35GeV_e->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_e = new TLine(1,1,146,1);
line_100_e->SetLineColor(kBlack);
line_100_e->SetLineWidth(1);




L1EM_eff_010_15GeV_e->Draw("A P");
line_100_e->Draw("same");

L1EM_eff_010_10GeV_e->Draw("P same");
L1EM_eff_010_25GeV_e->Draw("P same");
L1EM_eff_010_35GeV_e->Draw("P same");


fit_010_10GeV_e->Draw("same");
fit_010_15GeV_e->Draw("same");
fit_010_25GeV_e->Draw("same");
fit_010_35GeV_e->Draw("same");


l1em_photon_0_e->Draw("same");

new_lg_l1em_cent_e->AddEntry("L1EM_eff_010_10GeV_e","L1 RoI p_{T} > 10 GeV","pe");
new_lg_l1em_cent_e->AddEntry("L1EM_eff_010_15GeV_e","L1 RoI p_{T} > 15 GeV","pe");
new_lg_l1em_cent_e->AddEntry("L1EM_eff_010_25GeV_e","L1 RoI p_{T} > 25 GeV","pe");
new_lg_l1em_cent_e->AddEntry("L1EM_eff_010_35GeV_e","L1 RoI p_{T} > 35 GeV","pe");



new_lg_l1em_cent_e->Draw("same");



//------10-30% Centrality
//------Endcaps 

//------Creating Canvas & Legend
TCanvas *new_l1em_eff_centrality_1030_e = new TCanvas("new_l1em_eff_centrality_1030_e","new_l1em_eff_centrality_1030_e",600,500);
TLegend *new_lg_l1em_cent_1030_e = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_1030_15GeV_e = new TGraphAsymmErrors();
L1EM_eff_1030_15GeV_e->SetName("L1EM_eff_1030_15GeV_e");                                                                                                                           
L1EM_eff_1030_15GeV_e->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_1030_15GeV_e->SetMarkerStyle(24);
L1EM_eff_1030_15GeV_e->SetTitle("L1EM Efficiency");
L1EM_eff_1030_15GeV_e->SetLineColor(kRed);
L1EM_eff_1030_15GeV_e->SetMarkerColor(kRed);
L1EM_eff_1030_15GeV_e->BayesDivide(photon_w_L1EM_15GeV_1030_e,all_photons_1030_e);
L1EM_eff_1030_15GeV_e->SetMaximum(1.1);
L1EM_eff_1030_15GeV_e->SetMinimum(0);

TGraphAsymmErrors *L1EM_eff_1030_10GeV_e = new TGraphAsymmErrors();
L1EM_eff_1030_10GeV_e->SetName("L1EM_eff_1030_25GeV_e");
L1EM_eff_1030_10GeV_e->SetMarkerStyle(33);
L1EM_eff_1030_10GeV_e->SetMarkerColor(kCyan+1);
L1EM_eff_1030_10GeV_e->SetLineColor(kCyan+1);
L1EM_eff_1030_10GeV_e->BayesDivide(photon_w_L1EM_10GeV_1030_e,all_photons_1030_e);

TGraphAsymmErrors *L1EM_eff_1030_25GeV_e = new TGraphAsymmErrors();
L1EM_eff_1030_25GeV_e->SetName("L1EM_eff_1030_25GeV_e");
L1EM_eff_1030_25GeV_e->SetMarkerStyle(25);
L1EM_eff_1030_25GeV_e->SetMarkerColor(kBlue);
L1EM_eff_1030_25GeV_e->SetLineColor(kBlue);
L1EM_eff_1030_25GeV_e->BayesDivide(photon_w_L1EM_25GeV_1030_e,all_photons_1030_e);

TGraphAsymmErrors *L1EM_eff_1030_35GeV_e = new TGraphAsymmErrors();
L1EM_eff_1030_35GeV_e->SetName("L1EM_eff_1030_35GeV_e");
L1EM_eff_1030_35GeV_e->SetMarkerStyle(26);
L1EM_eff_1030_35GeV_e->SetMarkerColor(kMagenta);
L1EM_eff_1030_35GeV_e->SetLineColor(kMagenta);
L1EM_eff_1030_35GeV_e->BayesDivide(photon_w_L1EM_35GeV_1030_e,all_photons_1030_e);


TPaveText *l1em_photon_1030_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_1030_e->SetTextSize(0.03);
l1em_photon_1030_e->SetFillColor(0);
l1em_photon_1030_e->SetBorderSize(0);
l1em_photon_1030_e->SetShadowColor(0);
l1em_photon_1030_e->SetTextAlign(21);
l1em_photon_1030_e->AddText("");
l1em_photon_1030_e->AddText("Pb+Pb");
l1em_photon_1030_e->AddText("10-30%");
l1em_photon_1030_e->AddText("Tight Photons");
l1em_photon_1030_e->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_1030_e->AddText("DeltaR < 0.15");
l1em_photon_1030_e->AddText("1.56 < |#eta| < 2.37");
l1em_photon_1030_e->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_1030_10GeV_e = new TF1("fit_1030_10GeV_e", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_10GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_1030_e->GetXaxis()->GetNbins()+1));
TF1* fit_1030_15GeV_e = new TF1("fit_1030_15GeV_e", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_15GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_1030_e->GetXaxis()->GetNbins()+1));
TF1* fit_1030_25GeV_e = new TF1("fit_1030_25GeV_e","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_25GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_1030_e->GetXaxis()->GetNbins()+1));
TF1* fit_1030_35GeV_e = new TF1("fit_1030_35GeV_e","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_L1EM_35GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_1030_e->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_1030_10GeV_e->SetParameter(0,10.0);		
fit_1030_10GeV_e->SetParameter(1,2.0);

fit_1030_15GeV_e->SetParameter(0,15.0);		
fit_1030_15GeV_e->SetParameter(1,2.0);

fit_1030_25GeV_e->SetParameter(0,25.0);
fit_1030_25GeV_e->SetParameter(1,2.0);

fit_1030_35GeV_e->SetParameter(0,35.0);
fit_1030_35GeV_e->SetParameter(1,2.0);


L1EM_eff_1030_10GeV_e->Fit("fit_1030_10GeV_e","M Q N","",0,photon_w_L1EM_10GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_1030_e->GetXaxis()->GetNbins()+1));
L1EM_eff_1030_15GeV_e->Fit("fit_1030_15GeV_e","M Q N","",0,photon_w_L1EM_15GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_1030_e->GetXaxis()->GetNbins()+1));
L1EM_eff_1030_25GeV_e->Fit("fit_1030_25GeV_e", "M Q N","",0,photon_w_L1EM_25GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_1030_e->GetXaxis()->GetNbins()+1));
L1EM_eff_1030_35GeV_e->Fit("fit_1030_35GeV_e","M Q N","",0,photon_w_L1EM_35GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_1030_e->GetXaxis()->GetNbins()+1));


fit_1030_15GeV_e->SetLineStyle(2);
fit_1030_15GeV_e->SetLineColor(kCyan+1);
fit_1030_15GeV_e->SetMarkerColor(kCyan+1);


fit_1030_15GeV_e->SetLineStyle(2);
fit_1030_15GeV_e->SetLineColor(kRed);
fit_1030_15GeV_e->SetMarkerColor(kRed);

fit_1030_25GeV_e->SetLineStyle(2);
fit_1030_25GeV_e->SetLineColor(kBlue);
fit_1030_25GeV_e->SetMarkerColor(kBlue);

fit_1030_35GeV_e->SetLineStyle(2);
fit_1030_35GeV_e->SetLineColor(kMagenta);
fit_1030_35GeV_e->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_1_e = new TLine(1.0,1,147,1);
line_100_1_e->SetLineColor(kBlack);
line_100_1_e->SetLineWidth(1);




L1EM_eff_1030_15GeV_e->Draw("A P");
line_100_1_e->Draw("same");

L1EM_eff_1030_10GeV_e->Draw("P same");
L1EM_eff_1030_25GeV_e->Draw("P same");
L1EM_eff_1030_35GeV_e->Draw("P same");

fit_1030_10GeV_e->Draw("same");
fit_1030_15GeV_e->Draw("same");
fit_1030_25GeV_e->Draw("same");
fit_1030_35GeV_e->Draw("same");


l1em_photon_1030_e->Draw("same");

new_lg_l1em_cent_1030_e->AddEntry("L1EM_eff_1030_10GeV_e","L1 RoI p_{T} > 10 GeV","pe");
new_lg_l1em_cent_1030_e->AddEntry("L1EM_eff_1030_15GeV_e","L1 RoI p_{T} > 15 GeV","pe");
new_lg_l1em_cent_1030_e->AddEntry("L1EM_eff_1030_25GeV_e","L1 RoI p_{T} > 25 GeV","pe");
new_lg_l1em_cent_1030_e->AddEntry("L1EM_eff_1030_35GeV_e","L1 RoI p_{T} > 35 GeV","pe");



new_lg_l1em_cent_1030_e->Draw("same");

//------30-50% Centrality
//------Barrel & Endcaps

//------Creating Canvas & Legend
TCanvas *new_l1em_eff_centrality_3050_e = new TCanvas("new_l1em_eff_centrality_3050_e","new_l1em_eff_centrality_3050_e",600,500);
TLegend *new_lg_l1em_cent_3050_e = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_3050_15GeV_e = new TGraphAsymmErrors();
L1EM_eff_3050_15GeV_e->SetName("L1EM_eff_3050_15GeV_e");                                                                                                                           
L1EM_eff_3050_15GeV_e->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_3050_15GeV_e->SetMarkerStyle(24);
L1EM_eff_3050_15GeV_e->SetTitle("L1EM Efficiency");
L1EM_eff_3050_15GeV_e->SetLineColor(kRed);
L1EM_eff_3050_15GeV_e->SetMarkerColor(kRed);
L1EM_eff_3050_15GeV_e->BayesDivide(photon_w_L1EM_15GeV_3050_e,all_photons_3050_e);
L1EM_eff_3050_15GeV_e->SetMaximum(1.1);
L1EM_eff_3050_15GeV_e->SetMinimum(0);

TGraphAsymmErrors *L1EM_eff_3050_10GeV_e = new TGraphAsymmErrors();
L1EM_eff_3050_10GeV_e->SetName("L1EM_eff_3050_25GeV_e");
L1EM_eff_3050_10GeV_e->SetMarkerStyle(25);
L1EM_eff_3050_10GeV_e->SetMarkerColor(kBlue);
L1EM_eff_3050_10GeV_e->SetLineColor(kBlue);
L1EM_eff_3050_10GeV_e->BayesDivide(photon_w_L1EM_10GeV_3050_e,all_photons_3050_e);

TGraphAsymmErrors *L1EM_eff_3050_25GeV_e = new TGraphAsymmErrors();
L1EM_eff_3050_25GeV_e->SetName("L1EM_eff_3050_25GeV_e");
L1EM_eff_3050_25GeV_e->SetMarkerStyle(25);
L1EM_eff_3050_25GeV_e->SetMarkerColor(kBlue);
L1EM_eff_3050_25GeV_e->SetLineColor(kBlue);
L1EM_eff_3050_25GeV_e->BayesDivide(photon_w_L1EM_25GeV_3050_e,all_photons_3050_e);

TGraphAsymmErrors *L1EM_eff_3050_35GeV_e = new TGraphAsymmErrors();
L1EM_eff_3050_35GeV_e->SetName("L1EM_eff_3050_35GeV_e");
L1EM_eff_3050_35GeV_e->SetMarkerStyle(26);
L1EM_eff_3050_35GeV_e->SetMarkerColor(kMagenta);
L1EM_eff_3050_35GeV_e->SetLineColor(kMagenta);
L1EM_eff_3050_35GeV_e->BayesDivide(photon_w_L1EM_35GeV_3050_e,all_photons_3050_e);


TPaveText *l1em_photon_3050_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_3050_e->SetTextSize(0.03);
l1em_photon_3050_e->SetFillColor(0);
l1em_photon_3050_e->SetBorderSize(0);
l1em_photon_3050_e->SetShadowColor(0);
l1em_photon_3050_e->SetTextAlign(21);
l1em_photon_3050_e->AddText("");
l1em_photon_3050_e->AddText("Pb+Pb");
l1em_photon_3050_e->AddText("30-50%");
l1em_photon_3050_e->AddText("Tight Photons");
l1em_photon_3050_e->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_3050_e->AddText("DeltaR < 0.15");
l1em_photon_3050_e->AddText("1.56 < |#eta| < 2.37");
l1em_photon_3050_e->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_3050_10GeV_e = new TF1("fit_3050_10GeV_e", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_10GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_10GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_3050_e->GetXaxis()->GetNbins()+1));
TF1* fit_3050_15GeV_e = new TF1("fit_3050_15GeV_e", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_15GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_15GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_3050_e->GetXaxis()->GetNbins()+1));
TF1* fit_3050_25GeV_e = new TF1("fit_3050_25GeV_e","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_25GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_25GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_3050_e->GetXaxis()->GetNbins()+1));
TF1* fit_3050_35GeV_e = new TF1("fit_3050_35GeV_e","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_35GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_35GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_3050_e->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_3050_10GeV_e->SetParameter(0,10.0);		
fit_3050_10GeV_e->SetParameter(1,2.0);

fit_3050_15GeV_e->SetParameter(0,15.0);		
fit_3050_15GeV_e->SetParameter(1,2.0);

fit_3050_25GeV_e->SetParameter(0,25.0);
fit_3050_25GeV_e->SetParameter(1,2.0);

fit_3050_35GeV_e->SetParameter(0,35.0);
fit_3050_35GeV_e->SetParameter(1,2.0);


L1EM_eff_3050_10GeV_e->Fit("fit_3050_10GeV_e","M Q N","",photon_w_L1EM_10GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_10GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_3050_e->GetXaxis()->GetNbins()+1));
L1EM_eff_3050_15GeV_e->Fit("fit_3050_15GeV_e","M Q N","",photon_w_L1EM_15GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_15GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_3050_e->GetXaxis()->GetNbins()+1));
L1EM_eff_3050_25GeV_e->Fit("fit_3050_25GeV_e", "M Q N","",photon_w_L1EM_25GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_25GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_3050_e->GetXaxis()->GetNbins()+1));
L1EM_eff_3050_35GeV_e->Fit("fit_3050_35GeV_e","M Q N","",photon_w_L1EM_35GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_35GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_3050_e->GetXaxis()->GetNbins()+1));

fit_3050_10GeV_e->SetLineStyle(2);
fit_3050_10GeV_e->SetLineColor(kCyan+1);
fit_3050_10GeV_e->SetMarkerColor(kCyan+1);

fit_3050_15GeV_e->SetLineStyle(2);
fit_3050_15GeV_e->SetLineColor(kRed);
fit_3050_15GeV_e->SetMarkerColor(kRed);

fit_3050_25GeV_e->SetLineStyle(2);
fit_3050_25GeV_e->SetLineColor(kBlue);
fit_3050_25GeV_e->SetMarkerColor(kBlue);

fit_3050_35GeV_e->SetLineStyle(2);
fit_3050_35GeV_e->SetLineColor(kMagenta);
fit_3050_35GeV_e->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_2_e = new TLine(7.0,1,109,1);
line_100_2_e->SetLineColor(kBlack);
line_100_2_e->SetLineWidth(1);




L1EM_eff_3050_15GeV_e->Draw("A P");
line_100_2_e->Draw("same");

L1EM_eff_3050_10GeV_e->Draw("P same");
L1EM_eff_3050_25GeV_e->Draw("P same");
L1EM_eff_3050_35GeV_e->Draw("P same");

fit_3050_10GeV_e->Draw("same");
fit_3050_15GeV_e->Draw("same");
fit_3050_25GeV_e->Draw("same");
fit_3050_35GeV_e->Draw("same");


l1em_photon_3050_e->Draw("same");

new_lg_l1em_cent_3050_e->AddEntry("L1EM_eff_3050_10GeV_e","L1 RoI p_{T} > 15 GeV","pe");
new_lg_l1em_cent_3050_e->AddEntry("L1EM_eff_3050_15GeV_e","L1 RoI p_{T} > 15 GeV","pe");
new_lg_l1em_cent_3050_e->AddEntry("L1EM_eff_3050_25GeV_e","L1 RoI p_{T} > 25 GeV","pe");
new_lg_l1em_cent_3050_e->AddEntry("L1EM_eff_3050_35GeV_e","L1 RoI p_{T} > 35 GeV","pe");



new_lg_l1em_cent_3050_e->Draw("same");


//------50-80% Centrality
//------Barrel 

//------Creating Canvas & Legend
TCanvas *new_l1em_eff_centrality_5080_e = new TCanvas("new_l1em_eff_centrality_5080_e","new_l1em_eff_centrality_5080_e",600,500);
TLegend *new_lg_l1em_cent_5080_e = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *L1EM_eff_5080_15GeV_e = new TGraphAsymmErrors();
L1EM_eff_5080_15GeV_e->SetName("L1EM_eff_5080_15GeV_e");                                                                                                                           
L1EM_eff_5080_15GeV_e->GetXaxis()->SetTitle("p_{T} [GeV]");
L1EM_eff_5080_15GeV_e->SetMarkerStyle(24);
L1EM_eff_5080_15GeV_e->SetTitle("L1EM Efficiency");
L1EM_eff_5080_15GeV_e->SetLineColor(kRed);
L1EM_eff_5080_15GeV_e->SetMarkerColor(kRed);
L1EM_eff_5080_15GeV_e->BayesDivide(photon_w_L1EM_15GeV_5080_e,all_photons_5080_e);
L1EM_eff_5080_15GeV_e->SetMaximum(1.1);
L1EM_eff_5080_15GeV_e->SetMinimum(0);

TGraphAsymmErrors *L1EM_eff_5080_10GeV_e = new TGraphAsymmErrors();
L1EM_eff_5080_10GeV_e->SetName("L1EM_eff_5080_25GeV_e");
L1EM_eff_5080_10GeV_e->SetMarkerStyle(25);
L1EM_eff_5080_10GeV_e->SetMarkerColor(kBlue);
L1EM_eff_5080_10GeV_e->SetLineColor(kBlue);
L1EM_eff_5080_10GeV_e->BayesDivide(photon_w_L1EM_10GeV_5080_e,all_photons_5080_e);


TGraphAsymmErrors *L1EM_eff_5080_25GeV_e = new TGraphAsymmErrors();
L1EM_eff_5080_25GeV_e->SetName("L1EM_eff_5080_25GeV_e");
L1EM_eff_5080_25GeV_e->SetMarkerStyle(25);
L1EM_eff_5080_25GeV_e->SetMarkerColor(kBlue);
L1EM_eff_5080_25GeV_e->SetLineColor(kBlue);
L1EM_eff_5080_25GeV_e->BayesDivide(photon_w_L1EM_25GeV_5080_e,all_photons_5080_e);

TGraphAsymmErrors *L1EM_eff_5080_35GeV_e = new TGraphAsymmErrors();
L1EM_eff_5080_35GeV_e->SetName("L1EM_eff_5080_35GeV_e");
L1EM_eff_5080_35GeV_e->SetMarkerStyle(26);
L1EM_eff_5080_35GeV_e->SetMarkerColor(kMagenta);
L1EM_eff_5080_35GeV_e->SetLineColor(kMagenta);
L1EM_eff_5080_35GeV_e->BayesDivide(photon_w_L1EM_35GeV_5080_e,all_photons_5080_e);


TPaveText *l1em_photon_5080_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_5080_e->SetTextSize(0.03);
l1em_photon_5080_e->SetFillColor(0);
l1em_photon_5080_e->SetBorderSize(0);
l1em_photon_5080_e->SetShadowColor(0);
l1em_photon_5080_e->SetTextAlign(21);
l1em_photon_5080_e->AddText("");
l1em_photon_5080_e->AddText("Pb+Pb");
l1em_photon_5080_e->AddText("50-80%");
l1em_photon_5080_e->AddText("Tight Photons");
l1em_photon_5080_e->AddText("p^{#gamma}_{T} > 15 GeV");
l1em_photon_5080_e->AddText("DeltaR < 0.15");
l1em_photon_5080_e->AddText("1.56 < |#eta| < 2.37");
l1em_photon_5080_e->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_5080_10GeV_e = new TF1("fit_5080_10GeV_e", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_10GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_10GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_5080_e->GetXaxis()->GetNbins()+1));
TF1* fit_5080_15GeV_e = new TF1("fit_5080_15GeV_e", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_15GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_15GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_5080_e->GetXaxis()->GetNbins()+1));
TF1* fit_5080_25GeV_e = new TF1("fit_5080_25GeV_e","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_25GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_25GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_5080_e->GetXaxis()->GetNbins()+1));
TF1* fit_5080_35GeV_e = new TF1("fit_5080_35GeV_e","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_L1EM_35GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_35GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_5080_e->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_5080_10GeV_e->SetParameter(0,10.0);		
fit_5080_10GeV_e->SetParameter(1,2.0);

fit_5080_15GeV_e->SetParameter(0,15.0);		
fit_5080_15GeV_e->SetParameter(1,2.0);

fit_5080_25GeV_e->SetParameter(0,25.0);
fit_5080_25GeV_e->SetParameter(1,2.0);

fit_5080_35GeV_e->SetParameter(0,35.0);
fit_5080_35GeV_e->SetParameter(1,2.0);


L1EM_eff_5080_10GeV_e->Fit("fit_5080_10GeV_e","M Q N","",photon_w_L1EM_10GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_10GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_10GeV_5080_e->GetXaxis()->GetNbins()+1));
L1EM_eff_5080_15GeV_e->Fit("fit_5080_15GeV_e","M Q N","",photon_w_L1EM_15GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_15GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_15GeV_5080_e->GetXaxis()->GetNbins()+1));
L1EM_eff_5080_25GeV_e->Fit("fit_5080_25GeV_e", "M Q N","",photon_w_L1EM_25GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_25GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_25GeV_5080_e->GetXaxis()->GetNbins()+1));
L1EM_eff_5080_35GeV_e->Fit("fit_5080_35GeV_e","M Q N","",photon_w_L1EM_35GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_L1EM_35GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_L1EM_35GeV_5080_e->GetXaxis()->GetNbins()+1));

fit_5080_10GeV_e->SetLineStyle(2);
fit_5080_10GeV_e->SetLineColor(kRed);
fit_5080_10GeV_e->SetMarkerColor(kRed);

fit_5080_15GeV_e->SetLineStyle(2);
fit_5080_15GeV_e->SetLineColor(kRed);
fit_5080_15GeV_e->SetMarkerColor(kRed);

fit_5080_25GeV_e->SetLineStyle(2);
fit_5080_25GeV_e->SetLineColor(kBlue);
fit_5080_25GeV_e->SetMarkerColor(kBlue);

fit_5080_35GeV_e->SetLineStyle(2);
fit_5080_35GeV_e->SetLineColor(kMagenta);
fit_5080_35GeV_e->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_3_e = new TLine(9.0,1,87,1);
line_100_3_e->SetLineColor(kBlack);
line_100_3_e->SetLineWidth(1);




L1EM_eff_5080_15GeV_e->Draw("A P");
line_100_3_e->Draw("same");

L1EM_eff_5080_10GeV_e->Draw("P same");
L1EM_eff_5080_25GeV_e->Draw("P same");
L1EM_eff_5080_35GeV_e->Draw("P same");

fit_5080_10GeV_e->Draw("same");
fit_5080_15GeV_e->Draw("same");
fit_5080_25GeV_e->Draw("same");
fit_5080_35GeV_e->Draw("same");


l1em_photon_5080_e->Draw("same");

new_lg_l1em_cent_5080_e->AddEntry("L1EM_eff_5080_10GeV_e","HLT p_{T} > 10 GeV","pe");
new_lg_l1em_cent_5080_e->AddEntry("L1EM_eff_5080_15GeV_e","HLT p_{T} > 15 GeV","pe");
new_lg_l1em_cent_5080_e->AddEntry("L1EM_eff_5080_25GeV_e","HLT p_{T} > 25 GeV","pe");
new_lg_l1em_cent_5080_e->AddEntry("L1EM_eff_5080_35GeV_e","HLT p_{T} > 35 GeV","pe");



new_lg_l1em_cent_5080_e->Draw("same");

//-------Summary Plots
//-------Barrel 
//-------This is for the parameter 0
TCanvas *summary_new_e = new TCanvas("summary_new_e","summary_new_e",600,500);
TLegend *new_legend_e = new TLegend(0.6,0.5,0.85,0.6);

//-----0-10%
TH1D *summary_010_0_e = new TH1D("summary_010_0_e","",3,10,40);
summary_010_0_e->SetMarkerStyle(22);
summary_010_0_e->SetMarkerColor(kBlack);
summary_010_0_e->SetLineColor(kBlack);
summary_010_0_e->SetName("summary_010_0_e");
summary_010_0_e->GetXaxis()->SetTitle("L1 p_{T} [GeV]");
summary_010_0_e->GetYaxis()->SetTitle("Value at 50 Percent Efficiency [GeV]");

summary_010_0_e->SetBinContent(1,fit_010_15GeV_e->GetParameter(0));
summary_010_0_e->SetBinError(1,fit_010_15GeV_e->GetParError(0));
summary_010_0_e->SetBinContent(2,fit_010_25GeV_e->GetParameter(0));
summary_010_0_e->SetBinError(2,fit_010_25GeV_e->GetParError(0));
summary_010_0_e->SetBinContent(3,fit_010_35GeV_e->GetParameter(0));
summary_010_0_e->SetBinError(3,fit_010_35GeV_e->GetParError(0));

//------10-30%
TH1D *summary_1030_0_e = new TH1D("summary_1030_0_e","",3,10,40);
summary_1030_0_e->SetMarkerStyle(21);
summary_1030_0_e->SetMarkerColor(kRed);
summary_1030_0_e->SetLineColor(kRed);
summary_1030_0_e->SetName("summary_1030_0_e");

summary_1030_0_e->SetBinContent(1,fit_1030_15GeV_e->GetParameter(0));
summary_1030_0_e->SetBinError(1,fit_1030_15GeV_e->GetParError(0));
summary_1030_0_e->SetBinContent(2,fit_1030_25GeV_e->GetParameter(0));
summary_1030_0_e->SetBinError(2,fit_1030_25GeV_e->GetParError(0));
summary_1030_0_e->SetBinContent(3,fit_1030_35GeV_e->GetParameter(0));
summary_1030_0_e->SetBinError(3,fit_1030_35GeV_e->GetParError(0));

//------30-50%
TH1D *summary_3050_0_e = new TH1D("summary_3050_0_e","",3,10,40);
summary_3050_0_e->SetMarkerStyle(33);
summary_3050_0_e->SetMarkerColor(kMagenta);
summary_3050_0_e->SetLineColor(kMagenta);
summary_3050_0_e->SetName("summary_3050_0_e");

summary_3050_0_e->SetBinContent(1,fit_3050_15GeV_e->GetParameter(0));
summary_3050_0_e->SetBinError(1,fit_3050_15GeV_e->GetParError(0));
summary_3050_0_e->SetBinContent(2,fit_3050_25GeV_e->GetParameter(0));
summary_3050_0_e->SetBinError(2,fit_3050_25GeV_e->GetParError(0));
summary_3050_0_e->SetBinContent(3,fit_3050_35GeV_e->GetParameter(0));
summary_3050_0_e->SetBinError(3,fit_3050_35GeV_e->GetParError(0));

//------50-80%
TH1D *summary_5080_0_e = new TH1D("summary_5080_0_e","",3,10,40);
summary_5080_0_e->SetMarkerStyle(20);
summary_5080_0_e->SetMarkerColor(kBlue);
summary_5080_0_e->SetLineColor(kBlue);
summary_5080_0_e->SetName("summary_5080_0_e");

summary_5080_0_e->SetBinContent(1,fit_5080_15GeV_e->GetParameter(0));
summary_5080_0_e->SetBinError(1,fit_5080_15GeV_e->GetParError(0));
summary_5080_0_e->SetBinContent(2,fit_5080_25GeV_e->GetParameter(0));
summary_5080_0_e->SetBinError(2,fit_5080_25GeV_e->GetParError(0));
summary_5080_0_e->SetBinContent(3,fit_5080_35GeV_e->GetParameter(0));
summary_5080_0_e->SetBinError(3,fit_5080_35GeV_e->GetParError(0));


TPaveText *l1_summary_0_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1_summary_0_e->SetTextSize(0.03);
l1_summary_0_e->SetFillColor(0);
l1_summary_0_e->SetBorderSize(0);
l1_summary_0_e->SetShadowColor(0);
l1_summary_0_e->SetTextAlign(21);
l1_summary_0_e->AddText("");
l1_summary_0_e->AddText("Pb+Pb");
l1_summary_0_e->AddText("Tight Photons");
l1_summary_0_e->AddText("p^{#gamma}_{T} > 15 GeV");
l1_summary_0_e->AddText("HLT Loose Photons");
l1_summary_0_e->AddText("DeltaR < 0.15");
l1_summary_0_e->AddText("1.56 < |#eta| < 2.37");
l1_summary_0_e->SetTextSize(0.036);




summary_010_0_e->Draw();
summary_1030_0_e->Draw("same");
summary_3050_0_e->Draw("same");
summary_5080_0_e->Draw("same");

l1_summary_0_e->Draw("same");


new_legend_e->AddEntry("summary_010_0_e","0-10%","pe");
new_legend_e->AddEntry("summary_1030_0_e","10-30%","pe");
new_legend_e->AddEntry("summary_3050_0_e","30-50%","pe");
new_legend_e->AddEntry("summary_5080_0_e","50-80%","pe");
new_legend_e->Draw("same");

//-------This is for the parameter 1
TCanvas *summary_new_2_e = new TCanvas("summary_new_2_e","summary_new_2_e",600,500);
TLegend *new_legend_2_e = new TLegend(0.6,0.5,0.85,0.6);

//-----0-10%
TH1D *summary_010_1_e = new TH1D("summary_010_1_e","",3,10,40);
summary_010_1_e->SetMarkerStyle(22);
summary_010_1_e->SetMarkerColor(kBlack);
summary_010_1_e->SetLineColor(kBlack);
summary_010_1_e->SetName("summary_010_1_e");
summary_010_1_e->GetXaxis()->SetTitle("L1 p_{T} [GeV]");
summary_010_1_e->GetYaxis()->SetTitle("Width");

summary_010_1_e->SetBinContent(1,fit_010_15GeV_e->GetParameter(1));
summary_010_1_e->SetBinError(1,fit_010_15GeV_e->GetParError(1));
summary_010_1_e->SetBinContent(2,fit_010_25GeV_e->GetParameter(1));
summary_010_1_e->SetBinError(2,fit_010_25GeV_e->GetParError(1));
summary_010_1_e->SetBinContent(3,fit_010_35GeV_e->GetParameter(1));
summary_010_1_e->SetBinError(3,fit_010_35GeV_e->GetParError(1));

//------10-30%
TH1D *summary_1030_1_e = new TH1D("summary_1030_1_e","",3,10,40);
summary_1030_1_e->SetMarkerStyle(21);
summary_1030_1_e->SetMarkerColor(kRed);
summary_1030_1_e->SetLineColor(kRed);
summary_1030_1_e->SetName("summary_1030_1_e");

summary_1030_1_e->SetBinContent(1,fit_1030_15GeV_e->GetParameter(1));
summary_1030_1_e->SetBinError(1,fit_1030_15GeV_e->GetParError(1));
summary_1030_1_e->SetBinContent(2,fit_1030_25GeV_e->GetParameter(1));
summary_1030_1_e->SetBinError(2,fit_1030_25GeV_e->GetParError(1));
summary_1030_1_e->SetBinContent(3,fit_1030_35GeV_e->GetParameter(1));
summary_1030_1_e->SetBinError(3,fit_1030_35GeV_e->GetParError(1));

//------30-50%
TH1D *summary_3050_1_e = new TH1D("summary_3050_1_e","",3,10,40);
summary_3050_1_e->SetMarkerStyle(33);
summary_3050_1_e->SetMarkerColor(kMagenta);
summary_3050_1_e->SetLineColor(kMagenta);
summary_3050_1_e->SetName("summary_3050_1_e");

summary_3050_1_e->SetBinContent(1,fit_3050_15GeV_e->GetParameter(1));
summary_3050_1_e->SetBinError(1,fit_3050_15GeV_e->GetParError(1));
summary_3050_1_e->SetBinContent(2,fit_3050_25GeV_e->GetParameter(1));
summary_3050_1_e->SetBinError(2,fit_3050_25GeV_e->GetParError(1));
summary_3050_1_e->SetBinContent(3,fit_3050_35GeV_e->GetParameter(1));
summary_3050_1_e->SetBinError(3,fit_3050_35GeV_e->GetParError(1));

//------50-80%
TH1D *summary_5080_1_e = new TH1D("summary_5080_1_e","",3,10,40);
summary_5080_1_e->SetMarkerStyle(20);
summary_5080_1_e->SetMarkerColor(kBlue);
summary_5080_1_e->SetLineColor(kBlue);
summary_5080_1_e->SetName("summary_5080_1_e");

summary_5080_1_e->SetBinContent(1,fit_5080_15GeV_e->GetParameter(1));
summary_5080_1_e->SetBinError(1,fit_5080_15GeV_e->GetParError(1));
summary_5080_1_e->SetBinContent(2,fit_5080_25GeV_e->GetParameter(1));
summary_5080_1_e->SetBinError(2,fit_5080_25GeV_e->GetParError(1));
summary_5080_1_e->SetBinContent(3,fit_5080_35GeV_e->GetParameter(1));
summary_5080_1_e->SetBinError(3,fit_5080_35GeV_e->GetParError(1));


TPaveText *l1_summary_1_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1_summary_1_e->SetTextSize(0.03);
l1_summary_1_e->SetFillColor(0);
l1_summary_1_e->SetBorderSize(0);
l1_summary_1_e->SetShadowColor(0);
l1_summary_1_e->SetTextAlign(21);
l1_summary_1_e->AddText("");
l1_summary_1_e->AddText("Pb+Pb");
l1_summary_1_e->AddText("Tight Photons");
l1_summary_1_e->AddText("p^{#gamma}_{T} > 15 GeV");
l1_summary_1_e->AddText("HLT Loose Photons");
l1_summary_1_e->AddText("DeltaR < 0.15");
l1_summary_1_e->AddText("1.56 < |#eta| < 2.37");
l1_summary_1_e->SetTextSize(0.036);




summary_010_1_e->Draw();
summary_1030_1_e->Draw("same");
summary_3050_1_e->Draw("same");
summary_5080_1_e->Draw("same");

l1_summary_1_e->Draw("same");


new_legend_2_e->AddEntry("summary_010_1_e","0-10%","pe");
new_legend_2_e->AddEntry("summary_1030_1_e","10-30%","pe");
new_legend_2_e->AddEntry("summary_3050_1_e","30-50%","pe");
new_legend_2_e->AddEntry("summary_5080_1_e","50-80%","pe");
new_legend_2_e->Draw("same");

//--------------------------------------------------------------//
//-----------------------HLT Efficiency-------------------------//
//--------------------------------------------------------------//

//------HLT Efficiencies
//------0-10% Centrality
//------Barrel & Endcaps

//------Creating Canvas & Legend
TCanvas *new_hlt_eff_centrality = new TCanvas("new_hlt_eff_centrality","new_hlt_eff_centrality",600,500);
TLegend *new_lg_hlt_cent = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *HLT_eff_010_15GeV = new TGraphAsymmErrors();
HLT_eff_010_15GeV->SetName("HLT_eff_010_15GeV");                                                                                                                           
HLT_eff_010_15GeV->GetXaxis()->SetTitle("p_{T} [GeV]");
HLT_eff_010_15GeV->SetMarkerStyle(24);
HLT_eff_010_15GeV->SetTitle("HLT Efficiency");
HLT_eff_010_15GeV->SetLineColor(kRed);
HLT_eff_010_15GeV->SetMarkerColor(kRed);
HLT_eff_010_15GeV->BayesDivide(photon_w_HLT_15GeV_010,all_photons_010);
HLT_eff_010_15GeV->SetMaximum(1.1);
HLT_eff_010_15GeV->SetMinimum(0);

TGraphAsymmErrors *HLT_eff_010_25GeV = new TGraphAsymmErrors();
HLT_eff_010_25GeV ->SetName("HLT_eff_010_25GeV");
HLT_eff_010_25GeV->SetMarkerStyle(25);
HLT_eff_010_25GeV->SetMarkerColor(kBlue);
HLT_eff_010_25GeV->SetLineColor(kBlue);
HLT_eff_010_25GeV->BayesDivide(photon_w_HLT_25GeV_010,all_photons_010);

TGraphAsymmErrors *HLT_eff_010_35GeV = new TGraphAsymmErrors();
HLT_eff_010_35GeV->SetName("HLT_eff_010_35GeV");
HLT_eff_010_35GeV->SetMarkerStyle(26);
HLT_eff_010_35GeV->SetMarkerColor(kMagenta);
HLT_eff_010_35GeV->SetLineColor(kMagenta);
HLT_eff_010_35GeV->BayesDivide(photon_w_HLT_35GeV_010,all_photons_010);


TPaveText *hlt_photon_0 = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_0->SetTextSize(0.03);
hlt_photon_0->SetFillColor(0);
hlt_photon_0->SetBorderSize(0);
hlt_photon_0->SetShadowColor(0);
hlt_photon_0->SetTextAlign(21);
hlt_photon_0->AddText("");
hlt_photon_0->AddText("Pb+Pb");
hlt_photon_0->AddText("0-10%");
hlt_photon_0->AddText("Tight Photons");
hlt_photon_0->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_photon_0->AddText("HLT Loose Photons");
hlt_photon_0->AddText("DeltaR < 0.15");
hlt_photon_0->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_010_15GeV_hlt = new TF1("fit_010_15GeV_hlt", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_15GeV_010->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_010->GetXaxis()->GetNbins()+1));
TF1* fit_010_25GeV_hlt = new TF1("fit_010_25GeV_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_25GeV_010->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_010->GetXaxis()->GetNbins()+1));
TF1* fit_010_35GeV_hlt = new TF1("fit_010_35GeV_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_35GeV_010->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_010->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_010_15GeV_hlt->SetParameter(0,15.0);		
fit_010_15GeV_hlt->SetParameter(1,2.0);

fit_010_25GeV_hlt->SetParameter(0,25.0);
fit_010_25GeV_hlt->SetParameter(1,2.0);

fit_010_35GeV_hlt->SetParameter(0,35.0);
fit_010_35GeV_hlt->SetParameter(1,2.0);



HLT_eff_010_15GeV->Fit("fit_010_15GeV_hlt","M Q N","",0,photon_w_HLT_15GeV_010->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_010->GetXaxis()->GetNbins()+1));
HLT_eff_010_25GeV->Fit("fit_010_25GeV_hlt", "M Q N","",0,photon_w_HLT_25GeV_010->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_010->GetXaxis()->GetNbins()+1));
HLT_eff_010_35GeV->Fit("fit_010_35GeV_hlt","M Q N","",0,photon_w_HLT_35GeV_010->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_010->GetXaxis()->GetNbins()+1));


fit_010_15GeV_hlt->SetLineStyle(2);
fit_010_15GeV_hlt->SetLineColor(kRed);
fit_010_15GeV_hlt->SetMarkerColor(kRed);

fit_010_25GeV_hlt->SetLineStyle(2);
fit_010_25GeV_hlt->SetLineColor(kBlue);
fit_010_25GeV_hlt->SetMarkerColor(kBlue);

fit_010_35GeV_hlt->SetLineStyle(2);
fit_010_35GeV_hlt->SetLineColor(kMagenta);
fit_010_35GeV_hlt->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_hlt = new TLine(0.4,1,169,1);
line_100_hlt->SetLineColor(kBlack);
line_100_hlt->SetLineWidth(1);




HLT_eff_010_15GeV->Draw("A P");
line_100_hlt->Draw("same");

HLT_eff_010_25GeV->Draw("P same");
HLT_eff_010_35GeV->Draw("P same");


fit_010_15GeV_hlt->Draw("same");
fit_010_25GeV_hlt->Draw("same");
fit_010_35GeV_hlt->Draw("same");


hlt_photon_0->Draw("same");

new_lg_hlt_cent->AddEntry("HLT_eff_010_15GeV","HLT p_{T} > 15 GeV","pe");
new_lg_hlt_cent->AddEntry("HLT_eff_010_25GeV","HLT p_{T} > 25 GeV","pe");
new_lg_hlt_cent->AddEntry("HLT_eff_010_35GeV","HLT p_{T} > 35 GeV","pe");



new_lg_hlt_cent->Draw("same");



//------10-30% Centrality
//------Barrel & Endcaps

//------Creating Canvas & Legend
TCanvas *new_hlt_eff_centrality_1030 = new TCanvas("new_hlt_eff_centrality_1030","new_hlt_eff_centrality_1030",600,500);
TLegend *new_lg_hlt_cent_1030 = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *HLT_eff_1030_15GeV = new TGraphAsymmErrors();
HLT_eff_1030_15GeV->SetName("HLT_eff_1030_15GeV");                                                                                                                           
HLT_eff_1030_15GeV->GetXaxis()->SetTitle("p_{T} [GeV]");
HLT_eff_1030_15GeV->SetMarkerStyle(24);
HLT_eff_1030_15GeV->SetTitle("HLT Efficiency");
HLT_eff_1030_15GeV->SetLineColor(kRed);
HLT_eff_1030_15GeV->SetMarkerColor(kRed);
HLT_eff_1030_15GeV->BayesDivide(photon_w_HLT_15GeV_1030,all_photons_1030);
HLT_eff_1030_15GeV->SetMaximum(1.1);
HLT_eff_1030_15GeV->SetMinimum(0);

TGraphAsymmErrors *HLT_eff_1030_25GeV = new TGraphAsymmErrors();
HLT_eff_1030_25GeV->SetName("HLT_eff_1030_25GeV");
HLT_eff_1030_25GeV->SetMarkerStyle(25);
HLT_eff_1030_25GeV->SetMarkerColor(kBlue);
HLT_eff_1030_25GeV->SetLineColor(kBlue);
HLT_eff_1030_25GeV->BayesDivide(photon_w_HLT_25GeV_1030,all_photons_1030);

TGraphAsymmErrors *HLT_eff_1030_35GeV = new TGraphAsymmErrors();
HLT_eff_1030_35GeV->SetName("HLT_eff_1030_35GeV");
HLT_eff_1030_35GeV->SetMarkerStyle(26);
HLT_eff_1030_35GeV->SetMarkerColor(kMagenta);
HLT_eff_1030_35GeV->SetLineColor(kMagenta);
HLT_eff_1030_35GeV->BayesDivide(photon_w_HLT_35GeV_1030,all_photons_1030);


TPaveText *hlt_photon_1030 = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_1030->SetTextSize(0.03);
hlt_photon_1030->SetFillColor(0);
hlt_photon_1030->SetBorderSize(0);
hlt_photon_1030->SetShadowColor(0);
hlt_photon_1030->SetTextAlign(21);
hlt_photon_1030->AddText("");
hlt_photon_1030->AddText("Pb+Pb");
hlt_photon_1030->AddText("10-30%");
hlt_photon_1030->AddText("Tight Photons");
hlt_photon_1030->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_photon_1030->AddText("HLT Loose Photons");
hlt_photon_1030->AddText("DeltaR < 0.15");
hlt_photon_1030->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_1030_15GeV_hlt = new TF1("fit_1030_15GeV_hlt", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_15GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_1030->GetXaxis()->GetNbins()+1));
TF1* fit_1030_25GeV_hlt = new TF1("fit_1030_25GeV_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_25GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_1030->GetXaxis()->GetNbins()+1));
TF1* fit_1030_35GeV_hlt = new TF1("fit_1030_35GeV_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_35GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_1030->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_1030_15GeV_hlt->SetParameter(0,15.0);		
fit_1030_15GeV_hlt->SetParameter(1,2.0);

fit_1030_25GeV_hlt->SetParameter(0,25.0);
fit_1030_25GeV_hlt->SetParameter(1,2.0);

fit_1030_35GeV_hlt->SetParameter(0,35.0);
fit_1030_35GeV_hlt->SetParameter(1,2.0);



HLT_eff_1030_15GeV->Fit("fit_1030_15GeV_hlt","M Q N","",0,photon_w_HLT_15GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_1030->GetXaxis()->GetNbins()+1));
HLT_eff_1030_25GeV->Fit("fit_1030_25GeV_hlt", "M Q N","",0,photon_w_HLT_25GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_1030->GetXaxis()->GetNbins()+1));
HLT_eff_1030_35GeV->Fit("fit_1030_35GeV_hlt","M Q N","",0,photon_w_HLT_35GeV_1030->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_1030->GetXaxis()->GetNbins()+1));


fit_1030_15GeV_hlt->SetLineStyle(2);
fit_1030_15GeV_hlt->SetLineColor(kRed);
fit_1030_15GeV_hlt->SetMarkerColor(kRed);

fit_1030_25GeV_hlt->SetLineStyle(2);
fit_1030_25GeV_hlt->SetLineColor(kBlue);
fit_1030_25GeV->SetMarkerColor(kBlue);

fit_1030_35GeV_hlt->SetLineStyle(2);
fit_1030_35GeV_hlt->SetLineColor(kMagenta);
fit_1030_35GeV_hlt->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_1_hlt = new TLine(0.5,1,165,1);
line_100_1_hlt->SetLineColor(kBlack);
line_100_1_hlt->SetLineWidth(1);




HLT_eff_1030_15GeV->Draw("A P");
line_100_1_hlt->Draw("same");

HLT_eff_1030_25GeV->Draw("P same");
HLT_eff_1030_35GeV->Draw("P same");


fit_1030_15GeV_hlt->Draw("same");
fit_1030_25GeV_hlt->Draw("same");
fit_1030_35GeV_hlt->Draw("same");


hlt_photon_1030->Draw("same");

new_lg_hlt_cent_1030->AddEntry("HLT_eff_1030_15GeV","HLT p_{T} > 15 GeV","pe");
new_lg_hlt_cent_1030->AddEntry("HLT_eff_1030_25GeV","HLT p_{T} > 25 GeV","pe");
new_lg_hlt_cent_1030->AddEntry("HLT_eff_1030_35GeV","HLT p_{T} > 35 GeV","pe");



new_lg_hlt_cent_1030->Draw("same");

//------30-50% Centrality
//------Barrel & Endcaps

//------Creating Canvas & Legend
TCanvas *new_hlt_eff_centrality_3050 = new TCanvas("new_hlt_eff_centrality_3050","new_hlt_eff_centrality_3050",600,500);
TLegend *new_lg_hlt_cent_3050 = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *HLT_eff_3050_15GeV = new TGraphAsymmErrors();
HLT_eff_3050_15GeV->SetName("HLT_eff_3050_15GeV");                                                                                                                           
HLT_eff_3050_15GeV->GetXaxis()->SetTitle("p_{T} [GeV]");
HLT_eff_3050_15GeV->SetMarkerStyle(24);
HLT_eff_3050_15GeV->SetTitle("HLT Efficiency");
HLT_eff_3050_15GeV->SetLineColor(kRed);
HLT_eff_3050_15GeV->SetMarkerColor(kRed);
HLT_eff_3050_15GeV->BayesDivide(photon_w_HLT_15GeV_3050,all_photons_3050);
HLT_eff_3050_15GeV->SetMaximum(1.1);
HLT_eff_3050_15GeV->SetMinimum(0);

TGraphAsymmErrors *HLT_eff_3050_25GeV = new TGraphAsymmErrors();
HLT_eff_3050_25GeV->SetName("HLT_eff_3050_25GeV");
HLT_eff_3050_25GeV->SetMarkerStyle(25);
HLT_eff_3050_25GeV->SetMarkerColor(kBlue);
HLT_eff_3050_25GeV->SetLineColor(kBlue);
HLT_eff_3050_25GeV->BayesDivide(photon_w_HLT_25GeV_3050,all_photons_3050);

TGraphAsymmErrors *HLT_eff_3050_35GeV = new TGraphAsymmErrors();
HLT_eff_3050_35GeV->SetName("HLT_eff_3050_35GeV");
HLT_eff_3050_35GeV->SetMarkerStyle(26);
HLT_eff_3050_35GeV->SetMarkerColor(kMagenta);
HLT_eff_3050_35GeV->SetLineColor(kMagenta);
HLT_eff_3050_35GeV->BayesDivide(photon_w_HLT_35GeV_3050,all_photons_3050);


TPaveText *hlt_photon_3050 = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_3050->SetTextSize(0.03);
hlt_photon_3050->SetFillColor(0);
hlt_photon_3050->SetBorderSize(0);
hlt_photon_3050->SetShadowColor(0);
hlt_photon_3050->SetTextAlign(21);
hlt_photon_3050->AddText("");
hlt_photon_3050->AddText("Pb+Pb");
hlt_photon_3050->AddText("30-50%");
hlt_photon_3050->AddText("Tight Photons");
hlt_photon_3050->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_photon_3050->AddText("DeltaR < 0.15");
hlt_photon_3050->AddText("HLT Loose Photons");
hlt_photon_3050->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_3050_15GeV_hlt = new TF1("fit_3050_15GeV_hlt", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_15GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_15GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_3050->GetXaxis()->GetNbins()+1));
TF1* fit_3050_25GeV_hlt = new TF1("fit_3050_25GeV_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_25GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_25GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_3050->GetXaxis()->GetNbins()+1));
TF1* fit_3050_35GeV_hlt = new TF1("fit_3050_35GeV_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_35GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_35GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_3050->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_3050_15GeV_hlt->SetParameter(0,15.0);		
fit_3050_15GeV_hlt->SetParameter(1,2.0);

fit_3050_25GeV_hlt->SetParameter(0,25.0);
fit_3050_25GeV_hlt->SetParameter(1,2.0);

fit_3050_35GeV_hlt->SetParameter(0,35.0);
fit_3050_35GeV_hlt->SetParameter(1,2.0);



HLT_eff_3050_15GeV->Fit("fit_3050_15GeV_hlt","M Q N","",photon_w_HLT_15GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_15GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_3050->GetXaxis()->GetNbins()+1));
HLT_eff_3050_25GeV->Fit("fit_3050_25GeV_hlt", "M Q N","",photon_w_HLT_25GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_25GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_3050->GetXaxis()->GetNbins()+1));
HLT_eff_3050_35GeV->Fit("fit_3050_35GeV_hlt","M Q N","",photon_w_HLT_35GeV_3050->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_35GeV_3050->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_3050->GetXaxis()->GetNbins()+1));


fit_3050_15GeV_hlt->SetLineStyle(2);
fit_3050_15GeV_hlt->SetLineColor(kRed);
fit_3050_15GeV_hlt->SetMarkerColor(kRed);

fit_3050_25GeV_hlt->SetLineStyle(2);
fit_3050_25GeV_hlt->SetLineColor(kBlue);
fit_3050_25GeV_hlt->SetMarkerColor(kBlue);

fit_3050_35GeV_hlt->SetLineStyle(2);
fit_3050_35GeV_hlt->SetLineColor(kMagenta);
fit_3050_35GeV_hlt->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_2_hlt = new TLine(2.0,1,153,1);
line_100_2_hlt->SetLineColor(kBlack);
line_100_2_hlt->SetLineWidth(1);




HLT_eff_3050_15GeV->Draw("A P");
line_100_2_hlt->Draw("same");

HLT_eff_3050_25GeV->Draw("P same");
HLT_eff_3050_35GeV->Draw("P same");


fit_3050_15GeV_hlt->Draw("same");
fit_3050_25GeV_hlt->Draw("same");
fit_3050_35GeV_hlt->Draw("same");


hlt_photon_3050->Draw("same");

new_lg_hlt_cent_3050->AddEntry("HLT_eff_3050_15GeV","HLT p_{T} > 15 GeV","pe");
new_lg_hlt_cent_3050->AddEntry("HLT_eff_3050_25GeV","HLT p_{T} > 25 GeV","pe");
new_lg_hlt_cent_3050->AddEntry("HLT_eff_3050_35GeV","HLT p_{T} > 35 GeV","pe");



new_lg_hlt_cent_3050->Draw("same");


//------50-80% Centrality
//------Barrel & Endcaps

//------Creating Canvas & Legend
TCanvas *new__eff_centrality_5080_hlt = new TCanvas("new_hlt_eff_centrality_5080","new_hlt_eff_centrality_5080",600,500);
TLegend *new_lg_hlt_cent_5080 = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *HLT_eff_5080_15GeV = new TGraphAsymmErrors();
HLT_eff_5080_15GeV->SetName("HLT_eff_5080_15GeV");                                                                                                                           
HLT_eff_5080_15GeV->GetXaxis()->SetTitle("p_{T} [GeV]");
HLT_eff_5080_15GeV->SetMarkerStyle(24);
HLT_eff_5080_15GeV->SetTitle("HLT Efficiency");
HLT_eff_5080_15GeV->SetLineColor(kRed);
HLT_eff_5080_15GeV->SetMarkerColor(kRed);
HLT_eff_5080_15GeV->BayesDivide(photon_w_HLT_15GeV_5080,all_photons_5080);
HLT_eff_5080_15GeV->SetMaximum(1.1);
HLT_eff_5080_15GeV->SetMinimum(0);

TGraphAsymmErrors *HLT_eff_5080_25GeV = new TGraphAsymmErrors();
HLT_eff_5080_25GeV->SetName("HLT_eff_5080_25GeV");
HLT_eff_5080_25GeV->SetMarkerStyle(25);
HLT_eff_5080_25GeV->SetMarkerColor(kBlue);
HLT_eff_5080_25GeV->SetLineColor(kBlue);
HLT_eff_5080_25GeV->BayesDivide(photon_w_HLT_25GeV_5080,all_photons_5080);

TGraphAsymmErrors *HLT_eff_5080_35GeV = new TGraphAsymmErrors();
HLT_eff_5080_35GeV->SetName("HLT_eff_5080_35GeV");
HLT_eff_5080_35GeV->SetMarkerStyle(26);
HLT_eff_5080_35GeV->SetMarkerColor(kMagenta);
HLT_eff_5080_35GeV->SetLineColor(kMagenta);
HLT_eff_5080_35GeV->BayesDivide(photon_w_HLT_35GeV_5080,all_photons_5080);


TPaveText *hlt_photon_5080 = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_5080->SetTextSize(0.03);
hlt_photon_5080->SetFillColor(0);
hlt_photon_5080->SetBorderSize(0);
hlt_photon_5080->SetShadowColor(0);
hlt_photon_5080->SetTextAlign(21);
hlt_photon_5080->AddText("");
hlt_photon_5080->AddText("Pb+Pb");
hlt_photon_5080->AddText("50-80%");
hlt_photon_5080->AddText("Tight Photons");
hlt_photon_5080->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_photon_5080->AddText("DeltaR < 0.15");
hlt_photon_5080->AddText("HLT Loose Photons");
hlt_photon_5080->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_5080_15GeV_hlt = new TF1("fit_5080_15GeV_hlt", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_15GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_15GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_5080->GetXaxis()->GetNbins()+1));
TF1* fit_5080_25GeV_hlt = new TF1("fit_5080_25GeV_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_25GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_25GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_5080->GetXaxis()->GetNbins()+1));
TF1* fit_5080_35GeV_hlt = new TF1("fit_5080_35GeV_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_35GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_35GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_5080->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_5080_15GeV_hlt->SetParameter(0,15.0);		
fit_5080_15GeV_hlt->SetParameter(1,2.0);

fit_5080_25GeV_hlt->SetParameter(0,25.0);
fit_5080_25GeV_hlt->SetParameter(1,2.0);

fit_5080_35GeV_hlt->SetParameter(0,35.0);
fit_5080_35GeV_hlt->SetParameter(1,2.0);



HLT_eff_5080_15GeV->Fit("fit_5080_15GeV_hlt","M Q N","",photon_w_HLT_15GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_15GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_5080->GetXaxis()->GetNbins()+1));
HLT_eff_5080_25GeV->Fit("fit_5080_25GeV_hlt", "M Q N","",photon_w_HLT_25GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_25GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_5080->GetXaxis()->GetNbins()+1));
HLT_eff_5080_35GeV->Fit("fit_5080_35GeV_hlt","M Q N","",photon_w_HLT_35GeV_5080->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_35GeV_5080->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_5080->GetXaxis()->GetNbins()+1));


fit_5080_15GeV_hlt->SetLineStyle(2);
fit_5080_15GeV_hlt->SetLineColor(kRed);
fit_5080_15GeV_hlt->SetMarkerColor(kRed);

fit_5080_25GeV_hlt->SetLineStyle(2);
fit_5080_25GeV_hlt->SetLineColor(kBlue);
fit_5080_25GeV_hlt->SetMarkerColor(kBlue);

fit_5080_35GeV_hlt->SetLineStyle(2);
fit_5080_35GeV_hlt->SetLineColor(kMagenta);
fit_5080_35GeV_hlt->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_3_hlt = new TLine(2.0,1,200,1);
line_100_3_hlt->SetLineColor(kBlack);
line_100_3_hlt->SetLineWidth(1);




HLT_eff_5080_15GeV->Draw("A P");
line_100_3_hlt->Draw("same");

HLT_eff_5080_25GeV->Draw("P same");
HLT_eff_5080_35GeV->Draw("P same");


fit_5080_15GeV_hlt->Draw("same");
fit_5080_25GeV_hlt->Draw("same");
fit_5080_35GeV_hlt->Draw("same");


hlt_photon_5080->Draw("same");

new_lg_hlt_cent_5080->AddEntry("HLT_eff_5080_15GeV","HLT p_{T} > 15 GeV","pe");
new_lg_hlt_cent_5080->AddEntry("HLT_eff_5080_25GeV","HLT p_{T} > 25 GeV","pe");
new_lg_hlt_cent_5080->AddEntry("HLT_eff_5080_35GeV","HLT p_{T} > 35 GeV","pe");



new_lg_hlt_cent_5080->Draw("same");

//-------Summary Plots
//-------Barrel & Endcaps
//-------This is for the parameter 0
TCanvas *summary_new_hlt = new TCanvas("summary_new_hlt","summary_new_hlt",600,500);
TLegend *new_legend_hlt = new TLegend(0.6,0.5,0.85,0.6);

//-----0-10%
TH1D *summary_010_0_hlt = new TH1D("summary_010_0_hlt","",3,10,40);
summary_010_0_hlt->SetMarkerStyle(22);
summary_010_0_hlt->SetMarkerColor(kBlack);
summary_010_0_hlt->SetLineColor(kBlack);
summary_010_0_hlt->SetName("summary_010_0_hlt");
summary_010_0_hlt->GetXaxis()->SetTitle("HLT p_{T} [GeV]");
summary_010_0_hlt->GetYaxis()->SetTitle("Value at 50 Percent Efficiency [GeV]");

summary_010_0_hlt->SetBinContent(1,fit_010_15GeV_hlt->GetParameter(0));
summary_010_0_hlt->SetBinError(1,fit_010_15GeV_hlt->GetParError(0));
summary_010_0_hlt->SetBinContent(2,fit_010_25GeV_hlt->GetParameter(0));
summary_010_0_hlt->SetBinError(2,fit_010_25GeV_hlt->GetParError(0));
summary_010_0_hlt->SetBinContent(3,fit_010_35GeV_hlt->GetParameter(0));
summary_010_0_hlt->SetBinError(3,fit_010_35GeV_hlt->GetParError(0));

//------10-30%
TH1D *summary_1030_0_hlt = new TH1D("summary_1030_0_hlt","",3,10,40);
summary_1030_0_hlt->SetMarkerStyle(21);
summary_1030_0_hlt->SetMarkerColor(kRed);
summary_1030_0_hlt->SetLineColor(kRed);
summary_1030_0_hlt->SetName("summary_1030_0_hlt");

summary_1030_0_hlt->SetBinContent(1,fit_1030_15GeV_hlt->GetParameter(0));
summary_1030_0_hlt->SetBinError(1,fit_1030_15GeV_hlt->GetParError(0));
summary_1030_0_hlt->SetBinContent(2,fit_1030_25GeV_hlt->GetParameter(0));
summary_1030_0_hlt->SetBinError(2,fit_1030_25GeV_hlt->GetParError(0));
summary_1030_0_hlt->SetBinContent(3,fit_1030_35GeV_hlt->GetParameter(0));
summary_1030_0_hlt->SetBinError(3,fit_1030_35GeV_hlt->GetParError(0));

//------30-50%
TH1D *summary_3050_0_hlt = new TH1D("summary_3050_0_hlt","",3,10,40);
summary_3050_0_hlt->SetMarkerStyle(33);
summary_3050_0_hlt->SetMarkerColor(kMagenta);
summary_3050_0_hlt->SetLineColor(kMagenta);
summary_3050_0_hlt->SetName("summary_3050_0_hlt");

summary_3050_0_hlt->SetBinContent(1,fit_3050_15GeV_hlt->GetParameter(0));
summary_3050_0_hlt->SetBinError(1,fit_3050_15GeV_hlt->GetParError(0));
summary_3050_0_hlt->SetBinContent(2,fit_3050_25GeV_hlt->GetParameter(0));
summary_3050_0_hlt->SetBinError(2,fit_3050_25GeV_hlt->GetParError(0));
summary_3050_0_hlt->SetBinContent(3,fit_3050_35GeV_hlt->GetParameter(0));
summary_3050_0_hlt->SetBinError(3,fit_3050_35GeV_hlt->GetParError(0));

//------50-80%
TH1D *summary_5080_0_hlt = new TH1D("summary_5080_0_hlt","",3,10,40);
summary_5080_0_hlt->SetMarkerStyle(20);
summary_5080_0_hlt->SetMarkerColor(kBlue);
summary_5080_0_hlt->SetLineColor(kBlue);
summary_5080_0_hlt->SetName("summary_5080_0_hlt");

summary_5080_0_hlt->SetBinContent(1,fit_5080_15GeV_hlt->GetParameter(0));
summary_5080_0_hlt->SetBinError(1,fit_5080_15GeV_hlt->GetParError(0));
summary_5080_0_hlt->SetBinContent(2,fit_5080_25GeV_hlt->GetParameter(0));
summary_5080_0_hlt->SetBinError(2,fit_5080_25GeV_hlt->GetParError(0));
summary_5080_0_hlt->SetBinContent(3,fit_5080_35GeV_hlt->GetParameter(0));
summary_5080_0_hlt->SetBinError(3,fit_5080_35GeV_hlt->GetParError(0));


TPaveText *hlt_summary_0 = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_summary_0->SetTextSize(0.03);
hlt_summary_0->SetFillColor(0);
hlt_summary_0->SetBorderSize(0);
hlt_summary_0->SetShadowColor(0);
hlt_summary_0->SetTextAlign(21);
hlt_summary_0->AddText("");
hlt_summary_0->AddText("Pb+Pb");
hlt_summary_0->AddText("Tight Photons");
hlt_summary_0->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_summary_0->AddText("HLT Loose Photons");
hlt_summary_0->AddText("DeltaR < 0.15");
hlt_summary_0->SetTextSize(0.036);




summary_010_0_hlt->Draw();
summary_1030_0_hlt->Draw("same");
summary_3050_0_hlt->Draw("same");
summary_5080_0_hlt->Draw("same");

hlt_summary_0->Draw("same");


new_legend_hlt->AddEntry("summary_010_0_hlt","0-10%","pe");
new_legend_hlt->AddEntry("summary_1030_0_hlt","10-30%","pe");
new_legend_hlt->AddEntry("summary_3050_0_hlt","30-50%","pe");
new_legend_hlt->AddEntry("summary_5080_0_hlt","50-80%","pe");
new_legend_hlt->Draw("same");

//-------This is for the parameter 1
TCanvas *summary_new_2_hlt = new TCanvas("summary_new_2_hlt","summary_new_2_hlt",600,500);
TLegend *new_legend_2_hlt = new TLegend(0.6,0.5,0.85,0.6);

//-----0-10%
TH1D *summary_010_1_hlt = new TH1D("summary_010_1_hlt","",3,10,40);
summary_010_1_hlt->SetMarkerStyle(22);
summary_010_1_hlt->SetMarkerColor(kBlack);
summary_010_1_hlt->SetLineColor(kBlack);
summary_010_1_hlt->SetName("summary_010_1_hlt");
summary_010_1_hlt->GetXaxis()->SetTitle("HLT p_{T} [GeV]");
summary_010_1_hlt->GetYaxis()->SetTitle("Width");

summary_010_1_hlt->SetBinContent(1,fit_010_15GeV_hlt->GetParameter(1));
summary_010_1_hlt->SetBinError(1,fit_010_15GeV_hlt->GetParError(1));
summary_010_1_hlt->SetBinContent(2,fit_010_25GeV_hlt->GetParameter(1));
summary_010_1_hlt->SetBinError(2,fit_010_25GeV_hlt->GetParError(1));
summary_010_1_hlt->SetBinContent(3,fit_010_35GeV_hlt->GetParameter(1));
summary_010_1_hlt->SetBinError(3,fit_010_35GeV_hlt->GetParError(1));

//------10-30%
TH1D *summary_1030_1_hlt = new TH1D("summary_1030_1_hlt","",3,10,40);
summary_1030_1_hlt->SetMarkerStyle(21);
summary_1030_1_hlt->SetMarkerColor(kRed);
summary_1030_1_hlt->SetLineColor(kRed);
summary_1030_1_hlt->SetName("summary_1030_1_hlt");

summary_1030_1_hlt->SetBinContent(1,fit_1030_15GeV_hlt->GetParameter(1));
summary_1030_1_hlt->SetBinError(1,fit_1030_15GeV_hlt->GetParError(1));
summary_1030_1_hlt->SetBinContent(2,fit_1030_25GeV_hlt->GetParameter(1));
summary_1030_1_hlt->SetBinError(2,fit_1030_25GeV_hlt->GetParError(1));
summary_1030_1_hlt->SetBinContent(3,fit_1030_35GeV_hlt->GetParameter(1));
summary_1030_1_hlt->SetBinError(3,fit_1030_35GeV_hlt->GetParError(1));

//------30-50%
TH1D *summary_3050_1_hlt = new TH1D("summary_3050_1_hlt","",3,10,40);
summary_3050_1_hlt->SetMarkerStyle(33);
summary_3050_1_hlt->SetMarkerColor(kMagenta);
summary_3050_1_hlt->SetLineColor(kMagenta);
summary_3050_1_hlt->SetName("summary_3050_1_hlt");

summary_3050_1_hlt->SetBinContent(1,fit_3050_15GeV_hlt->GetParameter(1));
summary_3050_1_hlt->SetBinError(1,fit_3050_15GeV_hlt->GetParError(1));
summary_3050_1_hlt->SetBinContent(2,fit_3050_25GeV_hlt->GetParameter(1));
summary_3050_1_hlt->SetBinError(2,fit_3050_25GeV_hlt->GetParError(1));
summary_3050_1_hlt->SetBinContent(3,fit_3050_35GeV_hlt->GetParameter(1));
summary_3050_1_hlt->SetBinError(3,fit_3050_35GeV_hlt->GetParError(1));

//------50-80%
TH1D *summary_5080_1_hlt = new TH1D("summary_5080_1_hlt","",3,10,40);
summary_5080_1_hlt->SetMarkerStyle(20);
summary_5080_1_hlt->SetMarkerColor(kBlue);
summary_5080_1_hlt->SetLineColor(kBlue);
summary_5080_1_hlt->SetName("summary_5080_1_hlt");

summary_5080_1_hlt->SetBinContent(1,fit_5080_15GeV_hlt->GetParameter(1));
summary_5080_1_hlt->SetBinError(1,fit_5080_15GeV_hlt->GetParError(1));
summary_5080_1_hlt->SetBinContent(2,fit_5080_25GeV_hlt->GetParameter(1));
summary_5080_1_hlt->SetBinError(2,fit_5080_25GeV_hlt->GetParError(1));
summary_5080_1_hlt->SetBinContent(3,fit_5080_35GeV_hlt->GetParameter(1));
summary_5080_1_hlt->SetBinError(3,fit_5080_35GeV_hlt->GetParError(1));


TPaveText *hlt_summary_1 = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_summary_1->SetTextSize(0.03);
hlt_summary_1->SetFillColor(0);
hlt_summary_1->SetBorderSize(0);
hlt_summary_1->SetShadowColor(0);
hlt_summary_1->SetTextAlign(21);
hlt_summary_1->AddText("");
hlt_summary_1->AddText("Pb+Pb");
hlt_summary_1->AddText("Tight Photons");
hlt_summary_1->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_summary_1->AddText("HLT Loose Photons");
hlt_summary_1->AddText("DeltaR < 0.15");
hlt_summary_1->SetTextSize(0.036);




summary_010_1_hlt->Draw();
summary_1030_1_hlt->Draw("same");
summary_3050_1_hlt->Draw("same");
summary_5080_1_hlt->Draw("same");

hlt_summary_1->Draw("same");


new_legend_2_hlt->AddEntry("summary_010_1_hlt","0-10%","pe");
new_legend_2_hlt->AddEntry("summary_1030_1_hlt","10-30%","pe");
new_legend_2_hlt->AddEntry("summary_3050_1_hlt","30-50%","pe");
new_legend_2_hlt->AddEntry("summary_5080_1_hlt","50-80%","pe");
new_legend_2_hlt->Draw("same");

//------HLT Efficiencies
//------0-10% Centrality
//------Barrel

//------Creating Canvas & Legend
TCanvas *new_hlt_eff_centrality_b = new TCanvas("new_hlt_eff_centrality_b","new_hlt_eff_centrality_b",600,500);
TLegend *new_lg_hlt_cent_b = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *HLT_eff_010_15GeV_b = new TGraphAsymmErrors();
HLT_eff_010_15GeV_b->SetName("HLT_eff_010_15GeV_b");                                                                                                                           
HLT_eff_010_15GeV_b->GetXaxis()->SetTitle("p_{T} [GeV]");
HLT_eff_010_15GeV_b->SetMarkerStyle(24);
HLT_eff_010_15GeV_b->SetTitle("HLT Efficiency");
HLT_eff_010_15GeV_b->SetLineColor(kRed);
HLT_eff_010_15GeV_b->SetMarkerColor(kRed);
HLT_eff_010_15GeV_b->BayesDivide(photon_w_HLT_15GeV_010_b,all_photons_010_b);
HLT_eff_010_15GeV_b->SetMaximum(1.1);
HLT_eff_010_15GeV_b->SetMinimum(0);

TGraphAsymmErrors *HLT_eff_010_25GeV_b = new TGraphAsymmErrors();
HLT_eff_010_25GeV_b->SetName("HLT_eff_010_25GeV_b");
HLT_eff_010_25GeV_b->SetMarkerStyle(25);
HLT_eff_010_25GeV_b->SetMarkerColor(kBlue);
HLT_eff_010_25GeV_b->SetLineColor(kBlue);
HLT_eff_010_25GeV_b->BayesDivide(photon_w_HLT_25GeV_010_b,all_photons_010_b);

TGraphAsymmErrors *HLT_eff_010_35GeV_b = new TGraphAsymmErrors();
HLT_eff_010_35GeV_b->SetName("HLT_eff_010_35GeV_b");
HLT_eff_010_35GeV_b->SetMarkerStyle(26);
HLT_eff_010_35GeV_b->SetMarkerColor(kMagenta);
HLT_eff_010_35GeV_b->SetLineColor(kMagenta);
HLT_eff_010_35GeV_b->BayesDivide(photon_w_HLT_35GeV_010_b,all_photons_010_b);


TPaveText *hlt_photon_0_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_0_b->SetTextSize(0.03);
hlt_photon_0_b->SetFillColor(0);
hlt_photon_0_b->SetBorderSize(0);
hlt_photon_0_b->SetShadowColor(0);
hlt_photon_0_b->SetTextAlign(21);
hlt_photon_0_b->AddText("");
hlt_photon_0_b->AddText("Pb+Pb");
hlt_photon_0_b->AddText("0-10%");
hlt_photon_0_b->AddText("Tight Photons");
hlt_photon_0_b->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_photon_0_b->AddText("DeltaR < 0.15");
hlt_photon_0_b->AddText("HLT Loose Photon");
hlt_photon_0_b->AddText("|#eta| < 1.37");
hlt_photon_0_b->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_010_15GeV_b_hlt = new TF1("fit_010_15GeV_b_hlt", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_15GeV_010_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_010_b->GetXaxis()->GetNbins()+1));
TF1* fit_010_25GeV_b_hlt = new TF1("fit_010_25GeV_b_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_25GeV_010_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_010_b->GetXaxis()->GetNbins()+1));
TF1* fit_010_35GeV_b_hlt = new TF1("fit_010_35GeV_b_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_35GeV_010_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_010_b->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_010_15GeV_b_hlt->SetParameter(0,15.0);		
fit_010_15GeV_b_hlt->SetParameter(1,2.0);

fit_010_25GeV_b_hlt->SetParameter(0,25.0);
fit_010_25GeV_b_hlt->SetParameter(1,2.0);

fit_010_35GeV_b_hlt->SetParameter(0,35.0);
fit_010_35GeV_b_hlt->SetParameter(1,2.0);



HLT_eff_010_15GeV_b->Fit("fit_010_15GeV_b_hlt","M Q N","",0,photon_w_HLT_15GeV_010_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_010_b->GetXaxis()->GetNbins()+1));
HLT_eff_010_25GeV_b->Fit("fit_010_25GeV_b_hlt", "M Q N","",0,photon_w_HLT_25GeV_010_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_010_b->GetXaxis()->GetNbins()+1));
HLT_eff_010_35GeV_b->Fit("fit_010_35GeV_b_hlt","M Q N","",0,photon_w_HLT_35GeV_010_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_010_b->GetXaxis()->GetNbins()+1));


fit_010_15GeV_b_hlt->SetLineStyle(2);
fit_010_15GeV_b_hlt->SetLineColor(kRed);
fit_010_15GeV_b_hlt->SetMarkerColor(kRed);

fit_010_25GeV_b_hlt->SetLineStyle(2);
fit_010_25GeV_b_hlt->SetLineColor(kBlue);
fit_010_25GeV_b_hlt->SetMarkerColor(kBlue);

fit_010_35GeV_b_hlt->SetLineStyle(2);
fit_010_35GeV_b_hlt->SetLineColor(kMagenta);
fit_010_35GeV_b_hlt->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_b_hlt = new TLine(0.4,1,169,1);
line_100_b_hlt->SetLineColor(kBlack);
line_100_b_hlt->SetLineWidth(1);




HLT_eff_010_15GeV_b->Draw("A P");
line_100_b_hlt->Draw("same");

HLT_eff_010_25GeV_b->Draw("P same");
HLT_eff_010_35GeV_b->Draw("P same");


fit_010_15GeV_b_hlt->Draw("same");
fit_010_25GeV_b_hlt->Draw("same");
fit_010_35GeV_b_hlt->Draw("same");


hlt_photon_0_b->Draw("same");

new_lg_hlt_cent_b->AddEntry("HLT_eff_010_15GeV_b","HLT p_{T} > 15 GeV","pe");
new_lg_hlt_cent_b->AddEntry("HLT_eff_010_25GeV_b","HLT p_{T} > 25 GeV","pe");
new_lg_hlt_cent_b->AddEntry("HLT_eff_010_35GeV_b","HLT p_{T} > 35 GeV","pe");



new_lg_hlt_cent_b->Draw("same");



//------10-30% Centrality
//------Barrel 

//------Creating Canvas & Legend
TCanvas *new_hlt_eff_centrality_1030_b = new TCanvas("new_hlt_eff_centrality_1030_b","new_hlt_eff_centrality_1030_b",600,500);
TLegend *new_lg_hlt_cent_1030_b = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *HLT_eff_1030_15GeV_b = new TGraphAsymmErrors();
HLT_eff_1030_15GeV_b->SetName("HLT_eff_1030_15GeV_b");                                                                                                                           
HLT_eff_1030_15GeV_b->GetXaxis()->SetTitle("p_{T} [GeV]");
HLT_eff_1030_15GeV_b->SetMarkerStyle(24);
HLT_eff_1030_15GeV_b->SetTitle("HLT Efficiency");
HLT_eff_1030_15GeV_b->SetLineColor(kRed);
HLT_eff_1030_15GeV_b->SetMarkerColor(kRed);
HLT_eff_1030_15GeV_b->BayesDivide(photon_w_HLT_15GeV_1030_b,all_photons_1030_b);
HLT_eff_1030_15GeV_b->SetMaximum(1.1);
HLT_eff_1030_15GeV_b->SetMinimum(0);

TGraphAsymmErrors *HLT_eff_1030_25GeV_b = new TGraphAsymmErrors();
HLT_eff_1030_25GeV_b->SetName("HLT_eff_1030_25GeV_b");
HLT_eff_1030_25GeV_b->SetMarkerStyle(25);
HLT_eff_1030_25GeV_b->SetMarkerColor(kBlue);
HLT_eff_1030_25GeV_b->SetLineColor(kBlue);
HLT_eff_1030_25GeV_b->BayesDivide(photon_w_HLT_25GeV_1030_b,all_photons_1030_b);

TGraphAsymmErrors *HLT_eff_1030_35GeV_b = new TGraphAsymmErrors();
HLT_eff_1030_35GeV_b->SetName("HLT_eff_1030_35GeV_b");
HLT_eff_1030_35GeV_b->SetMarkerStyle(26);
HLT_eff_1030_35GeV_b->SetMarkerColor(kMagenta);
HLT_eff_1030_35GeV_b->SetLineColor(kMagenta);
HLT_eff_1030_35GeV_b->BayesDivide(photon_w_HLT_35GeV_1030_b,all_photons_1030_b);


TPaveText *hlt_photon_1030_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_1030_b->SetTextSize(0.03);
hlt_photon_1030_b->SetFillColor(0);
hlt_photon_1030_b->SetBorderSize(0);
hlt_photon_1030_b->SetShadowColor(0);
hlt_photon_1030_b->SetTextAlign(21);
hlt_photon_1030_b->AddText("");
hlt_photon_1030_b->AddText("Pb+Pb");
hlt_photon_1030_b->AddText("10-30%");
hlt_photon_1030_b->AddText("Tight Photons");
hlt_photon_1030_b->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_photon_1030_b->AddText("DeltaR < 0.15");
hlt_photon_1030_b->AddText("HLT Loose Photons");
hlt_photon_1030_b->AddText("|#eta| < 1.37");
hlt_photon_1030_b->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_1030_15GeV_b_hlt = new TF1("fit_1030_15GeV_b_hlt", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_15GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_1030_b->GetXaxis()->GetNbins()+1));
TF1* fit_1030_25GeV_b_hlt = new TF1("fit_1030_25GeV_b_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_25GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_1030_b->GetXaxis()->GetNbins()+1));
TF1* fit_1030_35GeV_b_hlt = new TF1("fit_1030_35GeV_b_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_35GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_1030_b->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_1030_15GeV_b_hlt->SetParameter(0,15.0);		
fit_1030_15GeV_b_hlt->SetParameter(1,2.0);

fit_1030_25GeV_b_hlt->SetParameter(0,25.0);
fit_1030_25GeV_b_hlt->SetParameter(1,2.0);

fit_1030_35GeV_b_hlt->SetParameter(0,35.0);
fit_1030_35GeV_b_hlt->SetParameter(1,2.0);



HLT_eff_1030_15GeV_b->Fit("fit_1030_15GeV_b_hlt","M Q N","",0,photon_w_HLT_15GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_1030_b->GetXaxis()->GetNbins()+1));
HLT_eff_1030_25GeV_b->Fit("fit_1030_25GeV_b_hlt", "M Q N","",0,photon_w_HLT_25GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_1030_b->GetXaxis()->GetNbins()+1));
HLT_eff_1030_35GeV_b->Fit("fit_1030_35GeV_b_hlt","M Q N","",0,photon_w_HLT_35GeV_1030_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_1030_b->GetXaxis()->GetNbins()+1));


fit_1030_15GeV_b_hlt->SetLineStyle(2);
fit_1030_15GeV_b_hlt->SetLineColor(kRed);
fit_1030_15GeV_b_hlt->SetMarkerColor(kRed);

fit_1030_25GeV_b_hlt->SetLineStyle(2);
fit_1030_25GeV_b_hlt->SetLineColor(kBlue);
fit_1030_25GeV_b_hlt->SetMarkerColor(kBlue);

fit_1030_35GeV_b_hlt->SetLineStyle(2);
fit_1030_35GeV_b_hlt->SetLineColor(kMagenta);
fit_1030_35GeV_b_hlt->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_1_b_hlt = new TLine(0.5,1,165,1);
line_100_1_b_hlt->SetLineColor(kBlack);
line_100_1_b_hlt->SetLineWidth(1);




HLT_eff_1030_15GeV_b->Draw("A P");
line_100_1_b_hlt->Draw("same");

HLT_eff_1030_25GeV_b->Draw("P same");
HLT_eff_1030_35GeV_b->Draw("P same");


fit_1030_15GeV_b_hlt->Draw("same");
fit_1030_25GeV_b_hlt->Draw("same");
fit_1030_35GeV_b_hlt->Draw("same");


hlt_photon_1030_b->Draw("same");

new_lg_hlt_cent_1030_b->AddEntry("HLT_eff_1030_15GeV_b","HLT p_{T} > 15 GeV","pe");
new_lg_hlt_cent_1030_b->AddEntry("HLT_eff_1030_25GeV_b","HLT p_{T} > 25 GeV","pe");
new_lg_hlt_cent_1030_b->AddEntry("HLT_eff_1030_35GeV_b","HLT p_{T} > 35 GeV","pe");



new_lg_hlt_cent_1030_b->Draw("same");

//------30-50% Centrality
//------Barrel & Endcaps

//------Creating Canvas & Legend
TCanvas *new_hlt_eff_centrality_3050_b = new TCanvas("new_hlt_eff_centrality_3050_b","new_hlt_eff_centrality_3050_b",600,500);
TLegend *new_lg_hlt_cent_3050_b = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *HLT_eff_3050_15GeV_b = new TGraphAsymmErrors();
HLT_eff_3050_15GeV_b->SetName("HLT_eff_3050_15GeV_b");                                                                                                                           
HLT_eff_3050_15GeV_b->GetXaxis()->SetTitle("p_{T} [GeV]");
HLT_eff_3050_15GeV_b->SetMarkerStyle(24);
HLT_eff_3050_15GeV_b->SetTitle("HLT Efficiency");
HLT_eff_3050_15GeV_b->SetLineColor(kRed);
HLT_eff_3050_15GeV_b->SetMarkerColor(kRed);
HLT_eff_3050_15GeV_b->BayesDivide(photon_w_HLT_15GeV_3050_b,all_photons_3050_b);
HLT_eff_3050_15GeV_b->SetMaximum(1.1);
HLT_eff_3050_15GeV_b->SetMinimum(0);

TGraphAsymmErrors *HLT_eff_3050_25GeV_b = new TGraphAsymmErrors();
HLT_eff_3050_25GeV_b->SetName("HLT_eff_3050_25GeV_b");
HLT_eff_3050_25GeV_b->SetMarkerStyle(25);
HLT_eff_3050_25GeV_b->SetMarkerColor(kBlue);
HLT_eff_3050_25GeV_b->SetLineColor(kBlue);
HLT_eff_3050_25GeV_b->BayesDivide(photon_w_HLT_25GeV_3050_b,all_photons_3050_b);

TGraphAsymmErrors *HLT_eff_3050_35GeV_b = new TGraphAsymmErrors();
HLT_eff_3050_35GeV_b->SetName("HLT_eff_3050_35GeV_b");
HLT_eff_3050_35GeV_b->SetMarkerStyle(26);
HLT_eff_3050_35GeV_b->SetMarkerColor(kMagenta);
HLT_eff_3050_35GeV_b->SetLineColor(kMagenta);
HLT_eff_3050_35GeV_b->BayesDivide(photon_w_HLT_35GeV_3050_b,all_photons_3050_b);


TPaveText *hlt_photon_3050_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_3050_b->SetTextSize(0.03);
hlt_photon_3050_b->SetFillColor(0);
hlt_photon_3050_b->SetBorderSize(0);
hlt_photon_3050_b->SetShadowColor(0);
hlt_photon_3050_b->SetTextAlign(21);
hlt_photon_3050_b->AddText("");
hlt_photon_3050_b->AddText("Pb+Pb");
hlt_photon_3050_b->AddText("30-50%");
hlt_photon_3050_b->AddText("Tight Photons");
hlt_photon_3050_b->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_photon_3050_b->AddText("DeltaR < 0.15");
hlt_photon_3050_b->AddText("HLT Loose Photons");
hlt_photon_3050_b->AddText("|#eta| < 1.37");
hlt_photon_3050_b->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_3050_15GeV_b_hlt = new TF1("fit_3050_15GeV_b_hlt", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_15GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_15GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_3050_b->GetXaxis()->GetNbins()+1));
TF1* fit_3050_25GeV_b_hlt = new TF1("fit_3050_25GeV_b_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_25GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_25GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_3050_b->GetXaxis()->GetNbins()+1));
TF1* fit_3050_35GeV_b_hlt = new TF1("fit_3050_35GeV_b_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_35GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_35GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_3050_b->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_3050_15GeV_b_hlt->SetParameter(0,15.0);		
fit_3050_15GeV_b_hlt->SetParameter(1,2.0);

fit_3050_25GeV_b_hlt->SetParameter(0,25.0);
fit_3050_25GeV_b_hlt->SetParameter(1,2.0);

fit_3050_35GeV_b_hlt->SetParameter(0,35.0);
fit_3050_35GeV_b_hlt->SetParameter(1,2.0);



HLT_eff_3050_15GeV_b->Fit("fit_3050_15GeV_b_hlt","M Q N","",photon_w_HLT_15GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_15GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_3050_b->GetXaxis()->GetNbins()+1));
HLT_eff_3050_25GeV_b->Fit("fit_3050_25GeV_b_hlt", "M Q N","",photon_w_HLT_25GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_25GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_3050_b->GetXaxis()->GetNbins()+1));
HLT_eff_3050_35GeV_b->Fit("fit_3050_35GeV_b_hlt","M Q N","",photon_w_HLT_35GeV_3050_b->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_35GeV_3050_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_3050_b->GetXaxis()->GetNbins()+1));


fit_3050_15GeV_b_hlt->SetLineStyle(2);
fit_3050_15GeV_b_hlt->SetLineColor(kRed);
fit_3050_15GeV_b_hlt->SetMarkerColor(kRed);

fit_3050_25GeV_b_hlt->SetLineStyle(2);
fit_3050_25GeV_b_hlt->SetLineColor(kBlue);
fit_3050_25GeV_b_hlt->SetMarkerColor(kBlue);

fit_3050_35GeV_b_hlt->SetLineStyle(2);
fit_3050_35GeV_b_hlt->SetLineColor(kMagenta);
fit_3050_35GeV_b_hlt->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_2_b_hlt = new TLine(2.0,1,153,1);
line_100_2_b_hlt->SetLineColor(kBlack);
line_100_2_b_hlt->SetLineWidth(1);




HLT_eff_3050_15GeV_b->Draw("A P");
line_100_2_b_hlt->Draw("same");

HLT_eff_3050_25GeV_b->Draw("P same");
HLT_eff_3050_35GeV_b->Draw("P same");


fit_3050_15GeV_b_hlt->Draw("same");
fit_3050_25GeV_b_hlt->Draw("same");
fit_3050_35GeV_b_hlt->Draw("same");


hlt_photon_3050_b->Draw("same");

new_lg_hlt_cent_3050_b->AddEntry("HLT_eff_3050_15GeV_b","HLT p_{T} > 15 GeV","pe");
new_lg_hlt_cent_3050_b->AddEntry("HLT_eff_3050_25GeV_b","HLT p_{T} > 25 GeV","pe");
new_lg_hlt_cent_3050_b->AddEntry("HLT_eff_3050_35GeV_b","HLT p_{T} > 35 GeV","pe");



new_lg_hlt_cent_3050_b->Draw("same");


//------50-80% Centrality
//------Barrel 

//------Creating Canvas & Legend
TCanvas *new_hlt_eff_centrality_5080_b = new TCanvas("new_hlt_eff_centrality_5080_b","new_hlt_eff_centrality_5080_b",600,500);
TLegend *new_lg_hlt_cent_5080_b = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *HLT_eff_5080_15GeV_b = new TGraphAsymmErrors();
HLT_eff_5080_15GeV_b->SetName("HLT_eff_5080_15GeV_b");                                                                                                                           
HLT_eff_5080_15GeV_b->GetXaxis()->SetTitle("p_{T} [GeV]");
HLT_eff_5080_15GeV_b->SetMarkerStyle(24);
HLT_eff_5080_15GeV_b->SetTitle("HLT Efficiency");
HLT_eff_5080_15GeV_b->SetLineColor(kRed);
HLT_eff_5080_15GeV_b->SetMarkerColor(kRed);
HLT_eff_5080_15GeV_b->BayesDivide(photon_w_HLT_15GeV_5080_b,all_photons_5080_b);
HLT_eff_5080_15GeV_b->SetMaximum(1.1);
HLT_eff_5080_15GeV_b->SetMinimum(0);

TGraphAsymmErrors *HLT_eff_5080_25GeV_b = new TGraphAsymmErrors();
HLT_eff_5080_25GeV_b->SetName("HLT_eff_5080_25GeV_b");
HLT_eff_5080_25GeV_b->SetMarkerStyle(25);
HLT_eff_5080_25GeV_b->SetMarkerColor(kBlue);
HLT_eff_5080_25GeV_b->SetLineColor(kBlue);
HLT_eff_5080_25GeV_b->BayesDivide(photon_w_HLT_25GeV_5080_b,all_photons_5080_b);

TGraphAsymmErrors *HLT_eff_5080_35GeV_b = new TGraphAsymmErrors();
HLT_eff_5080_35GeV_b->SetName("HLT_eff_5080_35GeV_b");
HLT_eff_5080_35GeV_b->SetMarkerStyle(26);
HLT_eff_5080_35GeV_b->SetMarkerColor(kMagenta);
HLT_eff_5080_35GeV_b->SetLineColor(kMagenta);
HLT_eff_5080_35GeV_b->BayesDivide(photon_w_HLT_35GeV_5080_b,all_photons_5080_b);


TPaveText *hlt_photon_5080_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_5080_b->SetTextSize(0.03);
hlt_photon_5080_b->SetFillColor(0);
hlt_photon_5080_b->SetBorderSize(0);
hlt_photon_5080_b->SetShadowColor(0);
hlt_photon_5080_b->SetTextAlign(21);
hlt_photon_5080_b->AddText("");
hlt_photon_5080_b->AddText("Pb+Pb");
hlt_photon_5080_b->AddText("50-80%");
hlt_photon_5080_b->AddText("Tight Photons");
hlt_photon_5080_b->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_photon_5080_b->AddText("DeltaR < 0.15");
hlt_photon_5080_b->AddText("HLT Loose Photons");
hlt_photon_5080_b->AddText("|#eta| < 1.37");
hlt_photon_5080_b->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_5080_15GeV_b_hlt = new TF1("fit_5080_15GeV_b_hlt", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_15GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_15GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_5080_b->GetXaxis()->GetNbins()+1));
TF1* fit_5080_25GeV_b_hlt = new TF1("fit_5080_25GeV_b_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_25GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_25GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_5080_b->GetXaxis()->GetNbins()+1));
TF1* fit_5080_35GeV_b_hlt = new TF1("fit_5080_35GeV_b_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_35GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_35GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_5080_b->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_5080_15GeV_b_hlt->SetParameter(0,15.0);		
fit_5080_15GeV_b_hlt->SetParameter(1,2.0);

fit_5080_25GeV_b_hlt->SetParameter(0,25.0);
fit_5080_25GeV_b_hlt->SetParameter(1,2.0);

fit_5080_35GeV_b_hlt->SetParameter(0,35.0);
fit_5080_35GeV_b_hlt->SetParameter(1,2.0);



HLT_eff_5080_15GeV_b->Fit("fit_5080_15GeV_b_hlt","M Q N","",photon_w_HLT_15GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_15GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_5080_b->GetXaxis()->GetNbins()+1));
HLT_eff_5080_25GeV_b->Fit("fit_5080_25GeV_b_hlt", "M Q N","",photon_w_HLT_25GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_25GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_5080_b->GetXaxis()->GetNbins()+1));
HLT_eff_5080_35GeV_b->Fit("fit_5080_35GeV_b_hlt","M Q N","",photon_w_HLT_35GeV_5080_b->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_35GeV_5080_b->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_5080_b->GetXaxis()->GetNbins()+1));


fit_5080_15GeV_b_hlt->SetLineStyle(2);
fit_5080_15GeV_b_hlt->SetLineColor(kRed);
fit_5080_15GeV_b_hlt->SetMarkerColor(kRed);

fit_5080_25GeV_b_hlt->SetLineStyle(2);
fit_5080_25GeV_b_hlt->SetLineColor(kBlue);
fit_5080_25GeV_b_hlt->SetMarkerColor(kBlue);

fit_5080_35GeV_b_hlt->SetLineStyle(2);
fit_5080_35GeV_b_hlt->SetLineColor(kMagenta);
fit_5080_35GeV_b_hlt->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_3_b_hlt = new TLine(2.0,1,200,1);
line_100_3_b_hlt->SetLineColor(kBlack);
line_100_3_b_hlt->SetLineWidth(1);




HLT_eff_5080_15GeV_b->Draw("A P");
line_100_3_b_hlt->Draw("same");

HLT_eff_5080_25GeV_b->Draw("P same");
HLT_eff_5080_35GeV_b->Draw("P same");


fit_5080_15GeV_b_hlt->Draw("same");
fit_5080_25GeV_b_hlt->Draw("same");
fit_5080_35GeV_b_hlt->Draw("same");


hlt_photon_5080_b->Draw("same");

new_lg_hlt_cent_5080_b->AddEntry("HLT_eff_5080_15GeV_b","HLT p_{T} > 15 GeV","pe");
new_lg_hlt_cent_5080_b->AddEntry("HLT_eff_5080_25GeV_b","HLT p_{T} > 25 GeV","pe");
new_lg_hlt_cent_5080_b->AddEntry("HLT_eff_5080_35GeV_b","HLT p_{T} > 35 GeV","pe");



new_lg_hlt_cent_5080_b->Draw("same");

//-------Summary Plots
//-------Barrel 
//-------This is for the parameter 0
TCanvas *summary_new_b_hlt = new TCanvas("summary_new_b_hlt","summary_new_b_hlt",600,500);
TLegend *new_legend_b_hlt = new TLegend(0.6,0.5,0.85,0.6);

//-----0-10%
TH1D *summary_010_0_b_hlt = new TH1D("summary_010_0_b_hlt","",3,10,40);
summary_010_0_b_hlt->SetMarkerStyle(22);
summary_010_0_b_hlt->SetMarkerColor(kBlack);
summary_010_0_b_hlt->SetLineColor(kBlack);
summary_010_0_b_hlt->SetName("summary_010_0_b_hlt");
summary_010_0_b_hlt->GetXaxis()->SetTitle("HLT p_{T} [GeV]");
summary_010_0_b_hlt->GetYaxis()->SetTitle("Value at 50 Percent Efficiency [GeV]");

summary_010_0_b_hlt->SetBinContent(1,fit_010_15GeV_b_hlt->GetParameter(0));
summary_010_0_b_hlt->SetBinError(1,fit_010_15GeV_b_hlt->GetParError(0));
summary_010_0_b_hlt->SetBinContent(2,fit_010_25GeV_b_hlt->GetParameter(0));
summary_010_0_b_hlt->SetBinError(2,fit_010_25GeV_b_hlt->GetParError(0));
summary_010_0_b_hlt->SetBinContent(3,fit_010_35GeV_b_hlt->GetParameter(0));
summary_010_0_b_hlt->SetBinError(3,fit_010_35GeV_b_hlt->GetParError(0));

//------10-30%
TH1D *summary_1030_0_b_hlt = new TH1D("summary_1030_0_b_hlt","",3,10,40);
summary_1030_0_b_hlt->SetMarkerStyle(21);
summary_1030_0_b_hlt->SetMarkerColor(kRed);
summary_1030_0_b_hlt->SetLineColor(kRed);
summary_1030_0_b_hlt->SetName("summary_1030_0_b_hlt");

summary_1030_0_b_hlt->SetBinContent(1,fit_1030_15GeV_b_hlt->GetParameter(0));
summary_1030_0_b_hlt->SetBinError(1,fit_1030_15GeV_b_hlt->GetParError(0));
summary_1030_0_b_hlt->SetBinContent(2,fit_1030_25GeV_b_hlt->GetParameter(0));
summary_1030_0_b_hlt->SetBinError(2,fit_1030_25GeV_b_hlt->GetParError(0));
summary_1030_0_b_hlt->SetBinContent(3,fit_1030_35GeV_b_hlt->GetParameter(0));
summary_1030_0_b_hlt->SetBinError(3,fit_1030_35GeV_b_hlt->GetParError(0));

//------30-50%
TH1D *summary_3050_0_b_hlt = new TH1D("summary_3050_0_b_hlt","",3,10,40);
summary_3050_0_b_hlt->SetMarkerStyle(33);
summary_3050_0_b_hlt->SetMarkerColor(kMagenta);
summary_3050_0_b_hlt->SetLineColor(kMagenta);
summary_3050_0_b_hlt->SetName("summary_3050_0_b_hlt");

summary_3050_0_b_hlt->SetBinContent(1,fit_3050_15GeV_b_hlt->GetParameter(0));
summary_3050_0_b_hlt->SetBinError(1,fit_3050_15GeV_b_hlt->GetParError(0));
summary_3050_0_b_hlt->SetBinContent(2,fit_3050_25GeV_b_hlt->GetParameter(0));
summary_3050_0_b_hlt->SetBinError(2,fit_3050_25GeV_b_hlt->GetParError(0));
summary_3050_0_b_hlt->SetBinContent(3,fit_3050_35GeV_b_hlt->GetParameter(0));
summary_3050_0_b_hlt->SetBinError(3,fit_3050_35GeV_b_hlt->GetParError(0));

//------50-80%
TH1D *summary_5080_0_b_hlt = new TH1D("summary_5080_0_b_hlt","",3,10,40);
summary_5080_0_b_hlt->SetMarkerStyle(20);
summary_5080_0_b_hlt->SetMarkerColor(kBlue);
summary_5080_0_b_hlt->SetLineColor(kBlue);
summary_5080_0_b_hlt->SetName("summary_5080_0_b_hlt");

summary_5080_0_b_hlt->SetBinContent(1,fit_5080_15GeV_b_hlt->GetParameter(0));
summary_5080_0_b_hlt->SetBinError(1,fit_5080_15GeV_b_hlt->GetParError(0));
summary_5080_0_b_hlt->SetBinContent(2,fit_5080_25GeV_b_hlt->GetParameter(0));
summary_5080_0_b_hlt->SetBinError(2,fit_5080_25GeV_b_hlt->GetParError(0));
summary_5080_0_b_hlt->SetBinContent(3,fit_5080_35GeV_b_hlt->GetParameter(0));
summary_5080_0_b_hlt->SetBinError(3,fit_5080_35GeV_b_hlt->GetParError(0));


TPaveText *hlt_summary_0_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_summary_0_b->SetTextSize(0.03);
hlt_summary_0_b->SetFillColor(0);
hlt_summary_0_b->SetBorderSize(0);
hlt_summary_0_b->SetShadowColor(0);
hlt_summary_0_b->SetTextAlign(21);
hlt_summary_0_b->AddText("");
hlt_summary_0_b->AddText("Pb+Pb");
hlt_summary_0_b->AddText("Tight Photons");
hlt_summary_0_b->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_summary_0_b->AddText("HLT Loose Photons");
hlt_summary_0_b->AddText("DeltaR < 0.15");
hlt_summary_0_b->AddText("HLT Loose Photons");
hlt_summary_0_b->AddText("|#eta| < 1.37");
hlt_summary_0_b->SetTextSize(0.036);




summary_010_0_b_hlt->Draw();
summary_1030_0_b_hlt->Draw("same");
summary_3050_0_b_hlt->Draw("same");
summary_5080_0_b_hlt->Draw("same");

hlt_summary_0_b->Draw("same");


new_legend_b_hlt->AddEntry("summary_010_0_b_hlt","0-10%","pe");
new_legend_b_hlt->AddEntry("summary_1030_0_b_hlt","10-30%","pe");
new_legend_b_hlt->AddEntry("summary_3050_0_b_hlt","30-50%","pe");
new_legend_b_hlt->AddEntry("summary_5080_0_b_hlt","50-80%","pe");
new_legend_b_hlt->Draw("same");

//-------This is for the parameter 1
TCanvas *summary_new_2_b_hlt = new TCanvas("summary_new_2_b_hlt","summary_new_2_b_hlt",600,500);
TLegend *new_legend_2_b_hlt = new TLegend(0.6,0.5,0.85,0.6);

//-----0-10%
TH1D *summary_010_1_b_hlt = new TH1D("summary_010_1_b_hlt","",3,10,40);
summary_010_1_b_hlt->SetMarkerStyle(22);
summary_010_1_b_hlt->SetMarkerColor(kBlack);
summary_010_1_b_hlt->SetLineColor(kBlack);
summary_010_1_b_hlt->SetName("summary_010_1_b_hlt");
summary_010_1_b_hlt->GetXaxis()->SetTitle("HLT p_{T} [GeV]");
summary_010_1_b_hlt->GetYaxis()->SetTitle("Width");

summary_010_1_b_hlt->SetBinContent(1,fit_010_15GeV_b_hlt->GetParameter(1));
summary_010_1_b_hlt->SetBinError(1,fit_010_15GeV_b_hlt->GetParError(1));
summary_010_1_b_hlt->SetBinContent(2,fit_010_25GeV_b_hlt->GetParameter(1));
summary_010_1_b_hlt->SetBinError(2,fit_010_25GeV_b_hlt->GetParError(1));
summary_010_1_b_hlt->SetBinContent(3,fit_010_35GeV_b_hlt->GetParameter(1));
summary_010_1_b_hlt->SetBinError(3,fit_010_35GeV_b_hlt->GetParError(1));

//------10-30%
TH1D *summary_1030_1_b_hlt = new TH1D("summary_1030_1_b_hlt","",3,10,40);
summary_1030_1_b_hlt->SetMarkerStyle(21);
summary_1030_1_b_hlt->SetMarkerColor(kRed);
summary_1030_1_b_hlt->SetLineColor(kRed);
summary_1030_1_b_hlt->SetName("summary_1030_1_b_hlt");

summary_1030_1_b_hlt->SetBinContent(1,fit_1030_15GeV_b_hlt->GetParameter(1));
summary_1030_1_b_hlt->SetBinError(1,fit_1030_15GeV_b_hlt->GetParError(1));
summary_1030_1_b_hlt->SetBinContent(2,fit_1030_25GeV_b_hlt->GetParameter(1));
summary_1030_1_b_hlt->SetBinError(2,fit_1030_25GeV_b_hlt->GetParError(1));
summary_1030_1_b_hlt->SetBinContent(3,fit_1030_35GeV_b_hlt->GetParameter(1));
summary_1030_1_b_hlt->SetBinError(3,fit_1030_35GeV_b_hlt->GetParError(1));

//------30-50%
TH1D *summary_3050_1_b_hlt = new TH1D("summary_3050_1_b_hlt","",3,10,40);
summary_3050_1_b_hlt->SetMarkerStyle(33);
summary_3050_1_b_hlt->SetMarkerColor(kMagenta);
summary_3050_1_b_hlt->SetLineColor(kMagenta);
summary_3050_1_b_hlt->SetName("summary_3050_1_b_hlt");

summary_3050_1_b_hlt->SetBinContent(1,fit_3050_15GeV_b_hlt->GetParameter(1));
summary_3050_1_b_hlt->SetBinError(1,fit_3050_15GeV_b_hlt->GetParError(1));
summary_3050_1_b_hlt->SetBinContent(2,fit_3050_25GeV_b_hlt->GetParameter(1));
summary_3050_1_b_hlt->SetBinError(2,fit_3050_25GeV_b_hlt->GetParError(1));
summary_3050_1_b_hlt->SetBinContent(3,fit_3050_35GeV_b_hlt->GetParameter(1));
summary_3050_1_b_hlt->SetBinError(3,fit_3050_35GeV_b_hlt->GetParError(1));

//------50-80%
TH1D *summary_5080_1_b_hlt = new TH1D("summary_5080_1_b_hlt","",3,10,40);
summary_5080_1_b_hlt->SetMarkerStyle(20);
summary_5080_1_b_hlt->SetMarkerColor(kBlue);
summary_5080_1_b_hlt->SetLineColor(kBlue);
summary_5080_1_b->SetName("summary_5080_1_b");

summary_5080_1_b_hlt->SetBinContent(1,fit_5080_15GeV_b_hlt->GetParameter(1));
summary_5080_1_b_hlt->SetBinError(1,fit_5080_15GeV_b_hlt->GetParError(1));
summary_5080_1_b_hlt->SetBinContent(2,fit_5080_25GeV_b_hlt->GetParameter(1));
summary_5080_1_b_hlt->SetBinError(2,fit_5080_25GeV_b_hlt->GetParError(1));
summary_5080_1_b_hlt->SetBinContent(3,fit_5080_35GeV_b_hlt->GetParameter(1));
summary_5080_1_b_hlt->SetBinError(3,fit_5080_35GeV_b_hlt->GetParError(1));


TPaveText *hlt_summary_1_b = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_summary_1_b->SetTextSize(0.03);
hlt_summary_1_b->SetFillColor(0);
hlt_summary_1_b->SetBorderSize(0);
hlt_summary_1_b->SetShadowColor(0);
hlt_summary_1_b->SetTextAlign(21);
hlt_summary_1_b->AddText("");
hlt_summary_1_b->AddText("Pb+Pb");
hlt_summary_1_b->AddText("Tight Photons");
hlt_summary_1_b->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_summary_1_b->AddText("HLT Loose Photons");
hlt_summary_1_b->AddText("DeltaR < 0.15");
hlt_summary_1_b->AddText("HTL Loose Photons");
hlt_summary_1_b->AddText("|#eta| < 1.37");
hlt_summary_1_b->SetTextSize(0.036);




summary_010_1_b_hlt->Draw();
summary_1030_1_b_hlt->Draw("same");
summary_3050_1_b_hlt->Draw("same");
summary_5080_1_b_hlt->Draw("same");

hlt_summary_1_b->Draw("same");


new_legend_2_b_hlt->AddEntry("summary_010_1_b_hlt","0-10%","pe");
new_legend_2_b_hlt->AddEntry("summary_1030_1_b_hlt","10-30%","pe");
new_legend_2_b_hlt->AddEntry("summary_3050_1_b_hlt","30-50%","pe");
new_legend_2_b_hlt->AddEntry("summary_5080_1_b_hlt","50-80%","pe");
new_legend_2_b_hlt->Draw("same");

//------L1 Efficiencies
//------0-10% Centrality
//------Endcaps

//------Creating Canvas & Legend
TCanvas *new_hlt_eff_centrality_e = new TCanvas("new_hlt_eff_centrality_e","new_hlt_eff_centrality_e",600,500);
TLegend *new_lg_hlt_cent_e = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *HLT_eff_010_15GeV_e = new TGraphAsymmErrors();
HLT_eff_010_15GeV_e->SetName("HLT_eff_010_15GeV_e");                                                                                                                           
HLT_eff_010_15GeV_e->GetXaxis()->SetTitle("p_{T} [GeV]");
HLT_eff_010_15GeV_e->SetMarkerStyle(24);
HLT_eff_010_15GeV_e->SetTitle("HLT Efficiency");
HLT_eff_010_15GeV_e->SetLineColor(kRed);
HLT_eff_010_15GeV_e->SetMarkerColor(kRed);
HLT_eff_010_15GeV_e->BayesDivide(photon_w_HLT_15GeV_010_e,all_photons_010_e);
HLT_eff_010_15GeV_e->SetMaximum(1.1);
HLT_eff_010_15GeV_e->SetMinimum(0);

TGraphAsymmErrors *HLT_eff_010_25GeV_e = new TGraphAsymmErrors();
HLT_eff_010_25GeV_e->SetName("L1EM_eff_010_25GeV_e");
HLT_eff_010_25GeV_e->SetMarkerStyle(25);
HLT_eff_010_25GeV_e->SetMarkerColor(kBlue);
HLT_eff_010_25GeV_e->SetLineColor(kBlue);
HLT_eff_010_25GeV_e->BayesDivide(photon_w_HLT_25GeV_010_e,all_photons_010_e);

TGraphAsymmErrors *HLT_eff_010_35GeV_e = new TGraphAsymmErrors();
HLT_eff_010_35GeV_e->SetName("HLT_eff_010_35GeV_e");
HLT_eff_010_35GeV_e->SetMarkerStyle(26);
HLT_eff_010_35GeV_e->SetMarkerColor(kMagenta);
HLT_eff_010_35GeV_e->SetLineColor(kMagenta);
HLT_eff_010_35GeV_e->BayesDivide(photon_w_HLT_35GeV_010_e,all_photons_010_e);


TPaveText *hlt_photon_0_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_0_e->SetTextSize(0.03);
hlt_photon_0_e->SetFillColor(0);
hlt_photon_0_e->SetBorderSize(0);
hlt_photon_0_e->SetShadowColor(0);
hlt_photon_0_e->SetTextAlign(21);
hlt_photon_0_e->AddText("");
hlt_photon_0_e->AddText("Pb+Pb");
hlt_photon_0_e->AddText("0-10%");
hlt_photon_0_e->AddText("Tight Photons");
hlt_photon_0_e->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_photon_0_e->AddText("DeltaR < 0.15");
hlt_photon_0_e->AddText("HLT Loose Photon");
hlt_photon_0_e->AddText("1.56 < |#eta| < 2.37");
hlt_photon_0_e->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_010_15GeV_e_hlt = new TF1("fit_010_15GeV_e_hlt", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_15GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_010_e->GetXaxis()->GetNbins()+1));
TF1* fit_010_25GeV_e_hlt = new TF1("fit_010_25GeV_e_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_25GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_010_e->GetXaxis()->GetNbins()+1));
TF1* fit_010_35GeV_e_hlt = new TF1("fit_010_35GeV_e_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_35GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_010_e->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_010_15GeV_e_hlt->SetParameter(0,15.0);		
fit_010_15GeV_e_hlt->SetParameter(1,2.0);

fit_010_25GeV_e_hlt->SetParameter(0,25.0);
fit_010_25GeV_e_hlt->SetParameter(1,2.0);

fit_010_35GeV_e_hlt->SetParameter(0,35.0);
fit_010_35GeV_e_hlt->SetParameter(1,2.0);



HLT_eff_010_15GeV_e->Fit("fit_010_15GeV_e_hlt","M Q N","",0,photon_w_HLT_15GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_010_e->GetXaxis()->GetNbins()+1));
HLT_eff_010_25GeV_e->Fit("fit_010_25GeV_e_hlt", "M Q N","",0,photon_w_HLT_25GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_010_e->GetXaxis()->GetNbins()+1));
HLT_eff_010_35GeV_e->Fit("fit_010_35GeV_e_hlt","M Q N","",0,photon_w_HLT_35GeV_010_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_010_e->GetXaxis()->GetNbins()+1));


fit_010_15GeV_e_hlt->SetLineStyle(2);
fit_010_15GeV_e_hlt->SetLineColor(kRed);
fit_010_15GeV_e_hlt->SetMarkerColor(kRed);

fit_010_25GeV_e_hlt->SetLineStyle(2);
fit_010_25GeV_e_hlt->SetLineColor(kBlue);
fit_010_25GeV_e_hlt->SetMarkerColor(kBlue);

fit_010_35GeV_e_hlt->SetLineStyle(2);
fit_010_35GeV_e_hlt->SetLineColor(kMagenta);
fit_010_35GeV_e_hlt->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_e_hlt = new TLine(1,1,146,1);
line_100_e_hlt->SetLineColor(kBlack);
line_100_e_hlt->SetLineWidth(1);




HLT_eff_010_15GeV_e->Draw("A P");
line_100_e_hlt->Draw("same");

HLT_eff_010_25GeV_e->Draw("P same");
HLT_eff_010_35GeV_e->Draw("P same");


fit_010_15GeV_e_hlt->Draw("same");
fit_010_25GeV_e_hlt->Draw("same");
fit_010_35GeV_e_hlt->Draw("same");


hlt_photon_0_e->Draw("same");

new_lg_hlt_cent_e->AddEntry("HLT_eff_010_15GeV_e","L1 RoI p_{T} > 15 GeV","pe");
new_lg_hlt_cent_e->AddEntry("HLT_eff_010_25GeV_e","L1 RoI p_{T} > 25 GeV","pe");
new_lg_hlt_cent_e->AddEntry("HLT_eff_010_35GeV_e","L1 RoI p_{T} > 35 GeV","pe");



new_lg_hlt_cent_e->Draw("same");



//------10-30% Centrality
//------Barrel 

//------Creating Canvas & Legend
TCanvas *new_hlt_eff_centrality_1030_e = new TCanvas("new_hlt_eff_centrality_1030_e","new_hlt_eff_centrality_1030_e",600,500);
TLegend *new_lg_hlt_cent_1030_e = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *HLT_eff_1030_15GeV_e = new TGraphAsymmErrors();
HLT_eff_1030_15GeV_e->SetName("HLT_eff_1030_15GeV_e");                                                                                                                           
HLT_eff_1030_15GeV_e->GetXaxis()->SetTitle("p_{T} [GeV]");
HLT_eff_1030_15GeV_e->SetMarkerStyle(24);
HLT_eff_1030_15GeV_e->SetTitle("HLT Efficiency");
HLT_eff_1030_15GeV_e->SetLineColor(kRed);
HLT_eff_1030_15GeV_e->SetMarkerColor(kRed);
HLT_eff_1030_15GeV_e->BayesDivide(photon_w_HLT_15GeV_1030_e,all_photons_1030_e);
HLT_eff_1030_15GeV_e->SetMaximum(1.1);
HLT_eff_1030_15GeV_e->SetMinimum(0);

TGraphAsymmErrors *HLT_eff_1030_25GeV_e = new TGraphAsymmErrors();
HLT_eff_1030_25GeV_e->SetName("L1EM_eff_1030_25GeV_e");
HLT_eff_1030_25GeV_e->SetMarkerStyle(25);
HLT_eff_1030_25GeV_e->SetMarkerColor(kBlue);
HLT_eff_1030_25GeV_e->SetLineColor(kBlue);
HLT_eff_1030_25GeV_e->BayesDivide(photon_w_HLT_25GeV_1030_e,all_photons_1030_e);

TGraphAsymmErrors *HLT_eff_1030_35GeV_e = new TGraphAsymmErrors();
HLT_eff_1030_35GeV_e->SetName("HLT_eff_1030_35GeV_e");
HLT_eff_1030_35GeV_e->SetMarkerStyle(26);
HLT_eff_1030_35GeV_e->SetMarkerColor(kMagenta);
HLT_eff_1030_35GeV_e->SetLineColor(kMagenta);
HLT_eff_1030_35GeV_e->BayesDivide(photon_w_HLT_35GeV_1030_e,all_photons_1030_e);


TPaveText *hlt_photon_1030_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_1030_e->SetTextSize(0.03);
hlt_photon_1030_e->SetFillColor(0);
hlt_photon_1030_e->SetBorderSize(0);
hlt_photon_1030_e->SetShadowColor(0);
hlt_photon_1030_e->SetTextAlign(21);
hlt_photon_1030_e->AddText("");
hlt_photon_1030_e->AddText("Pb+Pb");
hlt_photon_1030_e->AddText("10-30%");
hlt_photon_1030_e->AddText("Tight Photons");
hlt_photon_1030_e->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_photon_1030_e->AddText("DeltaR < 0.15");
hlt_photon_1030_e->AddText("HLT Loose Photons");
hlt_photon_1030_e->AddText("1.56 < |#eta| < 2.37");
hlt_photon_1030_e->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_1030_15GeV_e_hlt = new TF1("fit_1030_15GeV_e_hlt", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_15GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_1030_e->GetXaxis()->GetNbins()+1));
TF1* fit_1030_25GeV_e_hlt = new TF1("fit_1030_25GeV_e_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_25GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_1030_e->GetXaxis()->GetNbins()+1));
TF1* fit_1030_35GeV_e_hlt = new TF1("fit_1030_35GeV_e_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",0,photon_w_HLT_35GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_1030_e->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_1030_15GeV_e_hlt->SetParameter(0,15.0);		
fit_1030_15GeV_e_hlt->SetParameter(1,2.0);

fit_1030_25GeV_e_hlt->SetParameter(0,25.0);
fit_1030_25GeV_e_hlt->SetParameter(1,2.0);

fit_1030_35GeV_e_hlt->SetParameter(0,35.0);
fit_1030_35GeV_e_hlt->SetParameter(1,2.0);



HLT_eff_1030_15GeV_e->Fit("fit_1030_15GeV_e_hlt","M Q N","",0,photon_w_HLT_15GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_1030_e->GetXaxis()->GetNbins()+1));
HLT_eff_1030_25GeV_e->Fit("fit_1030_25GeV_e_hlt", "M Q N","",0,photon_w_HLT_25GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_1030_e->GetXaxis()->GetNbins()+1));
HLT_eff_1030_35GeV_e->Fit("fit_1030_35GeV_e_hlt","M Q N","",0,photon_w_HLT_35GeV_1030_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_1030_e->GetXaxis()->GetNbins()+1));


fit_1030_15GeV_e_hlt->SetLineStyle(2);
fit_1030_15GeV_e_hlt->SetLineColor(kRed);
fit_1030_15GeV_e_hlt->SetMarkerColor(kRed);

fit_1030_25GeV_e_hlt->SetLineStyle(2);
fit_1030_25GeV_e_hlt->SetLineColor(kBlue);
fit_1030_25GeV_e_hlt->SetMarkerColor(kBlue);

fit_1030_35GeV_e_hlt->SetLineStyle(2);
fit_1030_35GeV_e_hlt->SetLineColor(kMagenta);
fit_1030_35GeV_e_hlt->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_1_e_hlt = new TLine(1.0,1,147,1);
line_100_1_e_hlt->SetLineColor(kBlack);
line_100_1_e_hlt->SetLineWidth(1);




HLT_eff_1030_15GeV_e->Draw("A P");
line_100_1_e_hlt->Draw("same");

HLT_eff_1030_25GeV_e->Draw("P same");
HLT_eff_1030_35GeV_e->Draw("P same");


fit_1030_15GeV_e_hlt->Draw("same");
fit_1030_25GeV_e_hlt->Draw("same");
fit_1030_35GeV_e_hlt->Draw("same");


hlt_photon_1030_e->Draw("same");

new_lg_hlt_cent_1030_e->AddEntry("HLT_eff_1030_15GeV_e","HLT p_{T} > 15 GeV","pe");
new_lg_hlt_cent_1030_e->AddEntry("HLT_eff_1030_25GeV_e","HLT p_{T} > 25 GeV","pe");
new_lg_hlt_cent_1030_e->AddEntry("HLT_eff_1030_35GeV_e","HLT p_{T} > 35 GeV","pe");



new_lg_hlt_cent_1030_e->Draw("same");

//------30-50% Centrality
//------Barrel & Endcaps

//------Creating Canvas & Legend
TCanvas *new_hlt_eff_centrality_3050_e = new TCanvas("new_hlt_eff_centrality_3050_e","new_hlt_eff_centrality_3050_e",600,500);
TLegend *new_lg_hlt_cent_3050_e = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *HLT_eff_3050_15GeV_e = new TGraphAsymmErrors();
HLT_eff_3050_15GeV_e->SetName("HLT_eff_3050_15GeV_e");                                                                                                                           
HLT_eff_3050_15GeV_e->GetXaxis()->SetTitle("p_{T} [GeV]");
HLT_eff_3050_15GeV_e->SetMarkerStyle(24);
HLT_eff_3050_15GeV_e->SetTitle("HLT Efficiency");
HLT_eff_3050_15GeV_e->SetLineColor(kRed);
HLT_eff_3050_15GeV_e->SetMarkerColor(kRed);
HLT_eff_3050_15GeV_e->BayesDivide(photon_w_HLT_15GeV_3050_e,all_photons_3050_e);
HLT_eff_3050_15GeV_e->SetMaximum(1.1);
HLT_eff_3050_15GeV_e->SetMinimum(0);

TGraphAsymmErrors *HLT_eff_3050_25GeV_e = new TGraphAsymmErrors();
HLT_eff_3050_25GeV_e->SetName("HLT_eff_3050_25GeV_e");
HLT_eff_3050_25GeV_e->SetMarkerStyle(25);
HLT_eff_3050_25GeV_e->SetMarkerColor(kBlue);
HLT_eff_3050_25GeV_e->SetLineColor(kBlue);
HLT_eff_3050_25GeV_e->BayesDivide(photon_w_HLT_25GeV_3050_e,all_photons_3050_e);

TGraphAsymmErrors *HLT_eff_3050_35GeV_e = new TGraphAsymmErrors();
HLT_eff_3050_35GeV_e->SetName("HLT_eff_3050_35GeV_e");
HLT_eff_3050_35GeV_e->SetMarkerStyle(26);
HLT_eff_3050_35GeV_e->SetMarkerColor(kMagenta);
HLT_eff_3050_35GeV_e->SetLineColor(kMagenta);
HLT_eff_3050_35GeV_e->BayesDivide(photon_w_HLT_35GeV_3050_e,all_photons_3050_e);


TPaveText *hlt_photon_3050_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_3050_e->SetTextSize(0.03);
hlt_photon_3050_e->SetFillColor(0);
hlt_photon_3050_e->SetBorderSize(0);
hlt_photon_3050_e->SetShadowColor(0);
hlt_photon_3050_e->SetTextAlign(21);
hlt_photon_3050_e->AddText("");
hlt_photon_3050_e->AddText("Pb+Pb");
hlt_photon_3050_e->AddText("30-50%");
hlt_photon_3050_e->AddText("Tight Photons");
hlt_photon_3050_e->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_photon_3050_e->AddText("DeltaR < 0.15");
hlt_photon_3050_e->AddText("HLT Loose Photons");
hlt_photon_3050_e->AddText("1.56 < |#eta| < 2.37");
hlt_photon_3050_e->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_3050_15GeV_e_hlt = new TF1("fit_3050_15GeV_e_hlt", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_15GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_15GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_3050_e->GetXaxis()->GetNbins()+1));
TF1* fit_3050_25GeV_e_hlt = new TF1("fit_3050_25GeV_e_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_25GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_25GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_3050_e->GetXaxis()->GetNbins()+1));
TF1* fit_3050_35GeV_e_hlt = new TF1("fit_3050_35GeV_e_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_35GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_35GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_3050_e->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_3050_15GeV_e_hlt->SetParameter(0,15.0);		
fit_3050_15GeV_e_hlt->SetParameter(1,2.0);

fit_3050_25GeV_e_hlt->SetParameter(0,25.0);
fit_3050_25GeV_e_hlt->SetParameter(1,2.0);

fit_3050_35GeV_e_hlt->SetParameter(0,35.0);
fit_3050_35GeV_e_hlt->SetParameter(1,2.0);



HLT_eff_3050_15GeV_e->Fit("fit_3050_15GeV_e_hlt","M Q N","",photon_w_HLT_15GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_15GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_3050_e->GetXaxis()->GetNbins()+1));
HLT_eff_3050_25GeV_e->Fit("fit_3050_25GeV_e_hlt", "M Q N","",photon_w_HLT_25GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_25GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_3050_e->GetXaxis()->GetNbins()+1));
HLT_eff_3050_35GeV_e->Fit("fit_3050_35GeV_e_hlt","M Q N","",photon_w_HLT_35GeV_3050_e->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_35GeV_3050_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_3050_e->GetXaxis()->GetNbins()+1));


fit_3050_15GeV_e_hlt->SetLineStyle(2);
fit_3050_15GeV_e_hlt->SetLineColor(kRed);
fit_3050_15GeV_e_hlt->SetMarkerColor(kRed);

fit_3050_25GeV_e_hlt->SetLineStyle(2);
fit_3050_25GeV_e_hlt->SetLineColor(kBlue);
fit_3050_25GeV_e_hlt->SetMarkerColor(kBlue);

fit_3050_35GeV_e_hlt->SetLineStyle(2);
fit_3050_35GeV_e_hlt->SetLineColor(kMagenta);
fit_3050_35GeV_e_hlt->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_2_e_hlt = new TLine(7.0,1,109,1);
line_100_2_e_hlt->SetLineColor(kBlack);
line_100_2_e_hlt->SetLineWidth(1);




HLT_eff_3050_15GeV_e->Draw("A P");
line_100_2_e_hlt->Draw("same");

HLT_eff_3050_25GeV_e->Draw("P same");
HLT_eff_3050_35GeV_e->Draw("P same");


fit_3050_15GeV_e_hlt->Draw("same");
fit_3050_25GeV_e_hlt->Draw("same");
fit_3050_35GeV_e_hlt->Draw("same");


hlt_photon_3050_e->Draw("same");

new_lg_hlt_cent_3050_e->AddEntry("HLT_eff_3050_15GeV_e","L1 RoI p_{T} > 15 GeV","pe");
new_lg_hlt_cent_3050_e->AddEntry("HLT_eff_3050_25GeV_e","L1 RoI p_{T} > 25 GeV","pe");
new_lg_hlt_cent_3050_e->AddEntry("HLT_eff_3050_35GeV_e","L1 RoI p_{T} > 35 GeV","pe");



new_lg_hlt_cent_3050_e->Draw("same");


//------50-80% Centrality
//------Barrel 

//------Creating Canvas & Legend
TCanvas *new_hlt_eff_centrality_5080_e = new TCanvas("new_hlt_eff_centrality_5080_e","new_hlt_eff_centrality_5080_e",600,500);
TLegend *new_lg_hlt_cent_5080_e = new TLegend(0.6,0.5,0.85,0.6);

//------Creating Turn-On Curves
TGraphAsymmErrors *HLT_eff_5080_15GeV_e = new TGraphAsymmErrors();
HLT_eff_5080_15GeV_e->SetName("HLT_eff_5080_15GeV_e");                                                                                                                           
HLT_eff_5080_15GeV_e->GetXaxis()->SetTitle("p_{T} [GeV]");
HLT_eff_5080_15GeV_e->SetMarkerStyle(24);
HLT_eff_5080_15GeV_e->SetTitle("HLT Efficiency");
HLT_eff_5080_15GeV_e->SetLineColor(kRed);
HLT_eff_5080_15GeV_e->SetMarkerColor(kRed);
HLT_eff_5080_15GeV_e->BayesDivide(photon_w_HLT_15GeV_5080_e,all_photons_5080_e);
HLT_eff_5080_15GeV_e->SetMaximum(1.1);
HLT_eff_5080_15GeV_e->SetMinimum(0);

TGraphAsymmErrors *HLT_eff_5080_25GeV_e = new TGraphAsymmErrors();
HLT_eff_5080_25GeV_e->SetName("HLT_eff_5080_25GeV_e");
HLT_eff_5080_25GeV_e->SetMarkerStyle(25);
HLT_eff_5080_25GeV_e->SetMarkerColor(kBlue);
HLT_eff_5080_25GeV_e->SetLineColor(kBlue);
HLT_eff_5080_25GeV_e->BayesDivide(photon_w_HLT_25GeV_5080_e,all_photons_5080_e);

TGraphAsymmErrors *HLT_eff_5080_35GeV_e = new TGraphAsymmErrors();
HLT_eff_5080_35GeV_e->SetName("HLT_eff_5080_35GeV_e");
HLT_eff_5080_35GeV_e->SetMarkerStyle(26);
HLT_eff_5080_35GeV_e->SetMarkerColor(kMagenta);
HLT_eff_5080_35GeV_e->SetLineColor(kMagenta);
HLT_eff_5080_35GeV_e->BayesDivide(photon_w_HLT_35GeV_5080_e,all_photons_5080_e);


TPaveText *hlt_photon_5080_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_5080_e->SetTextSize(0.03);
hlt_photon_5080_e->SetFillColor(0);
hlt_photon_5080_e->SetBorderSize(0);
hlt_photon_5080_e->SetShadowColor(0);
hlt_photon_5080_e->SetTextAlign(21);
hlt_photon_5080_e->AddText("");
hlt_photon_5080_e->AddText("Pb+Pb");
hlt_photon_5080_e->AddText("50-80%");
hlt_photon_5080_e->AddText("Tight Photons");
hlt_photon_5080_e->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_photon_5080_e->AddText("DeltaR < 0.15");
hlt_photon_5080_e->AddText("HLT Loose Photons");
hlt_photon_5080_e->AddText("1.56 < |#eta| < 2.37");
hlt_photon_5080_e->SetTextSize(0.036);


//----Fit Functions
//----Declaring a error function fitter
//[0] will control posotion of the center of the curve
//[1] will control the width
TF1* fit_5080_15GeV_e_hlt = new TF1("fit_5080_15GeV_e_hlt", "0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_15GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_15GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_5080_e->GetXaxis()->GetNbins()+1));
TF1* fit_5080_25GeV_e_hlt = new TF1("fit_5080_25GeV_e_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_25GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_25GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_5080_e->GetXaxis()->GetNbins()+1));
TF1* fit_5080_35GeV_e_hlt = new TF1("fit_5080_35GeV_e_hlt","0.5*TMath::Erf((x - [0])/[1]) + 0.5",photon_w_HLT_35GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_35GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_5080_e->GetXaxis()->GetNbins()+1));

//----Guessing the fit parameters implemeted into our fit function
fit_5080_15GeV_e_hlt->SetParameter(0,15.0);		
fit_5080_15GeV_e_hlt->SetParameter(1,2.0);

fit_5080_25GeV_e_hlt->SetParameter(0,25.0);
fit_5080_25GeV_e_hlt->SetParameter(1,2.0);

fit_5080_35GeV_e_hlt->SetParameter(0,35.0);
fit_5080_35GeV_e_hlt->SetParameter(1,2.0);



HLT_eff_5080_15GeV_e->Fit("fit_5080_15GeV_e_hlt","M Q N","",photon_w_HLT_15GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_15GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_15GeV_5080_e->GetXaxis()->GetNbins()+1));
HLT_eff_5080_25GeV_e->Fit("fit_5080_25GeV_e_hlt", "M Q N","",photon_w_HLT_25GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_25GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_25GeV_5080_e->GetXaxis()->GetNbins()+1));
HLT_eff_5080_35GeV_e->Fit("fit_5080_35GeV_e_hlt","M Q N","",photon_w_HLT_35GeV_5080_e->GetXaxis()->GetBinLowEdge(1),photon_w_HLT_35GeV_5080_e->GetXaxis()->GetBinLowEdge(photon_w_HLT_35GeV_5080_e->GetXaxis()->GetNbins()+1));


fit_5080_15GeV_e_hlt->SetLineStyle(2);
fit_5080_15GeV_e_hlt->SetLineColor(kRed);
fit_5080_15GeV_e_hlt->SetMarkerColor(kRed);

fit_5080_25GeV_e_hlt->SetLineStyle(2);
fit_5080_25GeV_e_hlt->SetLineColor(kBlue);
fit_5080_25GeV_e_hlt->SetMarkerColor(kBlue);

fit_5080_35GeV_e_hlt->SetLineStyle(2);
fit_5080_35GeV_e_hlt->SetLineColor(kMagenta);
fit_5080_35GeV_e_hlt->SetMarkerColor(kMagenta);

//-----y=1 line
TLine *line_100_3_e_hlt = new TLine(9.0,1,87,1);
line_100_3_e_hlt->SetLineColor(kBlack);
line_100_3_e_hlt->SetLineWidth(1);




HLT_eff_5080_15GeV_e->Draw("A P");
line_100_3_e_hlt->Draw("same");

HLT_eff_5080_25GeV_e->Draw("P same");
HLT_eff_5080_35GeV_e->Draw("P same");


fit_5080_15GeV_e_hlt->Draw("same");
fit_5080_25GeV_e_hlt->Draw("same");
fit_5080_35GeV_e_hlt->Draw("same");


hlt_photon_5080_e->Draw("same");

new_lg_hlt_cent_5080_e->AddEntry("HLT_eff_5080_15GeV_e","HLT p_{T} > 15 GeV","pe");
new_lg_hlt_cent_5080_e->AddEntry("HLT_eff_5080_25GeV_e","HLT p_{T} > 25 GeV","pe");
new_lg_hlt_cent_5080_e->AddEntry("HLT_eff_5080_35GeV_e","HLT p_{T} > 35 GeV","pe");



new_lg_hlt_cent_5080_e->Draw("same");

//-------Summary Plots
//-------Barrel 
//-------This is for the parameter 0
TCanvas *summary_new_e_hlt = new TCanvas("summary_new_e_hlt","summary_new_e_hlt",600,500);
TLegend *new_legend_e_hlt = new TLegend(0.6,0.5,0.85,0.6);

//-----0-10%
TH1D *summary_010_0_e_hlt = new TH1D("summary_010_0_e_hlt","",3,10,40);
summary_010_0_e_hlt->SetMarkerStyle(22);
summary_010_0_e_hlt->SetMarkerColor(kBlack);
summary_010_0_e_hlt->SetLineColor(kBlack);
summary_010_0_e_hlt->SetName("summary_010_0_e_hlt");
summary_010_0_e_hlt->GetXaxis()->SetTitle("HLT p_{T} [GeV]");
summary_010_0_e_hlt->GetYaxis()->SetTitle("Value at 50 Percent Efficiency [GeV]");

summary_010_0_e_hlt->SetBinContent(1,fit_010_15GeV_e_hlt->GetParameter(0));
summary_010_0_e_hlt->SetBinError(1,fit_010_15GeV_e_hlt->GetParError(0));
summary_010_0_e_hlt->SetBinContent(2,fit_010_25GeV_e_hlt->GetParameter(0));
summary_010_0_e_hlt->SetBinError(2,fit_010_25GeV_e_hlt->GetParError(0));
summary_010_0_e_hlt->SetBinContent(3,fit_010_35GeV_e_hlt->GetParameter(0));
summary_010_0_e_hlt->SetBinError(3,fit_010_35GeV_e_hlt->GetParError(0));

//------10-30%
TH1D *summary_1030_0_e_hlt = new TH1D("summary_1030_0_e_hlt","",3,10,40);
summary_1030_0_e_hlt->SetMarkerStyle(21);
summary_1030_0_e_hlt->SetMarkerColor(kRed);
summary_1030_0_e_hlt->SetLineColor(kRed);
summary_1030_0_e_hlt->SetName("summary_1030_0_e_hlt");

summary_1030_0_e_hlt->SetBinContent(1,fit_1030_15GeV_e_hlt->GetParameter(0));
summary_1030_0_e_hlt->SetBinError(1,fit_1030_15GeV_e_hlt->GetParError(0));
summary_1030_0_e_hlt->SetBinContent(2,fit_1030_25GeV_e_hlt->GetParameter(0));
summary_1030_0_e_hlt->SetBinError(2,fit_1030_25GeV_e_hlt->GetParError(0));
summary_1030_0_e_hlt->SetBinContent(3,fit_1030_35GeV_e_hlt->GetParameter(0));
summary_1030_0_e_hlt->SetBinError(3,fit_1030_35GeV_e_hlt->GetParError(0));

//------30-50%
TH1D *summary_3050_0_e_hlt = new TH1D("summary_3050_0_e_hlt","",3,10,40);
summary_3050_0_e_hlt->SetMarkerStyle(33);
summary_3050_0_e_hlt->SetMarkerColor(kMagenta);
summary_3050_0_e_hlt->SetLineColor(kMagenta);
summary_3050_0_e_hlt->SetName("summary_3050_0_e_hlt");

summary_3050_0_e_hlt->SetBinContent(1,fit_3050_15GeV_e_hlt->GetParameter(0));
summary_3050_0_e_hlt->SetBinError(1,fit_3050_15GeV_e_hlt->GetParError(0));
summary_3050_0_e_hlt->SetBinContent(2,fit_3050_25GeV_e_hlt->GetParameter(0));
summary_3050_0_e_hlt->SetBinError(2,fit_3050_25GeV_e_hlt->GetParError(0));
summary_3050_0_e_hlt->SetBinContent(3,fit_3050_35GeV_e_hlt->GetParameter(0));
summary_3050_0_e_hlt->SetBinError(3,fit_3050_35GeV_e_hlt->GetParError(0));

//------50-80%
TH1D *summary_5080_0_e_hlt = new TH1D("summary_5080_0_e_hlt","",3,10,40);
summary_5080_0_e_hlt->SetMarkerStyle(20);
summary_5080_0_e_hlt->SetMarkerColor(kBlue);
summary_5080_0_e_hlt->SetLineColor(kBlue);
summary_5080_0_e_hlt->SetName("summary_5080_0_e_hlt");

summary_5080_0_e_hlt->SetBinContent(1,fit_5080_15GeV_e_hlt->GetParameter(0));
summary_5080_0_e_hlt->SetBinError(1,fit_5080_15GeV_e_hlt->GetParError(0));
summary_5080_0_e_hlt->SetBinContent(2,fit_5080_25GeV_e_hlt->GetParameter(0));
summary_5080_0_e_hlt->SetBinError(2,fit_5080_25GeV_e_hlt->GetParError(0));
summary_5080_0_e_hlt->SetBinContent(3,fit_5080_35GeV_e_hlt->GetParameter(0));
summary_5080_0_e_hlt->SetBinError(3,fit_5080_35GeV_e_hlt->GetParError(0));


TPaveText *hlt_summary_0_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_summary_0_e->SetTextSize(0.03);
hlt_summary_0_e->SetFillColor(0);
hlt_summary_0_e->SetBorderSize(0);
hlt_summary_0_e->SetShadowColor(0);
hlt_summary_0_e->SetTextAlign(21);
hlt_summary_0_e->AddText("");
hlt_summary_0_e->AddText("Pb+Pb");
hlt_summary_0_e->AddText("Tight Photons");
hlt_summary_0_e->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_summary_0_e->AddText("HLT Loose Photons");
hlt_summary_0_e->AddText("DeltaR < 0.15");
hlt_summary_0_e->AddText("1.56 < |#eta| < 2.37");
hlt_summary_0_e->SetTextSize(0.036);




summary_010_0_e_hlt->Draw();
summary_1030_0_e_hlt->Draw("same");
summary_3050_0_e_hlt->Draw("same");
summary_5080_0_e_hlt->Draw("same");

hlt_summary_0_e->Draw("same");


new_legend_e_hlt->AddEntry("summary_010_0_e_hlt","0-10%","pe");
new_legend_e_hlt->AddEntry("summary_1030_0_e_hlt","10-30%","pe");
new_legend_e_hlt->AddEntry("summary_3050_0_e_hlt","30-50%","pe");
new_legend_e_hlt->AddEntry("summary_5080_0_e_hlt","50-80%","pe");
new_legend_e_hlt->Draw("same");

//-------This is for the parameter 1
TCanvas *summary_new_2_e_hlt = new TCanvas("summary_new_2_e_hlt","summary_new_2_e_hlt",600,500);
TLegend *new_legend_2_e_hlt = new TLegend(0.6,0.5,0.85,0.6);

//-----0-10%
TH1D *summary_010_1_e_hlt = new TH1D("summary_010_1_e_hlt","",3,10,40);
summary_010_1_e_hlt->SetMarkerStyle(22);
summary_010_1_e_hlt->SetMarkerColor(kBlack);
summary_010_1_e_hlt->SetLineColor(kBlack);
summary_010_1_e_hlt->SetName("summary_010_1_e_hlt");
summary_010_1_e_hlt->GetXaxis()->SetTitle("HLT p_{T} [GeV]");
summary_010_1_e_hlt->GetYaxis()->SetTitle("Width");

summary_010_1_e_hlt->SetBinContent(1,fit_010_15GeV_e_hlt->GetParameter(1));
summary_010_1_e_hlt->SetBinError(1,fit_010_15GeV_e_hlt->GetParError(1));
summary_010_1_e_hlt->SetBinContent(2,fit_010_25GeV_e_hlt->GetParameter(1));
summary_010_1_e_hlt->SetBinError(2,fit_010_25GeV_e_hlt->GetParError(1));
summary_010_1_e_hlt->SetBinContent(3,fit_010_35GeV_e_hlt->GetParameter(1));
summary_010_1_e_hlt->SetBinError(3,fit_010_35GeV_e_hlt->GetParError(1));

//------10-30%
TH1D *summary_1030_1_e_hlt = new TH1D("summary_1030_1_e_hlt","",3,10,40);
summary_1030_1_e_hlt->SetMarkerStyle(21);
summary_1030_1_e_hlt->SetMarkerColor(kRed);
summary_1030_1_e_hlt->SetLineColor(kRed);
summary_1030_1_e_hlt->SetName("summary_1030_1_e_hlt");

summary_1030_1_e_hlt->SetBinContent(1,fit_1030_15GeV_e_hlt->GetParameter(1));
summary_1030_1_e_hlt->SetBinError(1,fit_1030_15GeV_e_hlt->GetParError(1));
summary_1030_1_e_hlt->SetBinContent(2,fit_1030_25GeV_e_hlt->GetParameter(1));
summary_1030_1_e_hlt->SetBinError(2,fit_1030_25GeV_e_hlt->GetParError(1));
summary_1030_1_e_hlt->SetBinContent(3,fit_1030_35GeV_e_hlt->GetParameter(1));
summary_1030_1_e_hlt->SetBinError(3,fit_1030_35GeV_e_hlt->GetParError(1));

//------30-50%
TH1D *summary_3050_1_e_hlt = new TH1D("summary_3050_1_e_hlt","",3,10,40);
summary_3050_1_e_hlt->SetMarkerStyle(33);
summary_3050_1_e_hlt->SetMarkerColor(kMagenta);
summary_3050_1_e_hlt->SetLineColor(kMagenta);
summary_3050_1_e_hlt->SetName("summary_3050_1_e_hlt");

summary_3050_1_e_hlt->SetBinContent(1,fit_3050_15GeV_e_hlt->GetParameter(1));
summary_3050_1_e_hlt->SetBinError(1,fit_3050_15GeV_e_hlt->GetParError(1));
summary_3050_1_e_hlt->SetBinContent(2,fit_3050_25GeV_e_hlt->GetParameter(1));
summary_3050_1_e_hlt->SetBinError(2,fit_3050_25GeV_e_hlt->GetParError(1));
summary_3050_1_e_hlt->SetBinContent(3,fit_3050_35GeV_e_hlt->GetParameter(1));
summary_3050_1_e_hlt->SetBinError(3,fit_3050_35GeV_e_hlt->GetParError(1));

//------50-80%
TH1D *summary_5080_1_e_hlt = new TH1D("summary_5080_1_e_hlt","",3,10,40);
summary_5080_1_e_hlt->SetMarkerStyle(20);
summary_5080_1_e_hlt->SetMarkerColor(kBlue);
summary_5080_1_e_hlt->SetLineColor(kBlue);
summary_5080_1_e_hlt->SetName("summary_5080_1_e_hlt");

summary_5080_1_e_hlt->SetBinContent(1,fit_5080_15GeV_e_hlt->GetParameter(1));
summary_5080_1_e_hlt->SetBinError(1,fit_5080_15GeV_e_hlt->GetParError(1));
summary_5080_1_e_hlt->SetBinContent(2,fit_5080_25GeV_e_hlt->GetParameter(1));
summary_5080_1_e_hlt->SetBinError(2,fit_5080_25GeV_e_hlt->GetParError(1));
summary_5080_1_e_hlt->SetBinContent(3,fit_5080_35GeV_e_hlt->GetParameter(1));
summary_5080_1_e_hlt->SetBinError(3,fit_5080_35GeV_e_hlt->GetParError(1));


TPaveText *hlt_summary_1_e = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_summary_1_e->SetTextSize(0.03);
hlt_summary_1_e->SetFillColor(0);
hlt_summary_1_e->SetBorderSize(0);
hlt_summary_1_e->SetShadowColor(0);
hlt_summary_1_e->SetTextAlign(21);
hlt_summary_1_e->AddText("");
hlt_summary_1_e->AddText("Pb+Pb");
hlt_summary_1_e->AddText("Tight Photons");
hlt_summary_1_e->AddText("p^{#gamma}_{T} > 15 GeV");
hlt_summary_1_e->AddText("HLT Loose Photons");
hlt_summary_1_e->AddText("DeltaR < 0.15");
hlt_summary_1_e->AddText("HLT Loose Photons");
hlt_summary_1_e->AddText("1.56 < |#eta| < 2.37");
hlt_summary_1_e->SetTextSize(0.036);




summary_010_1_e_hlt->Draw();
summary_1030_1_e_hlt->Draw("same");
summary_3050_1_e_hlt->Draw("same");
summary_5080_1_e_hlt->Draw("same");

hlt_summary_1_e->Draw("same");


new_legend_2_e_hlt->AddEntry("summary_010_1_e_hlt","0-10%","pe");
new_legend_2_e_hlt->AddEntry("summary_1030_1_e_hlt","10-30%","pe");
new_legend_2_e_hlt->AddEntry("summary_3050_1_e_hlt","30-50%","pe");
new_legend_2_e_hlt->AddEntry("summary_5080_1_e_hlt","50-80%","pe");
new_legend_2_e_hlt->Draw("same");

return 0;

}