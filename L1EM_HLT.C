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

using namespace std;

float pi_dos = 360;
float pi_mi = 180;
float original = 3.14159;

float delta_phi(float phi_uno,float phi_dos){                                                                                                                                   

  float final_phi = fmod(abs(phi_uno - phi_dos), 2*original);
  float phi = final_phi > original  ? 2*original - final_phi : final_phi;
  return phi;

}



int L1EM_HLT(){
  
  //-------Histograms----------//                                                                                                                                                
  TH1D *all_photons_barrel = new TH1D("all_photons_barrel","",40,0,200);               //----pT Distribution w/ L1 Trigger                                                         
  TH1D *all_photons_endcaps = new TH1D("all_photons_endcaps","",40,0,200);
  TH1D *all_photons= new TH1D("all_photons","",40,0,200);


  TH1D *all_photons_eta = new TH1D("all_photons_eta","",20,-4,4);
  TH1D *all_photons_phi = new TH1D("all_photons_phi","",20,-4,4);
  TH1D *all_photons_sumET = new TH1D("all_photons_sumET","",45,-0.5,5);

  TH1D *photons_wL1EM_15GeV = new TH1D("photons_wL1EM_15GeV","",40,0,200);          //----pT Distribution W/ Matched L1 Trigger
  TH1D *photons_wL1EM_15GeV_barrel = new TH1D("photons_wL1EM_15GeV_barrel","",40,0,200);
  TH1D *photons_wL1EM_15GeV_endcaps = new TH1D("photons_wL1EM_15GeV_endcaps","",40,0,200);

  TH1D *photons_wL1EM_25GeV = new TH1D("photons_wL1EM_25GeV","",40,0,200);
  TH1D *photons_wL1EM_25GeV_barrel = new TH1D("photons_wL1EM_25GeV_barrel","",40,0,200);
  TH1D *photons_wL1EM_25GeV_endcaps = new TH1D("photons_wL1EM_25GeV_endcaps","",40,0,200);


  TH1D *photons_wL1EM_35GeV = new TH1D("photons_wL1EM_35GeV","",40,0,200);
  TH1D *photons_wL1EM_35GeV_barrel = new TH1D("photons_wL1EM_35GeV_barrel","",40,0,200);
  TH1D *photons_wL1EM_35GeV_endcaps = new TH1D("photons_wL1EM_15GeV_endcaps","",40,0,200);


  TH1D *photons_wL1EM_eta_15GeV = new TH1D("photons_wL1EM_eta_15GeV","",20,-4,4);  
  TH1D *photons_wL1EM_eta_25GeV = new TH1D("photons_wL1EM_25GeV","",20,-4,4);
  TH1D *photons_wL1EM_eta_35GeV = new TH1D("photons_wL1EM_35GeV","",20,-4,4);


  TH1D *photons_wL1EM_phi_15GeV = new TH1D("photons_wL1EM_phi_15GeV","",20,-4,4);
  TH1D *photons_wL1EM_phi_25GeV = new TH1D("photons_wL1EM_phi_25GeV","",20,-4,4);
  TH1D *photons_wL1EM_phi_35GeV = new TH1D("photons_wL1EM_phi_35GeV","",20,-4,4);

  //----Photons That Fired HLT Triggers > 15 GeV
  TH1D *photons_w_HLT_15_barrel = new TH1D("photons_w_HLT_15_barrel","", 40, 0, 200);
  TH1D *photons_w_HLT_15_endcaps = new TH1D("photons_w_HLT_15_endcaps","", 40,0,200);
  TH1D *photons_w_HLT_15 = new TH1D("photons_w_HLT_15","",40,0,200);

  //----Photons That Fired HLT Triggers > 25 GeV
  TH1D *photons_w_HLT_25_barrel = new TH1D("photons_w_HLT_25_barrel","",40,0,200);
  TH1D *photons_w_HLT_25_endcaps = new TH1D("photons_w_HLT_25_endcaps","",40,0,200);
  TH1D *photons_w_HLT_25 = new TH1D("photons_w_HLT_25","",40,0,200);

  //----Photons THat Fired HLT Triggers > 35 GeV  
  TH1D *photons_w_HLT_35 = new TH1D("photons_w_HLT_35","",40,0,200);
  TH1D *photons_w_HLT_35_barrel = new TH1D("photons_w_HLT_35_barrel","",40,0,200);
  TH1D *photons_w_HLT_35_endcaps = new TH1D("photons_w_HLT_35_endcaps","",40,0,200);



  TH1D *nofire_photons_pt = new TH1D("nofire_photons_pt", "",20,0,500);
  TH1D *nofire_photons_eta = new TH1D("nofire_photons_eta", "",20,-4,4);
  TH1D *nofire_photons_phi = new TH1D("nofire_photons_phi", "",20,-3.5,3.5);

  TH1D *L1EM_num =new TH1D("L1EM_num", "",45,0,250);
  TH1D *HLT_num = new TH1D("HLT_num","",45,0,50);

  TH1D *nofire_sumET = new TH1D("nofire_sumET","",45,-0.5,5);
  TH1D *fire_sumET = new TH1D("fire_sumET","",45,-0.5,5);



  TH1D *hlt_fire_sumET_15 = new TH1D("hlt_fire_sumET_15","",45,-0.5,5);


  TH1D *hlt_fire_sumET_25 = new TH1D("hlt_fire_sumET_25","",45,-0.5,5);

  TH1D *hlt_fire_sumET_35_barrel = new TH1D("hlt_fire_sumET_35_barrel","",45,-0.5,5);
  TH1D *hlt_fire_sumET_35_endcaps = new TH1D("hlt_fire_sumET_35_endcaps","",45,-0.5,5);
  TH1D *hlt_fire_sumET_35 = new TH1D("hlt_fire_sumET_35","",45,-0.5,5);




  TH1D *hlt_nofire_sumET_15 = new TH1D("hlt_nofire_sumET_15","",45,-0.5,5);
  TH1D *hlt_nofire_sumET_25 = new TH1D("hlt_nofire_sumET_25","",45,-0.5,5);
  TH1D *hlt_nofire_sumET_35 = new TH1D("hlt_nofire_sumET_35","",45,-0.5,5);


  
  //-----HLT Plots-------//

  TH1D *hlt_pT_dis = new TH1D("hlt_pT_dis", "", 45, 0,250);
  TH1D *hlt_eta_dis = new TH1D("hlt_eta_dis", "", 45, 0, 3);
  TH1D *hlt_phi_dis = new TH1D("hlt_phi_dis", "", 45, -3.18,3.18);


  //---L1-------//
  TH1D *l1_pT_dis = new TH1D("l1_pT_dis", "", 45, 0,250);


  //-----2D Histograms-----//
  TH2D *pT_vs_L1EM = new TH2D("pT_vs_L1EM","",15,0,500,50,0,250); 

  TH2D *pT_vs_HLT_num = new TH2D("pT_vs_HLT_num","",15,0,500,50,0,45);


  //----Photons That Did Not Fire a HLT Trigger 
  TH1D *photons_nofire_HLT_pt = new TH1D("photons_nofire_HLT_pt","",45,0,500);
  TH1D *photons_nofire_HLT_eta = new TH1D("photons_nofire_HLT_eta","",45,-3,3);
  TH1D *photons_nofire_HLT_sumET = new TH1D("photons_nofire_HLT_sumET","",45,-0.5,5);
  
  //---How Many Events Pass These Triggers
  TH1D *HLT_g15_loose_ion_histo = new TH1D("HLT_g15_loose_ion_histo","",5,0,5);
  TH1D *HLT_g20_loose_ion_histo = new TH1D("HLT_g20_loose_ion_histo","",5,0,5);
  TH1D *HLT_g20_loose_histo = new TH1D("HLT_g20_loose_histo","",5,0,5);
  TH1D *HLT_g30_loose_ion_histo = new TH1D("HLT_g30_loose_ion_histo","",5,0,5);
  TH1D *HLT_g13_etcut_ion_histo = new TH1D("HLT_g13_etcut_ion_histo","",5,0,5);
  TH1D *HLT_g18_etcut_ion_histo = new TH1D("HLT_g18_etcut_ion_histo","",5,0,5);
  TH1D *HLT_g18_etcut_histo = new TH1D("HLT_g18_etcut_histo","", 5,0,5);
  TH1D *HLT_g28_etcut_ion_histo = new TH1D("HLT_g28_etcut_ion_histo","",5,0,5);
  TH1D *HLT_g50_loose_ion_histo = new TH1D("HLT_g50_loose_ion_histo","",5,0,5);
  TH1D *HLT_noalg_L1EM10_histo = new TH1D("HLT_noalg_L1EM10_histo","",5,0,5);
  TH1D *HLT_noalg_L1EM12_histo = new TH1D("HLT_noalg_L1EM12","",5,0,5);

  TH1D *HLT_noalg_eb_L1TE50_histo = new TH1D("HLT_noalg_eb_L1TE50","",5,0,5);
  
  //--------------------------------------------------------//
  //-----L1 & HLT Efficiencies Dependent onn Centrality-----//
  //--------------------------------------------------------//
  //**********************************L1 Histograms**********************************// 
  //----Barrel & EndCaps
  TH1D *all_photons_010 = new TH1D("all_photons_010","", 40,0,200); 
  TH1D *all_photons_1030 = new TH1D("all_photons_1030","",40,0,200);
  TH1D *all_photons_3050 = new TH1D("all_photons_3050","",40,0,200);
  TH1D *all_photons_5080 = new TH1D("all_photons_5080","",40,0,200);



  TH1D *photon_w_L1EM_15GeV_010 = new TH1D("photon_w_L1EM_15GeV_010","",40,0,200);
  TH1D *photon_w_L1EM_15GeV_1030 = new TH1D("photon_w_L1EM_15GeV_1030","",40,0,200);
  TH1D *photon_w_L1EM_15GeV_3050 = new TH1D("photon_w_L1EM_15GeV_3050","",40,0,200);
  TH1D *photon_w_L1EM_15GeV_5080 = new TH1D("photon_w_L1EM_15GeV_5080","",40,0,200);


  TH1D *photon_w_L1EM_25GeV_010 = new TH1D("photon_w_L1EM_25GeV_010","",40,0,200);
  TH1D *photon_w_L1EM_25GeV_1030 = new TH1D("photon_w_L1EM_25GeV_1030","",40,0,200);
  TH1D *photon_w_L1EM_25GeV_3050 = new TH1D("photon_w_L1EM_25GeV_3050","",40,0,200);
  TH1D *photon_w_L1EM_25GeV_5080 = new TH1D("photon_w_L1EM_25GeV_5080","",40,0,200);

  TH1D *photon_w_L1EM_35GeV_010 = new TH1D("photon_w_L1EM_35GeV_010","",40,0,200);
  TH1D *photon_w_L1EM_35GeV_1030 = new TH1D("photon_w_L1EM_35GeV_1030","",40,0,200);
  TH1D *photon_w_L1EM_35GeV_3050 = new TH1D("photon_w_L1EM_35GeV_3050","",40,0,200);
  TH1D *photon_w_L1EM_35GeV_5080 = new TH1D("photon_w_L1EM_35GeV_5080","",40,0,200);
  
  //----Barrel
  TH1D *all_photons_010_b = new TH1D("all_photons_010_b","",40,0,200); 
  TH1D *all_photons_1030_b = new TH1D("all_photons_1030_b","",40,0,200);
  TH1D *all_photons_3050_b = new TH1D("all_photons_3050_b","",40,0,200);
  TH1D *all_photons_5080_b = new TH1D("all_photons_5080_b","",40,0,200);

  
  TH1D *photon_w_L1EM_15GeV_010_b = new TH1D("photon_w_L1EM_15GeV_010_b","",40,0,200);
  TH1D *photon_w_L1EM_15GeV_1030_b = new TH1D("photon_w_L1EM_15GeV_1030_b","",40,0,200);
  TH1D *photon_w_L1EM_15GeV_3050_b = new TH1D("photon_w_L1EM_15GeV_3050_b","",40,0,200);
  TH1D *photon_w_L1EM_15GeV_5080_b = new TH1D("photon_w_L1EM_15GeV_5080_b","",40,0,200);

  TH1D *photon_w_L1EM_25GeV_010_b = new TH1D("photon_w_L1EM_25GeV_010_b","",40,0,200);
  TH1D *photon_w_L1EM_25GeV_1030_b = new TH1D("photon_w_L1EM_25GeV_1030_b","",40,0,200);
  TH1D *photon_w_L1EM_25GeV_3050_b = new TH1D("photon_w_L1EM_25GeV_3050_b","",40,0,200);
  TH1D *photon_w_L1EM_25GeV_5080_b = new TH1D("photon_w_L1EM_25GeV_5080_b","",40,0,200);

  TH1D *photon_w_L1EM_35GeV_010_b = new TH1D("photon_w_L1EM_35GeV_010_b","",40,0,200);
  TH1D *photon_w_L1EM_35GeV_1030_b = new TH1D("photon_w_L1EM_35GeV_1030_b","",40,0,200);
  TH1D *photon_w_L1EM_35GeV_3050_b = new TH1D("photon_w_L1EM_35GeV_3050_b","",40,0,200);
  TH1D *photon_w_L1EM_35GeV_5080_b = new TH1D("photon_w_L1EM_35GeV_5080_b","",40,0,200);
  

  //----Endcaps
  TH1D *all_photons_010_e = new TH1D("all_photons_010_e","",40,0,200); 
  TH1D *all_photons_1030_e = new TH1D("all_photons_1030_e","",40,0,200);
  TH1D *all_photons_3050_e = new TH1D("all_photons_3050_e","",40,0,200);
  TH1D *all_photons_5080_e = new TH1D("all_photons_5080_e","",40,0,200);


  TH1D *photon_w_L1EM_15GeV_010_e = new TH1D("photon_w_L1EM_15GeV_010_e","",40,0,200);
  TH1D *photon_w_L1EM_15GeV_1030_e = new TH1D("photon_w_L1EM_15GeV_1030_e","",40,0,200);
  TH1D *photon_w_L1EM_15GeV_3050_e = new TH1D("photon_w_L1EM_15GeV_3050_e","",40,0,200);
  TH1D *photon_w_L1EM_15GeV_5080_e = new TH1D("photon_w_L1EM_15GeV_5080_e","",40,0,200);

  TH1D *photon_w_L1EM_25GeV_010_e = new TH1D("photon_w_L1EM_25GeV_010_e","",40,0,200);
  TH1D *photon_w_L1EM_25GeV_1030_e = new TH1D("photon_w_L1EM_25GeV_1030_e","",40,0,200);
  TH1D *photon_w_L1EM_25GeV_3050_e = new TH1D("photon_w_L1EM_25GeV_3050_e","",40,0,200);
  TH1D *photon_w_L1EM_25GeV_5080_e = new TH1D("photon_w_L1EM_25GeV_5080_e","",40,0,200);

  TH1D *photon_w_L1EM_35GeV_010_e = new TH1D("photon_w_L1EM_35GeV_010_e","",40,0,200);
  TH1D *photon_w_L1EM_35GeV_1030_e = new TH1D("photon_w_L1EM_35GeV_1030_e","",40,0,200);
  TH1D *photon_w_L1EM_35GeV_3050_e = new TH1D("photon_w_L1EM_35GeV_3050_e","",40,0,200);
  TH1D *photon_w_L1EM_35GeV_5080_e = new TH1D("photon_w_L1EM_35GeV_5080_e","",40,0,200);

//***************************************************************************************//

//----Barrel & EndCaps
  TH1D *photon_w_HLT_15GeV_010 = new TH1D("photon_w_HLT_15GeV_010","",40,0,200);
  TH1D *photon_w_HLT_15GeV_1030 = new TH1D("photon_w_HLT_15GeV_1030","",40,0,200);
  TH1D *photon_w_HLT_15GeV_3050 = new TH1D("photon_w_HLT_15GeV_3050","",40,0,200);
  TH1D *photon_w_HLT_15GeV_5080 = new TH1D("photon_w_HLT_15GeV_5080","",40,0,200);

  TH1D *photon_w_HLT_25GeV_010 = new TH1D("photon_w_HLT_25GeV_010","",40,0,200);
  TH1D *photon_w_HLT_25GeV_1030 = new TH1D("photon_w_HLT_25GeV_1030","",40,0,200);
  TH1D *photon_w_HLT_25GeV_3050 = new TH1D("photon_w_HLT_25GeV_3050","",40,0,200);
  TH1D *photon_w_HLT_25GeV_5080 = new TH1D("photon_w_HLT_25GeV_5080","",40,0,200);

  TH1D *photon_w_HLT_35GeV_010 = new TH1D("photon_w_HLT_35GeV_010","",40,0,200);
  TH1D *photon_w_HLT_35GeV_1030 = new TH1D("photon_w_HLT_35GeV_1030","",40,0,200);
  TH1D *photon_w_HLT_35GeV_3050 = new TH1D("photon_w_HLT_35GeV_3050","",40,0,200);
  TH1D *photon_w_HLT_35GeV_5080 = new TH1D("photon_w_HLT_35GeV_5080","",40,0,200);
  
  //----Barrel
  TH1D *photon_w_HLT_15GeV_010_b = new TH1D("photon_w_HLT_15GeV_010_b","",40,0,200);
  TH1D *photon_w_HLT_15GeV_1030_b = new TH1D("photon_w_HLT_15GeV_1030_b","",40,0,200);
  TH1D *photon_w_HLT_15GeV_3050_b = new TH1D("photon_w_HLT_15GeV_3050_b","",40,0,200);
  TH1D *photon_w_HLT_15GeV_5080_b = new TH1D("photon_w_HLT_15GeV_5080_b","",40,0,200);

  TH1D *photon_w_HLT_25GeV_010_b = new TH1D("photon_w_HLT_25GeV_010_b","",40,0,200);
  TH1D *photon_w_HLT_25GeV_1030_b = new TH1D("photon_w_HLT_25GeV_1030_b","",40,0,200);
  TH1D *photon_w_HLT_25GeV_3050_b = new TH1D("photon_w_HLT_25GeV_3050_b","",40,0,200);
  TH1D *photon_w_HLT_25GeV_5080_b = new TH1D("photon_w_HLT_25GeV_5080_b","",40,0,200);

  TH1D *photon_w_HLT_35GeV_010_b = new TH1D("photon_w_HLT_35GeV_010_b","",40,0,200);
  TH1D *photon_w_HLT_35GeV_1030_b = new TH1D("photon_w_HLT_35GeV_1030_b","",40,0,200);
  TH1D *photon_w_HLT_35GeV_3050_b = new TH1D("photon_w_HLT_35GeV_3050_b","",40,0,200);
  TH1D *photon_w_HLT_35GeV_5080_b = new TH1D("photon_w_HLT_35GeV_5080_b","",40,0,200);
  

  //----Endcaps
  TH1D *photon_w_HLT_15GeV_010_e = new TH1D("photon_w_HLT_15GeV_010_e","",40,0,200);
  TH1D *photon_w_HLT_15GeV_1030_e = new TH1D("photon_w_HLT_15GeV_1030_e","",40,0,200);
  TH1D *photon_w_HLT_15GeV_3050_e = new TH1D("photon_w_HLT_15GeV_3050_e","",40,0,200);
  TH1D *photon_w_HLT_15GeV_5080_e = new TH1D("photon_w_HLT_15GeV_5080_e","",40,0,200);

  TH1D *photon_w_HLT_25GeV_010_e = new TH1D("photon_w_HLT_25GeV_010_e","",40,0,200);
  TH1D *photon_w_HLT_25GeV_1030_e = new TH1D("photon_w_HLT_25GeV_1030_e","",40,0,200);
  TH1D *photon_w_HLT_25GeV_3050_e = new TH1D("photon_w_HLT_25GeV_3050_e","",40,0,200);
  TH1D *photon_w_HLT_25GeV_5080_e = new TH1D("photon_w_HLT_25GeV_5080_e","",40,0,200);

  TH1D *photon_w_HLT_35GeV_010_e = new TH1D("photon_w_HLT_35GeV_010_e","",40,0,200);
  TH1D *photon_w_HLT_35GeV_1030_e = new TH1D("photon_w_HLT_35GeV_1030_e","",40,0,200);
  TH1D *photon_w_HLT_35GeV_3050_e = new TH1D("photon_w_HLT_35GeV_3050_e","",40,0,200);
  TH1D *photon_w_HLT_35GeV_5080_e = new TH1D("photon_w_HLT_35GeV_5080_e","",40,0,200);

  //Events & Run Number                                                                                                                                                          
  int event_num = 0;

  //FCal 
  float fcalC_et_data = 0;
  float fcalA_et_data = 0;

  //HLT 
  bool HLT_noalg_eb_L1TE50 = 0;
  bool HLT_noalg_eb_L1ZDC_A_C_VTE50 = 0;
  bool HLT_noalg_eb_L1MU4 = 0;

  //-------HLT Photon Triggers------//
  //---------PbPb 2018 Data---------//
  bool HLT_g15_loose_ion = false;
  bool HLT_g20_loose_ion = false;
  bool HLT_g20_loose = false;
  bool HLT_g30_loose_ion = false;
  bool HLT_g13_etcut_ion = false;
  bool HLT_g18_etcut_ion = false;
  bool HLT_g18_etcut = false;
  bool HLT_g28_etcut_ion = false;
  bool HLT_g50_loose_ion = false;
  bool HLT_noalg_L1EM10 = false;
  bool HLT_noalg_L1EM12 = false;

 std::vector<bool>* m_HLT_photon_loose_MC152 = NULL;

  std::vector<float>* b_HLT_pt = NULL;
  std::vector<float>* b_HLT_phi = NULL;
  std::vector<float>* b_HLT_eta = NULL;
  std::vector<float>* b_HLT_etaBE = NULL;

  int b_HLT_photons_n = 0;


  //---L1RO                                                                                                                                                                      
  int b_L1_n = 0;
  std::vector<float>* L1RO_pt = NULL;
  std::vector<float>* L1RO_eta = NULL;
  std::vector<float>* L1RO_phi = NULL;

  std::vector<bool>* L1RO_EM3 = NULL;
  std::vector<bool>* L1RO_EM5 = NULL;
  std::vector<bool>* L1RO_EM10 = NULL;
  std::vector<bool>* L1RO_EM12 = NULL;
  std::vector<bool>* L1RO_EM16 = NULL;
  std::vector<bool>* L1RO_EM20 = NULL;

  //---Offline Photons                                                                                                                                                           
  int photon_n = 0;
  std::vector<float>* photon_pt = NULL;
  std::vector<float>* photon_eta = NULL;
  std::vector<float>* photon_phi = NULL;
  std::vector<bool>* photon_loose = NULL;
  std::vector<bool>* photon_tight = NULL;


 
  //--------------------------------//
  //-------Loading Up Files---------//
  //--------------------------------//

  TChain *t = new TChain("analysis");

  t->Add("/Users/berenicegarcia/Desktop/ATLAS/Qual-Task/PbPb_Data/user.berenice.21817306._000001.ANALYSIS.root");
  t->Add("/Users/berenicegarcia/Desktop/ATLAS/Qual-Task/PbPb_Data/user.berenice.21817306._000002.ANALYSIS.root");
  t->Add("/Users/berenicegarcia/Desktop/ATLAS/Qual-Task/PbPb_Data/user.berenice.21817306._000003.ANALYSIS.root");


  /*Branches For Data Vectors*/

 /*Branches For Data Vectors*/

  //---FCAL
  t->SetBranchAddress("fcalA_et", &fcalA_et_data);
  t->SetBranchAddress("fcalC_et", &fcalC_et_data);

  



  //---HLT
  t->SetBranchAddress("HLT_noalg_eb_L1TE50", &HLT_noalg_eb_L1TE50);
  t->SetBranchAddress("HLT_noalg_eb_L1ZDC_A_C_VTE50", &HLT_noalg_eb_L1ZDC_A_C_VTE50);
  t->SetBranchAddress("HLT_noalg_eb_L1MU4", &HLT_noalg_eb_L1MU4);


  //----HLT Photon Triggers-----//
  t->SetBranchAddress("HLT_g15_loose_ion", &HLT_g15_loose_ion);
  t->SetBranchAddress("HLT_g20_loose_ion", &HLT_g20_loose_ion);
  t->SetBranchAddress("HLT_g20_loose", &HLT_g20_loose);
  t->SetBranchAddress("HLT_g30_loose_ion", &HLT_g30_loose_ion);
  t->SetBranchAddress("HLT_g13_etcut_ion", &HLT_g13_etcut_ion);
  t->SetBranchAddress("HLT_g18_etcut_ion", &HLT_g18_etcut_ion);
  t->SetBranchAddress("HLT_g18_etcut", &HLT_g18_etcut);
  t->SetBranchAddress("HLT_g28_etcut_ion", &HLT_g28_etcut_ion);
  t->SetBranchAddress("HLT_g50_loose_ion", &HLT_g50_loose_ion);
  t->SetBranchAddress("HLT_noalg_L1EM10", &HLT_noalg_L1EM10);
  t->SetBranchAddress("HLT_noalg_L1EM12", &HLT_noalg_L1EM12);

  t->SetBranchAddress("b_HLT_photons_n", &b_HLT_photons_n);
  t->SetBranchAddress("m_HLT_photon_loose_MC152", &m_HLT_photon_loose_MC152);

  t->SetBranchAddress("b_HLT_pt", &b_HLT_pt);
  t->SetBranchAddress("b_HLT_eta", &b_HLT_eta);
  t->SetBranchAddress("b_HLT_etaBE", &b_HLT_etaBE);
  t->SetBranchAddress("b_HLT_phi", &b_HLT_phi);



  //---L1RO                                                                                                                                                                      
  cout << "Setting Branch Addresses For L1ROs..." << endl;
  cout << endl;
  
  t->SetBranchAddress("b_L1_n", &b_L1_n);
  t->SetBranchAddress("L1RO_pt", &L1RO_pt);
  t->SetBranchAddress("L1RO_eta", &L1RO_eta);
  t->SetBranchAddress("L1RO_phi", &L1RO_phi);

  t->SetBranchAddress("L1RO_EM3", &L1RO_EM3);
  t->SetBranchAddress("L1RO_EM5", &L1RO_EM5);
  t->SetBranchAddress("L1RO_EM10", &L1RO_EM10);
  t->SetBranchAddress("L1RO_EM12", &L1RO_EM12);
  t->SetBranchAddress("L1RO_EM16", &L1RO_EM16);
  t->SetBranchAddress("L1RO_EM20", &L1RO_EM20);
  cout << "Done..." << endl;

  cout << "Setting Branch Address For Photons Offline.." << endl;

  t->SetBranchAddress("photon_n", &photon_n);
  t->SetBranchAddress("photon_pt", &photon_pt);
  t->SetBranchAddress("photon_eta", &photon_eta);
  t->SetBranchAddress("photon_phi", &photon_phi);
  t->SetBranchAddress("photon_loose", &photon_loose);
  t->SetBranchAddress("photon_tight", &photon_tight);

  cout << "Done .." << endl;

  event_num = t->GetEntries();



 for (int iEvent = 0; iEvent < event_num; ++iEvent){

  t->GetEntry(iEvent);
  

  if(HLT_g15_loose_ion){
    HLT_g15_loose_ion_histo->Fill(HLT_g15_loose_ion);
  }
  if(HLT_g20_loose_ion){
    HLT_g20_loose_ion_histo->Fill(HLT_g20_loose_ion);
  }
  if(HLT_g20_loose){
    HLT_g20_loose_histo->Fill(HLT_g20_loose);
  }
  if(HLT_g30_loose_ion){
    HLT_g30_loose_ion_histo->Fill(HLT_g30_loose_ion);
  }
  if(HLT_g13_etcut_ion){
    HLT_g13_etcut_ion_histo->Fill(HLT_g13_etcut_ion);
  }
  if(HLT_g18_etcut_ion){
    HLT_g18_etcut_ion_histo->Fill(HLT_g18_etcut_ion);
  }
  if(HLT_g18_etcut){
    HLT_g18_etcut_histo->Fill(HLT_g18_etcut);
  }
  if(HLT_g28_etcut_ion){
    HLT_g28_etcut_ion_histo->Fill(HLT_g28_etcut_ion);
  }
  if(HLT_g50_loose_ion){
    HLT_g50_loose_ion_histo->Fill(HLT_g50_loose_ion);
  }
  if(HLT_noalg_L1EM10){
    HLT_noalg_L1EM10_histo->Fill(HLT_noalg_L1EM10);
  }
  if(HLT_noalg_L1EM12){
    HLT_noalg_L1EM12_histo->Fill(HLT_noalg_L1EM12);
  }
  if(HLT_noalg_eb_L1TE50){
    HLT_noalg_eb_L1TE50_histo->Fill(HLT_noalg_eb_L1TE50);
  }



  //cout << "This is event:  " << iEvent << endl;

  //if(HLT_noalg_eb_L1ZDC_A_C_VTE50 !=1){continue;}   //----MB Trigger
  //if(HLT_noalg_eb_L1MU4 !=1){continue;}
  
  if(HLT_noalg_eb_L1TE50 !=1){continue;}  //---Requires there to be 10 GeV of energy trhoughout the whole Calorimeter (MB Trigger) 
  //if(HLT_noalg_L1EM10!=1){continue;}      //---Requiremment only at L1. This trigger will require to have 10 GeV L1 EM RoI and no ther requirement. 



 

  //if(SumEt_data_fire < 10){continue;}

  bool check = true;
  
  //---L1EM RoI Triggers
  bool check_trigger_15 = false;
  bool check_trigger_25 = false;
  bool check_trigger_35 = false;

  //---HLT Triggers
  bool check_trigger_hlt = false;
  bool check_trigger_hlt_25 = false;
  bool check_trigger_hlt_35 = false;
  bool check_trigger_hlt_50 = false;

 L1EM_num->Fill(b_L1_n);

  for (int iPhoton = 0; iPhoton < photon_n; ++iPhoton){
    
    if(photon_tight->at(iPhoton) != 1){continue;}            //------Tight Photons                                                                           
    //if(photon_loose->at(iPhoton) != 1){continue;}             //-------Loose Photons
    if(photon_pt->at(iPhoton) < 15){continue;}               //--------Photong > 20 GeV                                                                                                                                                            
    
    
    //---Endcaps
    bool state_4 = photon_eta->at(iPhoton) > 1.56;
    bool state_5 = photon_eta->at(iPhoton) < 2.37;
    bool state_6 = photon_eta->at(iPhoton) > -2.37;
    bool state_7 = photon_eta->at(iPhoton) < -1.56;
    
    //---Barrel
    bool state_8 = photon_eta->at(iPhoton) < 1.37;
    bool state_9 = photon_eta->at(iPhoton) > -1.37;
    
    
    //---SumEt
    float sumET_tot = fcalA_et_data + fcalC_et_data;

    if ( ((state_4 && state_5) || (state_6 && state_7) || (state_8 && state_9)) == 1 ) {
      
      all_photons_phi->Fill(photon_phi->at(iPhoton));
      all_photons_eta->Fill(photon_eta->at(iPhoton));
      all_photons->Fill(photon_pt->at(iPhoton));

      if(sumET_tot > 2989.31){//---0-10% Centrality
        all_photons_010->Fill(photon_pt->at(iPhoton));
      }else if((sumET_tot < 2989.31) && (sumET_tot > 1368.75)){ //-----10-30% Centrality
        all_photons_1030->Fill(photon_pt->at(iPhoton));
      }else if((sumET_tot < 1368.75) && (sumET_tot > 525.092)){//-----30-50% Centrality
        all_photons_3050->Fill(photon_pt->at(iPhoton));
      }else if((sumET_tot< 525.092) && (sumET_tot > 63.719)){//-----50-80% CEntrality
        all_photons_5080->Fill(photon_pt->at(iPhoton));
      }

    } //--Encaps && Barrel
    
    if((state_8 && state_9) == 1){
      all_photons_barrel->Fill(photon_pt->at(iPhoton));
      if(sumET_tot > 2989.31){//---0-10% Centrality
        all_photons_010_b->Fill(photon_pt->at(iPhoton));
      }else if((sumET_tot < 2989.31) && (sumET_tot > 1368.75)){ //-----10-30% Centrality
        all_photons_1030_b->Fill(photon_pt->at(iPhoton));
      }else if((sumET_tot < 1368.75) && (sumET_tot > 525.092)){//-----30-50% Centrality
        all_photons_3050_b->Fill(photon_pt->at(iPhoton));
      }else if((sumET_tot< 525.092) && (sumET_tot > 63.719)){//-----50-80% CEntrality
        all_photons_5080_b->Fill(photon_pt->at(iPhoton));
      }
    }//-----Photons that ONLY went through the barrel
    
    if ( ((state_4 && state_5) || (state_6 && state_7)) == 1 ) {
      all_photons_endcaps->Fill(photon_pt->at(iPhoton));
      if(sumET_tot > 2989.31){//---0-10% Centrality
        all_photons_010_e->Fill(photon_pt->at(iPhoton));
      }else if((sumET_tot < 2989.31) && (sumET_tot > 1368.75)){ //-----10-30% Centrality
        all_photons_1030_e->Fill(photon_pt->at(iPhoton));
      }else if((sumET_tot < 1368.75) && (sumET_tot > 525.092)){//-----30-50% Centrality
        all_photons_3050_e->Fill(photon_pt->at(iPhoton));
      }else if((sumET_tot< 525.092) && (sumET_tot > 63.719)){//-----50-80% CEntrality
        all_photons_5080_e->Fill(photon_pt->at(iPhoton));
      }

    }  //------Photons went through the endcaps ONLY
    
     


    //------Make sure that this photon went through the barrel or endcaps
    
    
    float SumEt_value = fcalC_et_data + fcalA_et_data;
    
    

    all_photons_sumET->Fill(SumEt_value*(0.001));
   
    pT_vs_L1EM->Fill(photon_pt->at(iPhoton),b_L1_n);

    L1EM_num->Fill(b_L1_n);


    for(int iTrig = 0; iTrig < b_L1_n; iTrig++){

      l1_pT_dis->Fill(L1RO_pt->at(iTrig));

      float deltaEta = abs(L1RO_eta->at(iTrig)-photon_eta->at(iPhoton));
      float deltaPhi = delta_phi(L1RO_phi->at(iTrig),photon_phi->at(iPhoton));
      float deltaR = sqrt(pow(deltaEta,2) + pow(deltaPhi,2));
    
      //cout << "This is L1 info: pt / eta / phi: " << L1RO_pt->at(iTrig) << " / " << L1RO_eta->at(iTrig) << " / " << L1RO_phi->at(iTrig) << endl;                                                                                                                                            
      
     
      if((L1RO_pt->at(iTrig) > 15) && (deltaR < 0.15)){
        check_trigger_15 = true; //---This Photon Matched W/ the L1EMX Trigger Fired
      }
      if((L1RO_pt->at(iTrig) > 25) && (deltaR < 0.15)){
        check_trigger_25 = true;
      }
      if((L1RO_pt->at(iTrig) > 35) && (deltaR < 0.15)){
        check_trigger_35 = true;
      }

      }// L1 Triger Loop        

  
    pT_vs_HLT_num->Fill(photon_pt->at(iPhoton),b_HLT_photons_n);
    HLT_num->Fill(b_HLT_photons_n);
    //cout << "This event is: " << iEvent << endl;
    //cout << "This is how many HLT Photons we have: " << b_HLT_photons_n << endl;
    //cout << "For event: " << iEvent << endl;
   
     for(int iHLT = 0; iHLT < b_HLT_photons_n; iHLT++){
      hlt_pT_dis->Fill(b_HLT_pt->at(iHLT));
      hlt_phi_dis->Fill(b_HLT_phi->at(iHLT));
      hlt_eta_dis->Fill(abs(b_HLT_eta->at(iHLT)));

      //cout << "Is this HLT photon loose: " << m_HLT_photon_loose_MC152->at(iHLT) << endl;
      
      //cout << "This is HLT info: pt / eta / phi / photon loose? : " << b_HLT_pt->at(iHLT) << " / " << b_HLT_eta->at(iHLT) << " / " << b_HLT_phi->at(iHLT) << " / " << m_HLT_photon_loose_MC152->at(iHLT) <<  endl;  

      if(m_HLT_photon_loose_MC152->at(iHLT) != 1){continue;} //----Make Sure We Look at loose HLT photons
    
      float deltaEta_hlt = abs(b_HLT_eta->at(iHLT)-photon_eta->at(iPhoton));
      float deltaPhi_hlt = delta_phi(b_HLT_phi->at(iHLT),photon_phi->at(iPhoton));
      float deltaR_hlt = sqrt(pow(deltaEta_hlt,2) + pow(deltaPhi_hlt,2));


                                                                                                                                            

      if((b_HLT_pt->at(iHLT) > 15) && (deltaR_hlt < 0.15)){
        check_trigger_hlt = true; //---This Photon Matched W/ the HLT Trigger Fired
        //cout << "This photon matched to a HLT > 15 GeV!! This is the pT: " << photon_pt->at(iPhoton) << endl;
      }
      if((b_HLT_pt->at(iHLT) > 25) && (deltaR_hlt < 0.15)){
        check_trigger_hlt_25 = true;
        //cout << "This photon matched to a HLT > 25 GeV!! This is the pT: " << photon_pt->at(iPhoton) << endl;
      }
      if((b_HLT_pt->at(iHLT) > 35) && (deltaR_hlt < 0.15)){
        check_trigger_hlt_35 = true;
        //cout << "This photon matched to a HLT > 35 GeV!! This is the pT: " << photon_pt->at(iPhoton) << endl;
      }



    }// HLT Loop


  //---Endcaps
  bool req_4 = photon_eta->at(iPhoton) > 1.56;
  bool req_5 = photon_eta->at(iPhoton) < 2.37;
  bool req_6 = photon_eta->at(iPhoton) > -2.37;
  bool req_7 = photon_eta->at(iPhoton) < -1.56;
    
  //---Barrel
  bool req_8 = photon_eta->at(iPhoton) < 1.37;
  bool req_9 = photon_eta->at(iPhoton) > -1.37;

  float SumEt_data_total = fcalC_et_data + fcalA_et_data;


  //-------Photon Matched W/ L1EM RoI Trigger  
    if(check_trigger_15){
    
    if ( ((req_4 && req_5) || (req_6 && req_7) || (req_8 && req_9)) == 1 ) {
        photons_wL1EM_15GeV->Fill(photon_pt->at(iPhoton));
        
        if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
          photon_w_L1EM_15GeV_010->Fill(photon_pt->at(iPhoton));
        }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
          photon_w_L1EM_15GeV_1030->Fill(photon_pt->at(iPhoton));
        }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
          photon_w_L1EM_15GeV_3050->Fill(photon_pt->at(iPhoton));
        }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
          photon_w_L1EM_15GeV_5080->Fill(photon_pt->at(iPhoton));
        }

     } //--Encaps && Barrel
    
    if((req_8 && req_9) == 1){
      photons_wL1EM_15GeV_barrel->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_L1EM_15GeV_010_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_L1EM_15GeV_1030_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_L1EM_15GeV_3050_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_L1EM_15GeV_5080_b->Fill(photon_pt->at(iPhoton));
      }


    }//-----Photons that ONLY went through the barrel
    
    if ( ((req_4 && req_5) || (state_6 && state_7)) == 1 ){
    photons_wL1EM_15GeV_endcaps->Fill(photon_pt->at(iPhoton));
    if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
      photon_w_L1EM_15GeV_010_e->Fill(photon_pt->at(iPhoton));
    }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
      photon_w_L1EM_15GeV_1030_e->Fill(photon_pt->at(iPhoton));
    }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
      photon_w_L1EM_15GeV_3050_e->Fill(photon_pt->at(iPhoton));
    }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
      photon_w_L1EM_15GeV_5080_e->Fill(photon_pt->at(iPhoton));
    }

    }//------Photons went through the endcaps ONLY

                                                                                             
    check_trigger_15 = false;

  }
  if(check_trigger_25){
    
    if ( ((req_4 && req_5) || (req_6 && req_7) || (req_8 && req_9)) == 1 ) {
      photons_wL1EM_25GeV->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_L1EM_25GeV_010->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
       photon_w_L1EM_25GeV_1030->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
       photon_w_L1EM_25GeV_3050->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
       photon_w_L1EM_25GeV_5080->Fill(photon_pt->at(iPhoton));
      }

    } //--Encaps && Barrel
    if((req_8 && req_9) == 1){
      photons_wL1EM_25GeV_barrel->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_L1EM_25GeV_010_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_L1EM_25GeV_1030_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
       photon_w_L1EM_25GeV_3050_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
       photon_w_L1EM_25GeV_5080_b->Fill(photon_pt->at(iPhoton));
      }
    }//-----Photons that ONLY went through the barrel
    
    if ( ((req_4 && req_5) || (state_6 && state_7)) == 1 ) {
      photons_wL1EM_25GeV_endcaps->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_L1EM_25GeV_010_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_L1EM_25GeV_1030_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_L1EM_25GeV_3050_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_L1EM_25GeV_5080_e->Fill(photon_pt->at(iPhoton));
      }
    } //------Photons went through the endcaps ONLY



    check_trigger_25 = false;
  }
    if(check_trigger_35){
    if ( ((req_4 && req_5) || (req_6 && req_7) || (req_8 && req_9)) == 1 ) {
      photons_wL1EM_35GeV->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_L1EM_35GeV_010->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_L1EM_35GeV_1030->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_L1EM_35GeV_3050->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_L1EM_35GeV_5080->Fill(photon_pt->at(iPhoton));
      }

    } //--Encaps && Barrel
    
    if((req_8 && req_9) == 1){
      photons_wL1EM_35GeV_barrel->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_L1EM_35GeV_010_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_L1EM_35GeV_1030_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_L1EM_35GeV_3050_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_L1EM_35GeV_5080_b->Fill(photon_pt->at(iPhoton));
      }

    }                 //-----Photons that ONLY went through the barrel
    if ( ((req_4 && req_5) || (state_6 && state_7)) == 1 ) {
      photons_wL1EM_35GeV_endcaps->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_L1EM_35GeV_010_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_L1EM_35GeV_1030_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_L1EM_35GeV_3050_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_L1EM_35GeV_5080_e->Fill(photon_pt->at(iPhoton));
      }

    }  //------Photons went through the endcaps ONLY


   
    check_trigger_35 = false;
  }

  //-------Photon Matched W/ HLT Trigger                                                                                                                             
  if(check_trigger_hlt){
    
    if ( ((req_4 && req_5) || (req_6 && req_7) || (req_8 && req_9)) == 1 ){
      photons_w_HLT_15->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_HLT_15GeV_010->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_HLT_15GeV_1030->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_HLT_15GeV_3050->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_HLT_15GeV_5080->Fill(photon_pt->at(iPhoton));
      }

    } //--Encaps && Barrel
    if((req_8 && req_9) == 1){
      photons_w_HLT_15_barrel->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_HLT_15GeV_010_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_HLT_15GeV_1030_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_HLT_15GeV_3050_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_HLT_15GeV_5080_b->Fill(photon_pt->at(iPhoton));
      }

    }                 //-----Photons that ONLY went through the barrel
    if ( ((req_4 && req_5) || (state_6 && state_7)) == 1 ) {
      photons_w_HLT_15_endcaps->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_HLT_15GeV_010_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_HLT_15GeV_1030_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_HLT_15GeV_3050_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_HLT_15GeV_5080_e->Fill(photon_pt->at(iPhoton));
      }
    }  //------Photons went through the endcaps ONLY

                                                                                             
    check_trigger_hlt = false;

  }
  if(check_trigger_hlt_25){
    
    if ( ((req_4 && req_5) || (req_6 && req_7) || (req_8 && req_9)) == 1 ){
      photons_w_HLT_25->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_HLT_25GeV_010->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_HLT_25GeV_1030->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_HLT_25GeV_3050->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_HLT_25GeV_5080->Fill(photon_pt->at(iPhoton));
      }

    } //--Encaps && Barrel
    if((req_8 && req_9) == 1){
      photons_w_HLT_25_barrel->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_HLT_25GeV_010_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_HLT_25GeV_1030_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_HLT_25GeV_3050_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_HLT_25GeV_5080_b->Fill(photon_pt->at(iPhoton));
      }
    }                 //-----Photons that ONLY went through the barrel
    if ( ((req_4 && req_5) || (state_6 && state_7)) == 1 ) {
      photons_w_HLT_25_endcaps->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_HLT_25GeV_010_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_HLT_25GeV_1030_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_HLT_25GeV_3050_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_HLT_25GeV_5080_e->Fill(photon_pt->at(iPhoton));
      }
    }  //------Photons went through the endcaps ONLY



    float num_sumET_2 = fcalC_et_data + fcalA_et_data;
    hlt_fire_sumET_25->Fill(num_sumET_2*(0.001));
   
    check_trigger_hlt_25 = false;
  }
  if(check_trigger_hlt_35!=1 && (photon_pt->at(iPhoton) > 50)){
    float sumET_num = fcalC_et_data + fcalA_et_data;
     
    photons_nofire_HLT_pt->Fill(photon_pt->at(iPhoton));
    photons_nofire_HLT_eta->Fill(photon_eta->at(iPhoton));
    photons_nofire_HLT_sumET->Fill(sumET_num*(0.001));
  }
  if(check_trigger_hlt_35){
    if ( ((req_4 && req_5) || (req_6 && req_7) || (req_8 && req_9)) == 1 ) {
      photons_w_HLT_35->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total > 2989.31){                                           //---0-10% Centrality
        photon_w_HLT_35GeV_010->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_HLT_35GeV_1030->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_HLT_35GeV_3050->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_HLT_35GeV_5080->Fill(photon_pt->at(iPhoton));
      }
    } //--Encaps && Barrel
    if((req_8 && req_9) == 1){
      photons_w_HLT_35_barrel->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total < 2989.31){                                           //---0-10% Centrality
        photon_w_HLT_35GeV_010_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_HLT_35GeV_1030_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_HLT_35GeV_3050_b->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_HLT_35GeV_5080_b->Fill(photon_pt->at(iPhoton));
      }
    }                 //-----Photons that ONLY went through the barrel
    if ( ((req_4 && req_5) || (state_6 && state_7)) == 1 ) {
      photons_w_HLT_35_endcaps->Fill(photon_pt->at(iPhoton));
      if(SumEt_data_total < 2989.31){                                           //---0-10% Centrality
        photon_w_HLT_35GeV_010_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 2989.31) && (SumEt_data_total > 1368.75)){ //-----10-30% Centrality
        photon_w_HLT_35GeV_1030_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 1368.75) && (SumEt_data_total > 525.092)){//-----30-50% Centrality
        photon_w_HLT_35GeV_3050_e->Fill(photon_pt->at(iPhoton));
      }else if((SumEt_data_total < 525.092) && (SumEt_data_total > 63.719)){//-----50-80% CEntrality
        photon_w_HLT_35GeV_5080_e->Fill(photon_pt->at(iPhoton));
      }
    }  //------Photons went through the endcaps ONLY




    float num_sumET_3 = fcalC_et_data + fcalA_et_data;
    hlt_fire_sumET_35->Fill(num_sumET_3*(0.001));
   
    check_trigger_hlt_35 = false;
  }


  if((check_trigger_hlt!= 1)&&(photon_pt->at(iPhoton) > 50)){
    float num_sumET_15 = fcalC_et_data + fcalA_et_data;
    hlt_nofire_sumET_15->Fill(num_sumET_15*(0.001));

  }
  if(check_trigger_hlt_25!=1 && (photon_pt->at(iPhoton) > 50)){
    float num_sumET_25 = fcalC_et_data + fcalA_et_data;
    hlt_nofire_sumET_25->Fill(num_sumET_25*(0.001));

  }
  if(check_trigger_hlt_35!=1 && (photon_pt->at(iPhoton) >50)){
    float num_sumET_35 = fcalC_et_data + fcalA_et_data;
    hlt_nofire_sumET_35->Fill(num_sumET_35*(0.001));

  }
  

  }// Photon Loop                                            


   

} //Event Loop




cout << "Number of events that fire HLT_g15_loose_ion : " << HLT_g15_loose_ion_histo->Integral() << endl;
cout << "Number of events that fire HLT_g20_loose_ion : " << HLT_g20_loose_ion_histo->Integral() << endl;
cout << "Number of events that fire HLT_g20_loose: " << HLT_g20_loose_histo->Integral() << endl;
cout << "Number of events that fire HLT_g30_loose_ion: " << HLT_g30_loose_ion_histo->Integral() << endl;
cout << "Number of events that fire HLT_g13_etcut_ion: " << HLT_g13_etcut_ion_histo->Integral() << endl;
cout << "Number of events that fire HLT_g18_etcut_ion: " << HLT_g18_etcut_ion_histo->Integral() << endl;
cout << "Number of events that fire HLT_g18_etcut: " << HLT_g18_etcut_histo->Integral() << endl;
cout << "Number of events that fire HLT_g28_etcut_ion: " << HLT_g28_etcut_ion_histo->Integral() << endl;
cout << "Number of events that fire HLT_g50_loose_ion: " << HLT_g50_loose_ion_histo->Integral() << endl;
cout << "Number of events that fire HLT_noalg_L1EM10: " << HLT_noalg_L1EM10_histo->Integral() << endl;
cout << "Number of events that fire HLT_noalg_L1EM12: " << HLT_noalg_L1EM12_histo->Integral() << endl << endl;

cout << "Number of events that fire HLT_noalg_L1EM10: " << HLT_noalg_L1EM10_histo->Integral() << endl;
cout << "Number of events that fire HLT_noalg_eb_L1TE50: " << HLT_noalg_eb_L1TE50_histo->Integral() << endl;

TFile *file_root = new TFile("/Users/berenicegarcia/Desktop/ATLAS/Qual-Task/histograms.root","RECREATE");
file_root->cd();

//---------------------------------------//
//----------ALl Offline Photons----------//
//---------------------------------------//
all_photons->Write("all_photons",TObject::kOverwrite);
all_photons_barrel->Write("all_photons_endcaps",TObject::kOverwrite);
all_photons_endcaps->Write("all_photons_barrel",TObject::kOverwrite);


//---------------------------------------//
//--------L1EM RoI pT Distributions------//
//---------------------------------------//
photons_wL1EM_15GeV->Write("photons_wL1EM_15GeV",TObject::kOverwrite);
photons_wL1EM_15GeV_barrel->Write("photons_wL1EM_15GeV_barrel", TObject::kOverwrite);
photons_wL1EM_15GeV_endcaps->Write("photons_wL1EM_15GeV_endcaps",TObject::kOverwrite);

photons_wL1EM_25GeV->Write("photons_wL1EM_25GeV",TObject::kOverwrite);
photons_wL1EM_25GeV_barrel->Write("photons_wL1EM_15GeV_barrel",TObject::kOverwrite);
photons_wL1EM_25GeV_endcaps->Write("photons_wL1EM_25GeV_endcaps", TObject::kOverwrite);

photons_wL1EM_35GeV->Write("photons_wL1EM_35GeV", TObject::kOverwrite);
photons_wL1EM_35GeV_barrel->Write("photons_wL1EM_35GeV_barrel", TObject::kOverwrite);
photons_wL1EM_35GeV_endcaps->Write("photons_wL1EM_35GeV_endcaps", TObject::kOverwrite);

//--------------------------------------//
//--------HLT pT Distribution-----------//
//--------------------------------------//
photons_w_HLT_15->Write("photons_w_HLT_15",TObject::kOverwrite);
photons_w_HLT_25->Write("photons_w_HLT_25", TObject::kOverwrite);
photons_w_HLT_35->Write("photons_w_HLT_35", TObject::kOverwrite);

photons_w_HLT_15_barrel->Write("photons_w_HLT_15_barrel",TObject::kOverwrite);
photons_w_HLT_25_barrel->Write("photons_w_HLT_25_barrel", TObject::kOverwrite);
photons_w_HLT_35_barrel->Write("photons_w_HLT_35_barrel", TObject::kOverwrite);

photons_w_HLT_15_endcaps->Write("photons_w_HLT_15_endcaps", TObject::kOverwrite);
photons_w_HLT_25_endcaps->Write("photons_w_HLT_25_endcaps", TObject::kOverwrite);
photons_w_HLT_35_endcaps->Write("photons_w_HLT_35_endcaps", TObject::kOverwrite);

//--------------------------------------------------//
//------Photons Failed TO Fire A HLT Trigger--------//
//--------------------------------------------------//

photons_nofire_HLT_pt->Write("photons_nofire_HLT_pt",TObject::kOverwrite);
photons_nofire_HLT_eta->Write("photons_nofire_HLT_eta", TObject::kOverwrite);
photons_nofire_HLT_sumET->Write("photons_nofire_HLT_sumET", TObject::kOverwrite);


//--------------------------------------------------------//
//-----L1 & HLT Efficiencies Dependent onn Centrality-----//
//--------------------------------------------------------//
//**********************************L1 Histograms**************************************// 

//----Barrel & EndCaps
//----All Photons
all_photons_010->Write("all_photons_010", TObject::kOverwrite); 
all_photons_1030->Write("all_photons_1030", TObject::kOverwrite);
all_photons_3050->Write("all_photons_3050", TObject::kOverwrite);
all_photons_5080->Write("all_photons_5080", TObject::kOverwrite);


photon_w_L1EM_15GeV_010->Write("photon_w_L1EM_15GeV_010", TObject::kOverwrite);
photon_w_L1EM_15GeV_1030->Write("photon_w_L1EM_15GeV_1030", TObject::kOverwrite);
photon_w_L1EM_15GeV_3050->Write("photon_w_L1EM_15GeV_3050", TObject::kOverwrite);
photon_w_L1EM_15GeV_5080->Write("photon_w_L1EM_15GeV_5080", TObject::kOverwrite);

photon_w_L1EM_25GeV_010->Write("photon_w_L1EM_25GeV_010", TObject::kOverwrite);
photon_w_L1EM_25GeV_1030->Write("photon_w_L1EM_25GeV_1030", TObject::kOverwrite);
photon_w_L1EM_25GeV_3050->Write("photon_w_L1EM_25GeV_3050", TObject::kOverwrite);
photon_w_L1EM_25GeV_5080->Write("photon_w_L1EM_25GeV_5080", TObject::kOverwrite);

photon_w_L1EM_35GeV_010->Write("photon_w_L1EM_35GeV_010", TObject::kOverwrite);
photon_w_L1EM_35GeV_1030->Write("photon_w_L1EM_35GeV_1030", TObject::kOverwrite);
photon_w_L1EM_35GeV_3050->Write("photon_w_L1EM_35GeV_3050", TObject::kOverwrite);
photon_w_L1EM_35GeV_5080->Write("photon_w_L1EM_35GeV_5080", TObject::kOverwrite);
  
//----Barrel
all_photons_010_b->Write("all_photons_010_b", TObject::kOverwrite); 
all_photons_1030_b->Write("all_photons_1030_b", TObject::kOverwrite);
all_photons_3050_b->Write("all_photons_3050_b", TObject::kOverwrite);
all_photons_5080_b->Write("all_photons_5080_b", TObject::kOverwrite);


photon_w_L1EM_15GeV_010_b->Write("photon_w_L1EM_15GeV_010_b", TObject::kOverwrite);
photon_w_L1EM_15GeV_1030_b->Write("photon_w_L1EM_15GeV_1030_b", TObject::kOverwrite);
photon_w_L1EM_15GeV_3050_b->Write("photon_w_L1EM_15GeV_3050_b",TObject::kOverwrite);
photon_w_L1EM_15GeV_5080_b->Write("photon_w_L1EM_15GeV_5080_b", TObject::kOverwrite);

photon_w_L1EM_25GeV_010_b->Write("photon_w_L1EM_25GeV_010_b", TObject::kOverwrite);
photon_w_L1EM_25GeV_1030_b->Write("photon_w_L1EM_25GeV_1030_b", TObject::kOverwrite);
photon_w_L1EM_25GeV_3050_b->Write("photon_w_L1EM_25GeV_3050_b", TObject::kOverwrite);
photon_w_L1EM_25GeV_5080_b->Write("photon_w_L1EM_25GeV_5080_b", TObject::kOverwrite);

photon_w_L1EM_35GeV_010_b->Write("photon_w_L1EM_35GeV_010_b", TObject::kOverwrite); 
photon_w_L1EM_35GeV_1030_b->Write("photon_w_L1EM_35GeV_1030_b", TObject::kOverwrite); 
photon_w_L1EM_35GeV_3050_b->Write("photon_w_L1EM_35GeV_3050_b", TObject::kOverwrite); 
photon_w_L1EM_35GeV_5080_b->Write("photon_w_L1EM_35GeV_5080_b", TObject::kOverwrite); 
  

//----Endcaps
all_photons_010_e->Write("all_photons_010_e", TObject::kOverwrite); 
all_photons_1030_e->Write("all_photons_1030_e", TObject::kOverwrite);
all_photons_3050_e->Write("all_photons_3050_e", TObject::kOverwrite);
all_photons_5080_e->Write("all_photons_5080_e", TObject::kOverwrite);


photon_w_L1EM_15GeV_010_e->Write("photon_w_L1EM_15GeV_010_e", TObject::kOverwrite); 
photon_w_L1EM_15GeV_1030_e->Write("photon_w_L1EM_15GeV_1030_e", TObject::kOverwrite); 
photon_w_L1EM_15GeV_3050_e->Write("photon_w_L1EM_15GeV_3050_e", TObject::kOverwrite); 
photon_w_L1EM_15GeV_5080_e->Write("photon_w_L1EM_15GeV_5080_e", TObject::kOverwrite); 

photon_w_L1EM_25GeV_010_e->Write("photon_w_L1EM_25GeV_010_e", TObject::kOverwrite);
photon_w_L1EM_25GeV_1030_e->Write("photon_w_L1EM_25GeV_1030_e", TObject::kOverwrite); 
photon_w_L1EM_25GeV_3050_e->Write("photon_w_L1EM_25GeV_3050_e", TObject::kOverwrite);
photon_w_L1EM_25GeV_5080_e->Write("photon_w_L1EM_25GeV_5080_e", TObject::kOverwrite);

photon_w_L1EM_35GeV_010_e->Write("photon_w_L1EM_35GeV_010_e", TObject::kOverwrite);
photon_w_L1EM_35GeV_1030_e->Write("photon_w_L1EM_35GeV_1030_e", TObject::kOverwrite);
photon_w_L1EM_35GeV_3050_e->Write("photon_w_L1EM_35GeV_3050_e", TObject::kOverwrite); 
photon_w_L1EM_35GeV_5080_e->Write("photon_w_L1EM_35GeV_5080_e", TObject::kOverwrite); 
//*************************************************************************************//

//---------------------------------------------------------------------------------//
//-------------------HLT Efficiencies Dependent on Centrality----------------------//
//---------------------------------------------------------------------------------//

//----Barrel & EndCaps
photon_w_HLT_15GeV_010->Write("photon_w_HLT_15GeV_010", TObject::kOverwrite);
photon_w_HLT_15GeV_1030->Write("photon_w_HLT_15GeV_1030", TObject::kOverwrite);
photon_w_HLT_15GeV_3050->Write("photon_w_HLT_15GeV_3050", TObject::kOverwrite);
photon_w_HLT_15GeV_5080->Write("photon_w_HLT_15GeV_5080", TObject::kOverwrite);

photon_w_HLT_25GeV_010->Write("photon_w_HLT_25GeV_010", TObject::kOverwrite);
photon_w_HLT_25GeV_1030->Write("photon_w_HLT_25GeV_1030", TObject::kOverwrite);
photon_w_HLT_25GeV_3050->Write("photon_w_HLT_25GeV_3050", TObject::kOverwrite);
photon_w_HLT_25GeV_5080->Write("photon_w_HLT_25GeV_5080", TObject::kOverwrite);

photon_w_HLT_35GeV_010->Write("photon_w_HLT_35GeV_010", TObject::kOverwrite);
photon_w_HLT_35GeV_1030->Write("photon_w_HLT_35GeV_1030", TObject::kOverwrite);
photon_w_HLT_35GeV_3050->Write("photon_w_HLT_35GeV_3050", TObject::kOverwrite);
photon_w_HLT_35GeV_5080->Write("photon_w_HLT_35GeV_5080", TObject::kOverwrite);
  
//----Barrel
photon_w_HLT_15GeV_010_b->Write("photon_w_HLT_15GeV_010_b", TObject::kOverwrite);
photon_w_HLT_15GeV_1030_b->Write("photon_w_HLT_15GeV_1030_b", TObject::kOverwrite); 
photon_w_HLT_15GeV_3050_b->Write("photon_w_HLT_15GeV_3050_b", TObject::kOverwrite);
photon_w_HLT_15GeV_5080_b->Write("photon_w_HLT_15GeV_5080_b", TObject::kOverwrite);

photon_w_HLT_25GeV_010_b->Write("photon_w_HLT_25GeV_010_b", TObject::kOverwrite);
photon_w_HLT_25GeV_1030_b->Write("photon_w_HLT_25GeV_1030_b", TObject::kOverwrite);
photon_w_HLT_25GeV_3050_b->Write("photon_w_HLT_25GeV_3050_b", TObject::kOverwrite);
photon_w_HLT_25GeV_5080_b->Write("photon_w_HLT_25GeV_5080_b", TObject::kOverwrite);

photon_w_HLT_35GeV_010_b->Write("photon_w_HLT_35GeV_010_b", TObject::kOverwrite);
photon_w_HLT_35GeV_1030_b->Write("photon_w_HLT_35GeV_1030_b", TObject::kOverwrite);
photon_w_HLT_35GeV_3050_b->Write("photon_w_HLT_35GeV_3050_b", TObject::kOverwrite);
photon_w_HLT_35GeV_5080_b->Write("photon_w_HLT_35GeV_5080_b", TObject::kOverwrite);
  

//----Endcaps
photon_w_HLT_15GeV_010_e->Write("photon_w_HLT_15GeV_010_e", TObject::kOverwrite);
photon_w_HLT_15GeV_1030_e->Write("photon_w_HLT_15GeV_1030_e", TObject::kOverwrite);
photon_w_HLT_15GeV_3050_e->Write("photon_w_HLT_15GeV_3050_e", TObject::kOverwrite);
photon_w_HLT_15GeV_5080_e->Write("photon_w_HLT_15GeV_5080_e", TObject::kOverwrite);

photon_w_HLT_25GeV_010_e->Write("photon_w_HLT_25GeV_010_e", TObject::kOverwrite);
photon_w_HLT_25GeV_1030_e->Write("photon_w_HLT_25GeV_1030_e", TObject::kOverwrite);
photon_w_HLT_25GeV_3050_e->Write("photon_w_HLT_25GeV_3050_e", TObject::kOverwrite);
photon_w_HLT_25GeV_5080_e->Write("photon_w_HLT_25GeV_5080_e", TObject::kOverwrite);

photon_w_HLT_35GeV_010_e->Write("photon_w_HLT_35GeV_010_e", TObject::kOverwrite);
photon_w_HLT_35GeV_1030_e->Write("photon_w_HLT_35GeV_1030_e", TObject::kOverwrite);
photon_w_HLT_35GeV_3050_e->Write("photon_w_HLT_35GeV_3050_e", TObject::kOverwrite);
photon_w_HLT_35GeV_5080_e->Write("photon_w_HLT_35GeV_5080_e", TObject::kOverwrite);
//**************************************************************************************//



  return 0;
}
