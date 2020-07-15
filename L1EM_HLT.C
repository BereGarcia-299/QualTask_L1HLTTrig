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

 float SumEt_data_fire = fcalC_et_data + fcalA_et_data;

 

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
    
    
    if ( ((state_4 && state_5) || (state_6 && state_7) || (state_8 && state_9)) == 1 ) {
      
      all_photons_phi->Fill(photon_phi->at(iPhoton));
      all_photons_eta->Fill(photon_eta->at(iPhoton));
      all_photons->Fill(photon_pt->at(iPhoton));

    } //--Encaps && Barrel
    if((state_8 && state_9) == 1){all_photons_barrel->Fill(photon_pt->at(iPhoton));}                 //-----Photons that ONLY went through the barrel
    if ( ((state_4 && state_5) || (state_6 && state_7)) == 1 ) {all_photons_endcaps->Fill(photon_pt->at(iPhoton));}  //------Photons went through the endcaps ONLY
    
     


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

  
  //-------Photon Matched W/ L1EM RoI Trigger  
    if(check_trigger_15){
    
    if ( ((req_4 && req_5) || (req_6 && req_7) || (req_8 && req_9)) == 1 ) {photons_wL1EM_15GeV->Fill(photon_pt->at(iPhoton));} //--Encaps && Barrel
    if((req_8 && req_9) == 1){photons_wL1EM_15GeV_barrel->Fill(photon_pt->at(iPhoton));}                 //-----Photons that ONLY went through the barrel
    if ( ((req_4 && req_5) || (state_6 && state_7)) == 1 ) {photons_wL1EM_15GeV_endcaps->Fill(photon_pt->at(iPhoton));}  //------Photons went through the endcaps ONLY

                                                                                             
    check_trigger_15 = false;

  }
  if(check_trigger_25){
    
    if ( ((req_4 && req_5) || (req_6 && req_7) || (req_8 && req_9)) == 1 ) {photons_wL1EM_25GeV->Fill(photon_pt->at(iPhoton));} //--Encaps && Barrel
    if((req_8 && req_9) == 1){photons_wL1EM_25GeV_barrel->Fill(photon_pt->at(iPhoton));}                 //-----Photons that ONLY went through the barrel
    if ( ((req_4 && req_5) || (state_6 && state_7)) == 1 ) {photons_wL1EM_25GeV_endcaps->Fill(photon_pt->at(iPhoton));}  //------Photons went through the endcaps ONLY



    check_trigger_25 = false;
  }
    if(check_trigger_35){
    if ( ((req_4 && req_5) || (req_6 && req_7) || (req_8 && req_9)) == 1 ) {photons_wL1EM_35GeV->Fill(photon_pt->at(iPhoton));} //--Encaps && Barrel
    if((req_8 && req_9) == 1){photons_wL1EM_35GeV_barrel->Fill(photon_pt->at(iPhoton));}                 //-----Photons that ONLY went through the barrel
    if ( ((req_4 && req_5) || (state_6 && state_7)) == 1 ) {photons_wL1EM_35GeV_endcaps->Fill(photon_pt->at(iPhoton));}  //------Photons went through the endcaps ONLY


   
    check_trigger_35 = false;
  }

  //-------Photon Matched W/ HLT Trigger                                                                                                                             
  if(check_trigger_hlt){
    
    if ( ((req_4 && req_5) || (req_6 && req_7) || (req_8 && req_9)) == 1 ) {photons_w_HLT_15->Fill(photon_pt->at(iPhoton));} //--Encaps && Barrel
    if((req_8 && req_9) == 1){photons_w_HLT_15_barrel->Fill(photon_pt->at(iPhoton));}                 //-----Photons that ONLY went through the barrel
    if ( ((req_4 && req_5) || (state_6 && state_7)) == 1 ) {photons_w_HLT_15_endcaps->Fill(photon_pt->at(iPhoton));}  //------Photons went through the endcaps ONLY

                                                                                             
    check_trigger_hlt = false;

  }
  if(check_trigger_hlt_25){
    
    if ( ((req_4 && req_5) || (req_6 && req_7) || (req_8 && req_9)) == 1 ) {photons_w_HLT_25->Fill(photon_pt->at(iPhoton));} //--Encaps && Barrel
    if((req_8 && req_9) == 1){photons_w_HLT_25_barrel->Fill(photon_pt->at(iPhoton));}                 //-----Photons that ONLY went through the barrel
    if ( ((req_4 && req_5) || (state_6 && state_7)) == 1 ) {photons_w_HLT_25_endcaps->Fill(photon_pt->at(iPhoton));}  //------Photons went through the endcaps ONLY



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
    if ( ((req_4 && req_5) || (req_6 && req_7) || (req_8 && req_9)) == 1 ) {photons_w_HLT_35->Fill(photon_pt->at(iPhoton));} //--Encaps && Barrel
    if((req_8 && req_9) == 1){photons_w_HLT_35_barrel->Fill(photon_pt->at(iPhoton));}                 //-----Photons that ONLY went through the barrel
    if ( ((req_4 && req_5) || (state_6 && state_7)) == 1 ) {photons_w_HLT_35_endcaps->Fill(photon_pt->at(iPhoton));}  //------Photons went through the endcaps ONLY




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

//----------------------------------------//
//---------------Histograms---------------//
//----------------------------------------//

//-----L1 pT Distribution------//
TCanvas *l1_pt = new TCanvas("l1_pt","l1_pt",600,500); 
TLegend *l1EMRoI_lg = new TLegend(0.6,0.5,0.85,0.6);

photons_wL1EM_15GeV->SetTitle("L1EM RoI p_{T} Distribution");
photons_wL1EM_15GeV->GetXaxis()->SetTitle("p_{T} [GeV]");
photons_wL1EM_15GeV->GetYaxis()->SetTitle("Counts");
photons_wL1EM_15GeV->SetMarkerStyle(20);
photons_wL1EM_15GeV->SetLineColor(kRed);
photons_wL1EM_15GeV->SetMarkerColor(kRed);
photons_wL1EM_15GeV->SetName("photons_wL1EM_15GeV");

photons_wL1EM_25GeV->SetMarkerStyle(22);
photons_wL1EM_25GeV->SetLineColor(kBlue);
photons_wL1EM_25GeV->SetMarkerColor(kBlue);
photons_wL1EM_25GeV->SetName("photons_wL1EM_25GeV");

photons_wL1EM_35GeV->SetMarkerStyle(21);
photons_wL1EM_35GeV->SetMarkerColor(kMagenta);
photons_wL1EM_35GeV->SetLineColor(kMagenta);
photons_wL1EM_35GeV->SetName("photons_wL1EM_35GeV");

all_photons->SetMarkerStyle(33);
all_photons->SetLineColor(kBlack);
all_photons->SetMarkerColor(kBlack);
all_photons->SetName("all_photons");

photons_wL1EM_15GeV->Draw();
photons_wL1EM_35GeV->Draw("same");

//all_photons->Draw("same HIST E1 P");
photons_wL1EM_25GeV->Draw("same");


//l1EMRoI_lg->AddEntry("all_photons","Offline Photons pT > 15 GeV", "pe");
l1EMRoI_lg->AddEntry("photons_wL1EM_15GeV","L1EM RoI > 15 GeV","l");
l1EMRoI_lg->AddEntry("photons_wL1EM_25GeV","L1EM RoI > 25 GeV", "l");
l1EMRoI_lg->AddEntry("photons_wL1EM_35GeV","L1EM RoI > 35 GeV", "l");
l1EMRoI_lg->SetTextSize(0.03);

TPaveText *l1_pT_text= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1_pT_text->SetTextSize(0.036);
l1_pT_text->SetFillColor(0);
l1_pT_text->SetBorderSize(0);
l1_pT_text->SetShadowColor(0);
l1_pT_text->SetTextAlign(21);
l1_pT_text->AddText("");
l1_pT_text->AddText("Pb+Pb");
l1_pT_text->AddText("Tight Photons");
l1_pT_text->AddText("Offline p_{T}^{#gamma} > 15 GeV");
l1_pT_text->AddText("DeltaR < 0.15");

l1_pT_text->Draw("same");

l1EMRoI_lg->Draw("same");

//----L1 pT Distribution
//----Barrel
TCanvas *l1_pt_barrel = new TCanvas("l1_pt_barrel","l1_pt_barrel",600,500); 
TLegend *l1EMRoI_lg_barrel = new TLegend(0.6,0.5,0.85,0.6);

photons_wL1EM_15GeV_barrel->SetTitle("L1EM RoI p_{T} Distribution");
photons_wL1EM_15GeV_barrel->GetXaxis()->SetTitle("p_{T} [GeV]");
photons_wL1EM_15GeV_barrel->GetYaxis()->SetTitle("Counts");
photons_wL1EM_15GeV_barrel->SetMarkerStyle(20);
photons_wL1EM_15GeV_barrel->SetLineColor(kRed);
photons_wL1EM_15GeV_barrel->SetMarkerColor(kRed);
photons_wL1EM_15GeV_barrel->SetName("photons_wL1EM_15GeV_barrel");

photons_wL1EM_25GeV_barrel->SetMarkerStyle(22);
photons_wL1EM_25GeV_barrel->SetLineColor(kBlue);
photons_wL1EM_25GeV_barrel->SetMarkerColor(kBlue);
photons_wL1EM_25GeV_barrel->SetName("photons_wL1EM_25GeV_barrel");

photons_wL1EM_35GeV_barrel->SetMarkerStyle(21);
photons_wL1EM_35GeV_barrel->SetMarkerColor(kMagenta);
photons_wL1EM_35GeV_barrel->SetLineColor(kMagenta);
photons_wL1EM_35GeV_barrel->SetName("photons_wL1EM_35GeV_barrel");

photons_wL1EM_15GeV_barrel->Draw();
photons_wL1EM_35GeV_barrel->Draw("same");


photons_wL1EM_25GeV_barrel->Draw("same");


//l1EMRoI_lg->AddEntry("all_photons","Offline Photons pT > 15 GeV", "pe");
l1EMRoI_lg_barrel->AddEntry("photons_wL1EM_15GeV_barrel","L1EM RoI > 15 GeV","l");
l1EMRoI_lg_barrel->AddEntry("photons_wL1EM_25GeV_barrel","L1EM RoI > 25 GeV", "l");
l1EMRoI_lg_barrel->AddEntry("photons_wL1EM_35GeV_barrel","L1EM RoI > 35 GeV", "l");
l1EMRoI_lg_barrel->SetTextSize(0.03);


TPaveText *l1_pT_text_barrel = new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1_pT_text_barrel->SetTextSize(0.036);
l1_pT_text_barrel->SetFillColor(0);
l1_pT_text_barrel->SetBorderSize(0);
l1_pT_text_barrel->SetShadowColor(0);
l1_pT_text_barrel->SetTextAlign(21);
l1_pT_text_barrel->AddText("");
l1_pT_text_barrel->AddText("Pb+Pb");
l1_pT_text_barrel->AddText("Tight Photons");
l1_pT_text_barrel->AddText("Offline p_{T}^{#gamma} > 15 GeV");
l1_pT_text_barrel->AddText("DeltaR < 0.15");
l1_pT_text_barrel->AddText("|#eta| < 1.37");

l1_pT_text_barrel->Draw("same");

l1EMRoI_lg_barrel->Draw("same");


//----Barrel
TCanvas *l1_pt_endcaps = new TCanvas("l1_pt_endcaps","l1_pt_endcaps",600,500); 
TLegend *l1EMRoI_lg_endcaps = new TLegend(0.6,0.5,0.85,0.6);

photons_wL1EM_15GeV_endcaps->SetTitle("L1EM RoI p_{T} Distribution");
photons_wL1EM_15GeV_endcaps->GetXaxis()->SetTitle("p_{T} [GeV]");
photons_wL1EM_15GeV_endcaps->GetYaxis()->SetTitle("Counts");
photons_wL1EM_15GeV_endcaps->SetMarkerStyle(20);
photons_wL1EM_15GeV_endcaps->SetLineColor(kRed);
photons_wL1EM_15GeV_endcaps->SetMarkerColor(kRed);
photons_wL1EM_15GeV_endcaps->SetName("photons_wL1EM_15GeV_endcaps");

photons_wL1EM_25GeV_endcaps->SetMarkerStyle(22);
photons_wL1EM_25GeV_endcaps->SetLineColor(kBlue);
photons_wL1EM_25GeV_endcaps->SetMarkerColor(kBlue);
photons_wL1EM_25GeV_endcaps->SetName("photons_wL1EM_25GeV_endcaps");

photons_wL1EM_35GeV_endcaps->SetMarkerStyle(21);
photons_wL1EM_35GeV_endcaps->SetMarkerColor(kMagenta);
photons_wL1EM_35GeV_endcaps->SetLineColor(kMagenta);
photons_wL1EM_35GeV_endcaps->SetName("photons_wL1EM_35GeV_endcaps");

photons_wL1EM_15GeV_endcaps->Draw();
photons_wL1EM_35GeV_endcaps->Draw("same");


photons_wL1EM_25GeV_endcaps->Draw("same");


//l1EMRoI_lg->AddEntry("all_photons","Offline Photons pT > 15 GeV", "pe");
l1EMRoI_lg_endcaps->AddEntry("photons_wL1EM_15GeV_endcaps","L1EM RoI > 15 GeV","l");
l1EMRoI_lg_endcaps->AddEntry("photons_wL1EM_25GeV_endcaps","L1EM RoI > 25 GeV", "l");
l1EMRoI_lg_endcaps->AddEntry("photons_wL1EM_35GeV_endcaps","L1EM RoI > 35 GeV", "l");
l1EMRoI_lg_endcaps->SetTextSize(0.03);


TPaveText *l1_pT_text_endcaps= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1_pT_text_endcaps->SetTextSize(0.036);
l1_pT_text_endcaps->SetFillColor(0);
l1_pT_text_endcaps->SetBorderSize(0);
l1_pT_text_endcaps->SetShadowColor(0);
l1_pT_text_endcaps->SetTextAlign(21);
l1_pT_text_endcaps->AddText("");
l1_pT_text_endcaps->AddText("Pb+Pb");
l1_pT_text_endcaps->AddText("Tight Photons");
l1_pT_text_endcaps->AddText("Offline p_{T}^{#gamma} > 15 GeV");
l1_pT_text_endcaps->AddText("DeltaR < 0.15");
l1_pT_text_endcaps->AddText("1.56 < |#eta| < 2.37");

l1_pT_text_endcaps->Draw("same");

l1EMRoI_lg_endcaps->SetBorderSize(0);
l1EMRoI_lg_endcaps->Draw("same");



//-----------------------------------------------//
//-------Photons TRough Barrel && Encap----------//
//-----------------------------------------------//
TCanvas *new_l1em_eff = new TCanvas("new_l1em_eff","new_l1em_eff",600,500);
TLegend *new_lg_l1em = new TLegend(0.6,0.5,0.85,0.6);

TH1D* L1EM_eff = new TH1D("L1EM_eff","",40,0,200);
L1EM_eff->SetName("L1EM_eff");                                                                                                                           
L1EM_eff->GetXaxis()->SetTitle("p_{T}");
L1EM_eff->SetMarkerStyle(20);
L1EM_eff->SetTitle("L1EM Efficiency");
L1EM_eff->SetLineColor(kRed);
L1EM_eff->SetMarkerColor(kRed);
L1EM_eff->Divide(photons_wL1EM_15GeV,all_photons,1,1,"B");
L1EM_eff->SetMaximum(1.1);

TH1D *L1EM_eff_25 = new TH1D("L1EM_eff_25","",40,0,200);
L1EM_eff_25->SetName("L1EM_eff_25");
L1EM_eff_25->SetMarkerStyle(22);
L1EM_eff_25->SetMarkerColor(kBlue);
L1EM_eff_25->SetLineColor(kBlue);
L1EM_eff_25->Divide(photons_wL1EM_25GeV,all_photons,1,1,"B");

TH1D *L1EM_eff_35 = new TH1D("L1EM_eff_35","",40,0,200);
L1EM_eff_35->SetName("L1EM_eff_35");
L1EM_eff_35->SetMarkerStyle(21);
L1EM_eff_35->SetMarkerColor(kMagenta);
L1EM_eff_35->SetLineColor(kMagenta);
L1EM_eff_35->Divide(photons_wL1EM_35GeV,all_photons,1,1,"B");

TPaveText *l1em_photon_txt= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_txt->SetTextSize(0.03);
l1em_photon_txt->SetFillColor(0);
l1em_photon_txt->SetBorderSize(0);
l1em_photon_txt->SetShadowColor(0);
l1em_photon_txt->SetTextAlign(21);
l1em_photon_txt->AddText("");
l1em_photon_txt->AddText("Pb+Pb");
l1em_photon_txt->AddText("Tight Photons");
l1em_photon_txt->AddText("p^{#gamma}_{T} > 15 GeV");
//phi_photon_2->AddText("HLT pT > 15 GeV");
l1em_photon_txt->AddText("DeltaR < 0.15");
//phi_photon_2->AddText("1.56 < |#eta| < 2.37");




L1EM_eff->Draw();
L1EM_eff_35->Draw("same");
L1EM_eff_25->Draw("same");
l1em_photon_txt->Draw("same");

new_lg_l1em->AddEntry("L1EM_eff","L1EM RoI > 15 GeV","pe");
new_lg_l1em->AddEntry("L1EM_eff_25","L1EM RoI > 25 GeV","pe");
new_lg_l1em->AddEntry("L1EM_eff_35","L1EM RoI > 35 GeV","pe");
new_lg_l1em->Draw("same");

//----L1EM Efficiency in Barrel
TCanvas *new_l1em_eff_barrel = new TCanvas("new_l1em_eff_barrel","new_l1em_eff_barrel",600,500);
TLegend *new_lg_l1em_barrel = new TLegend(0.6,0.5,0.85,0.6);

TH1D* L1EM_eff_barrel = new TH1D("L1EM_eff_barrel","",40,0,200);
L1EM_eff_barrel->SetName("L1EM_eff_barrel");                                                                                                                           
L1EM_eff_barrel->GetXaxis()->SetTitle("p_{T}");
L1EM_eff_barrel->SetMarkerStyle(20);
L1EM_eff_barrel->SetTitle("L1EM Efficiency");
L1EM_eff_barrel->SetLineColor(kRed);
L1EM_eff_barrel->SetMarkerColor(kRed);
L1EM_eff_barrel->Divide(photons_wL1EM_15GeV_barrel,all_photons_barrel,1,1,"B");
L1EM_eff_barrel->SetMaximum(1.1);

TH1D *L1EM_eff_25_barrel = new TH1D("L1EM_eff_25_barrel","",40,0,200);
L1EM_eff_25_barrel->SetName("L1EM_eff_25_barrel");
L1EM_eff_25_barrel->SetMarkerStyle(22);
L1EM_eff_25_barrel->SetMarkerColor(kBlue);
L1EM_eff_25_barrel->SetLineColor(kBlue);
L1EM_eff_25_barrel->Divide(photons_wL1EM_25GeV_barrel,all_photons_barrel,1,1,"B");

TH1D *L1EM_eff_35_barrel = new TH1D("L1EM_eff_35_barrel","",40,0,200);
L1EM_eff_35_barrel->SetName("L1EM_eff_35_barrel");
L1EM_eff_35_barrel->SetMarkerStyle(21);
L1EM_eff_35_barrel->SetMarkerColor(kMagenta);
L1EM_eff_35_barrel->SetLineColor(kMagenta);
L1EM_eff_35_barrel->Divide(photons_wL1EM_35GeV_barrel,all_photons_barrel,1,1,"B");

TPaveText *l1em_photon_txt_barrel= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_txt_barrel->SetTextSize(0.03);
l1em_photon_txt_barrel->SetFillColor(0);
l1em_photon_txt_barrel->SetBorderSize(0);
l1em_photon_txt_barrel->SetShadowColor(0);
l1em_photon_txt_barrel->SetTextAlign(21);
l1em_photon_txt_barrel->AddText("");
l1em_photon_txt_barrel->AddText("Pb+Pb");
l1em_photon_txt_barrel->AddText("Tight Photons");
l1em_photon_txt_barrel->AddText("p^{#gamma}_{T} > 15 GeV");
//phi_photon_2->AddText("HLT pT > 15 GeV");
l1em_photon_txt_barrel->AddText("DeltaR < 0.15");
l1em_photon_txt_barrel->AddText("|#eta| < 1.37");




L1EM_eff_barrel->Draw();
L1EM_eff_35_barrel->Draw("same");
L1EM_eff_25_barrel->Draw("same");
l1em_photon_txt_barrel->Draw("same");

new_lg_l1em_barrel->AddEntry("L1EM_eff_barrel","L1EM RoI > 15 GeV","pe");
new_lg_l1em_barrel->AddEntry("L1EM_eff_25_barrel","L1EM RoI > 25 GeV","pe");
new_lg_l1em_barrel->AddEntry("L1EM_eff_35_barrel","L1EM RoI > 35 GeV","pe");
new_lg_l1em_barrel->Draw("same");

//----L1EM Efficiency in Endcaps
TCanvas *new_l1em_eff_endcaps = new TCanvas("new_l1em_eff_endcaps","new_l1em_eff_endcaps",600,500);
TLegend *new_lg_l1em_endcaps = new TLegend(0.6,0.5,0.85,0.6);

TH1D* L1EM_eff_endcaps = new TH1D("L1EM_eff_endcaps","",40,0,200);
L1EM_eff_endcaps->SetName("L1EM_eff_endcaps");                                                                                                                           
L1EM_eff_endcaps->GetXaxis()->SetTitle("p_{T}");
L1EM_eff_endcaps->SetMarkerStyle(20);
L1EM_eff_endcaps->SetTitle("L1EM Efficiency");
L1EM_eff_endcaps->SetLineColor(kRed);
L1EM_eff_endcaps->SetMarkerColor(kRed);
L1EM_eff_endcaps->Divide(photons_wL1EM_15GeV_endcaps,all_photons_endcaps,1,1,"B");
L1EM_eff_endcaps->SetMaximum(1.1);

TH1D *L1EM_eff_25_endcaps = new TH1D("L1EM_eff_25_endcaps","",40,0,200);
L1EM_eff_25_endcaps->SetName("L1EM_eff_25_endcaps");
L1EM_eff_25_endcaps->SetMarkerStyle(22);
L1EM_eff_25_endcaps->SetMarkerColor(kBlue);
L1EM_eff_25_endcaps->SetLineColor(kBlue);
L1EM_eff_25_endcaps->Divide(photons_wL1EM_25GeV_endcaps,all_photons_endcaps,1,1,"B");

TH1D *L1EM_eff_35_endcaps = new TH1D("L1EM_eff_35_endcaps","",40,0,200);
L1EM_eff_35_endcaps->SetName("L1EM_eff_35_endcaps");
L1EM_eff_35_endcaps->SetMarkerStyle(21);
L1EM_eff_35_endcaps->SetMarkerColor(kMagenta);
L1EM_eff_35_endcaps->SetLineColor(kMagenta);
L1EM_eff_35_endcaps->Divide(photons_wL1EM_35GeV_endcaps,all_photons_endcaps,1,1,"B");

TPaveText *l1em_photon_txt_endcaps= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
l1em_photon_txt_endcaps->SetTextSize(0.03);
l1em_photon_txt_endcaps->SetFillColor(0);
l1em_photon_txt_endcaps->SetBorderSize(0);
l1em_photon_txt_endcaps->SetShadowColor(0);
l1em_photon_txt_endcaps->SetTextAlign(21);
l1em_photon_txt_endcaps->AddText("");
l1em_photon_txt_endcaps->AddText("Pb+Pb");
l1em_photon_txt_endcaps->AddText("Tight Photons");
l1em_photon_txt_endcaps->AddText("p^{#gamma}_{T} > 15 GeV");
//phi_photon_2->AddText("HLT pT > 15 GeV");
l1em_photon_txt_endcaps->AddText("DeltaR < 0.15");
l1em_photon_txt_endcaps->AddText("1.56 < |#eta| < 2.37");




L1EM_eff_endcaps->Draw();
L1EM_eff_35_endcaps->Draw("same");
L1EM_eff_25_endcaps->Draw("same");
l1em_photon_txt_endcaps->Draw("same");

new_lg_l1em_endcaps->AddEntry("L1EM_eff_endcaps","L1EM RoI > 15 GeV","pe");
new_lg_l1em_endcaps->AddEntry("L1EM_eff_25_endcaps","L1EM RoI > 25 GeV","pe");
new_lg_l1em_endcaps->AddEntry("L1EM_eff_35_endcaps","L1EM RoI > 35 GeV","pe");
new_lg_l1em_endcaps->Draw("same");






//-----HLT pT Distribution------//
TCanvas *hlt_pT = new TCanvas("hlt_pT","hlt_pT",600,500); 

hlt_pT_dis->SetTitle("HLT p_{T} Distribution");
hlt_pT_dis->GetXaxis()->SetTitle("p_{T}");
hlt_pT_dis->GetYaxis()->SetTitle("Counts");

hlt_pT_dis->Draw();

TPaveText *hlt_pT_text= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_pT_text->SetTextSize(0.03);
hlt_pT_text->SetFillColor(0);
hlt_pT_text->SetBorderSize(0);
hlt_pT_text->SetShadowColor(0);
hlt_pT_text->SetTextAlign(21);
hlt_pT_text->AddText("");
hlt_pT_text->AddText("Pb+Pb");
hlt_pT_text->AddText("Loose Photons");
hlt_pT_text->AddText("p_{T}^{#gamma} > 15 GeV");

hlt_pT_text->Draw("same");

//-----HLT Phi Distribution-----//

TCanvas *hlt_phi = new TCanvas("hlt_phi","hlt_phi",600,500);

hlt_phi_dis->SetTitle("HLT Phi Distribution");
hlt_phi_dis->GetXaxis()->SetTitle("#phi");
hlt_phi_dis->GetYaxis()->SetTitle("Counts");

hlt_phi_dis->Draw();

TPaveText *hlt_phi_text= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_phi_text->SetTextSize(0.03);
hlt_phi_text->SetFillColor(0);
hlt_phi_text->SetBorderSize(0);
hlt_phi_text->SetShadowColor(0);
hlt_phi_text->SetTextAlign(21);
hlt_phi_text->AddText("");
hlt_phi_text->AddText("Pb+Pb");
hlt_phi_text->AddText("Loose Photons");
hlt_phi_text->AddText("p_{T}^{#gamma} > 15 GeV");

hlt_phi_text->Draw("same");

//---HLT Eta Distribution----//
TCanvas *hlt_eta = new TCanvas("hlt_eta","hlt_eta",600,500);

hlt_eta_dis->SetTitle("HLT #eta Distribution");
hlt_eta_dis->GetXaxis()->SetTitle("#eta");
hlt_eta_dis->GetYaxis()->SetTitle("Counts");

hlt_eta_dis->Draw();

TPaveText *hlt_eta_text= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_eta_text->SetTextSize(0.03);
hlt_eta_text->SetFillColor(0);
hlt_eta_text->SetBorderSize(0);
hlt_eta_text->SetShadowColor(0);
hlt_eta_text->SetTextAlign(21);
hlt_eta_text->AddText("");
hlt_eta_text->AddText("Pb+Pb");
hlt_eta_text->AddText("Loose Photons");
hlt_eta_text->AddText("p_{T}^{#gamma} > 15 GeV");

hlt_eta_text->Draw();






//-----HLT
TCanvas *new_pT_HLT_0 = new TCanvas("new_pT_HLT_0","new_pT_HLT_0",600,500);
TLegend *new_pT_lg_hlt = new TLegend(0.6,0.5,0.85,0.6);

all_photons->SetTitle("Photon p_{T}");
all_photons->SetMarkerStyle(20);
all_photons->SetMarkerColor(kBlack);
all_photons->SetLineColor(kBlack);
all_photons->GetYaxis()->SetTitle("Counts");
all_photons->GetXaxis()->SetTitle("p_{T}");
all_photons->SetName("all_photons");

photons_w_HLT_15->SetMarkerStyle(20);
photons_w_HLT_15->SetLineColor(kRed);
photons_w_HLT_15->SetMarkerColor(kRed);
photons_w_HLT_15->SetName("photons_w_HLT_15");

photons_w_HLT_25->SetMarkerColor(kBlue);
photons_w_HLT_25->SetLineColor(kBlue);
photons_w_HLT_25->SetName("photons_w_HLT_25");


photons_w_HLT_35->SetLineColor(kMagenta);
photons_w_HLT_35->SetMarkerColor(kMagenta);
photons_w_HLT_35->SetName("photons_w_HLT_35");

all_photons->Draw();
photons_w_HLT_15->Draw("same");
photons_w_HLT_35->Draw("same");
photons_w_HLT_25->Draw("same");

TPaveText *pt_photon_1= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
pt_photon_1->SetTextSize(0.03);
pt_photon_1->SetFillColor(0);
pt_photon_1->SetBorderSize(0);
pt_photon_1->SetShadowColor(0);
pt_photon_1->SetTextAlign(21);
pt_photon_1->AddText("");
pt_photon_1->AddText("Pb+Pb");
pt_photon_1->AddText("Tight Photons");
pt_photon_1->AddText("HLT Loose Photons");
pt_photon_1->AddText("p^{#gamma}_{T} > 20 GeV");
//pt_photon_1->AddText("1.56 < |#eta| < 2.37");
pt_photon_1->AddText("dR < 0.15");

pt_photon_1->Draw("same");

new_pT_lg_hlt->AddEntry("all_photons","All Offline Photons","l");
new_pT_lg_hlt->AddEntry("photons_w_HLT_15", "HLT Photons p_{T} > 15 GeV","l");
new_pT_lg_hlt->AddEntry("photons_w_HLT_25","HLT Photons p_{T} > 25 GeV", "l");
new_pT_lg_hlt->AddEntry("photons_w_HLT_35","HLT Photons p_{T} > 35 GeV", "l");

new_pT_lg_hlt->Draw("same");


//photons_wL1EM->Draw("same");
//pt_photon_0->Draw("same");                            

//-------Photons Failed To Fire a HLT Trigger 
//-------Photon pT Distribution
TCanvas *new_photon_nofire_hlt_pt = new TCanvas("new_photon_nofire_hlt_pt","new_photon_nofire_hlt_pt",600,500);


photons_nofire_HLT_pt->SetTitle("Photon p_{T} Distribution");
photons_nofire_HLT_pt->SetMarkerStyle(20);
photons_nofire_HLT_pt->SetMarkerColor(kBlack);
photons_nofire_HLT_pt->SetLineColor(kBlack);
photons_nofire_HLT_pt->GetXaxis()->SetTitle("p_{T}");
photons_nofire_HLT_pt->GetYaxis()->SetTitle("Counts");

photons_nofire_HLT_pt->Draw();


//------Photon Eta Distribution
TCanvas *new_photon_nofire_hlt_eta = new TCanvas("new_photon_nofire_hlt_eta","new_photon_nofire_hlt_eta",600,500);

photons_nofire_HLT_eta->SetTitle("Photon #eta Distribution");
photons_nofire_HLT_eta->SetMarkerStyle(20);
photons_nofire_HLT_eta->SetMarkerColor(kBlack);
photons_nofire_HLT_eta->SetLineColor(kBlack);
photons_nofire_HLT_eta->GetXaxis()->SetTitle("#eta");
photons_nofire_HLT_eta->GetYaxis()->SetTitle("Counts");


photons_nofire_HLT_eta->Draw();

//-----Photon sumET Distribution 
TCanvas *new_photon_nofire_hlt_sumET = new TCanvas("new_photon_nofire_hlt_sumET","new_photon_nofire_hlt_sumET",600,500);

photons_nofire_HLT_sumET->SetTitle("#sum ET");
photons_nofire_HLT_sumET->SetMarkerColor(kBlack);
photons_nofire_HLT_sumET->SetLineColor(kBlack);
photons_nofire_HLT_sumET->SetMarkerStyle(20);
photons_nofire_HLT_sumET->GetYaxis()->SetTitle("Counts");
photons_nofire_HLT_sumET->GetXaxis()->SetTitle("#sum ET");

photons_nofire_HLT_sumET->Draw();



//-------eta
TCanvas *new_eta_0 = new TCanvas("new_eta_0","new_eta_0",600,500);
TLegend *new_eta_lg = new TLegend(0.6,0.5,0.85,0.6);

/**
all_photons_eta->SetTitle("Photon #eta");
all_photons_eta->SetMarkerStyle(20);
all_photons_eta->SetMarkerColor(kBlack);
all_photons_eta->SetLineColor(kBlack); 
all_photons_eta->GetYaxis()->SetTitle("Counts");
all_photons_eta->GetXaxis()->SetTitle("#eta");
all_photons_eta->SetName("all_photons_eta");
*/
/**
photons_wL1EM_eta->SetTitle("L1EM RoI #eta");
photons_wL1EM_eta->SetMarkerStyle(20);
photons_wL1EM_eta->SetMarkerColor(kRed);
photons_wL1EM_eta->SetLineColor(kRed);
photons_wL1EM_eta->SetName("photons_wL1EM_eta");
photons_wL1EM_eta->GetYaxis()->SetTitle("Counts");
photons_wL1EM_eta->GetXaxis()->SetTitle("#eta");
*/
TPaveText *eta_photon_0= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
eta_photon_0->SetTextSize(0.03);
eta_photon_0->SetFillColor(0);
eta_photon_0->SetBorderSize(0);
eta_photon_0->SetShadowColor(0);
eta_photon_0->SetTextAlign(21);
eta_photon_0->AddText("");
eta_photon_0->AddText("Pb+Pb");
eta_photon_0->AddText("Tight Photons");
eta_photon_0->AddText("p^{#gamma}_{T} > 20 GeV");
eta_photon_0->AddText("L1EM pT > 15");
eta_photon_0->AddText("DeltaR < 0.15");
eta_photon_0->AddText("1.56 < |#eta| < 2.37");
eta_photon_0->AddText("|#eta| < 1.37");


//all_photons_eta->Draw();
//photons_wL1EM_eta->Draw("same");
eta_photon_0->Draw("same");

new_eta_lg->AddEntry("all_photons_eta","All Offline Photons","pe");
//new_eta_lg->AddEntry("photons_wL1EM_eta", "Photons MAtched To An L1EM Triggers","pe");

//new_eta_lg->Draw("same");

//------phi
TCanvas *new_phi_0 = new TCanvas("new_phi_0","new_phi_0",600,500);
TLegend *new_phi_lg = new TLegend(0.6,0.5,0.85,0.6);
/**
all_photons_phi->SetTitle("Photon #phi");
all_photons_phi->SetMarkerStyle(20);
all_photons_phi->SetMarkerColor(kBlack);
all_photons_phi->SetLineColor(kBlack); 
all_photons_phi->GetYaxis()->SetTitle("Counts");
all_photons_phi->GetXaxis()->SetTitle("#phi");
all_photons_phi->SetName("all_photons_phi");
*/
/**
photons_wL1EM_phi->SetMarkerStyle(20);
photons_wL1EM_phi->SetMarkerColor(kRed);
photons_wL1EM_phi->SetLineColor(kRed);
photons_wL1EM_phi->SetName("photons_wL1EM_phi");
photons_wL1EM_phi->GetYaxis()->SetTitle("Counts");
photons_wL1EM_phi->GetXaxis()->SetTitle("#phi");
*/
TPaveText *phi_photon_0= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
phi_photon_0->SetTextSize(0.03);
phi_photon_0->SetFillColor(0);
phi_photon_0->SetBorderSize(0);
phi_photon_0->SetShadowColor(0);
phi_photon_0->SetTextAlign(21);
phi_photon_0->AddText("");
phi_photon_0->AddText("Pb+Pb");
phi_photon_0->AddText("Tight Photons");
phi_photon_0->AddText("p^{#gamma}_{T} > 15 GeV");
phi_photon_0->AddText(" 1.56 < |#eta| < 2.37");
phi_photon_0->AddText("|#eta| < 1.37");
phi_photon_0->AddText("L1EM RoI pT > 15 GeV");
phi_photon_0->AddText("DeltaR < 0.15");

//all_photons_phi->Draw();
//photons_wL1EM_phi->Draw("same");
//phi_photon_0->Draw("same");

new_phi_lg->AddEntry("all_photons_phi","All Offline Photons","pe");
//new_phi_lg->AddEntry("photons_wL1EM_phi", "Photons MAtched To An L1EM Triggers","pe");

//new_phi_lg->Draw("same");



//--------------------------//
//--------Ratio Plot--------//
//--------------------------//
/**
TCanvas *new_ratio_1 = new TCanvas("new_ratio_1","new_ratio_1",600,500);
TLegend *new_diff = new TLegend(0.6,0.5,0.85,0.6);

TH1D* L1EM_eff = new TH1D("L1EM_eff","",40,0,200);
L1EM_eff->GetYaxis()->SetRangeUser(0,2);
L1EM_eff->SetTitle("HLT Efficiency");
L1EM_eff->SetName("L1EM_eff");
L1EM_eff->GetYaxis()->SetRangeUser(0,2);                                                                                                                                     
L1EM_eff->GetXaxis()->SetTitle("p_{T}");
L1EM_eff->SetMarkerStyle(21);
L1EM_eff->SetLineColor(kBlack);
L1EM_eff->SetMarkerColor(kBlack);
L1EM_eff->Divide(photons_wL1EM,all_photons,1,1,"B");

//L1EM_eff->Draw();
*/
//-----------------------------------------------//
//-------Photons TRough Barrel && Encap----------//
//-----------------------------------------------//
TH1D* HLT_eff = new TH1D("HLT_eff","",40,0,200);
HLT_eff->SetName("HLT_eff");                                                                                                                           
HLT_eff->GetXaxis()->SetTitle("p_{T}");
HLT_eff->SetMarkerStyle(20);
HLT_eff->SetTitle("HLT Efficiency");
HLT_eff->SetLineColor(kRed);
HLT_eff->SetMarkerColor(kRed);
HLT_eff->Divide(photons_w_HLT_15,all_photons,1,1,"B");
HLT_eff->SetMaximum(1.1);

TH1D *HLT_eff_25 = new TH1D("HLT_eff","",40,0,200);
HLT_eff_25->SetName("HLT_eff_25");
HLT_eff_25->SetMarkerStyle(22);
HLT_eff_25->SetMarkerColor(kBlue);
HLT_eff_25->SetLineColor(kBlue);
HLT_eff_25->Divide(photons_w_HLT_25,all_photons,1,1,"B");

TH1D *HLT_eff_35 = new TH1D("HLT_eff_35","",40,0,200);
HLT_eff_35->SetName("HLT_eff_35");
HLT_eff_35->SetMarkerStyle(21);
HLT_eff_35->SetMarkerColor(kMagenta);
HLT_eff_35->SetLineColor(kMagenta);
HLT_eff_35->Divide(photons_w_HLT_35,all_photons,1,1,"B");

TPaveText *hlt_photon_txt= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_txt->SetTextSize(0.03);
hlt_photon_txt->SetFillColor(0);
hlt_photon_txt->SetBorderSize(0);
hlt_photon_txt->SetShadowColor(0);
hlt_photon_txt->SetTextAlign(21);
hlt_photon_txt->AddText("");
hlt_photon_txt->AddText("Pb+Pb");
hlt_photon_txt->AddText("Tight Photons");
hlt_photon_txt->AddText("HLT Loose Photons");
hlt_photon_txt->AddText("p^{#gamma}_{T} > 15 GeV");
//phi_photon_2->AddText("HLT pT > 15 GeV");
hlt_photon_txt->AddText("dR < 0.15");
//phi_photon_2->AddText("1.56 < |#eta| < 2.37");

HLT_eff->Draw();
HLT_eff_35->Draw("same");
HLT_eff_25->Draw("same");


//new_diff->AddEntry("L1EM_eff", "Photons fired L1 Triggers/All Photons", "pe");
//new_diff->AddEntry("HLT_eff", "HLT Photons p_{T} > 15 GeV", "pe");
//new_diff->AddEntry("HLT_eff_25","HLT Photons p_{T} > 25 GeV","pe");
//new_diff->AddEntry("HLT_eff_35","HLT Photons p_{T} > 35 GeV", "pe");
//new_diff->SetBorderSize(0);
//new_diff->Draw("same");
hlt_photon_txt->Draw("same");



//--------pT Distributions
TCanvas *new_ratio_pT = new TCanvas("new_ratio_pT","new_ratio_pT",600,500);
TLegend *new_diff_pT = new TLegend(0.6,0.5,0.85,0.6);

photons_w_HLT_15->SetTitle("Photons p_{T} Distribution");
photons_w_HLT_15->GetXaxis()->SetTitle("p_{T}");
photons_w_HLT_15->GetYaxis()->SetTitle("1/N");
photons_w_HLT_15->SetMarkerColor(kRed);
photons_w_HLT_15->SetLineColor(kRed);
photons_w_HLT_15->SetMarkerStyle(20);
//photons_w_HLT_15->Scale(1/photons_w_HLT_15->Integral());

photons_w_HLT_25->SetMarkerColor(kBlue);
photons_w_HLT_25->SetLineColor(kBlue);
photons_w_HLT_25->SetMarkerStyle(22);
//photons_w_HLT_25->Scale(1/photons_w_HLT_25->Integral());

photons_w_HLT_35->SetMarkerStyle(21);
photons_w_HLT_35->SetMarkerColor(kMagenta);
photons_w_HLT_35->SetLineColor(kMagenta);
//photons_w_HLT_35->Scale(1/photons_w_HLT_35->Integral());

photons_w_HLT_15->Draw();
photons_w_HLT_35->Draw("same");
photons_w_HLT_25->Draw("same");


TPaveText *hlt_photon_new= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_new->SetTextSize(0.03);
hlt_photon_new->SetFillColor(0);
hlt_photon_new->SetBorderSize(0);
hlt_photon_new->SetShadowColor(0);
hlt_photon_new->SetTextAlign(21);
hlt_photon_new->AddText("");
hlt_photon_new->AddText("Pb+Pb");
hlt_photon_new->AddText("Tight Photons");
hlt_photon_new->AddText("HLT Loose Photons");
hlt_photon_new->AddText("p^{#gamma}_{T} > 15 GeV");
//phi_photon_2->AddText("HLT pT > 15 GeV");
hlt_photon_new->AddText("dR < 0.15");




new_diff_pT->AddEntry("photons_w_HLT_15","HLT Photons p_{T} > 15 GeV","pe");
new_diff_pT->AddEntry("photons_w_HLT_25","HLT Photons p_{T} > 25 GeV", "pe");
new_diff_pT->AddEntry("photons_w_HLT_35","HLT Photons p_{T} > 35 GeV", "pe");

new_diff_pT->Draw("same");
hlt_photon_new->Draw("same");


//----pT Distribution Through Barrel

TCanvas *new_ratio_pT_barrel = new TCanvas("new_ratio_pT_barrel","new_ratio_pT_barrel",600,500);
TLegend *new_diff_pT_barrel = new TLegend(0.6,0.5,0.85,0.6);

photons_w_HLT_15_barrel->SetTitle("Photons p_{T} Distribution");
photons_w_HLT_15_barrel->GetXaxis()->SetTitle("p_{T}");
photons_w_HLT_15_barrel->GetYaxis()->SetTitle("Counts");
photons_w_HLT_15_barrel->SetMarkerColor(kRed);
photons_w_HLT_15_barrel->SetLineColor(kRed);
photons_w_HLT_15_barrel->SetMarkerStyle(20);
//photons_w_HLT_15->Scale(1/photons_w_HLT_15->Integral());

photons_w_HLT_25_barrel->SetMarkerColor(kBlue);
photons_w_HLT_25_barrel->SetLineColor(kBlue);
photons_w_HLT_25_barrel->SetMarkerStyle(22);
//photons_w_HLT_25->Scale(1/photons_w_HLT_25->Integral());

photons_w_HLT_35_barrel->SetMarkerStyle(21);
photons_w_HLT_35_barrel->SetMarkerColor(kMagenta);
photons_w_HLT_35_barrel->SetLineColor(kMagenta);
//photons_w_HLT_35->Scale(1/photons_w_HLT_35->Integral());

photons_w_HLT_15_barrel->Draw();
photons_w_HLT_35_barrel->Draw("same");
photons_w_HLT_25_barrel->Draw("same");


TPaveText *hlt_photon_new_barrel= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_new_barrel->SetTextSize(0.03);
hlt_photon_new_barrel->SetFillColor(0);
hlt_photon_new_barrel->SetBorderSize(0);
hlt_photon_new_barrel->SetShadowColor(0);
hlt_photon_new_barrel->SetTextAlign(21);
hlt_photon_new_barrel->AddText("");
hlt_photon_new_barrel->AddText("Pb+Pb");
hlt_photon_new_barrel->AddText("Tight Photons");
hlt_photon_new_barrel->AddText("HLT Loose Photons");
hlt_photon_new_barrel->AddText("p^{#gamma}_{T} > 15 GeV");
//phi_photon_2->AddText("HLT pT > 15 GeV");
hlt_photon_new_barrel->AddText("dR < 0.15");
hlt_photon_new_barrel->AddText("|#eta| < 1.37");




new_diff_pT_barrel->AddEntry("photons_w_HLT_15_barrel","HLT Photons p_{T} > 15 GeV","pe");
new_diff_pT_barrel->AddEntry("photons_w_HLT_25_barrel","HLT Photons p_{T} > 25 GeV", "pe");
new_diff_pT_barrel->AddEntry("photons_w_HLT_35_barrel","HLT Photons p_{T} > 35 GeV", "pe");

new_diff_pT_barrel->Draw("same");
hlt_photon_new_barrel->Draw("same");

//----pT Distribution Through Endcaps

TCanvas *new_ratio_pT_endcaps = new TCanvas("new_ratio_pT_endcaps","new_ratio_pT_endcaps",600,500);
TLegend *new_diff_pT_endcaps = new TLegend(0.6,0.5,0.85,0.6);

photons_w_HLT_15_endcaps->SetTitle("Photons p_{T} Distribution");
photons_w_HLT_15_endcaps->GetXaxis()->SetTitle("p_{T}");
photons_w_HLT_15_endcaps->GetYaxis()->SetTitle("Counts");
photons_w_HLT_15_endcaps->SetMarkerColor(kRed);
photons_w_HLT_15_endcaps->SetLineColor(kRed);
photons_w_HLT_15_endcaps->SetMarkerStyle(20);
//photons_w_HLT_15->Scale(1/photons_w_HLT_15->Integral());

photons_w_HLT_25_endcaps->SetMarkerColor(kBlue);
photons_w_HLT_25_endcaps->SetLineColor(kBlue);
photons_w_HLT_25_endcaps->SetMarkerStyle(22);
//photons_w_HLT_25->Scale(1/photons_w_HLT_25->Integral());

photons_w_HLT_35_endcaps->SetMarkerStyle(21);
photons_w_HLT_35_endcaps->SetMarkerColor(kMagenta);
photons_w_HLT_35_endcaps->SetLineColor(kMagenta);
//photons_w_HLT_35->Scale(1/photons_w_HLT_35->Integral());

photons_w_HLT_15_endcaps->Draw();
photons_w_HLT_35_endcaps->Draw("same");
photons_w_HLT_25_endcaps->Draw("same");


TPaveText *hlt_photon_new_endcaps= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_new_endcaps->SetTextSize(0.036);
hlt_photon_new_endcaps->SetFillColor(0);
hlt_photon_new_endcaps->SetBorderSize(0);
hlt_photon_new_endcaps->SetShadowColor(0);
hlt_photon_new_endcaps->SetTextAlign(21);
hlt_photon_new_endcaps->AddText("");
hlt_photon_new_endcaps->AddText("Pb+Pb");
hlt_photon_new_endcaps->AddText("Tight Photons");
hlt_photon_new_endcaps->AddText("HLT Loose Photons");
hlt_photon_new_endcaps->AddText("p^{#gamma}_{T} > 15 GeV");
//phi_photon_2->AddText("HLT pT > 15 GeV");
hlt_photon_new_endcaps->AddText("dR < 0.15");
hlt_photon_new_endcaps->AddText("1.56 < |#eta| < 2.37");




new_diff_pT_endcaps->AddEntry("photons_w_HLT_15_endcaps","HLT Photons p_{T} > 15 GeV","pe");
new_diff_pT_endcaps->AddEntry("photons_w_HLT_25_endcaps","HLT Photons p_{T} > 25 GeV", "pe");
new_diff_pT_endcaps->AddEntry("photons_w_HLT_35_endcaps","HLT Photons p_{T} > 35 GeV", "pe");
new_diff_pT_endcaps->SetTextSize(0.3);

new_diff_pT_endcaps->Draw("same");
hlt_photon_new_endcaps->Draw("same");



//-------------------------------------------//
//------Photons Trhough the Encaps ONLY------//
//-------------------------------------------//
TCanvas *HLT_eff_endcaps = new TCanvas("HLT_eff_endcaps","HLT_eff_endcaps",600,500);
TLegend *new_diff_1 = new TLegend(0.6,0.5,0.85,0.6);

TH1D* HLT_eff_15_endcaps = new TH1D("HLT_eff_15_endcaps","",40,0,200);
HLT_eff_15_endcaps->SetName("HLT_eff_15_endcaps");                                                                                                                           
HLT_eff_15_endcaps->GetXaxis()->SetTitle("p_{T}");
HLT_eff_15_endcaps->SetMarkerStyle(20);
HLT_eff_15_endcaps->SetTitle("HLT Efficiency");
HLT_eff_15_endcaps->SetLineColor(kRed);
HLT_eff_15_endcaps->SetMarkerColor(kRed);
HLT_eff_15_endcaps->Divide(photons_w_HLT_15_endcaps,all_photons_endcaps,1,1,"B");
HLT_eff_15_endcaps->SetMaximum(1.1);

TH1D *HLT_eff_25_endcaps = new TH1D("HLT_eff_25_endcaps","",40,0,200);
HLT_eff_25_endcaps->SetName("HLT_eff_25_endcaps");
HLT_eff_25_endcaps->SetMarkerStyle(22);
HLT_eff_25_endcaps->SetMarkerColor(kBlue);
HLT_eff_25_endcaps->SetLineColor(kBlue);
HLT_eff_25_endcaps->Divide(photons_w_HLT_25_endcaps,all_photons_endcaps,1,1,"B");

TH1D *HLT_eff_35_endcaps = new TH1D("HLT_eff_35_endcaps","",40,0,200);
HLT_eff_35_endcaps->SetName("HLT_eff_35_endcaps");
HLT_eff_35_endcaps->SetMarkerStyle(21);
HLT_eff_35_endcaps->SetMarkerColor(kMagenta);
HLT_eff_35_endcaps->SetLineColor(kMagenta);
HLT_eff_35_endcaps->Divide(photons_w_HLT_35_endcaps,all_photons_endcaps,1,1,"B");


TPaveText *hlt_photon_1= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_1->SetTextSize(0.036);
hlt_photon_1->SetFillColor(0);
hlt_photon_1->SetBorderSize(0);
hlt_photon_1->SetShadowColor(0);
hlt_photon_1->SetTextAlign(21);
hlt_photon_1->AddText("");
hlt_photon_1->AddText("Pb+Pb");
hlt_photon_1->AddText("Tight Photons");
hlt_photon_1->AddText("HLT Loose Photons");
hlt_photon_1->AddText("p^{#gamma}_{T} > 15 GeV");
//phi_photon_2->AddText("HLT pT > 15 GeV");
hlt_photon_1->AddText("dR < 0.15");
hlt_photon_1->AddText("1.56 < |#eta| < 2.37");

HLT_eff_15_endcaps->Draw();
HLT_eff_35_endcaps->Draw("same");
HLT_eff_25_endcaps->Draw("same");


//new_diff->AddEntry("L1EM_eff", "Photons fired L1 Triggers/All Photons", "pe");
new_diff_1->AddEntry("HLT_eff_15_endcaps", "HLT Photons p_{T} > 15 GeV", "pe");
new_diff_1->AddEntry("HLT_eff_25_endcaps","HLT Photons p_{T} > 25 GeV","pe");
new_diff_1->AddEntry("HLT_eff_35_endcaps","HLT Photons p_{T} > 35 GeV", "pe");
new_diff_1->SetTextSize(0.03);
new_diff_1->SetBorderSize(0);
new_diff_1->Draw("same");
hlt_photon_1->Draw("same");



//-------------------------------------------//
//------Photons Trhough the Barrel ONLY------//
//-------------------------------------------//
TCanvas *HLT_eff_barrel = new TCanvas("HLT_eff_barrel","HLT_eff_barrel",600,500);
TLegend *new_diff_2 = new TLegend(0.6,0.5,0.85,0.6);


TH1D* HLT_eff_15_barrel = new TH1D("HLT_eff_15_barrel","",40,0,200);
HLT_eff_15_barrel->SetName("HLT_eff_15_barrel");                                                                                                                           
HLT_eff_15_barrel->GetXaxis()->SetTitle("p_{T}");
HLT_eff_15_barrel->SetMarkerStyle(20);
HLT_eff_15_barrel->SetTitle("HLT Efficiency");
HLT_eff_15_barrel->SetLineColor(kRed);
HLT_eff_15_barrel->SetMarkerColor(kRed);
HLT_eff_15_barrel->Divide(photons_w_HLT_15_barrel,all_photons_barrel,1,1,"B");
HLT_eff_15_barrel->SetMaximum(1.1);

TH1D *HLT_eff_25_barrel = new TH1D("HLT_eff_25_barrel","",40,0,200);
HLT_eff_25_barrel->SetName("HLT_eff_25_barrel");
HLT_eff_25_barrel->SetMarkerStyle(22);
HLT_eff_25_barrel->SetMarkerColor(kBlue);
HLT_eff_25_barrel->SetLineColor(kBlue);
HLT_eff_25_barrel->Divide(photons_w_HLT_25_barrel,all_photons_barrel,1,1,"B");

TH1D *HLT_eff_35_barrel = new TH1D("HLT_eff_35_barrel","",40,0,200);
HLT_eff_35_barrel->SetName("HLT_eff_35_barrel");
HLT_eff_35_barrel->SetMarkerStyle(21);
HLT_eff_35_barrel->SetMarkerColor(kMagenta);
HLT_eff_35_barrel->SetLineColor(kMagenta);
HLT_eff_35_barrel->Divide(photons_w_HLT_35_barrel,all_photons_barrel,1,1,"B");


TPaveText *hlt_photon_2= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
hlt_photon_2->SetTextSize(0.036);
hlt_photon_2->SetFillColor(0);
hlt_photon_2->SetBorderSize(0);
hlt_photon_2->SetShadowColor(0);
hlt_photon_2->SetTextAlign(21);
hlt_photon_2->AddText("");
hlt_photon_2->AddText("Pb+Pb");
hlt_photon_2->AddText("Tight Photons");
hlt_photon_2->AddText("HLT Loose Photons");
hlt_photon_2->AddText("p^{#gamma}_{T} > 15 GeV");
//phi_photon_2->AddText("HLT pT > 15 GeV");
hlt_photon_2->AddText("dR < 0.15");
hlt_photon_2->AddText("|#eta| < 1.37");

HLT_eff_15_barrel->Draw();
HLT_eff_35_barrel->Draw("same");
HLT_eff_25_barrel->Draw("same");


//new_diff->AddEntry("L1EM_eff", "Photons fired L1 Triggers/All Photons", "pe");
new_diff_2->AddEntry("HLT_eff_15_barrel", "HLT Photons p_{T} > 15 GeV", "pe");
new_diff_2->AddEntry("HLT_eff_25_barrel","HLT Photons p_{T} > 25 GeV","pe");
new_diff_2->AddEntry("HLT_eff_35_barrel","HLT Photons p_{T} > 35 GeV", "pe");
new_diff_2->SetTextSize(0.03);
new_diff_2->SetBorderSize(0);
new_diff_2->Draw("same");
hlt_photon_2->Draw("same");

/**
TCanvas *new_pT_diff = new TCanvas("new_pT_diff", "new_pT_diff",600,500);

photons_wL1EM->SetMarkerStyle(21);
photons_wL1EM->SetLineColor(kBlack);
photons_wL1EM->SetMarkerColor(kBlack);
photons_wL1EM->SetTitle("Photon p_{T} Distribution (Triggered)");
photons_wL1EM->GetXaxis()->SetTitle("p_{T}");
photons_wL1EM->GetYaxis()->SetTitle("1/N");

photons_wL1EM->Scale(1/photons_wL1EM->Integral());


photons_wL1EM->Draw();

TCanvas *new_pT_new = new TCanvas("new_pT_pT", "new_pT_new",600,500);

photons_w_HLT->SetMarkerColor(kBlue);
photons_w_HLT->SetLineColor(kBlue);
photons_w_HLT->SetMarkerStyle(20);
photons_w_HLT->Scale(1/photons_w_HLT->Integral());
photons_w_HLT->SetTitle("Photon p_{T} Distribution (Triggered HLT)");
photons_w_HLT->GetXaxis()->SetTitle("p_{T}");
photons_w_HLT->GetYaxis()->SetTitle("1/N");


photons_w_HLT->Draw("same");


*/

/**
photons_wL1EM->Sumw2();
all_photons->Sumw2();

TH1D *ratio_test = (TH1D*) photons_wL1EM->Clone("photons_wL1EM");    
ratio_test->Sumw2();
ratio_test->SetName("ratio_test");
ratio_test->Divide(all_photons);
ratio_test->SetMarkerStyle(20);
ratio_test->SetMarkerColor(kBlack);
ratio_test->SetLineColor(kBlack);
ratio_test->SetTitle("L1EM Efficiency");
ratio_test->GetXaxis()->SetTitle("p_{T}");

ratio_test->Draw();


*/





TCanvas *new_L1EM_num = new TCanvas("new_L1EM_num","new_L1EM_num",600,500);

L1EM_num->SetTitle("Number of L1EM Triggers / Event");
L1EM_num->SetMarkerStyle(20);
L1EM_num->SetMarkerColor(kBlack);
L1EM_num->SetLineColor(kBlack); 
L1EM_num->GetYaxis()->SetTitle("Counts");
L1EM_num->GetXaxis()->SetTitle("Total Number of L1EM Triggers Fired");
L1EM_num->SetName("L1EM_num");

L1EM_num->Draw();

TCanvas *new_HLT_num_1 = new TCanvas("new_HLT_num_1","new_HLT_num_1",600,500);
//-----Total HLT Triggers Fired
HLT_num->SetTitle("Number of HLT Triggers / Event");
HLT_num->SetMarkerColor(kBlack);
HLT_num->SetLineColor(kBlack);
HLT_num->SetMarkerStyle(20);
HLT_num->GetYaxis()->SetTitle("Counts");
HLT_num->GetXaxis()->SetTitle("Total Number of HLT Triggers Fired");
HLT_num->SetName("HLT_num");

HLT_num->Draw();

//--------photon pT VS Total L1EM Triggers Fired
TCanvas *new_pt_L1EM = new TCanvas("new_pt_L1EM","new_pt_L1EM",600,500);
TLegend *new_pt_L1EM_lg = new TLegend(0.6,0.5,0.85,0.6);

pT_vs_L1EM->SetTitle("Total Number of L1EM Triggers VS Photon pT");
pT_vs_L1EM->GetXaxis()->SetTitle("Photon p_{T}");
pT_vs_L1EM->GetYaxis()->SetTitle("Total Number of L1EM Triggers Fired");
pT_vs_L1EM->SetName("pT_vs_L1EM");

pT_vs_L1EM->Draw("colz");

//-----photon pT VS Total HLT Triggers Fired
TCanvas *new_HLT_num = new TCanvas("new_HLT_num","new_HLT_num",600,500);

pT_vs_HLT_num->SetTitle("Total Number of HLT Triggers Fired VS Photon pT");
pT_vs_HLT_num->GetXaxis()->SetTitle("Photon p_{T}");
pT_vs_HLT_num->GetYaxis()->SetTitle("Total Number of HLT Triggers Fired");
pT_vs_HLT_num->SetName("pT_vs_HLT_num");

pT_vs_HLT_num->Draw("colz");

/**
//-----------------------------------------------------//
//---------0 Triggers Fired For These Photons----------//
//-----------------------------------------------------//

//--------pT
TCanvas *new_nofire = new TCanvas("new_nofire","new_nofire",600,500);

nofire_photons_pt->SetTitle("Photon p_{T}");
nofire_photons_pt->SetMarkerStyle(20);
nofire_photons_pt->SetMarkerColor(kBlack);
nofire_photons_pt->SetLineColor(kBlack); 
nofire_photons_pt->GetYaxis()->SetTitle("Counts");
nofire_photons_pt->GetXaxis()->SetTitle("p_{T}");
nofire_photons_pt->SetName("nofire_photons_pt");


TPaveText *nofire_txt= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
nofire_txt->SetTextSize(0.03);
nofire_txt->SetFillColor(0);
nofire_txt->SetBorderSize(0);
nofire_txt->SetShadowColor(0);
nofire_txt->SetTextAlign(21);
nofire_txt->AddText("");
nofire_txt->AddText("Pb+Pb");
nofire_txt->AddText("Loose Photons");
nofire_txt->AddText("p^{#gamma}_{T} > 60 GeV");
nofire_txt->AddText("L1EMRO > 30 GeV");
nofire_txt->AddText("dR < 0.15");

nofire_photons_pt->Draw();
nofire_txt->Draw("same");

//--------eta
TCanvas *new_nofire_eta = new TCanvas("new_nofire_eta","new_nofire_eta",600,500);

nofire_photons_eta->SetTitle("Photon #eta");
nofire_photons_eta->SetMarkerStyle(20);
nofire_photons_eta->SetMarkerColor(kBlack);
nofire_photons_eta->SetLineColor(kBlack); 
nofire_photons_eta->GetYaxis()->SetTitle("Counts");
nofire_photons_eta->GetXaxis()->SetTitle("#eta");
nofire_photons_eta->SetName("nofire_photons_eta");


TPaveText *nofire_tx= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
nofire_tx->SetTextSize(0.03);
nofire_tx->SetFillColor(0);
nofire_tx->SetBorderSize(0);
nofire_tx->SetShadowColor(0);
nofire_tx->SetTextAlign(21);
nofire_tx->AddText("");
nofire_tx->AddText("Pb+Pb");
nofire_tx->AddText("Loose Photons");
nofire_tx->AddText("p^{#gamma}_{T} > 60 GeV");
nofire_tx->AddText("L1EMRO > 30 GeV");
nofire_tx->AddText("dR < 0.15");

nofire_photons_eta->Draw();
nofire_tx->Draw("same");

//--------phi
TCanvas *new_nofire_phi = new TCanvas("new_nofire_phi","new_nofire_phi",600,500);

nofire_photons_phi->SetTitle("Photon #phi");
nofire_photons_phi->SetMarkerStyle(20);
nofire_photons_phi->SetMarkerColor(kBlack);
nofire_photons_phi->SetLineColor(kBlack); 
nofire_photons_phi->GetYaxis()->SetTitle("Counts");
nofire_photons_phi->GetXaxis()->SetTitle("#phi");
nofire_photons_phi->SetName("nofire_photons_phi");


TPaveText *nofire_t= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
nofire_t->SetTextSize(0.03);
nofire_t->SetFillColor(0);
nofire_t->SetBorderSize(0);
nofire_t->SetShadowColor(0);
nofire_t->SetTextAlign(21);
nofire_t->AddText("");
nofire_t->AddText("Pb+Pb");
nofire_t->AddText("Loose Photons");
nofire_t->AddText("p^{#gamma}_{T} > 60 GeV");
nofire_t->AddText("L1EMRO > 30 GeV");
nofire_t->AddText("dR < 0.15");

nofire_photons_phi->Draw();
nofire_t->Draw("same");

*/

/**
//-----SumET
//---photons with no L1EM Trigger Fired
TCanvas *new_nofire_sumET = new TCanvas("new_nofire_sumET","new_nofire_sumET",600,500);
nofire_sumET->SetTitle("SumET");
nofire_sumET->SetMarkerStyle(20);
nofire_sumET->SetMarkerColor(kRed);
nofire_sumET->SetLineColor(kRed); 
nofire_sumET->GetYaxis()->SetTitle("Counts");
nofire_sumET->GetXaxis()->SetTitle("SumET [TeV]");
nofire_sumET->SetName("nofire_sumET");

TPaveText *nofire_sumET_txt= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
nofire_sumET_txt->SetTextSize(0.03);
nofire_sumET_txt->SetFillColor(0);
nofire_sumET_txt->SetBorderSize(0);
nofire_sumET_txt->SetShadowColor(0);
nofire_sumET_txt->SetTextAlign(21);
nofire_sumET_txt->AddText("");
nofire_sumET_txt->AddText("Pb+Pb");
nofire_sumET_txt->AddText("Loose Photons");
nofire_sumET_txt->AddText("p^{#gamma}_{T} > 20 GeV");
nofire_sumET_txt->AddText("L1EMRO > 30 GeV");
nofire_sumET_txt->AddText("dR < 0.15");

//nofire_sumET->Draw();
//nofire_sumET_txt->Draw("same");

//---photons with L1EM Trigger Fired
//TCanvas *new_fire_sumET = new TCanvas("new_nofire_sumET","new_nofire_sumET",600,500);
TLegend *new_sumET_lg = new TLegend(0.6,0.5,0.85,0.6);

all_photons_sumET->SetTitle("SumET");
all_photons_sumET->SetMarkerStyle(20);
all_photons_sumET->SetMarkerColor(kBlack);
all_photons_sumET->SetLineColor(kBlack); 
all_photons_sumET->GetYaxis()->SetTitle("Counts");
all_photons_sumET->GetXaxis()->SetTitle("SumET [TeV]");
all_photons_sumET->SetName("all_photons_sumET");

fire_sumET->SetMarkerStyle(20);
fire_sumET->SetMarkerColor(kBlack);
fire_sumET->SetLineColor(kBlack);


TPaveText *fire_sumET_txt= new TPaveText(0.57,0.68,0.77,0.88,"NDC");
fire_sumET_txt->SetTextSize(0.03);
fire_sumET_txt->SetFillColor(0);
fire_sumET_txt->SetBorderSize(0);
fire_sumET_txt->SetShadowColor(0);
fire_sumET_txt->SetTextAlign(21);
fire_sumET_txt->AddText("");
fire_sumET_txt->AddText("Pb+Pb");
fire_sumET_txt->AddText("Loose Photons");
fire_sumET_txt->AddText("p^{#gamma}_{T} > 20 GeV");
fire_sumET_txt->AddText("L1EMRO > 30 GeV");
fire_sumET_txt->AddText("dR < 0.15");




//all_photons_sumET->Draw();
nofire_sumET->Draw();
fire_sumET->Draw("same");



fire_sumET_txt->Draw("same");

new_sumET_lg->AddEntry("fire_sumET","SumET for Photons That Fired An L1EMX Trigger","l");
new_sumET_lg->AddEntry("nofire_sumET", "SumET for Photons w/ no fired L1 EM RoI Triggers","l");

new_sumET_lg->Draw("same");
*/
//----HLT
//---HLT Photons pT > 15 GeV
TCanvas *hlt_nofire_sumET_15L = new TCanvas("hlt_nofire_sumET_15L","hlt_nofire_sumET_15L",600,500);
TLegend *hlt_fire_sumET_lg_15 = new TLegend(0.6,0.5,0.85,0.6);
TPaveText *hlt_fire_sumET_15_txt= new TPaveText(0.57,0.68,0.77,0.88,"NDC");

hlt_fire_sumET_15->SetTitle("SumET");
hlt_fire_sumET_15->SetMarkerStyle(20);
hlt_fire_sumET_15->SetMarkerColor(kRed);
hlt_fire_sumET_15->SetLineColor(kRed); 
hlt_fire_sumET_15->GetYaxis()->SetTitle("Counts");
hlt_fire_sumET_15->GetXaxis()->SetTitle("SumET [TeV]");
hlt_fire_sumET_15->SetName("hlt_fire_sumET_15");

hlt_nofire_sumET_15->SetMarkerColor(kBlack);
hlt_nofire_sumET_15->SetLineColor(kBlack);

hlt_nofire_sumET_15->SetName("hlt_nofire_sumET_15");

hlt_fire_sumET_15->Draw();
hlt_nofire_sumET_15->Draw("same");


hlt_fire_sumET_lg_15->AddEntry("hlt_fire_sumET_15","SumET for Photons That Fired HLT Trigger","l");
hlt_fire_sumET_lg_15->AddEntry("hlt_nofire_sumET_15","SumET for Photons w/ no fired HLT Triggers");

hlt_fire_sumET_lg_15->Draw("same");


hlt_fire_sumET_15_txt->SetTextSize(0.03);
hlt_fire_sumET_15_txt->SetFillColor(0);
hlt_fire_sumET_15_txt->SetBorderSize(0);
hlt_fire_sumET_15_txt->SetShadowColor(0);
hlt_fire_sumET_15_txt->SetTextAlign(21);
hlt_fire_sumET_15_txt->AddText("HLT Photon p_{T} > 15 GeV");
hlt_fire_sumET_15_txt->Draw("same");


//---HLT Photons pT > 25 GeV
TCanvas *hlt_nofire_sumET_25L = new TCanvas("hlt_nofire_sumET_25L","hlt_nofire_sumET_25L",600,500);
TLegend *hlt_fire_sumET_lg_25 = new TLegend(0.6,0.5,0.85,0.6);

TPaveText *hlt_fire_sumET_25_txt= new TPaveText(0.57,0.68,0.77,0.88,"NDC");

hlt_fire_sumET_25->SetTitle("SumET");
hlt_fire_sumET_25->SetMarkerStyle(20);
hlt_fire_sumET_25->SetMarkerColor(kRed);
hlt_fire_sumET_25->SetLineColor(kRed); 
hlt_fire_sumET_25->GetYaxis()->SetTitle("Counts");
hlt_fire_sumET_25->GetXaxis()->SetTitle("SumET [TeV]");
hlt_fire_sumET_25->SetName("hlt_fire_sumET_25");

hlt_nofire_sumET_25->SetMarkerColor(kBlack);
hlt_nofire_sumET_25->SetLineColor(kBlack);

hlt_nofire_sumET_25->SetName("hlt_nofire_sumET_25");

hlt_fire_sumET_25->Draw();
hlt_nofire_sumET_25->Draw("same");

hlt_fire_sumET_lg_25->AddEntry("hlt_fire_sumET_25","SumET for Photons THat Fired HLT Trigger", "l");
hlt_fire_sumET_lg_25->AddEntry("hlt_nofire_sumET_25","SumET for Photons w/ not fire HLT Triggers","l");


hlt_fire_sumET_lg_25->Draw("same");


hlt_fire_sumET_25_txt->SetTextSize(0.03);
hlt_fire_sumET_25_txt->SetFillColor(0);
hlt_fire_sumET_25_txt->SetBorderSize(0);
hlt_fire_sumET_25_txt->SetShadowColor(0);
hlt_fire_sumET_25_txt->SetTextAlign(21);
hlt_fire_sumET_25_txt->AddText("HLT Photon p_{T} > 25 GeV");
hlt_fire_sumET_25_txt->Draw("same");

//----HLT_Photons pT> 35 GeV
TCanvas *hlt_nofire_sumET_35L = new TCanvas("hlt_nofire_sumET_35L","hlt_nofire_sumET_35L",600,500);
TLegend *hlt_fire_sumET_lg_35 = new TLegend(0.6,0.5,0.85,0.6);

TPaveText *hlt_fire_sumET_35_txt= new TPaveText(0.57,0.68,0.77,0.88,"NDC");

hlt_fire_sumET_35->SetTitle("SumET");
hlt_fire_sumET_35->SetMarkerStyle(20);
hlt_fire_sumET_35->SetMarkerColor(kRed);
hlt_fire_sumET_35->SetLineColor(kRed); 
hlt_fire_sumET_35->GetYaxis()->SetTitle("Counts");
hlt_fire_sumET_35->GetXaxis()->SetTitle("SumET [TeV]");
hlt_fire_sumET_35->SetName("hlt_fire_sumET_35");


hlt_nofire_sumET_35->SetMarkerColor(kBlack);
hlt_nofire_sumET_35->SetLineColor(kBlack);

hlt_nofire_sumET_35->SetName("hlt_nofire_sumET_35");

hlt_fire_sumET_35->Draw();
hlt_nofire_sumET_35->Draw("same");

hlt_fire_sumET_lg_35->AddEntry("hlt_fire_sumET_35","SumET for Photons THat Fired HLT Trigger", "l");
hlt_fire_sumET_lg_35->AddEntry("hlt_nofire_sumET_35","SumET for Photons w/ not fire HLT Triggers","l");


hlt_fire_sumET_lg_35->Draw("same");

hlt_fire_sumET_35_txt->SetTextSize(0.03);
hlt_fire_sumET_35_txt->SetFillColor(0);
hlt_fire_sumET_35_txt->SetBorderSize(0);
hlt_fire_sumET_35_txt->SetShadowColor(0);
hlt_fire_sumET_35_txt->SetTextAlign(21);
hlt_fire_sumET_35_txt->AddText("HLT Photon p_{T} > 35 GeV");
hlt_fire_sumET_35_txt->Draw("same");

  return 0;
}
