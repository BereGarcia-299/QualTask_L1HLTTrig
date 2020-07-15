#include <assert.h>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <tgmath.h>
#include <iostream>
#include <TChain.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"
#include <TAttMarker.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TPad.h>
#include <TGaxis.h> 
#include <TObject.h>

using namespace std;

float pi_dos = 360;
float pi_mi = 180;
float original = 3.14159;

float delta_phi(float phi_uno,float phi_dos){                                                                                                                                   

  float final_phi = fmod(abs(phi_uno - phi_dos), 2*original);
  float phi = final_phi > original  ? 2*original - final_phi : final_phi;
  return phi;

}


int practice_hlt(){

  TChain *t_data = new TChain("analysis");

  TH1D *photon_phi = new TH1D("photon_phi","",20,0,3.14);

  //t_data->Add("/Users/berenicegarcia/Desktop/ATLAS/Photon+JetInSitu/data17_5TeV.root");

  t_data->Add("/Users/berenicegarcia/Desktop/ATLAS/Qual-Task/PbPb_Data18/user.berenice.21685833._000001.ANALYSIS.root");

  std::vector<float> *b_HLT_pt;
  std::vector<float> *b_HLT_eta;
  std::vector<float> *b_HLT_etaBE;
  std::vector<float> *b_HLT_phi;

  int b_HLT_photons_n = 0;

  t_data->SetBranchAddress("b_HLT_pt", &b_HLT_pt);
  t_data->SetBranchAddress("b_HLT_eta", &b_HLT_eta);
  t_data->SetBranchAddress("b_HLT_etaBE", &b_HLT_etaBE);
  t_data->SetBranchAddress("b_HLT_phi", &b_HLT_phi);
  t_data->SetBranchAddress("b_HLT_photons_n", &b_HLT_photons_n);


    int event_num = t_data->GetEntries();


 	for (int iEvent = 0; iEvent < event_num; ++iEvent){

 		for (int iphi = 0; iphi < b_HLT_photons_n; ++iphi)
 		{
 			photon_phi->Fill(abs(b_HLT_phi->at(iphi)));
 		}


 	}

 	photon_phi->Draw();

	return 0;
}