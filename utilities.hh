
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TApplication.h"
#include "TRandom2.h"
#include "TMath.h"

#include "ConfigFile.hh"

using namespace std;

struct photon
{
   float time;
   float energy;
   float x_hit;
   float y_hit;
   int creator_proc;
};

double SiliconReflectivity(const double &energy);
double PDE(const double &energy);
double PDEHPK_MPPC(const double &energy);
void Fill(std::vector<photon> &Ph_event, std::vector<float>* Ph_time, /*std::vector<float>* Ph_energy, std::vector<float>* Ph_x_hit, 
									std::vector<float>* Ph_y_hit,*/ std::vector<int>* Ph_creator_proc,int cherenkov);
void Position_selection(std::vector<photon> &Ph_event,double SiPMLength);
void PDERandomDeleting(TRandom2 &gen,std::vector<photon> &Ph_event);
void SPTR (TRandom2 &gen,std::vector<photon> &Ph_event);
void Dark_current(TRandom2 &gen,std::vector<photon> &Ph_event,const double rate);
void Digitize(TRandom2 &gen, std::vector<photon> &Ph_event, float *V,float *t,const int Nsample); 
float LeadEdgeDiscr(float *V,float *t,const double threeshold);
