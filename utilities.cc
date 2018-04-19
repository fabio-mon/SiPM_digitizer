
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
#include "utilities.hh"
using namespace std;

double SiliconReflectivity(const double &energy)
{
   return 0.33249 - 0.0469204*energy + 0.028069*energy*energy + 
          0.0856852 * exp(-(energy-3.36542)*(energy-3.36542) / (2*0.183577*0.183577)) +
          0.05148 * exp(-(energy-4.39853)*(energy-4.39853) / (2*0.14405*0.14405));
}

double PDE(const double &energy)
{
   //cout<<"energy = "<<energy<<endl;
   //cout<<"PDE = "<<0.39784 * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512))<<endl;
   //cout<<"Reflec = "<<SiliconReflectivity(energy)<<endl;
   //cout<<"TOT PDE = "<<0.39784/0.48/SiliconReflectivity(energy) * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512))<<endl<<endl;
   return 0.39784/0.7*0.25/*/SiliconReflectivity(energy) */ * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512));
}

double PDEHPK_MPPC(const double &energy)
{
  return 0.37043 * exp(-(energy - 2.69767)*(energy - 2.69767) / (2 * 0.692431 * 0.692431));
}

void Fill(std::vector<photon> &Ph_event, std::vector<float>* Ph_time, /*std::vector<float>* Ph_energy, std::vector<float>* Ph_x_hit, 
std::vector<float>* Ph_y_hit,*/ std::vector<int>* Ph_creator_proc,int cherenkov)
{
   //cout<<"\n>>In function Fill"<<endl;
   struct photon newPhoton;
   unsigned glob_it = 0;
   for(std::vector<Float_t>::iterator it = Ph_time->begin() ; it != Ph_time->end(); ++it)
   {
      if( Ph_time->at(glob_it) > 0.)
      {
         //cout<<Ph_time->at(glob_it)<<"\t"<<Ph_energy->at(glob_it)<<"\t"<<Ph_x_hit->at(glob_it)<<"\t"<<Ph_y_hit->at(glob_it)<<"\t"<<Ph_creator_proc->at(glob_it)<<endl;
         if(cherenkov == 0)//only scint
         {
            if(Ph_creator_proc->at(glob_it) == 0)
            {
               newPhoton.time = Ph_time->at(glob_it)+10;//Set an offset of 5 ns is useful for the baseline subtraction
               newPhoton.energy = 0.;//Ph_energy->at(glob_it);
               newPhoton.x_hit =0.;// Ph_x_hit->at(glob_it);
               newPhoton.y_hit = 0.;//Ph_y_hit->at(glob_it);
               newPhoton.creator_proc = Ph_creator_proc->at(glob_it);
               Ph_event.push_back(newPhoton);
            }
         }
         else
            if(cherenkov == 2)//only cherenkov
            {
               if(Ph_creator_proc->at(glob_it) == 1)
               {
                  newPhoton.time = Ph_time->at(glob_it)+10;//Set an offset of 5 ns is useful for the baseline subtraction
                  newPhoton.energy = 0.;//Ph_energy->at(glob_it);
                  newPhoton.x_hit =0.;// Ph_x_hit->at(glob_it);
                  newPhoton.y_hit = 0.;//Ph_y_hit->at(glob_it);
                  newPhoton.creator_proc = Ph_creator_proc->at(glob_it);
                  Ph_event.push_back(newPhoton);
               }
            }
            else //cherenkov & scint
            {
               newPhoton.time = Ph_time->at(glob_it)+10;//Set an offset of 5 ns is useful for the baseline subtraction
               newPhoton.energy = 0.;//Ph_energy->at(glob_it);
               newPhoton.x_hit =0.;// Ph_x_hit->at(glob_it);
               newPhoton.y_hit = 0.;//Ph_y_hit->at(glob_it);
               newPhoton.creator_proc = Ph_creator_proc->at(glob_it);
               Ph_event.push_back(newPhoton);
            }
      }
      glob_it++;
   }
}

void Position_selection(std::vector<photon> &Ph_event,double SiPMLength)
{
   cout<<"\n>>In function Position_selection"<<endl;
   cout<<"# Ph in = "<<Ph_event.size()<<endl;
   const double CellSize = 0.02;
   const double CellDist = 0.1;

   const double NcellX = SiPMLength / CellDist;
   double CellXCenter;
   double CellYCenter;
   
   for(std::vector<photon>::iterator it = Ph_event.begin() ; it != Ph_event.end(); ++it)
   {
      CellXCenter=CellDist/2.;
      CellYCenter=CellDist/2.;
      //cout<<"x_hit = "<< it->x_hit + SiPMLength/2.<<"\ny_hit = "<< it->y_hit + SiPMLength/2. <<endl;
      //cout<<"Scanning X stripes"<<endl;
      
      while(fabs(CellXCenter - (it->x_hit + SiPMLength/2.))>(CellSize/2.) && CellXCenter<SiPMLength)
      {
         //cout<<CellXCenter<<endl;
         CellXCenter += CellDist;
      }
      if(CellXCenter>SiPMLength)
      {
         Ph_event.erase(it);
         it--;
         //cout<<"No strips found"<<endl;
      }
      else
      {
         //cout<<"Strips matched at center "<<CellXCenter<<endl;
         while(fabs(CellYCenter - (it->y_hit + SiPMLength/2.))>(CellSize/2.) && CellYCenter<SiPMLength)
            CellYCenter += CellDist;
         if(CellYCenter>SiPMLength)
         {  
            Ph_event.erase(it);
            it--;
         }
      }
   }
   cout<<"# Ph fin = "<<Ph_event.size()<<endl;
}

void PDERandomDeleting(TRandom2 &gen,std::vector<photon> &Ph_event)
{
   //cout<<"\n>>In function RandomDeleting"<<endl;
   //cout<<"# Ph in = "<<Ph_event.size()<<endl;
   for(std::vector<photon>::iterator it = Ph_event.begin() ; it != Ph_event.end(); ++it)
   {
      double rndm = gen.Uniform(0.,1.);
      //cout<< "PDE("<<it->energy<<") = "<<PDE(it->energy*1.E+06)<<endl;
      if(rndm>0.3)
      {
         //cout<<"deleted"<<endl;
         Ph_event.erase(it);
         it--;
      }
   }
   //cout<<"# Ph fin = "<<Ph_event.size()<<endl;
}


void SPTR (TRandom2 &gen,std::vector<photon> &Ph_event)
{
   //cout<<"\n>>In function SPTR"<<endl;
   for(std::vector<photon>::iterator it = Ph_event.begin() ; it != Ph_event.end(); ++it)
     it->time += gen.Gaus(0.,0.066);   
}

void Dark_current(TRandom2 &gen,std::vector<photon> &Ph_event,const double rate)
{
     //cout<<"In function dark current"<<endl;
     const double time_interval = 40+100;
     const double mean_dark_ph = time_interval*rate;
     int N_dark_ph = gen.PoissonD(mean_dark_ph);
     double ph_t;
     
     //cout<<"generate"<<N_dark_ph<<" dark current photons"<<endl;
     struct photon newPhoton;
     for(int i=0;i<N_dark_ph;i++)
     {
         ph_t = gen.Uniform(-100.,time_interval-100);
         //cout<<ph_t<<endl;
         newPhoton.time = ph_t;
         newPhoton.energy = 0.;//Ph_energy->at(glob_it);
         newPhoton.x_hit =0.;// Ph_x_hit->at(glob_it);
         newPhoton.y_hit = 0.;//Ph_y_hit->at(glob_it);
         newPhoton.creator_proc = 3;
         Ph_event.push_back(newPhoton);
     }
}

void Digitize(TRandom2 &gen, std::vector<photon> &Ph_event, float *V,float *t,const int Nsample) 
{
  //cout<<"\nIn function Digitize"<<endl;
  const double A_mean=1.;
  const double t_rise = 0.2;
  const double t_decay = 10.;
  const double t_sample = 0.01;  

  double A = gen.Gaus(A_mean,A_mean/10.);
  const double b = t_decay/t_rise;
  const double denom = pow(b,1./(1.-b))+pow(b,1./(1./b - 1.));
  double T;
  int i=0;

  //Initialize V and set the time sampling values
  for(i=0;i<Nsample;i++)
  {
     t[i]=i * t_sample;
     V[i]=0;
  }
  for(std::vector<photon>::iterator it = Ph_event.begin() ; it != Ph_event.end(); ++it)
  {
    T=0.;
    //cout<<"\nadding photon arriving at "<<it->time<<" ns"<<endl;
    for (int i=0;i<Nsample;i++)
    {
       if(T <= it->time)
       {
          //cout<<T<<"<="<<it->time<<endl;  
          V[i] += 0.;
       }
       else
       {
          V[i] += A * ( exp(-(T - it->time)/t_decay) - exp(-(T - it->time)/t_rise) ) / denom; 
          //cout<<"V("<<T<<") = "<<V[i]<<endl; 
       }
       T += t_sample;
    }   
  }     
}

float LeadEdgeDiscr(float *V,float *t,const double threeshold)
{
   float time = -1;
   for(int i=0;i<4000;i++)
   {
      //cout<<"V("<<t[i]<<") = "<<V[i]<<endl;
      if(V[i]>threeshold && time==-1)
      {
         time = t[i];
         //cout<<"time over threshold = "<<time<<endl;
      }
   }
   return time;
}
