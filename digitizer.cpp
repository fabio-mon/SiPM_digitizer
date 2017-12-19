
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
std::vector<float>* Ph_y_hit,*/ std::vector<int>* Ph_creator_proc,const bool cherenkov)
{
   //cout<<"\n>>In function Fill"<<endl;
   struct photon newPhoton;
   unsigned glob_it = 0;
   for(std::vector<Float_t>::iterator it = Ph_time->begin() ; it != Ph_time->end(); ++it)
   {
      if( Ph_time->at(glob_it) > 0.)
      {
         //cout<<Ph_time->at(glob_it)<<"\t"<<Ph_energy->at(glob_it)<<"\t"<<Ph_x_hit->at(glob_it)<<"\t"<<Ph_y_hit->at(glob_it)<<"\t"<<Ph_creator_proc->at(glob_it)<<endl;
         if(cherenkov == false)
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
/*
TCanvas* Draw(std::vector<photon> &Ph_event)
{
   TCanvas* c = new TCanvas();
   TH2F* Ph_hit = new TH2F("Photon hit position","Photon hit position",100,-6.,6.,100,-6.,6.);
   for(std::vector<photon>::iterator it = Ph_event.begin() ; it != Ph_event.end(); ++it)
      Ph_hit->Fill(it->x_hit,it->y_hit);
   Ph_hit->Draw("COLZ");         
   return c;
}
*/
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

float TimeOverThresh(float *V,float *t,const double threeshold)
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

int main(int argc, char** argv)
{
   if(argc ==1 || argc>=3)
   {
      cout<<"ERROR: unvalid number of input parameters\n";
      exit(EXIT_FAILURE);
   }
     
   ConfigFile config(argv[1]);

   TString filename;
   if (config.keyExists("filename"))
     filename = config.read<string> ("filename");
   else
   {
      cout<<"ERROR: filename not found in configfile\n";
      exit(EXIT_FAILURE);
   }

   double rate;
   if (config.keyExists("dcr"))
     rate = config.read<double> ("dcr");
   else
   {
      cout<<"ERROR: dcr not found in configfile\n";
      exit(EXIT_FAILURE);
   }

   bool cherenkov = true;
   if (config.keyExists("cherenkov"))
     cherenkov = config.read<bool> ("cherenkov");

  //define and initialize the random generator
  TRandom2 gen;
  gen.SetSeed();

  //define the variables useful to read the input tree
  Int_t eventNb;
  Float_t mu_x_hit;
  Float_t mu_y_hit;
  Float_t mu_Edep;
  Float_t mu_TrackLength;
  Int_t Ph_tot_Nb;
  Int_t Ph_detected_Nb;
  Int_t Ph_lost_Nb;
  Int_t Ph_labs_Nb;
  std::vector<Float_t> * Ph_time = 0;
  std::vector<Float_t> * Ph_energy = 0;
  std::vector<Int_t> * Ph_creator_proc = 0;
  std::vector<Float_t> * Ph_x_hit = 0;
  std::vector<Float_t> * Ph_y_hit = 0;
  std::vector<Int_t> * Ph_SiPM_number = 0;

  //The struct photon contains all the useful information of each photon 
  std::vector<photon> Ph_event;

  //Open the input file and set the tree branch
  TFile* file = TFile::Open(filename.Data());
  TTree* ntu = (TTree*)file->Get("ntu");
  ntu->SetBranchAddress("eventNb",&eventNb);
  ntu->SetBranchAddress("mu_x_hit",&mu_x_hit);
  ntu->SetBranchAddress("mu_y_hit",&mu_y_hit);
  ntu->SetBranchAddress("mu_Edep",&mu_Edep);
  ntu->SetBranchAddress("mu_TrackLength",&mu_TrackLength);
  ntu->SetBranchAddress("Ph_tot_Nb",&Ph_tot_Nb);
  ntu->SetBranchAddress("Ph_detected_Nb",&Ph_detected_Nb);
  ntu->SetBranchAddress("Ph_lost_Nb",&Ph_lost_Nb);
  ntu->SetBranchAddress("Ph_labs_Nb",&Ph_labs_Nb);
  ntu->SetBranchAddress("Ph_time",&Ph_time);
  //ntu->SetBranchAddress("Ph_energy",&Ph_energy);
  ntu->SetBranchAddress("Ph_creator_proc",&Ph_creator_proc);
  //ntu->SetBranchAddress("Ph_x_hit",&Ph_x_hit);
  //ntu->SetBranchAddress("Ph_y_hit",&Ph_y_hit);
  //ntu->SetBranchAddress("Ph_SiPM_number",&Ph_SiPM_number);

/*
  TBranch* bPh_time = ntu->GetBranch("Ph_time");
  TBranch* bPh_energy = ntu->GetBranch("Ph_energy");
  TBranch* bPh_creator_proc = ntu->GetBranch("Ph_creator_proc");
  TBranch* bPh_x_hit = ntu->GetBranch("Ph_x_hit");
  TBranch* bPh_y_hit = ntu->GetBranch("Ph_y_hit");
  TBranch* bmu_x_hit = ntu->GetBranch("mu_x_hit");
  TBranch* bmu_y_hit = ntu->GetBranch("mu_y_hit");

  bPh_time->SetAddress(&Ph_time);
  bPh_energy->SetAddress(&Ph_energy);
  bPh_creator_proc->SetAddress(&Ph_creator_proc);
  bPh_x_hit->SetAddress(&Ph_x_hit);
  bPh_y_hit->SetAddress(&Ph_y_hit);
  bmu_x_hit->SetAddress(&mu_x_hit);
  bmu_y_hit->SetAddress(&mu_y_hit);
*/
  //Define some variables for the signal sampling 
  const int Nsample = 4000;
  Float_t V[Nsample];
  Float_t t[Nsample];
  Float_t LDE2;
  Float_t LDE5;
  Float_t LDE10;
  Float_t LDE20;
  Float_t LDE50;
  Float_t LDE100;
  Float_t PH2;
  Float_t PH5;
  Float_t PH10;
  Float_t PH20;
  Float_t PH50;
  Float_t PH100;
  Float_t AMP_MAX;
  Float_t baseline;


//Define the output Tree
  TFile *outfile;
  TString outfilename(filename);

  if(cherenkov == true)
     outfilename.Insert(outfilename.Last('/')+1,Form("digi_DCR%fMHz_",rate*1000));
  else
     outfilename.Insert(outfilename.Last('/')+1,Form("digi_NOcherenkov_DCR%fMHz_",rate*1000));
  
  outfilename.Remove(0,outfilename.Last('/')+1);
  outfile=new TFile(outfilename.Data(),"RECREATE");
  TTree *outtree=new TTree("digi","Digitized SiPM signal");
  outtree->Branch("eventNb",&eventNb,"eventNb/I");
  outtree->Branch("mu_x_hit",&mu_x_hit,"mu_x_hit/F");
  outtree->Branch("mu_y_hit",&mu_y_hit,"mu_y_hit/F");
  outtree->Branch("mu_Edep",&mu_Edep,"mu_Edep/F");
  outtree->Branch("mu_TrackLength",&mu_TrackLength,"mu_TrackLength/F");
  outtree->Branch("Ph_tot_Nb",&Ph_tot_Nb,"Ph_tot_Nb/I");
  outtree->Branch("Ph_detected_Nb",&Ph_detected_Nb,"Ph_detected_Nb/I");
  outtree->Branch("Ph_lost_Nb",&Ph_lost_Nb,"Ph_lost_Nb/I");
  outtree->Branch("Ph_labs_Nb",&Ph_labs_Nb,"Ph_labs_Nb/I");
  outtree->Branch("LDE2",&LDE2,"LDE2/F");
  outtree->Branch("LDE5",&LDE5,"LDE5/F");
  outtree->Branch("LDE10",&LDE10,"LDE10/F");
  outtree->Branch("LDE20",&LDE20,"LDE20/F");
  outtree->Branch("LDE50",&LDE50,"LDE50/F");
  outtree->Branch("LDE100",&LDE100,"LDE100/F");

  outtree->Branch("PH2",&PH2,"PH2/F");
  outtree->Branch("PH5",&PH5,"PH5/F");
  outtree->Branch("PH10",&PH10,"PH10/F");
  outtree->Branch("PH20",&PH20,"PH20/F");
  outtree->Branch("PH50",&PH50,"PH50/F");
  outtree->Branch("PH100",&PH100,"PH100/F");

  outtree->Branch("AMP_MAX",&AMP_MAX,"AMP_MAX/F");
  outtree->Branch("baseline",&baseline,"baseline/F");

  outtree->Branch("WFval",&V,("WFval[" + TString::Itoa(Nsample,10) +"]/F").Data());
  outtree->Branch("WFtime",&t,("WFtime[" + TString::Itoa(Nsample,10) +"]/F").Data());

  //outtree->Branch("mu_x_hit",&mu_x_hit,"mu_x_hit/F");
  //outtree->Branch("mu_y_hit",&mu_y_hit,"mu_y_hit/F");

  TH1F *histo_time = new TH1F("Photon arrival time","Photon arrival time",100,0.,20.);
  TH1F *histo_timePDE = new TH1F("Photon PDE selected arrival time","Photon PDE selected arrival time",100,0.,20.);
  cout<<"Reading ntu"<<endl;
  Long64_t nentries = ntu->GetEntries();
  for (Long64_t i_entry=0;i_entry<nentries;i_entry++)
  {
      //if(i_entry%5==0)
        cout<<"Reading entry "<<i_entry<< "\r" << std::flush;
/*
      bPh_time->GetEntry(i_entry);
      bPh_energy->GetEntry(i_entry);
      bPh_x_hit->GetEntry(i_entry);
      bPh_y_hit->GetEntry(i_entry);
      bmu_x_hit->GetEntry(i_entry);
      bmu_y_hit->GetEntry(i_entry);
      bPh_creator_proc->GetEntry(i_entry);
      */
      ntu->GetEntry(i_entry);
      //if(i_entry%100==0)
         //for(std::vector<Float_t>::iterator it = Ph_time->begin() ; it != Ph_time->end(); ++it)
           // if(*it>0.)
              // histo_time->Fill(*it);
      Fill(Ph_event, Ph_time, /*Ph_energy, Ph_x_hit, Ph_y_hit, */Ph_creator_proc,cherenkov); 
      //PDERandomDeleting(gen,Ph_event);//Delete randomly photons depending on the SiPM PDE (except the FILLL FACTOR that is already accounted)
      if(Ph_event.size()>=2)
        PH2 = Ph_event.at(1).time;
      else
        PH2 = -1;
      if(Ph_event.size()>=5)
        PH5 = Ph_event.at(4).time;
      else
        PH5 = -1;
      if(Ph_event.size()>=10)
        PH10 = Ph_event.at(9).time;
      else
        PH10 = -1;
      if(Ph_event.size()>=20)
        PH20 = Ph_event.at(19).time;
      else
        PH20 = -1;
      if(Ph_event.size()>=50)
        PH50 = Ph_event.at(49).time;
      else
        PH50 = -1;
      if(Ph_event.size()>=100)
        PH100 = Ph_event.at(99).time;
      else
        PH100 = -1;
         
      //Position_selection(Ph_event,SiPMLength);//Delete photon that doesn't hit SiPM active cells   

      //if(i_entry%100==0)
         /*for(std::vector<photon>::iterator it = Ph_event.begin() ; it != Ph_event.end(); ++it)
            histo_timePDE->Fill(it->time);*/

      SPTR(gen,Ph_event);//Modify randomly the photon arrival time because of the SiPM time jitter
      Dark_current(gen,Ph_event,rate);//Add dark current rate
      Digitize(gen,Ph_event,V,t,Nsample);//Get the sampled voltage signal (a.u.)
      baseline = TMath::Mean(1000, V);
      AMP_MAX = TMath::MaxElement(4000,V)-baseline;//Get the signal max amplitude
      LDE2 = TimeOverThresh(V,t,2.-0.5+baseline);//Get the time over threshold
      LDE5 = TimeOverThresh(V,t,5.-0.5+baseline);//Get the time over threshold
      LDE10 = TimeOverThresh(V,t,10.-0.5+baseline);//Get the time over threshold
      LDE20 = TimeOverThresh(V,t,20.-0.5+baseline);//Get the time over threshold
      LDE50 = TimeOverThresh(V,t,50.-0.5+baseline);//Get the time over threshold
      LDE100 = TimeOverThresh(V,t,100.-0.5+baseline);//Get the time over threshold
      /*if(i_entry%1000==0)
      {
         c1=Draw(Ph_event);
         c1->Print((TString::Itoa(i_entry,10) + ".pdf").Data());
         c2->cd();
         TGraph *WF = new TGraph(1000,t,V);
         WF->Draw("APL");
      }*/

      outtree->Fill();//Fill the output ntuple
      Ph_event.clear();
  }
  
  //TCanvas *c2 = new TCanvas();
  //c2->Divide(2,1);
  //c2->cd();
  //histo_time->Scale(1./nentries*100.);
  //histo_time->Draw();
  //histo_timePDE->Scale(1./nentries*100.);
  //histo_timePDE->Draw("SAME");
  //c2->BuildLegend();
  //c2->cd(2);
  //histo_time->GetCumulative()->Draw();
  //histo_timePDE->GetCumulative()->Draw("SAME");
  //file->Close();
  outtree->Write();
  //outtree->Print();
  outfile->Close();
}
