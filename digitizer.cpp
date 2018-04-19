
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

   int cherenkov = 1; // 0->only scintillation , 1->cherenkov+scintillation , 2->only cherenkov 
   if (config.keyExists("cherenkov"))
     cherenkov = config.read<int> ("cherenkov");
   if(cherenkov<0 || cherenkov>2)
   {
      cout<<"ERROR: cherenkov in configfile has an invalid value\n";
      exit(EXIT_FAILURE);
   }

   long LoopOnEntry = -1;
   if (config.keyExists("LoopOnEntry"))
     LoopOnEntry = config.read<long> ("LoopOnEntry");

   long Nloop = -1;
   if (config.keyExists("Nloop"))
     Nloop = config.read<long> ("Nloop");


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
  if(!file)
  {
     cout<<"ERROR: file "<<filename<<"NOT exist"<<endl;
     exit(EXIT_FAILURE);     
  }
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

  if(cherenkov == 1)
     outfilename.Insert(outfilename.Last('/')+1,Form("digi_DCR%fMHz_",rate*1000));
  else
     if(cherenkov == 0)
        outfilename.Insert(outfilename.Last('/')+1,Form("digi_ONLYSCINT_DCR%fMHz_",rate*1000));
     else
        if(cherenkov == 2)
           outfilename.Insert(outfilename.Last('/')+1,Form("digi_ONLYCHERENKOV_DCR%fMHz_",rate*1000));
  
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
  if(Nloop>nentries && LoopOnEntry==-1)
  {
    cout<<"ERROR: Nloop>nentries, check configfile"<<endl;
    exit(EXIT_FAILURE);
  }
  else
    if(Nloop==-1)
      Nloop=nentries;
  for (Long64_t i_entry=0;i_entry<Nloop;i_entry++)
  {
      if(i_entry%5==0)
        cout<<"Loop "<<i_entry<< "\r" << std::flush;
      if(LoopOnEntry==-1) 
	ntu->GetEntry(i_entry);
      else
	ntu->GetEntry(LoopOnEntry);
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
      LDE2 = LeadEdgeDiscr(V,t,2.-0.5+baseline);//Get the time over threshold
      LDE5 = LeadEdgeDiscr(V,t,5.-0.5+baseline);//Get the time over threshold
      LDE10 = LeadEdgeDiscr(V,t,10.-0.5+baseline);//Get the time over threshold
      LDE20 = LeadEdgeDiscr(V,t,20.-0.5+baseline);//Get the time over threshold
      LDE50 = LeadEdgeDiscr(V,t,50.-0.5+baseline);//Get the time over threshold
      LDE100 = LeadEdgeDiscr(V,t,100.-0.5+baseline);//Get the time over threshold
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
