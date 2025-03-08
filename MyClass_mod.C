#define MyClass_mod_cxx
#include "MyClass_mod.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cstring>
#include <vector>
#include <cmath>
#include <iostream>

void PlotHistograms(THStack *hstack, const char* canvasName, const char* title, const char* saveName, TH1D* hParticleHist[5][7], int histIndex, int startIdx, int endIdx) {
    TCanvas *canvas = new TCanvas(canvasName, title, 800, 600);
    canvas->SetGrid();
    
    
    hParticleHist[0][histIndex]->SetLineColor(kRed);
    hParticleHist[1][histIndex]->SetLineColor(kBlue);
    hParticleHist[2][histIndex]->SetLineColor(kGreen);
    hParticleHist[3][histIndex]->SetLineColor(kMagenta);
    hParticleHist[4][histIndex]->SetLineColor(kOrange);
    
    
    for (int i = startIdx; i < endIdx; i++) {
    hstack->Add(hParticleHist[i][histIndex]);
    }
    
    hstack->Draw("nostack");
    
    
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hParticleHist[0][histIndex], "Proton", "l");
    legend->AddEntry(hParticleHist[1][histIndex], "Neutron", "l");
    legend->AddEntry(hParticleHist[2][histIndex], "+Pion", "l");
    legend->AddEntry(hParticleHist[3][histIndex], "Pi0", "l");
    legend->AddEntry(hParticleHist[4][histIndex], "-Pion", "l");
    legend->Draw();
    
    hstack->GetXaxis()->SetTitle("GeV/c");
    hstack->GetYaxis()->SetTitle("Entries");
    
    canvas->Update();
    canvas->SaveAs(saveName);
    }
    
    
void AutomateHistograms() {
    THStack *hstack = new THStack();
    THStack *pxstack = new THStack();
    THStack *pystack = new THStack();
    THStack *pzstack = new THStack();
    THStack *vxstack = new THStack();
    THStack *vystack = new THStack();
    THStack *vzstack = new THStack();
    
    TH1D* hParticleHist[5][7]; 
    
    
    PlotHistograms(hstack, "Total Momentum", "Momentum Histograms", "Momentum_histograms.png", hParticleHist, 0, 0, 5);
    PlotHistograms(pxstack, "Px", "Px Histograms", "Px.png", hParticleHist, 1, 0, 5);
    PlotHistograms(pystack, "Py", "Py Histograms", "Py.png", hParticleHist, 2, 0, 5);
    PlotHistograms(pzstack, "Pz", "Pz Histograms", "Pz.png", hParticleHist, 3, 0, 5);
    PlotHistograms(vxstack, "Vx", "Vx Histograms", "Vx.png", hParticleHist, 4, 0, 5);
    PlotHistograms(vystack, "Vy", "Vy Histograms", "Vy.png", hParticleHist, 5, 0, 5);
    PlotHistograms(vzstack, "Vz", "Vz Histograms", "Vz.png", hParticleHist, 6, 0, 5);
    }

void MyClass_mod::Loop() {
  
  TH1D* hwgt = new TH1D("hwgt", "wgt", 100, -1, 2);
  TH1D* hEin = new TH1D("hEin", "Ein", 100, -1, 2);
  TH1D* hnparts = new TH1D("hnparts", "", 100, -1, 50);
  TH1D* hPDG = new TH1D("hPDG", "PDG", 100, -1000, 2500);
  TH1D* hinter = new TH1D("hinter", "", 100, -3000, 250000);
  TH1D* hhist =new TH1D("hhist", "", 100,-1, 2); 

  const int Npart = 5;
  const char* particles[Npart] = {"proton", "neutron", "pip", "pi0", "pim"};
  const char* properties[] = {"Tmom", "px", "py", "pz", "vx", "vy", "vz"};

    //Energy histograms
  const int ENint = 5;
  const char* Ecint[ENint] = {"proton","neutron","pip","pi0","pim"};
  int ENbins = 1002;
  double Emine = 0;
  double Emaxe = 2;
  TH1D* hE[ENint];
  for(int i=0;i<ENint;i++){
    hE[i] = new TH1D(Form("hE_%s", Ecint[i]),"",ENbins,Emine,Emaxe);
  }


    //New histograms:
  const int NintXsec = 6;
  const char* cintXsec[NintXsec] = {"rea","abs","pip","pim","pi0","2pi"};    
  int NbinsXsec   =    1002;
  double mineXsec = -0.0005;
  double maxeXsec = 1.0015;
  TH1D* hxsec[NintXsec];
  for(int i=0;i<NintXsec;i++){
    hxsec[i] = new TH1D(Form("hxsec_%s",cintXsec[i]),"",NbinsXsec,mineXsec,maxeXsec);
  }

  const int Ncomponents = 7;
  TH1D* hParticleHist[Npart][Ncomponents];
  for(int i =0; i < Npart; i++) {
    for (int j = 0; j < Ncomponents; j++) {   
      if (strcmp(properties[j], "vx") == 0 || strcmp(properties[j], "vy") == 0 || strcmp(properties[j], "vz") == 0) {
        hParticleHist[i][j] = new TH1D(Form("h%s_%s",particles[i], properties[j]), "", 894, -40, 40);
        } else {
            hParticleHist[i][j] = new TH1D(Form("h%s_%s",particles[i], properties[j]), "", 894, -1, 1);         
        }
     }      
  }

  if (fChain == 0) return;

  TH1F *hProtonMomentum = new TH1F("hProtonMomentum", "Proton Momentum in Pion Absorption Events", 100, 0, 2);

  std::vector<std::pair<double, double>> protonData;  //(momentum, angle)

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if(inter[0] == -999){
	    if(nparts>1)std::cout<<"=> Warning (inter[0] == -999) events with nparts = "<<nparts<<std::endl;
	    continue;
	    }

	    //first, filling the reaction (total) cross section
	  hxsec[0]->Fill(Ein,wgt);

    //
	  int Npip = 0; //number of pi+
	  int Npim = 0; //number of pi-
	  int Npi0 = 0; //number of pi0

    //
    double highestEnergy = -1.0;
    double highestMomentum = 0.0;
    double highestAngle = 0.0;

	  for (int i = 0; i < nparts; i++) {   
      hwgt->Fill(wgt);
      hEin->Fill(Ein);
      hnparts->Fill(nparts);
      hPDG->Fill(pdg[i]);
      hinter->Fill(inter[i]);
      hhist->Fill(hist[i]);
       
      double totalmomentum = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
       
      int idy = -1;
      if (pdg[i] == 2212 && E[i] > 0.010) idy = 0;   //Proton
      else if (pdg[i] == 2112 && E[i] > 0.010) idy = 1;  //Neutron
      else if (pdg[i] == 211) idy = 2;   //+Pion
      else if (pdg[i] == 111) idy = 3;   //Pi0
      else if (pdg[i] == -211) idy = 4;   //-pion
       
      if (idy >= 0) {
      hE[idy]->Fill(E[i], wgt);
      hParticleHist[idy][0]->Fill(totalmomentum, wgt);
      hParticleHist[idy][1]->Fill(px[i], wgt);
      hParticleHist[idy][2]->Fill(py[i], wgt);
      hParticleHist[idy][3]->Fill(pz[i], wgt);
      hParticleHist[idy][4]->Fill(vx[i], wgt);
      hParticleHist[idy][5]->Fill(vy[i], wgt);
      hParticleHist[idy][6]->Fill(vz[i], wgt);

      }

      //counting pions
  	  if(pdg[i] ==  211) ++Npip;
  	  if(pdg[i] == -211) ++Npim; 
  	  if(pdg[i] ==  111) ++Npi0; 

      if (pdg[i] == 2212 && E[i] > highestEnergy) {
        highestEnergy = E[i];
        highestMomentum = totalmomentum;
        highestAngle = acos(pz[i] / highestMomentum);
      }
    }

    if (highestEnergy > 0 && Npip == 0 && Npim == 0 && Npi0 == 0) {
      hProtonMomentum->Fill(highestMomentum);
      //protonData.push_back(std::make_pair(highestMomentum,highestAngle))
    }
	    //Finding the index of the cross-section type
	    int idx = -1;
	    if(Npip==0 && Npim==0 && Npi0==0) idx = 1;
	    if(Npip==1 && Npim==0 && Npi0==0) idx = 2;
	    if(Npip==0 && Npim==1 && Npi0==0) idx = 3;
	    if(Npip==0 && Npim==0 && Npi0==1) idx = 4;
	    if( (Npip + Npim + Npi0) > 1    ) idx = 5;

	    if(idx>0)hxsec[idx]->Fill(Ein,wgt);
    
    
  }
  
  AutomateHistograms();
}

    





