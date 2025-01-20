#define MyClass_mod_cxx
#include "MyClass_mod.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cstring>

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
  int ENbins = 1000;
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
        hParticleHist[i][j] = new TH1D(Form("h%s_%s",particles[i], properties[j]), "", 750, -40, 40);
        } else {
            hParticleHist[i][j] = new TH1D(Form("h%s_%s",particles[i], properties[j]), "", 750, -1, 1);         
        }
     }      
  }

  if (fChain == 0) return;

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
	
	  int Npip = 0; //number of pi+
	  int Npim = 0; //number of pi-
	  int Npi0 = 0; //number of pi0

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

	    //Finding the index of the cross-section type
	    int idx = -1;
	    if(Npip==0 && Npim==0 && Npi0==0) idx = 1;
	    if(Npip==1 && Npim==0 && Npi0==0) idx = 2;
	    if(Npip==0 && Npim==1 && Npi0==0) idx = 3;
	    if(Npip==0 && Npim==0 && Npi0==1) idx = 4;
	    if( (Npip + Npim + Npi0) > 1    ) idx = 5;

	    if(idx>0)hxsec[idx]->Fill(Ein,wgt);
	 
    }
  }


  //Momentum histograms
  THStack *hstack = new THStack();

  TH1D* Mhist = new TH1D("Mhist", "Momentum Histograms", 100, 0, 1);
  
  Mhist->Draw();

  TCanvas *c1 = new TCanvas("c1", "Total Momentum", 800, 600);
  
  c1->SetGrid();

  hParticleHist[0][0]->SetLineColor(kRed);
  hParticleHist[1][0]->SetLineColor(kBlue);
  hParticleHist[2][0]->SetLineColor(kGreen);
  hParticleHist[3][0]->SetLineColor(kMagenta);
  hParticleHist[4][0]->SetLineColor(kOrange);

  hstack->Add(hParticleHist[0][0]);
  hstack->Add(hParticleHist[1][0]);
  hstack->Add(hParticleHist[2][0]);
  hstack->Add(hParticleHist[3][0]);
  hstack->Add(hParticleHist[4][0]);
  hstack->Draw("nostack");

  TLegend* legend1 = new TLegend(0.7,0.7,0.9,0.9);
  legend1->AddEntry(hParticleHist[0][0], "Proton","l");
  legend1->AddEntry(hParticleHist[1][0], "Neutron", "l");
  legend1->AddEntry(hParticleHist[2][0], "+Pion","l");
  legend1->AddEntry(hParticleHist[3][0], "Pi0","l");
  legend1->AddEntry(hParticleHist[4][0], "-Pion","l");
  legend1->Draw();

  hstack->GetXaxis()->SetTitle("GeV/c");
  hstack->GetYaxis()->SetTitle("Entries");

  c1->Update();

  c1->SaveAs("Momentum_histograms.png");

  //Px
  THStack *pxstack = new THStack();
  
  TH1D* pxhist = new TH1D("Pxhist", "Px Histograms", 100, -1, 1); 
  
  pxhist->Draw();  

  TCanvas *c2 = new TCanvas("c2", "Px", 800, 600);

  c2->SetGrid();

  hParticleHist[0][1]->SetLineColor(kRed);
  hParticleHist[1][1]->SetLineColor(kBlue);
  hParticleHist[2][1]->SetLineColor(kGreen);
  hParticleHist[3][1]->SetLineColor(kMagenta);
  hParticleHist[4][1]->SetLineColor(kOrange);

  pxstack->Add(hParticleHist[0][1]);
  pxstack->Add(hParticleHist[1][1]);
  pxstack->Add(hParticleHist[2][1]);
  pxstack->Add(hParticleHist[3][1]);
  pxstack->Add(hParticleHist[4][1]);
  pxstack->Draw("nostack");

  TLegend* legend2 = new TLegend(0.7,0.7,0.9,0.9);
  legend2->AddEntry(hParticleHist[0][1], "Proton","l");
  legend2->AddEntry(hParticleHist[1][1], "Neutron", "l");
  legend2->AddEntry(hParticleHist[2][1], "+Pion","l");
  legend2->AddEntry(hParticleHist[3][1], "Pi0","l");
  legend2->AddEntry(hParticleHist[4][1], "-Pion","l");
  legend2->Draw();

  pxstack->GetXaxis()->SetTitle("GeV/c");
  pxstack->GetYaxis()->SetTitle("Entries");

  c2->Update();

  c2->SaveAs("Px.png");
  
  //py
  THStack *pystack = new THStack();
  
  TH1D* pyhist = new TH1D("Pyhist", "Py Histograms", 100, -1, 1); 
  
  pyhist->Draw();

  TCanvas *c3 = new TCanvas("c3", "Py", 800, 600);

  c3->SetGrid();

  hParticleHist[0][2]->SetLineColor(kRed);
  hParticleHist[1][2]->SetLineColor(kBlue);
  hParticleHist[2][2]->SetLineColor(kGreen);
  hParticleHist[3][2]->SetLineColor(kMagenta);
  hParticleHist[4][2]->SetLineColor(kOrange);

  pystack->Add(hParticleHist[0][2]);
  pystack->Add(hParticleHist[1][2]);
  pystack->Add(hParticleHist[2][2]);
  pystack->Add(hParticleHist[3][2]);
  pystack->Add(hParticleHist[4][2]);
  pystack->Draw("nostack");

  TLegend* legend3 = new TLegend(0.7,0.7,0.9,0.9);
  legend3->AddEntry(hParticleHist[0][2], "Proton","l");
  legend3->AddEntry(hParticleHist[1][2], "Neutron", "l");
  legend3->AddEntry(hParticleHist[2][2], "+Pion","l");
  legend3->AddEntry(hParticleHist[3][2], "Pi0","l");
  legend3->AddEntry(hParticleHist[4][2], "-Pion","l");
  legend3->Draw();

  pystack->GetXaxis()->SetTitle("GeV/c");
  pystack->GetYaxis()->SetTitle("Entries");

  c3->Update();

  c3->SaveAs("Py.png");
  
  //Pz
  THStack *pzstack = new THStack();
  
  TH1D* pzhist = new TH1D("Pzhist", "Pz Histograms", 100, -1, 1); 
  
  pzhist->Draw();

  TCanvas *c4 = new TCanvas("c4", "Pz", 800, 600);

  c4->SetGrid();

  hParticleHist[0][3]->SetLineColor(kRed);
  hParticleHist[1][3]->SetLineColor(kBlue);
  hParticleHist[2][3]->SetLineColor(kGreen);
  hParticleHist[3][3]->SetLineColor(kMagenta);
  hParticleHist[4][3]->SetLineColor(kOrange);

  pzstack->Add(hParticleHist[0][3]);
  pzstack->Add(hParticleHist[1][3]);
  pzstack->Add(hParticleHist[2][3]);
  pzstack->Add(hParticleHist[3][3]);
  pzstack->Add(hParticleHist[4][3]);
  pzstack->Draw("nostack");

  TLegend* legend4 = new TLegend(0.7,0.7,0.9,0.9);
  legend4->AddEntry(hParticleHist[0][1], "Proton","l");
  legend4->AddEntry(hParticleHist[1][1], "Neutron", "l");
  legend4->AddEntry(hParticleHist[2][1], "+Pion","l");
  legend4->AddEntry(hParticleHist[3][1], "Pi0","l");
  legend4->AddEntry(hParticleHist[4][1], "-Pion","l");
  legend4->Draw();

  pzstack->GetXaxis()->SetTitle("GeV/c");
  pzstack->GetYaxis()->SetTitle("Entries");

  c4->Update();

  c4->SaveAs("Pz.png");
  
  //Vx
  THStack *vxstack = new THStack();
  
  TH1D* vxhist = new TH1D("Vxhist", "Vx Histograms", 100, -1, 1); 
  
  vxhist->Draw();

  TCanvas *c5 = new TCanvas("c5", "Vx", 800, 600);

  c5->SetGrid();

  hParticleHist[0][4]->SetLineColor(kRed);
  hParticleHist[1][4]->SetLineColor(kBlue);
  hParticleHist[2][4]->SetLineColor(kGreen);
  hParticleHist[3][4]->SetLineColor(kMagenta);
  hParticleHist[4][4]->SetLineColor(kOrange);

  vxstack->Add(hParticleHist[0][4]);
  vxstack->Add(hParticleHist[1][4]);
  vxstack->Add(hParticleHist[2][4]);
  vxstack->Add(hParticleHist[3][4]);
  vxstack->Add(hParticleHist[4][4]);
  vxstack->Draw("nostack");

  TLegend* legend5 = new TLegend(0.7,0.7,0.9,0.9);
  legend5->AddEntry(hParticleHist[0][4], "Proton","l");
  legend5->AddEntry(hParticleHist[1][4], "Neutron", "l");
  legend5->AddEntry(hParticleHist[2][4], "+Pion","l");
  legend5->AddEntry(hParticleHist[3][4], "Pi0","l");
  legend5->AddEntry(hParticleHist[4][4], "-Pion","l");
  legend5->Draw();

  vxstack->GetXaxis()->SetTitle("");
  vxstack->GetYaxis()->SetTitle("Entries");

  c5->Update();

  c5->SaveAs("Vx.png");
  
  //Vy
  THStack *vystack = new THStack();
  
  TH1D* vyhist = new TH1D("Vyhist", "Vy Histograms", 100, -1, 1); 
  
  vyhist->Draw();

  TCanvas *c6 = new TCanvas("c6", "Vy", 800, 600);

  c6->SetGrid();

  hParticleHist[0][5]->SetLineColor(kRed);
  hParticleHist[1][5]->SetLineColor(kBlue);
  hParticleHist[2][5]->SetLineColor(kGreen);
  hParticleHist[3][5]->SetLineColor(kMagenta);
  hParticleHist[4][5]->SetLineColor(kOrange);

  vystack->Add(hParticleHist[0][5]);
  vystack->Add(hParticleHist[1][5]);
  vystack->Add(hParticleHist[2][5]);
  vystack->Add(hParticleHist[3][5]);
  vystack->Add(hParticleHist[4][5]);
  vystack->Draw("nostack");

  TLegend* legend6 = new TLegend(0.7,0.7,0.9,0.9);
  legend6->AddEntry(hParticleHist[0][5], "Proton","l");
  legend6->AddEntry(hParticleHist[1][5], "Neutron", "l");
  legend6->AddEntry(hParticleHist[2][5], "+Pion","l");
  legend6->AddEntry(hParticleHist[3][5], "Pi0","l");
  legend6->AddEntry(hParticleHist[4][5], "-Pion","l");
  legend6->Draw();

  vystack->GetXaxis()->SetTitle("");
  vystack->GetYaxis()->SetTitle("Entries");

  c6->Update();

  c6->SaveAs("Vy.png");
  
  //Vz
  THStack *vzstack = new THStack();
  
  TH1D* vzhist = new TH1D("Vzhist", "Vz Histograms", 100, -1, 1); 
  
  vzhist->Draw();

  TCanvas *c7 = new TCanvas("c7", "Vz", 800, 600);

  c7->SetGrid();

  hParticleHist[0][6]->SetLineColor(kRed);
  hParticleHist[1][6]->SetLineColor(kBlue);
  hParticleHist[2][6]->SetLineColor(kGreen);
  hParticleHist[3][6]->SetLineColor(kMagenta);
  hParticleHist[4][6]->SetLineColor(kOrange);

  vzstack->Add(hParticleHist[0][6]);
  vzstack->Add(hParticleHist[1][6]);
  vzstack->Add(hParticleHist[2][6]);
  vzstack->Add(hParticleHist[3][6]);
  vzstack->Add(hParticleHist[4][6]);
  vzstack->Draw("nostack");

  TLegend* legend7 = new TLegend(0.7,0.7,0.9,0.9);
  legend7->AddEntry(hParticleHist[0][6], "Proton","l");
  legend7->AddEntry(hParticleHist[1][6], "Neutron", "l");
  legend7->AddEntry(hParticleHist[2][6], "+Pion","l");
  legend7->AddEntry(hParticleHist[3][6], "Pi0","l");
  legend7->AddEntry(hParticleHist[4][6], "-Pion","l");
  legend7->Draw();

  vzstack->GetXaxis()->SetTitle("");
  vzstack->GetYaxis()->SetTitle("Entries");

  c7->Update();

  c7->SaveAs("Vz.png");
  
  //Vx pions
  THStack *vx2stack = new THStack();
  
  TH1D* vx2hist = new TH1D("Vx2hist", "Vx2 Histograms", 100, -1, 1); 
  
  vx2hist->Draw();

  TCanvas *c8 = new TCanvas("c8", "Vx2", 800, 600);

  c8->SetGrid();

  //hParticleHist[0][6]->SetLineColor(kRed);
  //hParticleHist[1][6]->SetLineColor(kBlue);
  hParticleHist[2][4]->SetLineColor(kGreen);
  hParticleHist[3][4]->SetLineColor(kMagenta);
  hParticleHist[4][4]->SetLineColor(kOrange);

  //vzstack->Add(hParticleHist[0][6]);
  //vzstack->Add(hParticleHist[1][6]);
  vx2stack->Add(hParticleHist[2][4]);
  vx2stack->Add(hParticleHist[3][4]);
  vx2stack->Add(hParticleHist[4][4]);
  vx2stack->Draw("nostack");

  TLegend* legend8 = new TLegend(0.7,0.7,0.9,0.9);
  //legend8->AddEntry(hParticleHist[0][4], "Proton","l");
  //legend8->AddEntry(hParticleHist[1][4], "Neutron", "l");
  legend8->AddEntry(hParticleHist[2][4], "+Pion","l");
  legend8->AddEntry(hParticleHist[3][4], "Pi0","l");
  legend8->AddEntry(hParticleHist[4][4], "-Pion","l");
  legend8->Draw();

  vx2stack->GetXaxis()->SetTitle("");
  vx2stack->GetYaxis()->SetTitle("Entries");

  c8->Update();

  c8->SaveAs("Vx2.png");
  
  //Vy pions
  THStack *vy2stack = new THStack();
  
  TH1D* vy2hist = new TH1D("Vy2hist", "Vy2 Histograms", 100, -1, 1); 
  
  vy2hist->Draw();

  TCanvas *c9 = new TCanvas("c9", "Vy2", 800, 600);

  c9->SetGrid();

  //hParticleHist[0][6]->SetLineColor(kRed);
  //hParticleHist[1][6]->SetLineColor(kBlue);
  hParticleHist[2][5]->SetLineColor(kGreen);
  hParticleHist[3][5]->SetLineColor(kMagenta);
  hParticleHist[4][5]->SetLineColor(kOrange);

  //vzstack->Add(hParticleHist[0][6]);
  //vzstack->Add(hParticleHist[1][6]);
  vy2stack->Add(hParticleHist[2][5]);
  vy2stack->Add(hParticleHist[3][5]);
  vy2stack->Add(hParticleHist[4][5]);
  vy2stack->Draw("nostack");

  TLegend* legend9 = new TLegend(0.7,0.7,0.9,0.9);
  //legend9->AddEntry(hParticleHist[0][5], "Proton","l");
  //legend9->AddEntry(hParticleHist[1][5], "Neutron", "l");
  legend9->AddEntry(hParticleHist[2][5], "+Pion","l");
  legend9->AddEntry(hParticleHist[3][5], "Pi0","l");
  legend9->AddEntry(hParticleHist[4][5], "-Pion","l");
  legend9->Draw();

  vy2stack->GetXaxis()->SetTitle("");
  vy2stack->GetYaxis()->SetTitle("Entries");

  c9->Update();

  c9->SaveAs("Vy2.png"); 
  
  //Vz pions only
  THStack *vz2stack = new THStack();
  
  TH1D* vz2hist = new TH1D("Vz2hist", "Vz2 Histograms", 100, -1, 1); 
  
  vz2hist->Draw();

  TCanvas *c10 = new TCanvas("c10", "Vz2", 800, 600);

  c10->SetGrid();

  //hParticleHist[0][6]->SetLineColor(kRed);
  //hParticleHist[1][6]->SetLineColor(kBlue);
  hParticleHist[2][6]->SetLineColor(kGreen);
  hParticleHist[3][6]->SetLineColor(kMagenta);
  hParticleHist[4][6]->SetLineColor(kOrange);

  //vzstack->Add(hParticleHist[0][6]);
  //vzstack->Add(hParticleHist[1][6]);
  vz2stack->Add(hParticleHist[2][6]);
  vz2stack->Add(hParticleHist[3][6]);
  vz2stack->Add(hParticleHist[4][6]);
  vz2stack->Draw("nostack");

  TLegend* legend10 = new TLegend(0.7,0.7,0.9,0.9);
  //legend9->AddEntry(hParticleHist[0][4], "Proton","l");
  //legend9->AddEntry(hParticleHist[1][4], "Neutron", "l");
  legend10->AddEntry(hParticleHist[2][6], "+Pion","l");
  legend10->AddEntry(hParticleHist[3][6], "Pi0","l");
  legend10->AddEntry(hParticleHist[4][6], "-Pion","l");
  legend10->Draw();

  vz2stack->GetXaxis()->SetTitle("");
  vz2stack->GetYaxis()->SetTitle("Entries");

  c10->Update();

  c10->SaveAs("Vz2.png");  


  //Storing cross-section wanted histograms:
  fOut->mkdir("xsec");
  fOut->cd("xsec");
    
  for(int i=0;i<NintXsec;i++){
    hxsec[i]->Scale(1./double(fNtrees));
    hxsec[i]->Write();
  }

  fOut->cd();    
  fOut->mkdir("Energy");
  fOut->cd("Energy");

  for(int i=0;i<ENint;i++){
    hE[i]->Scale(1./double(fNtrees));
    hE[i]->Write();
  }

  fOut->cd();

  fOut->mkdir("Components");
  fOut->cd("Components");
  for(int i =0; i < Npart; i++) {
      for (int j = 0; j < Ncomponents; j++){
      hParticleHist[i][j]->Write();
    }
  }

  fOut->cd();

  hwgt->Write();
  hEin->Write();
  hnparts->Write();
  hPDG->Write();
  hinter->Write();
  hhist->Write();
    
    


  fOut->Close();
  
}

    





