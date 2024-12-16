#define MyClass_mod_cxx
#include "MyClass_mod.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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
  double Emine = -0.0005;
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
  for(int i =0; i < Npart; ++i) {
    for (int j = 0; j < Ncomponents; ++j) {
      hParticleHist[i][j] = new TH1D(Form("h%s_%s",particles[i], properties[j]), "", 750, -10, 10);
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
      if (pdg[i] == 2212 && E[i] > 10.0) idy = 0;   //Proton
      else if (pdg[i] == 2112) idy = 1;  //Neutron
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

  //momentum histograms
  TH1D* Mhist = new TH1D("Mhist", "Momentum Histograms", 100, -10, 10);

  TCanvas *c1 = new TCanvas("c1", "Total Momentum", 800, 600);
  
  Mhist->Draw();
  
  c1->SetGrid();

  hParticleHist[0][0]->SetLineColor(kRed);
  hParticleHist[1][0]->SetLineColor(kBlue);
  hParticleHist[2][0]->SetLineColor(kGreen);
  hParticleHist[3][0]->SetLineColor(kMagenta);
  hParticleHist[4][0]->SetLineColor(kOrange);

  hParticleHist[0][0]->Draw("HIST");
  hParticleHist[1][0]->Draw("HIST SAME");
  hParticleHist[2][0]->Draw("HIST SAME");
  hParticleHist[3][0]->Draw("HIST SAME");
  hParticleHist[4][0]->Draw("HIST SAME");

  TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->AddEntry(hParticleHist[0][0], "Proton","l");
  legend->AddEntry(hParticleHist[1][0], "Neutron", "l");
  legend->AddEntry(hParticleHist[2][0], "+Pion","l");
  legend->AddEntry(hParticleHist[3][0], "Pi0","l");
  legend->AddEntry(hParticleHist[4][0], "-Pion","l");
  legend->Draw();

  Mhist->GetXaxis()->SetTitle("GeV/c");
  Mhist->GetYaxis()->SetTitle("Entries");

  c1->Update();

  c1->SaveAs("Momentum_histograms.png");

  //Px
  TH1D* Pxhist = new TH1D("Pxhist", "Px Histograms", 100, -10, 10); 

  TCanvas *c2 = new TCanvas("c2", "Px", 800, 600);

  Pxhist->Draw();

  c2->SetGrid();

  hParticleHist[0][1]->SetLineColor(kRed);
  hParticleHist[1][1]->SetLineColor(kBlue);
  hParticleHist[2][1]->SetLineColor(kGreen);
  hParticleHist[3][1]->SetLineColor(kMagenta);
  hParticleHist[4][1]->SetLineColor(kOrange);

  hParticleHist[0][1]->Draw("HIST");
  hParticleHist[1][1]->Draw("HIST SAME");
  hParticleHist[2][1]->Draw("HIST SAME");
  hParticleHist[3][1]->Draw("HIST SAME");
  hParticleHist[4][1]->Draw("HIST SAME");

  TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->AddEntry(hParticleHist[0][1], "Proton","l");
  legend->AddEntry(hParticleHist[1][1], "Neutron", "l");
  legend->AddEntry(hParticleHist[2][1], "+Pion","l");
  legend->AddEntry(hParticleHist[3][1], "Pi0","l");
  legend->AddEntry(hParticleHist[4][1], "-Pion","l");
  legend->Draw();

  Pxhist->GetXaxis()->SetTitle("GeV/c");
  Pxhist->GetYaxis()->SetTitle("Entries");

  c2->Update();

  c2->SaveAs("Px.png");


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
  for(int i =0; i < Npart; ++i) {
      for (int j = 0; j < Ncomponents; ++j){
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

    




