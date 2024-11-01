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

    //Proton
    TH1D* hE_proton = new TH1D("hE_proton", "Proton Energy", 100, 0.75, 1.25);
    TH1D* hTmom_proton = new TH1D("hTmom_proton", "Total proton momentum", 100, 0, 1.25);
    TH1D* hpx_proton = new TH1D("hpx_proton", "Proton Momentum px",100,-1,1);
    TH1D* hpy_proton = new TH1D("hpy_proton", "Proton Momentum py",100,-1,1);
    TH1D* hpz_proton = new TH1D("hpz_proton", "Proton Momentum pz",100,-2,2);
    TH1D* hvx_proton = new TH1D("hvx_proton", "Proton Vertex x",100,-10,10);
    TH1D* hvy_proton = new TH1D("hvy_proton", "Proton Vertex y",100,-10,10);
    TH1D* hvz_proton = new TH1D("hvz_proton", "Proton Vertex z",100,-10,10);
    
    
    //Neutron
    TH1D* hE_neutron = new TH1D ("hE_neutron","Neutron Energy",100, 0.5, 2);
    TH1D* hTmom_neutron = new TH1D("hTmom_neutron", "Total neutron momentum", 100, 0, 2);
    TH1D* hpx_neutron = new TH1D("hpx_neutron", "Neutron Momentum px",100,-1,1);
    TH1D* hpy_neutron = new TH1D("hpy_neutron", "Neutron Momentum py",100,-1,1);
    TH1D* hpz_neutron = new TH1D("hpz_neutron", "Neutron Momentum pz",100,-2,2);
    TH1D* hvx_neutron = new TH1D("hvx_neutron", "Neutron Vertex x",100,-2,2);
    TH1D* hvy_neutron = new TH1D("hvy_neutron", "Neutron Vertex y",100,-2,2);
    TH1D* hvz_neutron = new TH1D("hvz_neutron", "Neutron Vertex z",100,-2,2);
    
    
    //+pion
    TH1D* hE_pip = new TH1D ("hE_pip","Pion+ Energy",100, 0, 1.5);
    TH1D* hTmom_pip = new TH1D("hTmom_pip", "Total Pion+ momentum", 100, 0,1.5);
    TH1D* hpx_pip  = new TH1D("hpx_pip", "Pion+ Momentum px",100,-1,1);
    TH1D* hpy_pip  = new TH1D("hpy_pip", "Pion+ Momentum py",100,-1,1);
    TH1D* hpz_pip  = new TH1D("hpz_pip", "Pion+ Momentum pz",100,-2,2);
    TH1D* hvx_pip  = new TH1D("hvx_pip", "Pion+ Vertex x",100,-10,10);
    TH1D* hvy_pip = new TH1D("hvy_pip", "Pion+ Vertex y",100,-10,10);
    TH1D* hvz_pip  = new TH1D("hvz_pip", "Pion+ Vertex z",100,-10,10);
    
    //-pion
    TH1D* hE_pim = new TH1D("hE_pim", "Pion- Energy", 100, 0, 1);
    TH1D* hTmom_pim = new TH1D("hTmom_pim", "Total Pion- momentum", 100, 0, 2);
    TH1D* hpx_pim = new TH1D("hpx_pim", "Pion- Momentum px", 100, -2, 2);
    TH1D* hpy_pim = new TH1D("hpy_pim", "Pion- Momentum py", 100, -2, 2);
    TH1D* hpz_pim = new TH1D("hpz_pim", "Pion- Momentum pz", 100, -2, 2);
    TH1D* hvx_pim = new TH1D("hvx_pim", "Pion- Vertex vx", 100, -10, 10);
    TH1D* hvy_pim = new TH1D("hvy_pim", "Pion- Vertex vy", 100, -10, 10);
    TH1D* hvz_pim = new TH1D("hvz_pim", "Pion- Vertex vz", 100, -10, 10);
    

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        for (int i = 0; i < nparts; i++) {
            
            hwgt->Fill(wgt);
            hEin->Fill(Ein);
            hnparts->Fill(nparts);
            hPDG->Fill(pdg[i]);
            hinter->Fill(inter[i]);
            hhist->Fill(hist[i]);
            
            
            //Proton
            if (pdg[i] == 2212) {
            hE_proton->Fill(E[i], wgt);
            double totalmomentum_proton = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
            hTmom_proton->Fill(totalmomentum_proton, wgt);
            hpx_proton->Fill(px[i], wgt);
            hpy_proton->Fill(py[i], wgt);
            hpz_proton->Fill(pz[i], wgt);
            hvx_proton->Fill(vx[i], wgt);
            hvy_proton->Fill(vy[i], wgt);
            hvz_proton->Fill(vz[i], wgt);
            
            }
            
            //Neutron
            if (pdg[i] == 2112) {
            hE_neutron->Fill(E[i], wgt);
            double totalmomentum_neutron = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
            hTmom_neutron->Fill(totalmomentum_neutron, wgt);
            hpx_neutron->Fill(px[i], wgt);
            hpy_neutron->Fill(py[i], wgt);
            hpz_neutron->Fill(pz[i], wgt);
            hvx_neutron->Fill(pz[i], wgt);
            hvy_neutron->Fill(vx[i], wgt);
            hvz_neutron->Fill(vy[i], wgt);
            
            }
            
            
            //+pion
            if (pdg[i] == 211) {
            
            hE_pip->Fill(E[i], wgt);
            double totalmomentum_pip = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
            hTmom_pip->Fill(totalmomentum_pip, wgt);
            hpx_pip->Fill(px[i], wgt);
            hpy_pip->Fill(py[i], wgt);
            hpz_pip->Fill(pz[i], wgt);
            hvx_pip->Fill(vx[i], wgt);
            hvy_pip->Fill(vy[i], wgt);
            hvz_pip->Fill(vz[i], wgt);
            
            }
            
            //-pion
            if (pdg[i] == -211) {
            hE_pim->Fill(E[i], wgt);
            double totalmomentum_pim = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
            hTmom_pim->Fill(totalmomentum_pim, wgt);
            hpx_pim->Fill(px[i], wgt);
            hpy_pim->Fill(py[i], wgt);
            hpz_pim->Fill(pz[i], wgt);
            hvx_pim->Fill(vx[i], wgt);
            hvy_pim->Fill(vy[i], wgt);
            hvz_pim->Fill(vz[i], wgt);
            
            }
            
        }
    }
    
    

    hwgt->Write();
    hEin->Write();
    hnparts->Write();
    hPDG->Write();
    hinter->Write();
    hhist->Write();
    
    
    
    
    //Proton
    
    
    hE_proton->Write();
    hTmom_proton->Write();
    hpx_proton->Write();
    hpy_proton->Write();
    hpz_proton->Write();
    hvx_proton->Write();
    hvy_proton->Write();
    hvz_proton->Write();
    
    //Neutron
    hE_neutron->Write();
    hTmom_neutron->Write();
    hpx_neutron->Write();
    hpy_neutron->Write();
    hpz_neutron->Write();
    hvx_neutron->Write();
    hvy_neutron->Write();
    hvz_neutron->Write();
    
    
    //+pion
    hE_pip->Write();
    hTmom_pip->Write();
    hpx_pip->Write();
    hpy_pip->Write();
    hpz_pip->Write();
    hvx_pip->Write();
    hvy_pip->Write();
    hvz_pip->Write();
    
    //-pion
    hE_pim->Write();
    hTmom_pim->Write();
    hpx_pim->Write();
    hpy_pim->Write();
    hpz_pim->Write();
    hvx_pim->Write();
    hvy_pim->Write();
    hvz_pim->Write();
    
    fOut->Close();


}


