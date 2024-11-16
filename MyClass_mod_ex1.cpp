#define MyClass_mod_cxx
#include "MyClass_mod.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MyClass_mod::Loop() {
    
    // Basic histograms
    TH1D* hwgt = new TH1D("hwgt", "wgt", 100, -1, 2);
    TH1D* hEin = new TH1D("hEin", "Ein", 100, -1, 2);
    TH1D* hnparts = new TH1D("hnparts", "", 100, -1, 50);
    TH1D* hPDG = new TH1D("hPDG", "PDG", 100, -1000, 2500);
    TH1D* hinter = new TH1D("hinter", "", 100, -3000, 250000);
    TH1D* hhist = new TH1D("hhist", "", 100, -1, 2); 

    // Create histograms for different particle types
    std::map<std::string, std::vector<TH1D*>> particleHistograms;
    std::vector<std::string> particles = {"proton", "neutron", "pip", "pim"};

    // For each particle type, initialize a set of histograms
    for (auto& particle : particles) {
        particleHistograms[particle].push_back(new TH1D(("hE_" + particle).c_str(), (particle + " Energy").c_str(), 100, 0.5, 2));
        particleHistograms[particle].push_back(new TH1D(("hTmom_" + particle).c_str(), ("Total " + particle + " momentum").c_str(), 100, 0, 2));
        particleHistograms[particle].push_back(new TH1D(("hpx_" + particle).c_str(), (particle + " Momentum px").c_str(), 100, -2, 2));
        particleHistograms[particle].push_back(new TH1D(("hpy_" + particle).c_str(), (particle + " Momentum py").c_str(), 100, -2, 2));
        particleHistograms[particle].push_back(new TH1D(("hpz_" + particle).c_str(), (particle + " Momentum pz").c_str(), 100, -2, 2));
        particleHistograms[particle].push_back(new TH1D(("hvx_" + particle).c_str(), (particle + " Vertex x").c_str(), 100, -10, 10));
        particleHistograms[particle].push_back(new TH1D(("hvy_" + particle).c_str(), (particle + " Vertex y").c_str(), 100, -10, 10));
        particleHistograms[particle].push_back(new TH1D(("hvz_" + particle).c_str(), (particle + " Vertex z").c_str(), 100, -10, 10));
    }

    // New histograms for cross-section types
    const int Nint = 6;
    const char* cint[Nint] = {"rea", "abs", "pip", "pim", "pi0", "2pi"};
    int Nbins = 1002;
    double mine = -0.0005;
    double maxe = 1.0015;
    TH1D* hxsec[Nint];
    for (int i = 0; i < Nint; i++) {
        hxsec[i] = new TH1D(Form("hxsec_%s", cint[i]), "", Nbins, mine, maxe);
    }

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        if (inter[0] == -999) {
            if (nparts > 1) std::cout << "=> Warning (inter[0] == -999) events with nparts = " << nparts << std::endl;
            continue;
        }

        // Fill cross-section histograms
        hxsec[0]->Fill(Ein, wgt);

        int Npip = 0, Npim = 0, Npi0 = 0;
        for (int i = 0; i < nparts; i++) {
            // Count pions
            if (pdg[i] == 211) ++Npip;
            if (pdg[i] == -211) ++Npim;
            if (pdg[i] == 111) ++Npi0;

            // Fill basic histograms
            hwgt->Fill(wgt);
            hEin->Fill(Ein);
            hnparts->Fill(nparts);
            hPDG->Fill(pdg[i]);
            hinter->Fill(inter[i]);
            hhist->Fill(hist[i]);

            // Loop through each particle and fill the appropriate histograms
            if (pdg[i] == 2212) { // Proton
                FillParticleHistograms(particleHistograms["proton"], E[i], px[i], py[i], pz[i], vx[i], vy[i], vz[i], wgt);
            } 
            else if (pdg[i] == 2112) { // Neutron
                FillParticleHistograms(particleHistograms["neutron"], E[i], px[i], py[i], pz[i], vx[i], vy[i], vz[i], wgt);
            } 
            else if (pdg[i] == 211) { // Pion+
                FillParticleHistograms(particleHistograms["pip"], E[i], px[i], py[i], pz[i], vx[i], vy[i], vz[i], wgt);
            } 
            else if (pdg[i] == -211) { // Pion-
                FillParticleHistograms(particleHistograms["pim"], E[i], px[i], py[i], pz[i], vx[i], vy[i], vz[i], wgt);
            }
        }

        // Determine cross-section type
        int idx = -1;
        if (Npip == 0 && Npim == 0 && Npi0 == 0) idx = 1;
        if (Npip == 1 && Npim == 0 && Npi0 == 0) idx = 2;
        if (Npip == 0 && Npim == 1 && Npi0 == 0) idx = 3;
        if (Npip == 0 && Npim == 0 && Npi0 == 1) idx = 4;
        if ((Npip + Npim + Npi0) > 1) idx = 5;
        
        if (idx > 0) hxsec[idx]->Fill(Ein, wgt);
    }

    // Write histograms
    fOut->cd();
    fOut->mkdir("xsec");
    fOut->cd("xsec");
    for (int i = 0; i < Nint; i++) {
        hxsec[i]->Scale(1. / double(fNtrees));
        hxsec[i]->Write();
    }

    fOut->cd();
    hwgt->Write();
    hEin->Write();
    hnparts->Write();
    hPDG->Write();
    hinter->Write();
    hhist->Write();

    // Write particle histograms
    for (auto& particle : particles) {
        for (auto& hist : particleHistograms[particle]) {
            hist->Write();
        }
    }

    // Display energy distribution for protons
    TCanvas* c1 = new TCanvas("c1", "Energy distribution", 800, 600);
    particleHistograms["proton"][0]->SetTitle("Energy distribution");
    particleHistograms["proton"][0]->GetXaxis()->SetTitle("Energy (GeV)");
    particleHistograms["proton"][0]->GetYaxis()->SetTitle("Entries");

    fOut->Close();
}

// Function to fill histograms for a given particle
void MyClass_mod::FillParticleHistograms(std::vector<TH1D*>& histograms, double E, double px, double py, double pz, double vx, double vy, double vz, double wgt) {
    double totalmomentum = sqrt(px * px + py * py + pz * pz);
    histograms[0]->Fill(E, wgt); // Energy histogram
    histograms[1]->Fill(totalmomentum, wgt); // Total momentum histogram
    histograms[2]->Fill(px, wgt); // Momentum px histogram
    histograms[3]->Fill(py, wgt); // Momentum py histogram
    histograms[4]->

