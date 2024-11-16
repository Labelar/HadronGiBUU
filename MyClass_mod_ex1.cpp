#define MyClass_mod_cxx
#include "MyClass_mod.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>

void MyClass_mod::Loop() {

    // Declare histograms in vectors to simplify management
    std::vector<TH1D*> hwgt, hEin, hnparts, hPDG, hinter, hhist;
    std::vector<TH1D*> hE_proton, hTmom_proton, hpx_proton, hpy_proton, hpz_proton, hvx_proton, hvy_proton, hvz_proton;
    std::vector<TH1D*> hE_neutron, hTmom_neutron, hpx_neutron, hpy_neutron, hpz_neutron, hvx_neutron, hvy_neutron, hvz_neutron;
    std::vector<TH1D*> hE_pip, hTmom_pip, hpx_pip, hpy_pip, hpz_pip, hvx_pip, hvy_pip, hvz_pip;
    std::vector<TH1D*> hE_pim, hTmom_pim, hpx_pim, hpy_pim, hpz_pim, hvx_pim, hvy_pim, hvz_pim;

    // New histograms for cross-section
    const int Nint = 6;
    const char* cint[Nint] = {"rea","abs","pip","pim","pi0","2pi"};
    int Nbins = 1002;
    double mine = -0.0005;
    double maxe = 1.0015;
    std::vector<TH1D*> hxsec(Nint);

    for (int i = 0; i < Nint; ++i) {
        hxsec[i] = new TH1D(Form("hxsec_%s", cint[i]), "", Nbins, mine, maxe);
    }

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    // Store the branches into local variables using SetBranchAddress
    std::vector<int> pdg(nparts), inter(nparts), hist(nparts), nparts(nentries);
    std::vector<double> px(nparts), py(nparts), pz(nparts), E(nparts), vx(nparts), vy(nparts), vz(nparts);

    fChain->SetBranchAddress("pdg", &pdg[0]);
    fChain->SetBranchAddress("inter", &inter[0]);
    fChain->SetBranchAddress("hist", &hist[0]);
    fChain->SetBranchAddress("nparts", &nparts[0]);
    fChain->SetBranchAddress("px", &px[0]);
    fChain->SetBranchAddress("py", &py[0]);
    fChain->SetBranchAddress("pz", &pz[0]);
    fChain->SetBranchAddress("E", &E[0]);
    fChain->SetBranchAddress("vx", &vx[0]);
    fChain->SetBranchAddress("vy", &vy[0]);
    fChain->SetBranchAddress("vz", &vz[0]);

    for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        // Skip events with inter[0] == -999
        if (inter[0] == -999) {
            if (nparts > 1) std::cout << "=> Warning (inter[0] == -999) events with nparts = " << nparts << std::endl;
            continue;
        }

        // Fill total cross-section histogram
        hxsec[0]->Fill(Ein, wgt);

        int Npip = 0, Npim = 0, Npi0 = 0;

        // Process each particle in the event
        for (int i = 0; i < nparts; ++i) {
            // Count pions
            if (pdg[i] == 211) ++Npip;
            if (pdg[i] == -211) ++Npim;
            if (pdg[i] == 111) ++Npi0;

            // Fill histograms with particle properties
            hwgt.push_back(new TH1D("hwgt", "wgt", 100, -1, 2)); // Add once
            hEin.push_back(new TH1D("hEin", "Ein", 100, -1, 2));
            hnparts.push_back(new TH1D("hnparts", "", 100, -1, 50));
            hPDG.push_back(new TH1D("hPDG", "PDG", 100, -1000, 2500));
            hinter.push_back(new TH1D("hinter", "", 100, -3000, 250000));
            hhist.push_back(new TH1D("hhist", "", 100, -1, 2)); 

            // Proton
            if (pdg[i] == 2212) {
                double totalmomentum_proton = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
                hE_proton.push_back(new TH1D("hE_proton", "Proton Energy", 100, 0.75, 1.25));  // Create histos once
                hE_proton.back()->Fill(E[i], wgt);
                hTmom_proton.push_back(new TH1D("hTmom_proton", "Total proton momentum", 100, 0, 1.25));
                hTmom_proton.back()->Fill(totalmomentum_proton, wgt);
                // Continue for other proton histograms...
            }

            // Neutron, Pion+, Pion- histograms follow the same pattern

        }

        // Cross-section assignment
        int idx = -1;
        if (Npip == 0 && Npim == 0 && Npi0 == 0) idx = 1;
        if (Npim == 1 && Npip == 0) idx = 2;
        if (Npim == 0 && Npip == 1) idx = 3;
        if (Npi0 == 1) idx = 4;
        if (Npip + Npim + Npi0 > 1) idx = 5;
        if (idx > 0) hxsec[idx]->Fill(Ein, wgt);
    }

    // Write histograms to file
    fOut->cd();
    fOut->mkdir("xsec");
    fOut->cd("xsec");
    for (auto& hist : hxsec) {
        hist->Scale(1. / double(fNtrees));  // Normalize cross-section histograms
        hist->Write();
    }

    fOut->cd();
    // Write other histograms
    for (auto& hist : hwgt) hist->Write();
    // Continue for all histograms...

    fOut->Close();
}
