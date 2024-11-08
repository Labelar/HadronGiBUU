//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct  4 18:34:42 2024 by ROOT version 6.28/00
// from TTree records/records
// found on file: scan_lowpionabsorption_000_std_pipAr40_base0001MeV.root
//////////////////////////////////////////////////////////

#ifndef MyClass_mod_h
#define MyClass_mod_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

const int maxnparts = 2000;

// Header file for the classes stored in the TTree if any.

class MyClass_mod {
public :
   TFile          *fOut;     //!output file 
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
  Int_t            fNtrees;  //! number of trees
  
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         wgt;
   Float_t         Ein;
   Int_t           nparts;
   Int_t           pdg[maxnparts];   //[nparts]
   Int_t           inter[maxnparts];   //[nparts]
   Int_t           hist[maxnparts];   //[nparts]
   Float_t         E[maxnparts];   //[nparts]
   Float_t         px[maxnparts];   //[nparts]
   Float_t         py[maxnparts];   //[nparts]
   Float_t         pz[maxnparts];   //[nparts]
   Float_t         vx[maxnparts];   //[nparts]
   Float_t         vy[maxnparts];   //[nparts]
   Float_t         vz[maxnparts];   //[nparts]

   // List of branches
   TBranch        *b_wgt;   //!
   TBranch        *b_Ein;   //!
   TBranch        *b_nparts;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_inter;   //!
   TBranch        *b_hist;   //!
   TBranch        *b_E;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!

  MyClass_mod(const char* cfiles,const char* outroot);
   virtual ~MyClass_mod();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyClass_mod_cxx
MyClass_mod::MyClass_mod(const char* cfiles,const char* outroot){

  TChain* evts   = new TChain("records");
  fNtrees = 0;
  std::ifstream fIn(cfiles);
  std::string line;
  while(std::getline(fIn, line)){
    if( line.find(".root") == std::string::npos) continue;
    std::cout<<"=> Input file: "<<line<<std::endl;
    evts->Add(line.c_str());
    fNtrees++;
  }

  Init(evts);

  //Output file

  fOut = new TFile(outroot,"recreate");
}

MyClass_mod::~MyClass_mod()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass_mod::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass_mod::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass_mod::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("wgt", &wgt, &b_wgt);
   fChain->SetBranchAddress("Ein", &Ein, &b_Ein);
   fChain->SetBranchAddress("nparts", &nparts, &b_nparts);
   fChain->SetBranchAddress("pdg", pdg, &b_pdg);
   fChain->SetBranchAddress("inter", inter, &b_inter);
   fChain->SetBranchAddress("hist", hist, &b_hist);
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   Notify();
}

Bool_t MyClass_mod::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass_mod::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass_mod::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_mod_cxx
