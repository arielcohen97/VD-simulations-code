//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug  9 22:18:01 2023 by ROOT version 6.26/07
// from TTree FlashMatchTree/FlashMatchTree
// found on file: ana_hist.root
//////////////////////////////////////////////////////////

#ifndef Sigbackground_h
#define Sigbackground_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class Sigbackground {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           EventID;
   Float_t         TrueX;
   Float_t         TrueY;
   Float_t         TrueZ;
   Float_t         TrueT;
   Float_t         DetectedT;
   Float_t         TrueE;
   Float_t         Edep;
   Float_t         EdepSimE;
   Int_t           TruePDG;
   Int_t           NFlashes;
   vector<int>     *FlashIDVector;
   vector<float>   *YCenterVector;
   vector<float>   *ZCenterVector;
   vector<float>   *YWidthVector;
   vector<float>   *ZWidthVector;
   vector<float>   *TimeVector;
   vector<float>   *TimeWidthVector;
   vector<float>   *TimeDiffVector;
   vector<float>   *TotalPEVector;
   Float_t         SumPE;
   Int_t           NOpDets;
   vector<int>     *NHitOpDetVector;
   vector<float>   *Purity;
   vector<float>   *Distance;
   vector<float>   *RecoXVector;
   vector<short>   *DAQch_ID;
   vector<vector<float> > *onda_ID;
   vector<int>     *OpHittime;
   vector<int>     *OpHitChannel;
   vector<float>   *OpHitPE;

   // List of branches
   TBranch        *b_EventID;   //!
   TBranch        *b_TrueX;   //!
   TBranch        *b_TrueY;   //!
   TBranch        *b_TrueZ;   //!
   TBranch        *b_TrueT;   //!
   TBranch        *b_DetectedT;   //!
   TBranch        *b_TrueE;   //!
   TBranch        *b_Edep;   //!
   TBranch        *b_EdepSimE;   //!
   TBranch        *b_TruePDG;   //!
   TBranch        *b_NFlashes;   //!
   TBranch        *b_FlashIDVector;   //!
   TBranch        *b_YCenterVector;   //!
   TBranch        *b_ZCenterVector;   //!
   TBranch        *b_YWidthVector;   //!
   TBranch        *b_ZWidthVector;   //!
   TBranch        *b_TimeVector;   //!
   TBranch        *b_TimeWidthVector;   //!
   TBranch        *b_TimeDiffVector;   //!
   TBranch        *b_TotalPEVector;   //!
   TBranch        *b_SumPE;   //!
   TBranch        *b_NOpDets;   //!
   TBranch        *b_NHitOpDetVector;   //!
   TBranch        *b_Purity;   //!
   TBranch        *b_Distance;   //!
   TBranch        *b_RecoXVector;   //!
   TBranch        *b_DAQch_ID;   //!
   TBranch        *b_onda_ID;   //!
   TBranch        *b_OpHittime;   //!
   TBranch        *b_OpHitChannel;   //!
   TBranch        *b_OpHitPE;   //!

   Sigbackground(TTree *tree=0);
   virtual ~Sigbackground();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Sigbackground_cxx
Sigbackground::Sigbackground(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("25MeV_signal_only.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("25MeV_signal_only.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("25MeV_signal_only.root:/vdflashmatch");
      dir->GetObject("FlashMatchTree",tree);

   }
   Init(tree);
}

Sigbackground::~Sigbackground()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Sigbackground::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Sigbackground::LoadTree(Long64_t entry)
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

void Sigbackground::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   FlashIDVector = 0;
   YCenterVector = 0;
   ZCenterVector = 0;
   YWidthVector = 0;
   ZWidthVector = 0;
   TimeVector = 0;
   TimeWidthVector = 0;
   TimeDiffVector = 0;
   TotalPEVector = 0;
   NHitOpDetVector = 0;
   Purity = 0;
   Distance = 0;
   RecoXVector = 0;
   DAQch_ID = 0;
   onda_ID = 0;
   OpHittime = 0;
   OpHitChannel = 0;
   OpHitPE = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventID", &EventID, &b_EventID);
   fChain->SetBranchAddress("TrueX", &TrueX, &b_TrueX);
   fChain->SetBranchAddress("TrueY", &TrueY, &b_TrueY);
   fChain->SetBranchAddress("TrueZ", &TrueZ, &b_TrueZ);
   fChain->SetBranchAddress("TrueT", &TrueT, &b_TrueT);
   fChain->SetBranchAddress("DetectedT", &DetectedT, &b_DetectedT);
   fChain->SetBranchAddress("TrueE", &TrueE, &b_TrueE);
   fChain->SetBranchAddress("Edep", &Edep, &b_Edep);
   fChain->SetBranchAddress("EdepSimE", &EdepSimE, &b_EdepSimE);
   fChain->SetBranchAddress("TruePDG", &TruePDG, &b_TruePDG);
   fChain->SetBranchAddress("NFlashes", &NFlashes, &b_NFlashes);
   fChain->SetBranchAddress("FlashIDVector", &FlashIDVector, &b_FlashIDVector);
   fChain->SetBranchAddress("YCenterVector", &YCenterVector, &b_YCenterVector);
   fChain->SetBranchAddress("ZCenterVector", &ZCenterVector, &b_ZCenterVector);
   fChain->SetBranchAddress("YWidthVector", &YWidthVector, &b_YWidthVector);
   fChain->SetBranchAddress("ZWidthVector", &ZWidthVector, &b_ZWidthVector);
   fChain->SetBranchAddress("TimeVector", &TimeVector, &b_TimeVector);
   fChain->SetBranchAddress("TimeWidthVector", &TimeWidthVector, &b_TimeWidthVector);
   fChain->SetBranchAddress("TimeDiffVector", &TimeDiffVector, &b_TimeDiffVector);
   fChain->SetBranchAddress("TotalPEVector", &TotalPEVector, &b_TotalPEVector);
   fChain->SetBranchAddress("SumPE", &SumPE, &b_SumPE);
   fChain->SetBranchAddress("NOpDets", &NOpDets, &b_NOpDets);
   fChain->SetBranchAddress("NHitOpDetVector", &NHitOpDetVector, &b_NHitOpDetVector);
   fChain->SetBranchAddress("Purity", &Purity, &b_Purity);
   fChain->SetBranchAddress("Distance", &Distance, &b_Distance);
   fChain->SetBranchAddress("RecoXVector", &RecoXVector, &b_RecoXVector);
   fChain->SetBranchAddress("DAQch_ID", &DAQch_ID, &b_DAQch_ID);
   fChain->SetBranchAddress("onda_ID", &onda_ID, &b_onda_ID);
   fChain->SetBranchAddress("OpHittime", &OpHittime, &b_OpHittime);
   fChain->SetBranchAddress("OpHitChannel", &OpHitChannel, &b_OpHitChannel);
   fChain->SetBranchAddress("OpHitPE", &OpHitPE, &b_OpHitPE);
   Notify();
}

Bool_t Sigbackground::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Sigbackground::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Sigbackground::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Sigbackground_cxx
