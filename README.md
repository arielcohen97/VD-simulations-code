# VD-simulations-code
Codes and files for VD simulations DUNE
#include "TString.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <string>


using namespace std;
vector<TString> list_files(string folder){
        vector<TString> filelist;
        TSystemDirectory dir(folder.c_str(), folder.c_str());
        TList *files = dir.GetListOfFiles();
        if(files){
                TSystemFile *file;
                TString fname;
                TIter next(files);
                while ((file=(TSystemFile*)next())){
                        fname = file->GetName();
                        if (!file->IsDirectory() && fname.EndsWith(".root")){
                                filelist.push_back(folder + fname.Data());
                        }
                }
        }

        return filelist;
}


void for_ariel(){
        const char* folder = "/home/acohen/Documents/ejercicios_root/VD_simulations/august_run_energy/hist_all_energ_norm/";
        const char* dirname = "vdflashmatch";
        const char* tree_name = "FlashMatchTree";

        vector<TString> filenames = list_files(folder);
//        std::sort(filenames.begin(), filenames.end());

        uint counter = 0;
        vector<TH1F*> hists;
        vector<Color_t> color = {kRed, kBlack, kGreen, kMagenta, kCyan};
//        vector<Int_t> color = {30, 31, 32, 34} //, 35};
        vector<TString> legend = {"16-20 MeV", "26-30 MeV", "5-10 MeV", "21-25 MeV", "11-15 MeV"}; //"16-20 MeV", "21-25 MeV", "26-30 MeV"};
        Int_t aux = 0;

        //Loading data in hists
        for(const TString& fname: filenames){
                TFile *f = new TFile(fname);
                TDirectory *dir = f->Get<TDirectory>(dirname);
                TTree *t = dir->Get<TTree>(tree_name);

                TString hname = Form("h%d", counter);

                TH1F *h = new TH1F(hname, hname, 100, -350, 350);
                hists.push_back(h);
                t->Draw(Form("RecoXVector>>%s", hname.Data()), "",  "goff");

                counter++;
        }

        //Rescaling
        for(TH1F *h : hists){
                h->Scale(1.0/h->GetMaximum());
        }

        //Drawing

        TCanvas *c = new TCanvas();
        TString drawOpt = "hist";
        for(Int_t i = 0; i < hists.size(); i++){
               hists[i]->GetXaxis()->SetTitle("Reconstruced X");
               hists[i]->GetYaxis()->SetTitle("Entries");
               //hists[i]->SetTitle("Reconstructed X 5-30 MeV");
               hists[i]->SetTitle(legend[i]);
               hists[i]->SetLineColor(color[i]);
               hists[i]->Draw(drawOpt);
               if (drawOpt == "hist"){
                       drawOpt += " same";
               }
       }
	c->BuildLegend();
//        auto legend = new TLegend();
//        for(TH1F *h : hists){
//                h->GetXaxis()->SetTitle("Reconstruced X");
//                h->GetYaxis()->SetTitle("Entries");
//                h->SetTitle("Reconstructed X 5-30 MeV");
//                h->Draw(drawOpt);
//                if (drawOpt == "hist"){
//                        drawOpt += " same";
//                }
        //}


}
