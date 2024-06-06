#define PE_calc_cxx
#include "PE_calc.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TFile.h>
#include <TGraph.h>
#include <map>
#include <TLegend.h>
#include <string>
#include <numeric>
#include <TH1F.h>

enum Plane
{
   kC,
   kRW,
   kLW,
   kNoPlane
};

struct Flash_info
{
   Plane plane;
   float w_sumX = 0;
   float w_sumY = 0;
   float w_sumZ = 0;
   float sum_PE = 0;
   float oldw_sumX = 0;
   float oldw_sumY = 0;
   float oldw_sumZ = 0;
   float old_sumweight = 0;
   float min_dist_C = 1000000;
   float flash_time = 0;
   float flash_PEs = 0;
   int channel_counter = 0;
   int hit_counter = 0;
   std::vector<Float_t> hit_vectorX;
   std::vector<Float_t> hit_vectorY;
   std::vector<Float_t> hit_vectorZ;
   std::vector<Float_t> weight_hit_vector;
   std::vector<std::vector<Float_t>> flash_channels;
   float meanX = 0;
   float meanY = 0;
   float meanZ = 0;
   float errX = 0;
   float errY = 0;
   float errZ = 0;
};

float rms_calc(std::vector<float> coord, std::vector<float> weight, float mean)
{
   float value = 0;
   float sum_weights = 0;
   for (int k = 0; k < coord.size(); k++)
   {
      value += pow((coord.at(k) - mean), 2) * weight.at(k);
      sum_weights += weight.at(k);
   }
   if (sum_weights > 0)
   {
      value = pow(value / sum_weights, 0.5);
   }
   else
   {
      value = 0;
   }
   return value;
}

float combined_mean_calc(std::vector<float> coord, std::vector<float> weight)
{
   float value = 0;
   float sum_weights = 0;
   for (int i = 0; i < coord.size(); i++)
   {
      value += coord.at(i) * weight.at(i);
      sum_weights += weight.at(i);
   }
   if (sum_weights != 0)
   {
      value = value / sum_weights;
   }
   else
   {

      value = 0;
   }
   return value;
}

float calculate_PEs(float pe_vec, float dist_vec, float sin_t)
{
   float N_PEs = 0;
   float number = 0;
   float area = 60 * 60;

   number = 4 * TMath::Pi() * dist_vec * pe_vec;

   number /= area * TMath::Sin(sin_t);
   // std::cout << "Chann PEs " << pe_vec << std::endl;
   // std::cout << "Channel position " << dist_vec << std::endl;
   // // std::cout << "Theta " << sin_t << std::endl;
   // std::cout << "#PEs " << number << std::endl;

   N_PEs = number;

   return N_PEs;
}

// Groups the flashes in clusters with nearby Z and time means.
std::vector<std::vector<Flash_info>> coincidence_mean(std::vector<Flash_info> flashes)
{
   //    return {flashes};
   std::vector<std::vector<Flash_info>> clusters;
   while (flashes.size() != 0)
   {
      Flash_info new_flash = flashes[0];
      flashes.erase(flashes.begin());
      bool found_cluster = false;
      for (auto &cluster : clusters)
      {
         if (abs(cluster[0].flash_time - new_flash.flash_time) <= 1 && abs(cluster[0].meanZ - new_flash.meanZ) <= 400)
         {
            if ((cluster[0].plane == kC && new_flash.plane != kC) || (cluster[0].plane != kC && new_flash.plane == kC))
            {
               cluster.push_back(new_flash);
               found_cluster = true;
               break;
            }
         }
      }
      if (!found_cluster)
      {
         clusters.push_back({new_flash});
      }
   }
   return clusters;
}

Plane getPlane(float chX, float chY, float chZ)
{
   if (chX < -320)
   {
      return kC;
   }
   if (chY < -700)
   {
      return kLW;
   }
   if (chY > 700)
   {
      return kRW;
   }
   return kNoPlane;
}

std::vector<float> global_PE_vec;
std::vector<float> trueE_vec;
std::vector<float> event_PE_vec;
float flash_counter = 0;
float flash_counter_cut = 0;
float scale;
std::vector<Float_t> channel_vector(3);

std::vector<std::string> filenames = {
    "bckg_only_v35.root",
    "opslicer_time_360_400_400.root"};

void run_multiple(std::vector<std::string> filenames)
{
   std::map<std::string, TH1F *> hists;

   TFile *ofile = new TFile("output_hists.root", "RECREATE");
   auto legend = new TLegend();
   int aux = 0;
   for (std::string const &fname : filenames)
   {
      TFile *f = new TFile(fname.c_str());
      TTree *tree;
      TDirectory *dir = (TDirectory *)f->Get("vdflashmatch");
      dir->GetObject("FlashMatchTree", tree);
      if (!tree)
      {
         std::cout << "File " << fname << "does not exist!" << std::endl;
         abort();
      }

      TH1F *hist;
      PE_calc calculator(tree);
      calculator.Loop(&hist);
      hists[fname] = hist;

      ofile->cd();
      if (aux == 0)
      {
         //  hist -> Scale(scale*30);
         hist->SetLineColor(kRed);
         legend->AddEntry(hist, Form("Mean X-ARAPUCA hit for background, cathode: 7.745"));
         //  legend -> AddEntry(hist, Form("#Mean amount of PEs, signal: 2.949e+04"));
      }
      if (aux == 1)
      {
         hist->Scale(5);
         legend->AddEntry(hist, Form("Mean X-ARAPUCA hit for signal, cathode: 16.34"));
         // legend -> AddEntry(hist, Form("#Mean amount of PEs, background: 4.385e+04"));
      }
      // std::cout << ofile << std::endl;
      hist->Write(fname.c_str());
      aux++;
      f->Close();
   }

   ofile->Close();

   TCanvas *c = new TCanvas;
   bool isFirst = true;
   for (auto const &[fname, hist] : hists)
   {
      const char *drawOpt = isFirst ? "HIST" : "SAME HIST";
      isFirst = false;
      //   cout << "Yoooo" << endl;
      hist->Draw(drawOpt);
   }
   legend->Draw();
}

void PE_calc::Loop(TH1F **ret_hist)
{

   //   In a ROOT session, you can do:
   //      root> .L PE_calc.C
   //      root> PE_calc t
   //      root> t.GetEntry(12); // Fill t data members with entry number 12
   //      root> t.Show();       // Show values of entry 12
   //      root> t.Show(16);     // Read and show values of entry 16
   //      root> t.Loop();       // Loop on all entries
   //

   //     This is the loop skeleton where:
   //    jentry is the global entry number in the chain
   //    ientry is the entry number in the current Tree
   //  Note that the argument to GetEntry must be:
   //    jentry for TChain::GetEntry
   //    ientry for TTree::GetEntry and TBranch::GetEntry
   //
   //       To read only selected branches, Insert statements like:
   // METHOD1:
   //    fChain->SetBranchStatus("*",0);  // disable all branches
   //    fChain->SetBranchStatus("branchname",1);  // activate branchname
   // METHOD2: replace line
   //    fChain->GetEntry(jentry);       //read all branches
   // by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0)
      return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   TH1F *hist = new TH1F("hist", "PEs per flash", 100, 0, 20);
   *ret_hist = hist;
   hist->SetDirectory(0);
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      std::vector<float> xmean_vec;
      std::vector<float> ymean_vec;
      std::vector<float> zmean_vec;
      std::vector<float> xweight_vec;
      std::vector<float> yweight_vec;
      std::vector<float> zweight_vec;
      std::vector<Float_t> hit_vectorXR;
      std::vector<Float_t> hit_vectorXL;
      std::vector<Float_t> hit_vectorZC;
      std::vector<Float_t> hit_vectorZR;
      std::vector<Float_t> hit_vectorZL;
      std::vector<Float_t> hit_vectorYC;
      std::vector<Float_t> hit_vectorYback;
      std::vector<Float_t> hit_vectorYfront;
      std::vector<Float_t> hit_vectorXback;
      std::vector<Float_t> hit_vectorXfront;
      std::vector<Float_t> weight_hit_vectorXR;
      std::vector<Float_t> weight_hit_vectorXL;
      std::vector<Float_t> weight_hit_vectorZC;
      std::vector<Float_t> weight_hit_vectorYC;
      std::vector<Float_t> weight_hit_vectorZR;
      std::vector<Float_t> weight_hit_vectorZL;
      std::vector<Float_t> weight_hit_vectorXback;
      std::vector<Float_t> weight_hit_vectorXfront;
      std::vector<Float_t> weight_hit_vectorYback;
      std::vector<Float_t> weight_hit_vectorYfront;
      std::vector<Float_t> RMS_vectorX;
      std::vector<Float_t> RMS_vectorY;
      std::vector<Float_t> RMS_vectorZ;
      std::vector<Float_t> TrueY_vector;
      std::vector<Float_t> TrueX_vector;
      std::vector<Float_t> TrueZ_vector;
      std::vector<Float_t> vector_newresy;
      std::vector<Float_t> bckg_X;
      std::vector<Float_t> bckg_Y;
      float counter = 0;
      float flash_mean_PEs = 0;
      scale = 1. / nentries;
      std::map<int, Flash_info> flashes;
      for (int j = 0; j < OpHitChannel_Flash->size(); j++)
      {
         if (flashes.count(OpFlashNumber->at(j)) == 0)
         {
            flashes[OpFlashNumber->at(j)] = Flash_info();
            flashes[OpFlashNumber->at(j)].plane = getPlane(Channel_X->at(j), Channel_Y->at(j), Channel_Z->at(j));
         }
         Flash_info &flash = flashes[OpFlashNumber->at(j)];
         channel_vector[0] = Channel_X->at(j);
         channel_vector[1] = Channel_Y->at(j);
         channel_vector[2] = Channel_Z->at(j);
         if (flash.flash_channels.size() == 0)
         {
            flash.flash_channels.push_back(channel_vector);
            flash.channel_counter++;
         }
         flash.flash_time = TimeVector->at(OpFlashNumber->at(j));
         flash.flash_PEs = TotalPEVector->at(OpFlashNumber->at(j));
         flash.w_sumX += Channel_X->at(j) * OpHitPE_Flash->at(j);
         flash.w_sumY += Channel_Y->at(j) * OpHitPE_Flash->at(j);
         flash.w_sumZ += Channel_Z->at(j) * OpHitPE_Flash->at(j);
         flash.sum_PE += OpHitPE_Flash->at(j);
         flash.hit_vectorX.push_back(Channel_X->at(j));
         flash.hit_vectorY.push_back(Channel_Y->at(j));
         flash.hit_vectorZ.push_back(Channel_Z->at(j));
         flash.weight_hit_vector.push_back(OpHitPE_Flash->at(j));
         int match_counter = 0;
         int nchannels = flash.flash_channels.size();
         for (int q = 0; q < nchannels; q++)
         {
            if ((channel_vector[0] == flash.flash_channels.at(q)[0]) && (channel_vector[1] == flash.flash_channels.at(q)[1]) && (channel_vector[2] == flash.flash_channels.at(q)[2]))
            {
               flash.hit_counter++;
               break;
            }
            // if ((channel_vector[0] != flash.flash_channels.at(q)[0]) || (channel_vector[1] != flash.flash_channels.at(q)[1]) || (channel_vector[2] != flash.flash_channels.at(q)[2])){
            else
            {
               match_counter++;
            }
         }
         if (match_counter == nchannels)
         {
            flash.flash_channels.push_back(channel_vector);
            flash.channel_counter++;
         }

         //          flash.oldw_sumX += Channel_X->at(j)*OpHitPE_Flash->at(j);
         //          flash.oldw_sumY += Channel_Y->at(j)*OpHitPE_Flash->at(j);
         //          flash.oldw_sumZ += Channel_Z->at(j)*OpHitPE_Flash->at(j);
         //          flash.old_sumweight += OpHitPE_Flash->at(j);
      }

      // Compute means and errors for all flashes
      for (auto &[hitchannel, flash] : flashes)
      {
         //          if (hitchannel > 0){
         //             break;
         //          }
         //          oldw_wmeanX = flash.oldw_sumX/flash.old_sumweight;
         //          oldw_wmeanY = flash.oldw_sumY/flash.old_sumweight;
         //          oldw_wmeanZ = flash.oldw_sumZ/flash.old_sumweight;

         flash.meanX = flash.w_sumX / flash.sum_PE;
         flash.errX = rms_calc(flash.hit_vectorX, flash.weight_hit_vector, flash.meanX) / pow(flash.weight_hit_vector.size(), 0.5);

         flash.meanY = flash.w_sumY / flash.sum_PE;
         flash.errY = rms_calc(flash.hit_vectorY, flash.weight_hit_vector, flash.meanY) / pow(flash.weight_hit_vector.size(), 0.5);

         flash.meanZ = flash.w_sumZ / flash.sum_PE;
         flash.errZ = rms_calc(flash.hit_vectorZ, flash.weight_hit_vector, flash.meanZ) / pow(flash.weight_hit_vector.size(), 0.5);
         if (flash.plane == kC)
         {
            hist->Fill(flash.channel_counter);
         }
      }

      // Need to sort the flashes by "importance" (PE) desc
      std::vector<Flash_info> flashes_vec;

      for (auto const &[key, value] : flashes)
      {
         flashes_vec.push_back(value);
      }

      std::sort(flashes_vec.begin(), flashes_vec.end(), [](const Flash_info &a, const Flash_info &b)
                { return a.sum_PE > b.sum_PE; });

      std::vector<std::vector<Flash_info>> flash_clusters = coincidence_mean(flashes_vec);
      //       std::cout << "Nb clusters -> " << flash_clusters.size() << std::endl;

      for (auto const &cluster : flash_clusters)
      {
         flash_counter += cluster.size();
         // Need to filter the clusters here (ex: count number of planes)
         //  if(cluster.size() < 2 && cluster[0].flash_PEs < 215){ //can add && #PEs > mean to not filter all 1 flashes ( && cluster[0].flash_PEs < 215)
         //  counter++;
         //  continue;
         // }
         flash_counter_cut += cluster.size();

         std::vector<float> xmean_vec, ymean_vec, zmean_vec;
         std::vector<float> xweight_vec, yweight_vec, zweight_vec;
         float theta_C = 0;
         float theta_W = 0;
         float res_2 = 0;
         float avg_N_pes = 0;
         std::vector<float> npe_C;
         std::vector<float> npe_W;
         std::vector<float> N_pe_vec;

         for (Flash_info const &flash : cluster)
         {
            if (flash.plane == kC)
            {
               ymean_vec.push_back(flash.meanY);
               zmean_vec.push_back(flash.meanZ);
               yweight_vec.push_back(flash.w_sumY);
               zweight_vec.push_back(flash.w_sumZ);
            }
            else
            {
               xmean_vec.push_back(flash.meanX);
               zmean_vec.push_back(flash.meanZ);
               xweight_vec.push_back(flash.w_sumX);
               zweight_vec.push_back(flash.w_sumZ);
            }
         }

         float combined_x = combined_mean_calc(xmean_vec, xweight_vec);
         float combined_y = combined_mean_calc(ymean_vec, yweight_vec);
         float combined_z = combined_mean_calc(zmean_vec, zweight_vec);

         if (combined_x != 0 && combined_y != 0 && combined_z != 0)
         { // can change which coordinate I use to allow single flashes in specific planes (x for cathode and y for walls)
            for (Flash_info const &flash : cluster)
            {
               for (int k = 0; k < flash.weight_hit_vector.size(); k++)
               {
                  res_2 = pow(flash.hit_vectorX.at(k) - combined_x, 2) + pow(flash.hit_vectorY.at(k) - combined_y, 2) + pow(flash.hit_vectorZ.at(k) - combined_z, 2);
                  if (flash.plane == kC)
                  {
                     theta_C = TMath::ASin(abs(flash.hit_vectorX.at(k) - combined_x) / TMath::Sqrt(res_2));
                     N_pe_vec.push_back(calculate_PEs(flash.weight_hit_vector.at(k), res_2, theta_C));
                     npe_C.push_back(calculate_PEs(flash.weight_hit_vector.at(k), res_2, theta_C));
                  }
                  else
                  {
                     theta_W = TMath::ASin(abs(flash.hit_vectorY.at(k) - combined_y) / TMath::Sqrt(res_2));
                     N_pe_vec.push_back(calculate_PEs(flash.weight_hit_vector.at(k), res_2, theta_W));
                     npe_W.push_back(calculate_PEs(flash.weight_hit_vector.at(k), res_2, theta_W));
                  }
               }
            }
            avg_N_pes = std::accumulate(N_pe_vec.begin(), N_pe_vec.end(), 0) / N_pe_vec.size();

            global_PE_vec.push_back(avg_N_pes);
            trueE_vec.push_back(TrueE);
            //    hist->Fill(avg_N_pes);
         }

         float newres_x = 0;
         float newres_y = 0;
         float newres_z = 0;

         float combined_errx = 0;
         float combined_erry = 0;
         float combined_errz = 0;
         if (xmean_vec.size() > 1)
         {
            combined_errx = rms_calc(xmean_vec, xweight_vec, combined_x) / pow(xmean_vec.size(), 0.5);
         }
         else if (xmean_vec.size() == 1)
         {
            combined_errx = xweight_vec.at(0);
         }
         if (ymean_vec.size() > 1)
         {
            combined_erry = rms_calc(ymean_vec, yweight_vec, combined_y) / pow(ymean_vec.size(), 0.5);
         }
         else if (ymean_vec.size() == 1)
         {
            combined_erry = yweight_vec.at(0);
         }
         if (zmean_vec.size() > 1)
         {
            combined_errz = rms_calc(zmean_vec, zweight_vec, combined_z) / pow(zmean_vec.size(), 0.5);
         }
         else if (zmean_vec.size() == 1)
         {
            combined_errz = zweight_vec.at(0);
         }

         if (combined_errx > 0 && combined_x != 0)
         {
            newres_x = abs(combined_x - TrueX);
            bckg_X.push_back(combined_x);
            // hist -> Fill(newres_x);
         }
         if (combined_erry > 0 && combined_y != 0)
         {
            newres_y = abs(combined_y - TrueY);
            bckg_Y.push_back(combined_y);
            //    vector_newresy.push_back(newres_y);
            //    hist1 -> Fill(newres_y);
            RMS_vectorY.push_back(rms_calc(ymean_vec, yweight_vec, combined_y));
            TrueY_vector.push_back(TrueY);
            TrueX_vector.push_back(TrueX);
         }
         else
         {
            //    hist2D -> Fill(TrueY, TrueZ);
            //    hist2D2 -> Fill(TrueX, TrueY);
            //                   hist6 -> Fill(flash.sumC + flash.sumL +flash.sumR + flash.sumback + flash.sumfront);
         }
         if (combined_errz > 0 && combined_z != 0)
         {
            newres_z = abs(combined_z - TrueZ);
            //    hist2 -> Fill(newres_z);
            RMS_vectorZ.push_back(rms_calc(zmean_vec, zweight_vec, combined_z));
            TrueZ_vector.push_back(TrueZ);
         }
         else
         {
            // hist2D1 -> Fill(TrueY, TrueZ);
            // hist2D3 -> Fill(TrueX, TrueY);
         }
      }
      // if (global_PE_vec.size() != 0){
      // flash_mean_PEs = std::accumulate(global_PE_vec.begin(), global_PE_vec.end(), 0)/global_PE_vec.size();
      // event_PE_vec.push_back(flash_mean_PEs);
      // trueE_vec.push_back(TrueE);
      // }
      // if (counter == flash_clusters.size()){
      //    hist7 -> Fill(TrueE*1000);
      //    hist2D6 -> Fill(TrueX, TrueY);
      // }

      // if (Cut(ientry) < 0) continue;
      //       std::map<int, Flash_info> flashes1;
      //       for (int j = 0; j<OpHitPE_Flash->size(); j++){
      //          if(flashes1.count(OpFlashNumber->at(j)) == 0){
      //             flashes1[OpFlashNumber->at(j)] = Flash_info();
      //          }
      //          Flash_info &flash = flashes1[OpFlashNumber->at(j)];
      //          flash.flash_PEs += OpHitPE_Flash->at(j);

      // }
      //       for(const auto& [hitchannel, flash] : flashes1){
      //          hist -> Fill(flash.flash_PEs);
      //       }
   }
   hist->GetXaxis()->SetTitle("#X-ARAPUCA hit");
   hist->GetYaxis()->SetTitle("Entries");
   hist->SetTitle("#X-ARAPUCA per flash on the cathode, signal and background");
   //    hist -> Draw();
}
