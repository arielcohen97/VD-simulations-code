#define RMS_calculation_cxx
#include "RMS_calculation.h"
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
#include <TMath.h>
#include <math.h>

#define PI 3.1415926535897932384626
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
   float DPE_ratio = 0;
   float charge_calcX = 0;
   float mean_time_charge = 0;
   std::vector<Float_t> hit_vectorX;
   std::vector<Float_t> hit_vectorY;
   std::vector<Float_t> hit_vectorZ;
   std::vector<Float_t> weight_hit_vector;
   std::vector<std::vector<Float_t>> flash_channels;
   std::vector<Float_t> charge_recox;
   std::vector<Float_t> charge_time_vector;
   float meanX = 0;
   float meanY = 0;
   float meanZ = 0;
   float errX = 0;
   float errY = 0;
   float errZ = 0;
};

struct ChargeHit
{
   float time;
   float charge;
   float plane;
   float angle;
   float wireCenterX;
   float wireCenterY;
   float wireCenterZ;
};

float
rms_calc(std::vector<float> coord, std::vector<float> weight, float mean)
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
               // if(new_flash.plane != kC && new_flash.DPE_ratio < 35){
               // if (new_flash.plane != kC && new_flash.channel_counter < 5 && new_flash.hit_counter < 2)
               // if (new_flash.plane == kC && new_flash.channel_counter < 8 && new_flash.hit_counter < 3)
               // {
               //    {
               cluster.push_back(new_flash);
               found_cluster = true;
               break;
            }
            // }
            // }
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

void find_intersection(const ChargeHit a, const ChargeHit b, float& hitY, float& hitZ){
   
}

std::vector<std::pair<ChargeHit, ChargeHit>> match_hits(
   const std::vector<ChargeHit> &vec1,
   const std::vector<ChargeHit> &vec2,
   float dt){

   std::vector<std::pair<ChargeHit, ChargeHit>> matched_pairs;

   int idx1 = 0;
   int idx2 = 0;

   while(idx1 < vec1.size() && idx2 < vec2.size()){
      const ChargeHit &hit1 = vec1[idx1];
      const ChargeHit &hit2 = vec2[idx2];

      if(abs(hit1.time - hit2.time) < dt){
         matched_pairs.push_back({hit1, hit2});
         idx1++;
         idx2++;
      }

      else{
         if(hit1.time < hit2.time){
            idx1++;
         }
         else{
            idx2++; 
         }
      }
   }

   return matched_pairs;
}

std::vector<std::string> filenames = {
    "signal_only_charge.root"
    //  "bckg_v35_charge.root"
};

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

      TH1F *hist_loop;
      RMS_calculation calculator(tree);
      calculator.Loop(&hist_loop);
      hists[fname] = hist_loop;

      ofile->cd();
      if (aux == 0)
      {
         // hist_loop->Scale(scale * 30);
         hist_loop->SetLineColor(kRed);
         // legend->AddEntry(hist_loop, Form("Mean ADC charge, signal: 144.20"));
         legend->AddEntry(hist_loop, Form("#Mean amount of PEs, signal: 4.726e+04"));
      }
      if (aux == 1)
      {
         // hist_loop->Scale(9);
         // legend->AddEntry(hist_loop, Form("Mean ADC charge, background: 44.81"));
         legend->AddEntry(hist_loop, Form("#Mean amount of PEs, background: 4.524e+04"));
      }
      // std::cout << ofile << std::endl;
      hist_loop->Write(fname.c_str());
      aux++;
      f->Close();
   }

   ofile->Close();

   // TCanvas *c = new TCanvas;
   bool isFirst = true;
   for (auto const &[fname, hist_loop] : hists)
   {
      const char *drawOpt = isFirst ? "HIST" : "SAME HIST";
      isFirst = false;
      // hist_loop->Draw(drawOpt);
   }
   // legend->Draw();
}

std::vector<float> global_PE_vec;
std::vector<float> trueE_vec;
std::vector<float> event_PE_vec;
float flash_counter = 0;
float flash_counter_cut = 0;
float scale;
std::vector<Float_t> channel_vector(3);

void RMS_calculation::Loop(TH1F **ret_hist)

{
   TFile *f = new TFile("RMS_calculation.root", "RECREATE");

   if (fChain == 0)
      return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   TH1F *hist = new TH1F("hist", "Resolution for 2D reconstruction", 100, 0, 500);
   TH1F *hist1 = new TH1F("hist1", "Resolution for 2D reconstruction", 100, 0, 500);
   TH1F *hist2 = new TH1F("hist2", "Resolution for 2D reconstruction", 100, 0, 500);
   TH1F *hist3 = new TH1F("hist3", "Resolution for 3D reconstruction", 100, 0, 0);
   TH1F *hist4 = new TH1F("hist4", "Resolution for 3D reconstruction", 100, 0, 0);
   TH1F *hist5 = new TH1F("hist5", "Resolution for 3D reconstruction", 100, 0, 0);
   TH1F *hist6 = new TH1F("hist6", "Error 0 points #PEs", 100, 0, 0);
   TH1F *hist7 = new TH1F("hist7", "TrueE missed events", 200, 0, 200000);
   TH1F *hist8 = new TH1F("hist8", "X-ARAPUCA counter, cathode", 100, 0, 0);
   TH1F *hist9 = new TH1F("hist9", "X-ARAPUCA counter, wall", 100, 0, 0);
   TH1F *hist10 = new TH1F("hist10", "X-ARAPUCA counter, wall", 100, 0, 0);
   TH1F *hist_loop = new TH1F("hist_loop", "X-ARAPUCA counter, wall", 100, 0, 0);
   *ret_hist = hist_loop;
   TH2F *hist2D = new TH2F("hist2D", "True yz", 100, 0, 0, 100, 0, 2200);
   TH2F *hist2D1 = new TH2F("hist2D1", "True yz", 100, 0, 0, 100, 0, 2200);
   TH2F *hist2D2 = new TH2F("hist2D2", "True xy", 100, 0, 0, 100, 0, 2200);
   TH2F *hist2D3 = new TH2F("hist2D3", "True xy", 100, 0, 0, 100, 0, 2200);
   TH2F *hist2D4 = new TH2F("hist2D4", "Reco yz background", 100, 0, 0, 100, 0, 2200);
   TH2F *hist2D5 = new TH2F("hist2D5", "Reco xz background", 100, 0, 0, 100, 0, 2200);
   TH2F *hist2D6 = new TH2F("hist2D", "True yz", 100, 0, 0, 100, 0, 0);
   hist->SetDirectory(0);
   hist1->SetDirectory(0);
   hist2->SetDirectory(0);
   hist3->SetDirectory(0);
   hist4->SetDirectory(0);
   hist5->SetDirectory(0);
   hist6->SetDirectory(0);
   hist7->SetDirectory(0);
   hist8->SetDirectory(0);
   hist9->SetDirectory(0);
   hist10->SetDirectory(0);
   hist_loop->SetDirectory(0);
   hist2D->SetDirectory(0);
   hist2D1->SetDirectory(0);
   hist2D2->SetDirectory(0);
   hist2D3->SetDirectory(0);
   hist2D4->SetDirectory(0);
   hist2D5->SetDirectory(0);
   hist2D6->SetDirectory(0);
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      // if (Cut(ientry) < 0) continue;

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
      float vdrift = 0.16; // cm/us
      scale = 1. / nentries;
      // std::cout << scale << std::endl;
      // for (int k=0; k < Charge_adc->size(); k++){
      //    hist_loop -> Fill(Charge_adc->at(k));
      // }

      std::map<int, Flash_info> flashes;
      std::map<int, std::vector<ChargeHit>> charges;

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
         hist10->Fill(OpHittime_Flash->at(j) - flash.flash_time);
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

      for (int k = 0; k < Charge_time->size(); k++)
      {
         ChargeHit hit;
         hit.time = Charge_time->at(k);
         hit.charge = Charge_adc->at(k);
         hit.plane = Wire_plane->at(k);
         hit.wireCenterX = Wire_Xcenter->at(k);
         hit.wireCenterY = Wire_Ycenter->at(k);
         hit.wireCenterZ = Wire_Zcenter->at(k);
         hit.angle = Wire_angle->at(k);

         if(charges.count(hit.plane) == 0){
            charges[hit.plane] = {};
         }
         charges[hit.plane].push_back(hit);
      }

      // Ordering charge hits in time
      for(auto &[plane, vec] : charges)
      {
         std::sort(vec.begin(), vec.end(), [](const ChargeHit &a, const ChargeHit &b)
                { return a.time > b.time; });
      }

      //Making the associations

      std::vector<std::pair<ChargeHit, ChargeHit>> matched_hits = match_hits(charges[0], charges[1], 10);

      // 

      // Compute means and errors for all flashes
      for (auto &[hitchannel, flash] : flashes)
      {
         flash.meanX = flash.w_sumX / flash.sum_PE;
         flash.errX = rms_calc(flash.hit_vectorX, flash.weight_hit_vector, flash.meanX) / pow(flash.weight_hit_vector.size(), 0.5);

         flash.meanY = flash.w_sumY / flash.sum_PE;
         flash.errY = rms_calc(flash.hit_vectorY, flash.weight_hit_vector, flash.meanY) / pow(flash.weight_hit_vector.size(), 0.5);

         flash.meanZ = flash.w_sumZ / flash.sum_PE;
         flash.errZ = rms_calc(flash.hit_vectorZ, flash.weight_hit_vector, flash.meanZ) / pow(flash.weight_hit_vector.size(), 0.5);
         flash.DPE_ratio = flash.flash_PEs / flash.channel_counter;

         for (const auto& mypair : matched_hits){
            // std::cout << "Plane angle 1: " << mypair.first.angle << " Plane angle 2: " << mypair.second.angle << std::endl;
            float pair_z = (mypair.second.wireCenterY - mypair.first.wireCenterY);
            pair_z += (TMath::Tan(PI/2. - mypair.first.angle)*mypair.first.wireCenterZ - TMath::Tan(PI/2. - mypair.second.angle)*mypair.second.wireCenterZ);
            pair_z /= TMath::Tan(PI/2. - mypair.first.angle) - TMath::Tan(PI/2. - mypair.second.angle);
            float pair_y = TMath::Tan(PI/2. - mypair.second.angle)*(pair_z- mypair.second.wireCenterZ) + (mypair.second.wireCenterY);
            // std::cout << "Pair Y: " << pair_y << " Pair Z: " << pair_z << " True Z: " <<  TrueZ << " TrueY: " << TrueY << std::endl;
            if (abs(flash.meanY - pair_y) < 100 && abs(flash.meanZ - pair_z) < 100){
               flash.mean_time_charge = (mypair.first.time + mypair.second.time)/2;
               flash.charge_calcX = (abs(flash.mean_time_charge - flash.flash_time)) * vdrift;
               flash.charge_recox.push_back(flash.charge_calcX);
               hist2D6->Fill(TrueX, flash.charge_calcX);
            }
         }
         
         // if (flash.plane == kC)
         // {
         //    // hist_loop->Fill(flash.DPE_ratio);
         // }
         // else
         // {
         //    hist9->Fill(flash.DPE_ratio);
         // }
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
         //  if(cluster.size() < 2 && cluster[0].flash_PEs < 217){ //can add && #PEs > mean to not filter all 1 flashes ( && cluster[0].flash_PEs < 215)
         //  counter++;
         //  continue;
         //  }
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
               hist2D5->Fill(flash.meanZ, flash.meanY);
            }
            else
            {
               xmean_vec.push_back(flash.meanX);
               zmean_vec.push_back(flash.meanZ);
               xweight_vec.push_back(flash.w_sumX);
               zweight_vec.push_back(flash.w_sumZ);
               hist2D4->Fill(flash.meanZ, flash.meanX);
            }
         }

         float combined_x = combined_mean_calc(xmean_vec, xweight_vec);
         float combined_y = combined_mean_calc(ymean_vec, yweight_vec);
         float combined_z = combined_mean_calc(zmean_vec, zweight_vec);

         if (cluster.size() == 1 && combined_x == 0)
         {
            for (int k = 0; k < cluster[0].weight_hit_vector.size(); k++)
            {
               res_2 = pow(cluster[0].hit_vectorX.at(k) - (cluster[0].charge_calcX - 327), 2) + pow(cluster[0].hit_vectorY.at(k) - cluster[0].meanY, 2) + pow(cluster[0].hit_vectorZ.at(k) - cluster[0].meanZ, 2);
               if (cluster[0].plane == kC)
               {
                  theta_C = TMath::ASin(abs(cluster[0].hit_vectorX.at(k) - (cluster[0].charge_calcX - 327)) / TMath::Sqrt(res_2));
                  N_pe_vec.push_back(calculate_PEs(cluster[0].weight_hit_vector.at(k), res_2, theta_C));
                  npe_C.push_back(calculate_PEs(cluster[0].weight_hit_vector.at(k), res_2, theta_C));
               }
            }

            avg_N_pes = std::accumulate(N_pe_vec.begin(), N_pe_vec.end(), 0) / N_pe_vec.size();

            global_PE_vec.push_back(avg_N_pes);
            trueE_vec.push_back(TrueE);
            // hist2D6->Fill(TrueE * 1000, avg_N_pes);
            if (TrueE * 1000 > 9 && TrueE * 1000 < 11)
            {
               // hist7->Fill(avg_N_pes);
            }
         }

         if (combined_x != 0 && combined_y != 0 && combined_z != 0)
         { // can change which coordinate I use to allow single flashes in specific planes (x for cathode and y for walls)
            for (Flash_info const &flash : cluster)
            {
               for (int k = 0; k < flash.weight_hit_vector.size(); k++)
               {
                  res_2 = pow(flash.hit_vectorX.at(k) - (cluster[0].charge_calcX - 327), 2) + pow(flash.hit_vectorY.at(k) - combined_y, 2) + pow(flash.hit_vectorZ.at(k) - combined_z, 2);
                  if (flash.plane == kC)
                  {
                     theta_C = TMath::ASin(abs(flash.hit_vectorX.at(k) - (cluster[0].charge_calcX - 327)) / TMath::Sqrt(res_2));
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
            // hist2D6->Fill(TrueE * 1000, avg_N_pes);
            if (TrueE * 1000 > 68 && TrueE * 1000 < 70)
            {
               hist7->Fill(avg_N_pes);
            }
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
            hist->Fill(newres_x);
         }
         if (combined_erry > 0 && combined_y != 0)
         {
            newres_y = abs(combined_y - TrueY);
            bckg_Y.push_back(combined_y);
            //    vector_newresy.push_back(newres_y);
            hist1->Fill(newres_y);
            RMS_vectorY.push_back(rms_calc(ymean_vec, yweight_vec, combined_y));
            TrueY_vector.push_back(TrueY);
            TrueX_vector.push_back(TrueX);
         }
         else
         {
            hist2D->Fill(TrueY, TrueZ);
            hist2D2->Fill(TrueX, TrueY);
            //                   hist6 -> Fill(flash.sumC + flash.sumL +flash.sumR + flash.sumback + flash.sumfront);
         }
         if (combined_errz > 0 && combined_z != 0)
         {
            newres_z = abs(combined_z - TrueZ);
            hist2->Fill(newres_z);
            RMS_vectorZ.push_back(rms_calc(zmean_vec, zweight_vec, combined_z));
            TrueZ_vector.push_back(TrueZ);
         }
         else
         {
            hist2D1->Fill(TrueY, TrueZ);
            hist2D3->Fill(TrueX, TrueY);
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
      // std::cout << "NEXT EVENT!!!!!!" << std::endl;
   }
   // std::cout << "Signal flash rate per event: " << " " << flash_counter/nentries << std::endl;
   // std::cout << "Signal flash rate per event after cut: " << " " << flash_counter_cut/nentries << std::endl;
   // Plot new 2D reco

   auto c = new TCanvas("histogram", "True yz, erry = 0");

   // TGraph *gr4 = new TGraph(trueE_vec.size(), &trueE_vec[0], &global_PE_vec[0]);
   // gr4 -> SetTitle("Reconstructed #PEs vs true neutrino energy");
   // gr4 -> GetXaxis() -> SetTitle("True energy (MeV)");
   // gr4 -> GetYaxis() -> SetTitle("NPEs");
   // gr4 -> Draw("AP");

   // hist2 -> GetXaxis()->SetTitle("Resolution (cm)");
   // hist2 -> SetTitle("Resolution for 2D reconstruction, largest flash");
   // hist2 -> Draw();
   // hist1 -> SetLineColor(kRed);
   // hist1 -> GetXaxis()->SetTitle("Resolution (cm)");
   // hist1 -> Draw("SAME");
   // hist -> SetLineColor(kGreen);
   // hist -> GetXaxis()->SetTitle("Resolution (cm)");
   // hist -> Draw("SAME");
   // hist2D4 -> GetXaxis() -> SetTitle("Reco Z");
   // hist2D4 -> GetYaxis() -> SetTitle("Reco X");
   // hist2D4 -> SetTitle("Reconstructed signal flashes, XZ plane, flash + PE + PE/Ara ratio cut");
   // hist2D4 -> Draw("COLZ");
   hist2D6->GetXaxis()->SetTitle("#True X (cm)");
   hist2D6->GetYaxis()->SetTitle("RecoX (cm)");
   hist2D6->SetTitle("Reconstructed drift position with charge vs true neutrino position, signal");
   hist2D6->Draw("COLZ");
   // hist8 -> SetTitle("Background, cathode");
   // hist8 -> GetXaxis()->SetTitle("Repeated hits");
   // hist8 -> Draw();
   // hist10->SetTitle("Flash time profile, background");
   // hist10->GetXaxis()->SetTitle("(Hit time - T_centroid) [us]");
   // hist10->Draw();

   // auto legend1 = new TLegend();
   // legend1 -> AddEntry(hist, Form("New mean X:  %5.2f (cm)", hist->GetMean())); //Form("New mean: " %5.2f, std::to_string(hist->GetMean()).c_str())
   // legend1 -> AddEntry(hist1, Form("New mean Y:  %5.2f (cm)", hist1->GetMean()));
   // legend1 -> AddEntry(hist2, Form("New mean Z:  %5.2f (cm)", hist2->GetMean()));
   // legend1 -> Draw();

   // Plot old 3D reco

   // auto h = new TCanvas("histogram3", "True yz, errz = 0");
   //       hist5 -> GetXaxis()->SetTitle("Resolution (cm)");
   //       hist5 -> Draw();
   //       hist3 -> SetLineColor(kGreen);
   //       hist3 -> GetXaxis()->SetTitle("Resolution (cm)");
   //       hist3 -> Draw("SAME");
   //       hist4 -> Draw("SAME");
   //       hist4 -> GetXaxis()->SetTitle("Resolution (cm)");
   //       hist4 -> SetLineColor(kRed);
   //    hist2D1 -> Draw("COLZ");
   //    TGraph *gr3 = new TGraph(vector_newresy.size(), &TrueY_vector[0], &vector_newresy[0]);
   //    gr3 -> SetTitle("Resolution Y vs true neutrino y");
   //    gr3 -> GetXaxis() -> SetTitle("True Y (cm)");
   //    gr3 -> GetYaxis() -> SetTitle("Resolution Y (cm)");
   //    gr3 -> Draw("AP");
   // hist2D5 -> GetXaxis() -> SetTitle("Reco Z");
   // hist2D5 -> GetYaxis() -> SetTitle("Reco Y");
   // hist2D5 -> SetTitle("Reconstructed signal flashes, YZ plane, flash + PE + PE/Ara ratio cut");
   // hist2D5 -> Draw("COLZ");
   // hist7->Scale(scale);
   // hist7->SetTitle("#NPEs per drift window, background");
   // hist7->GetXaxis()->SetTitle("#PEs");
   // hist7->Draw("HIST");
   // hist9 -> SetTitle("Background, cathode");
   // hist9 -> GetXaxis()->SetTitle("Repeated hits");
   // hist9 -> Draw();

   // auto legend2 = new TLegend();
   // legend2 -> AddEntry(hist7, Form("Avg flashes per drift window:  %5.2f", hist7->Integral())); //Form("New mean: " %5.2f, std::to_string(hist->GetMean()).c_str());
   // legend2 -> Draw();

   // auto h2 = new TCanvas("histogram2", "True xy, erry = 0");
   //    //hist2D2 -> Draw("COLZ");
   //    TGraph *gr1 = new TGraph(RMS_vectorY.size(), &TrueY_vector[0], &RMS_vectorY[0]);
   //    gr1 -> SetTitle("RMS y vs true neutrino y");
   //    gr1 -> GetXaxis() -> SetTitle("True Y (cm)");
   //    gr1 -> GetYaxis() -> SetTitle("RMS Y");
   //    gr1 -> Draw("AP");
   // hist_loop->GetXaxis()->SetTitle("#PEs");
   // hist_loop->GetYaxis()->SetTitle("Entries");
   // hist_loop->SetTitle("#NPEs per drift window, signal and background");

   //    auto h1 = new TCanvas("histogram1", "True xy, errz = 0");
   //    hist2D3 -> Draw("COLZ");
   // hist6 -> Draw();
   //    TGraph *gr2 = new TGraph(RMS_vectorY.size(), &TrueX_vector[0], &RMS_vectorY[0]);
   //    gr2 -> SetTitle("RMS y vs true neutrino x");
   //    gr2 -> GetXaxis() -> SetTitle("True X (cm)");
   //    gr2 -> GetYaxis() -> SetTitle("RMS Y");
   //    gr2 -> Draw("AP");

   //    auto legend1 = new TLegend();
   //    legend1 -> AddEntry(hist3, Form("New mean X:  %5.2f (cm)", hist3->GetMean())); //Form("New mean: " %5.2f, std::to_string(hist->GetMean()).c_str())
   //    legend1 -> AddEntry(hist4, Form("New mean Y:  %5.2f (cm)", hist4->GetMean()));
   //    legend1 -> AddEntry(hist5, Form("New mean Z:  %5.2f (cm)", hist5->GetMean()));
   //    legend1 -> Draw();
   //    auto gr = new TGraph(n, x, y);

   f->Close();
}
