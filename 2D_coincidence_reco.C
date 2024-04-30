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

enum Plane{
  kC,
  kRW,
  kLW,
  kNoPlane
};



struct Flash_info {
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
      std::vector<Float_t> hit_vectorX;
      std::vector<Float_t> hit_vectorY;
      std::vector<Float_t> hit_vectorZ;
      std::vector<Float_t> weight_hit_vector;
      float meanX = 0;
      float meanY = 0;
      float meanZ = 0;
      float errX = 0;
      float errY = 0;
      float errZ = 0;


};

float rms_calc(std::vector<float> coord, std::vector<float> weight, float mean){
   float value = 0;
   float sum_weights = 0;
   for(int k=0; k<coord.size();k++){
      value += pow((coord.at(k) - mean), 2)*weight.at(k);
      sum_weights += weight.at(k);
   }
   if (sum_weights > 0){
      value = pow(value/sum_weights, 0.5);
   }
   else{
      value = 0;
   }
   return value;

}

float combined_mean_calc(std::vector<float> coord, std::vector<float> weight){
   float value = 0;
   float sum_weights = 0;
   for (int i=0; i<coord.size();i++){
      value += coord.at(i)*weight.at(i);
      sum_weights += weight.at(i);
   }
   if (sum_weights > 0){
      value = value/sum_weights;
   }
   else{

      value = 0;
   }
   return value;

}

//Groups the flashes in clusters with nearby Z and time means.
std::vector<std::vector<Flash_info>> coincidence_mean(std::vector<Flash_info> flashes){
//    return {flashes};
   std::vector<std::vector<Flash_info>> clusters;
   while(flashes.size() != 0){
      Flash_info new_flash = flashes[0];
      flashes.erase(flashes.begin());
      bool found_cluster = false;
      for(auto& cluster : clusters){
         if (abs(cluster[0].flash_time - new_flash.flash_time) <= 1 && abs(cluster[0].meanZ - new_flash.meanZ) <= 400){
            cluster.push_back(new_flash);
            found_cluster = true;
            break;
         }
      }
      if(!found_cluster){
         clusters.push_back({new_flash});
      }
   }
   return clusters;
}

Plane getPlane(float chX, float chY, float chZ){
   if (chX < -320){
      return kC;
   }
   if (chY < -700){
      return kLW;
   }
   if(chY > 700){
      return kRW;
   }
   return kNoPlane;
}



void RMS_calculation::Loop()



{
   TFile *f=new TFile("RMS_calculation.root","RECREATE");
//   In a ROOT session, you can do:
//      root> .L RMS_calc.C
//      root> RMS_calc t
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
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


   Long64_t nbytes = 0, nb = 0;
   TH1F* hist = new TH1F("hist","Resolution for 2D reconstruction", 100, 0, 500);
   TH1F* hist1 = new TH1F("hist1","Resolution for 2D reconstruction", 100, 0, 500);
   TH1F* hist2 = new TH1F("hist2","Resolution for 2D reconstruction", 100, 0, 500);
   TH1F* hist3 = new TH1F("hist3","Resolution for 3D reconstruction", 100, 0, 0);
   TH1F* hist4 = new TH1F("hist4","Resolution for 3D reconstruction", 100, 0, 0);
   TH1F* hist5 = new TH1F("hist5","Resolution for 3D reconstruction", 100, 0, 0);
   TH1F* hist6 = new TH1F("hist6","Error 0 points #PEs", 100, 0, 0);
   TH2F* hist2D = new TH2F("hist2D", "True yz", 100, 0, 0, 100, 0, 2200);
   TH2F* hist2D1 = new TH2F("hist2D1", "True yz", 100, 0, 0, 100, 0, 2200);
   TH2F* hist2D2 = new TH2F("hist2D2", "True xy", 100, 0, 0, 100, 0, 2200);
   TH2F* hist2D3 = new TH2F("hist2D3", "True xy", 100, 0, 0, 100, 0, 2200);
   TH2F* hist2D4 = new TH2F("hist2D4", "Reco yz background", 100, 0, 0, 100, 0, 2200);
   TH2F* hist2D5 = new TH2F("hist2D5", "Reco xz background", 100, 0, 0, 100, 0, 2200);
   hist -> SetDirectory(0);
   hist1 -> SetDirectory(0);
   hist2 -> SetDirectory(0);
   hist3 -> SetDirectory(0);
   hist4 -> SetDirectory(0);
   hist5 -> SetDirectory(0);
   hist6 -> SetDirectory(0);
   hist2D -> SetDirectory(0);
   hist2D1 -> SetDirectory(0);
   hist2D2 -> SetDirectory(0);
   hist2D3 -> SetDirectory(0);
   hist2D4 -> SetDirectory(0);
   hist2D5 -> SetDirectory(0);
   for (Long64_t jentry=0; jentry<1;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
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


      std::map<int, Flash_info> flashes;
      for (int j = 0; j<OpHitChannel_Flash->size(); j++){
         if(flashes.count(OpFlashNumber->at(j)) == 0){
            flashes[OpFlashNumber->at(j)] = Flash_info();
            flashes[OpFlashNumber->at(j)].plane = getPlane(Channel_X->at(j), Channel_Y->at(j), Channel_Z->at(j));
         }
         Flash_info &flash = flashes[OpFlashNumber->at(j)];

         flash.flash_time = TimeVector->at(OpFlashNumber->at(j));
         flash.w_sumX += Channel_X->at(j)*OpHitPE_Flash->at(j);
         flash.w_sumY += Channel_Y->at(j)*OpHitPE_Flash->at(j);
         flash.w_sumZ += Channel_Z->at(j)*OpHitPE_Flash->at(j);
         flash.sum_PE += OpHitPE_Flash->at(j);
         flash.hit_vectorX.push_back(Channel_X->at(j));
         flash.hit_vectorY.push_back(Channel_Y->at(j));
         flash.hit_vectorZ.push_back(Channel_Z->at(j));
         flash.weight_hit_vector.push_back(OpHitPE_Flash->at(j));

//          flash.oldw_sumX += Channel_X->at(j)*OpHitPE_Flash->at(j);
//          flash.oldw_sumY += Channel_Y->at(j)*OpHitPE_Flash->at(j);
//          flash.oldw_sumZ += Channel_Z->at(j)*OpHitPE_Flash->at(j);
//          flash.old_sumweight += OpHitPE_Flash->at(j);

      }

      //Compute means and errors for all flashes
      for(auto& [hitchannel, flash] : flashes){
//          if (hitchannel > 0){
//             break;
//          }
//          oldw_wmeanX = flash.oldw_sumX/flash.old_sumweight;
//          oldw_wmeanY = flash.oldw_sumY/flash.old_sumweight;
//          oldw_wmeanZ = flash.oldw_sumZ/flash.old_sumweight;

            flash.meanX = flash.w_sumX/flash.sum_PE;
            flash.errX = rms_calc(flash.hit_vectorX, flash.weight_hit_vector, flash.meanX)/pow(flash.weight_hit_vector.size(), 0.5);

            flash.meanY = flash.w_sumY/flash.sum_PE;
            flash.errY = rms_calc(flash.hit_vectorY, flash.weight_hit_vector, flash.meanY)/pow(flash.weight_hit_vector.size(), 0.5);

            flash.meanZ = flash.w_sumZ/flash.sum_PE;
            flash.errZ = rms_calc(flash.hit_vectorZ, flash.weight_hit_vector, flash.meanZ)/pow(flash.weight_hit_vector.size(), 0.5);



      }

      //Need to sort the flashes by "importance" (PE) desc
      std::vector<Flash_info> flashes_vec;

      for(auto const& [key, value] : flashes){
         flashes_vec.push_back(value);
      }

      std::sort(flashes_vec.begin(), flashes_vec.end(), [](const Flash_info &a, const Flash_info &b){
         return a.sum_PE > b.sum_PE;
      });

      std::vector<std::vector<Flash_info>> flash_clusters = coincidence_mean(flashes_vec);
//       std::cout << "Nb clusters -> " << flash_clusters.size() << std::endl;

      for(auto const& cluster : flash_clusters){
         //Need to filter the clusters here (ex: count number of planes)
         if(cluster.size() < 2){ //can add && #PEs > mean to not filter all 1 flashes  
            continue;
         }


         std::vector<float> xmean_vec, ymean_vec, zmean_vec;
         std::vector<float> xweight_vec, yweight_vec, zweight_vec;

         for(Flash_info const& flash: cluster){
            if(flash.plane == kC){
               ymean_vec.push_back(flash.meanY);
               zmean_vec.push_back(flash.meanZ);
               yweight_vec.push_back(flash.w_sumY);
               zweight_vec.push_back(flash.w_sumZ);
               hist2D5 -> Fill(flash.meanY, flash.meanZ);
            }
            else{
               xmean_vec.push_back(flash.meanX);
               zmean_vec.push_back(flash.meanZ);
               xweight_vec.push_back(flash.w_sumX);
               zweight_vec.push_back(flash.w_sumZ);
               hist2D4 -> Fill(flash.meanX, flash.meanZ);
            }
         }

         float combined_x = combined_mean_calc(xmean_vec, xweight_vec);
         float combined_y = combined_mean_calc(ymean_vec, yweight_vec);
         float combined_z = combined_mean_calc(zmean_vec, zweight_vec);

         // if (combined_y == 0){
         //    hist2D4 -> Fill(combined_x, combined_z);
         // }
         // if(combined_y != 0){
         //    hist2D5 -> Fill(combined_y, combined_z);
         // }

         float newres_x = 0;
         float newres_y = 0;
         float newres_z = 0;

         float combined_errx = 0;
         float combined_erry = 0;
         float combined_errz = 0;
         if (xmean_vec.size() > 1){
            combined_errx = rms_calc(xmean_vec, xweight_vec, combined_x)/pow(xmean_vec.size(), 0.5);
         }
         else if(xmean_vec.size() == 1){
            combined_errx = xweight_vec.at(0);

         }
         if (ymean_vec.size() > 1){
            combined_erry = rms_calc(ymean_vec, yweight_vec, combined_y)/pow(ymean_vec.size(), 0.5);
         }
         else if(ymean_vec.size() == 1){
            combined_erry = yweight_vec.at(0);

         }
         if (zmean_vec.size() > 1){
            combined_errz = rms_calc(zmean_vec, zweight_vec, combined_z)/pow(zmean_vec.size(), 0.5);
         }
         else if(zmean_vec.size() == 1){
            combined_errz = zweight_vec.at(0);

         }

         if(combined_errx > 0 && combined_x != 0){
            newres_x = abs(combined_x - TrueX);
            bckg_X.push_back(combined_x);
//             std::cout << newres_x << " " << TrueX << " " << combined_x << std::endl;
            hist -> Fill(newres_x);
         }
         if(combined_erry > 0 && combined_y != 0){
               newres_y = abs(combined_y - TrueY);
               bckg_Y.push_back(combined_y);
            //    vector_newresy.push_back(newres_y);
               hist1 -> Fill(newres_y);
               RMS_vectorY.push_back(rms_calc(ymean_vec, yweight_vec, combined_y));
               TrueY_vector.push_back(TrueY);
               TrueX_vector.push_back(TrueX);

            }
            else{
               hist2D -> Fill(TrueY, TrueZ);
               hist2D2 -> Fill(TrueX, TrueY);
//                   hist6 -> Fill(flash.sumC + flash.sumL +flash.sumR + flash.sumback + flash.sumfront);

         }
         if(combined_errz > 0 && combined_z != 0){
               newres_z = abs(combined_z - TrueZ);
               hist2 -> Fill(newres_z);
               RMS_vectorZ.push_back(rms_calc(zmean_vec, zweight_vec, combined_z));
               TrueZ_vector.push_back(TrueZ);
            }
         else{
            hist2D1 -> Fill(TrueY, TrueZ);
            hist2D3 -> Fill(TrueX, TrueY);

         }

      }
   }

   //Plot new 2D reco
   auto c = new TCanvas("histogram", "True yz, erry = 0");
//    TGraph *gr4 = new TGraph(vector_newresy.size(), &vector_newresy[0], &RMS_vectorY[0]);
//    gr4 -> SetTitle("RMS Y vs Resolution Y");
//    gr4 -> GetXaxis() -> SetTitle("Resolution Y (cm)");
//    gr4 -> GetYaxis() -> SetTitle("RMS Y");
//    gr4 -> Draw("AP");

      // hist2 -> GetXaxis()->SetTitle("Resolution (cm)");
      // hist2 -> SetTitle("Resolution for 2D reconstruction, largest flash");
      // hist2 -> Draw();
      // hist1 -> SetLineColor(kRed);
      // hist1 -> GetXaxis()->SetTitle("Resolution (cm)");
      // hist1 -> Draw("SAME");
      // hist -> SetLineColor(kGreen);
      // hist -> GetXaxis()->SetTitle("Resolution (cm)");
      // hist -> Draw("SAME");
         hist2D4 -> GetXaxis() -> SetTitle("Reco X");
         hist2D4 -> GetYaxis() -> SetTitle("Reco Z");
         hist2D4 -> SetTitle("Reconstructed background events, XZ plane");
         hist2D4 -> Draw("COLZ");

   // auto legend1 = new TLegend();
   // legend1 -> AddEntry(hist, Form("New mean X:  %5.2f (cm)", hist->GetMean())); //Form("New mean: " %5.2f, std::to_string(hist->GetMean()).c_str())
   // legend1 -> AddEntry(hist1, Form("New mean Y:  %5.2f (cm)", hist1->GetMean()));
   // legend1 -> AddEntry(hist2, Form("New mean Z:  %5.2f (cm)", hist2->GetMean()));
   // legend1 -> Draw();

   //Plot old 3D reco

   auto h = new TCanvas("histogram3", "True yz, errz = 0");
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
         hist2D5 -> GetXaxis() -> SetTitle("Reco Y");
         hist2D5 -> GetYaxis() -> SetTitle("Reco Z");
         hist2D5 -> SetTitle("Reconstructed background events, YZ plane");
         hist2D5 -> Draw("COLZ");


//    auto h2 = new TCanvas("histogram2", "True xy, erry = 0");
//    //hist2D2 -> Draw("COLZ");
//    TGraph *gr1 = new TGraph(RMS_vectorY.size(), &TrueY_vector[0], &RMS_vectorY[0]);
//    gr1 -> SetTitle("RMS y vs true neutrino y");
//    gr1 -> GetXaxis() -> SetTitle("True Y (cm)");
//    gr1 -> GetYaxis() -> SetTitle("RMS Y");
//    gr1 -> Draw("AP");

//    auto h1 = new TCanvas("histogram1", "True xy, errz = 0");
//    hist2D3 -> Draw("COLZ");
   //hist6 -> Draw();
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
