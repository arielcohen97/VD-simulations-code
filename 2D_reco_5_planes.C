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

float oldw_wmeanX;
float oldw_wmeanY;
float oldw_wmeanZ;
float w_meanXR;
float w_meanXL;
float w_meanYC;
float w_meanZC;
float w_meanZR;
float w_meanZL;
float w_meanYback;
float w_meanYfront;
float w_meanXback;
float w_meanXfront;
float w_errXR;
float w_errXL;
float w_errYC;
float w_errZC;
float w_errZR;
float w_errZL;
float w_errYback;
float w_errYfront;
float w_errXback;
float w_errXfront;
float w_sum;
float sum;
float res_C;
float res_R_wall;
float res_L_wall;
float res_back_wall;
float res_front_wall;
float res_3D;
float newres_x;
float newres_y;
float newres_z;
float oldres_x;
float oldres_y;
float oldres_z;
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


struct Flash_info {
      float w_sumXR = 0;
      float w_sumXL = 0;
      float w_sumYC = 0;
      float w_sumZC = 0;
      float w_sumZL = 0;
      float w_sumZR = 0;
      float w_sumXback = 0;
      float w_sumXfront = 0;
      float w_sumYback = 0;
      float w_sumYfront = 0;
      float sumC = 0;
      float sumR = 0;
      float sumL = 0;
      float sumback = 0;
      float sumfront = 0;
      float oldw_sumX = 0;
      float oldw_sumY = 0;
      float oldw_sumZ = 0;
      float old_sumweight = 0;
      float min_dist_C = 1000000;
      std::vector<Float_t> hit_vectorXR;
      std::vector<Float_t> hit_vectorXL;
      std::vector<Float_t> hit_vectorZC;
      std::vector<Float_t> hit_vectorZR;
      std::vector<Float_t> hit_vectorZL;
      std::vector<Float_t> hit_vectorYC;
      std::vector<Float_t> hit_vectorXback;
      std::vector<Float_t> hit_vectorXfront;
      std::vector<Float_t> hit_vectorYback;
      std::vector<Float_t> hit_vectorYfront;
      std::vector<Float_t> weight_hit_vectorXR;
      std::vector<Float_t> weight_hit_vectorXL;
      std::vector<Float_t> weight_hit_vectorZC;
      std::vector<Float_t> weight_hit_vectorZR;
      std::vector<Float_t> weight_hit_vectorZL;
      std::vector<Float_t> weight_hit_vectorYC;
      std::vector<Float_t> weight_hit_vectorXback;
      std::vector<Float_t> weight_hit_vectorXfront;
      std::vector<Float_t> weight_hit_vectorYback;
      std::vector<Float_t> weight_hit_vectorYfront;

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
   TH1F* hist = new TH1F("hist","Resolution for 2D reconstruction", 100, 0, 0);
   TH1F* hist1 = new TH1F("hist1","Resolution for 2D reconstruction", 100, 0, 0);
   TH1F* hist2 = new TH1F("hist2","Resolution for 2D reconstruction", 100, 0, 0);
   TH1F* hist3 = new TH1F("hist3","Resolution for 3D reconstruction", 100, 0, 0);
   TH1F* hist4 = new TH1F("hist4","Resolution for 3D reconstruction", 100, 0, 0);
   TH1F* hist5 = new TH1F("hist5","Resolution for 3D reconstruction", 100, 0, 0);
   TH1F* hist6 = new TH1F("hist6","Error 0 points #PEs", 100, 0, 0);
   TH2F* hist2D = new TH2F("hist2D", "True yz, erry = 0", 100, 0, 0, 100, 0, 2200);
   TH2F* hist2D1 = new TH2F("hist2D1", "True yz, errz = 0", 100, 0, 0, 100, 0, 2200);
   TH2F* hist2D2 = new TH2F("hist2D2", "True xy, erry = 0", 100, 0, 0, 100, 0, 2200);
   TH2F* hist2D3 = new TH2F("hist2D3", "True xy, errz = 0", 100, 0, 0, 100, 0, 2200);
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
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      std::map<int, Flash_info> flashes;
      for (int j = 0; j<OpHitChannel_Flash->size(); j++){
         if(flashes.count(OpFlashNumber->at(j)) == 0){
            flashes[OpFlashNumber->at(j)] = Flash_info();
         }
         Flash_info &flash = flashes[OpFlashNumber->at(j)];


         if (Channel_X->at(j) < -320){ //here I would add something like && abs(flash.chY - Channel_Y->at(j)) < distZ or something to account for the same distancing in Z
//             flash.w_sumX += Channel_X->at(j)*OpHitPE_Flash->at(j);
            flash.w_sumYC += Channel_Y->at(j)*OpHitPE_Flash->at(j);
            flash.w_sumZC += Channel_Z->at(j)*OpHitPE_Flash->at(j);
            flash.sumC += OpHitPE_Flash->at(j);
            flash.hit_vectorYC.push_back(Channel_Y->at(j));
            flash.hit_vectorZC.push_back(Channel_Z->at(j));
            flash.weight_hit_vectorZC.push_back(OpHitPE_Flash->at(j));
            flash.weight_hit_vectorYC.push_back(OpHitPE_Flash->at(j));
            for (int k=j+1; k<OpHitChannel_Flash->size(); k++){
            float norm = pow(pow(Channel_Y->at(j) - Channel_Y->at(k),2) + pow(Channel_Z->at(j) - Channel_Z->at(k),2), 0.5);
            if (norm < flash.min_dist_C && norm > 0){
               flash.min_dist_C = norm;
               }
         }
         }

         //Left wall reconstruction (x,z)
         if (Channel_Y->at(j) < -700){
            flash.w_sumXL += Channel_X->at(j)*OpHitPE_Flash->at(j);
//             flash.w_sumY += Channel_Y->at(j)*OpHitPE_Flash->at(j);
            flash.w_sumZL += Channel_Z->at(j)*OpHitPE_Flash->at(j);
            flash.sumL += OpHitPE_Flash->at(j);
            flash.hit_vectorXL.push_back(Channel_X->at(j));
            flash.hit_vectorZL.push_back(Channel_Z->at(j));
            flash.weight_hit_vectorXL.push_back(OpHitPE_Flash->at(j));
            flash.weight_hit_vectorZL.push_back(OpHitPE_Flash->at(j));
//             for (int k=j+1; k<OpHitChannel_Flash->size(); k++){
//             float norm = pow(pow(Channel_Y->at(j) - Channel_Y->at(k),2) + pow(Channel_Z->at(j) - Channel_Z->at(k),2), 0.5);
//             if (norm < flash.min_dist_C && norm > 0){
//                flash.min_dist_C = norm;
//                }
//          }
         }

         //Right wall reconstruction (x,z)
         if (Channel_Y->at(j) > 700){ //here I would add something like && abs(flash.chY - Channel_Y->at(j)) < distZ or something to account for the same distancing in Z
            flash.w_sumXR += Channel_X->at(j)*OpHitPE_Flash->at(j);
//             flash.w_sumY += Channel_Y->at(j)*OpHitPE_Flash->at(j);
            flash.w_sumZR += Channel_Z->at(j)*OpHitPE_Flash->at(j);
            flash.sumR += OpHitPE_Flash->at(j);
            flash.hit_vectorXR.push_back(Channel_X->at(j));
            flash.hit_vectorZR.push_back(Channel_Z->at(j));
            flash.weight_hit_vectorZR.push_back(OpHitPE_Flash->at(j));
            flash.weight_hit_vectorXR.push_back(OpHitPE_Flash->at(j));
//             for (int k=j+1; k<OpHitChannel_Flash->size(); k++){
//             float norm = pow(pow(Channel_Y->at(j) - Channel_Y->at(k),2) + pow(Channel_Z->at(j) - Channel_Z->at(k),2), 0.5);
//             if (norm < flash.min_dist_C && norm > 0){
//                flash.min_dist_C = norm;
//                }
//          }
         }
         //Back wall reco (x,y)
         if (Channel_Z->at(j) < 0){ //here I would add something like && abs(flash.chY - Channel_Y->at(j)) < distZ or something to account for the same distancing in Z
            flash.w_sumXback += Channel_X->at(j)*OpHitPE_Flash->at(j);
//             flash.w_sumY += Channel_Y->at(j)*OpHitPE_Flash->at(j);
            flash.w_sumYback += Channel_Y->at(j)*OpHitPE_Flash->at(j);
            flash.sumback += OpHitPE_Flash->at(j);
            flash.hit_vectorXback.push_back(Channel_X->at(j));
            flash.hit_vectorYback.push_back(Channel_Y->at(j));
            flash.weight_hit_vectorYback.push_back(OpHitPE_Flash->at(j));
            flash.weight_hit_vectorXback.push_back(OpHitPE_Flash->at(j));
//             for (int k=j+1; k<OpHitChannel_Flash->size(); k++){
//             float norm = pow(pow(Channel_Y->at(j) - Channel_Y->at(k),2) + pow(Channel_Z->at(j) - Channel_Z->at(k),2), 0.5);
//             if (norm < flash.min_dist_C && norm > 0){
//                flash.min_dist_C = norm;
//                }
//          }
         }
         //Front wall reco (x,y)
         if (Channel_Z->at(j) > 2100){ //here I would add something like && abs(flash.chY - Channel_Y->at(j)) < distZ or something to account for the same distancing in Z
            flash.w_sumXfront += Channel_X->at(j)*OpHitPE_Flash->at(j);
//             flash.w_sumY += Channel_Y->at(j)*OpHitPE_Flash->at(j);
            flash.w_sumYfront += Channel_Y->at(j)*OpHitPE_Flash->at(j);
            flash.sumfront += OpHitPE_Flash->at(j);
            flash.hit_vectorXfront.push_back(Channel_X->at(j));
            flash.hit_vectorYfront.push_back(Channel_Y->at(j));
            flash.weight_hit_vectorYfront.push_back(OpHitPE_Flash->at(j));
            flash.weight_hit_vectorXfront.push_back(OpHitPE_Flash->at(j));
//             for (int k=j+1; k<OpHitChannel_Flash->size(); k++){
//             float norm = pow(pow(Channel_Y->at(j) - Channel_Y->at(k),2) + pow(Channel_Z->at(j) - Channel_Z->at(k),2), 0.5);
//             if (norm < flash.min_dist_C && norm > 0){
//                flash.min_dist_C = norm;
//                }
//          }
         }
         flash.oldw_sumX += Channel_X->at(j)*OpHitPE_Flash->at(j);
         flash.oldw_sumY += Channel_Y->at(j)*OpHitPE_Flash->at(j);
         flash.oldw_sumZ += Channel_Z->at(j)*OpHitPE_Flash->at(j);
         flash.old_sumweight += OpHitPE_Flash->at(j);



         }
      for(const auto& [hitchannel, flash] : flashes){
//          if (hitchannel > 0){
//             break;
//          }
         oldw_wmeanX = flash.oldw_sumX/flash.old_sumweight;
         oldw_wmeanY = flash.oldw_sumY/flash.old_sumweight;
         oldw_wmeanZ = flash.oldw_sumZ/flash.old_sumweight;
         w_meanXR = 0;
         w_errXR = 0;
         if (flash.hit_vectorXR.size() > 0){
            w_meanXR = flash.w_sumXR/flash.sumR;
            w_errXR = rms_calc(flash.hit_vectorXR, flash.weight_hit_vectorXR, w_meanXR)/pow(flash.weight_hit_vectorXR.size(), 0.5);
         }
         w_meanXL = 0;
         w_errXL = 0;
         if (flash.hit_vectorXL.size() > 0){
            w_meanXL = flash.w_sumXL/flash.sumL;
            w_errXL = rms_calc(flash.hit_vectorXL, flash.weight_hit_vectorXL, w_meanXL)/pow(flash.weight_hit_vectorXL.size(), 0.5);
         }
         w_meanYC = 0;
         w_errYC = 0;
         if(flash.hit_vectorYC.size() > 0){
            w_meanYC = flash.w_sumYC/flash.sumC;
            w_errYC = rms_calc(flash.hit_vectorYC, flash.weight_hit_vectorYC, w_meanYC)/pow(flash.weight_hit_vectorYC.size(), 0.5);
         }
         w_meanZC = 0;
         w_errZC = 0;
         if (flash.hit_vectorZC.size() > 0){
            w_meanZC = flash.w_sumZC/flash.sumC;
            w_errZC = rms_calc(flash.hit_vectorZC, flash.weight_hit_vectorZC, w_meanZC)/pow(flash.weight_hit_vectorZC.size(), 0.5);
         }
         w_meanZR = 0;
         w_errZR = 0;
         if (flash.hit_vectorZR.size() > 0){
            w_meanZR = flash.w_sumZR/flash.sumR;
            w_errZR = rms_calc(flash.hit_vectorZR, flash.weight_hit_vectorZR, w_meanZR)/pow(flash.weight_hit_vectorZR.size(), 0.5);
         }
         w_meanZL = 0;
         w_errZL = 0;
         if (flash.hit_vectorZL.size() > 0){
            w_meanZL = flash.w_sumZL/flash.sumL;
            w_errZL = rms_calc(flash.hit_vectorZL, flash.weight_hit_vectorZL, w_meanZL)/pow(flash.weight_hit_vectorZL.size(), 0.5);
         }
         w_meanYback = 0;
         w_errYback = 0;
         if (flash.hit_vectorYback.size() > 0){
            w_meanYback = flash.w_sumYback/flash.sumback;
            w_errYback = rms_calc(flash.hit_vectorYback, flash.weight_hit_vectorYback, w_meanYback)/pow(flash.weight_hit_vectorYback.size(), 0.5);
         }
         w_meanYfront = 0;
         w_errYfront = 0;
         if (flash.hit_vectorYfront.size() > 0){
            w_meanYfront = flash.w_sumYfront/flash.sumfront;
            w_errYfront = rms_calc(flash.hit_vectorYfront, flash.weight_hit_vectorYfront, w_meanYfront)/pow(flash.weight_hit_vectorYfront.size(), 0.5);
         }
         w_meanXback = 0;
         w_errXback = 0;
         if (flash.hit_vectorXback.size() > 0){
            w_meanXback = flash.w_sumback/flash.sumback;
            w_errXback = rms_calc(flash.hit_vectorXback, flash.weight_hit_vectorXback, w_meanXback)/pow(flash.weight_hit_vectorXback.size(), 0.5);
         }
         w_meanXfront = 0;
         w_errXfront = 0;
         if (flash.hit_vectorXfront.size() > 0){
            w_meanXfront = flash.w_sumXfront/flash.sumfront;
            w_errXfront = rms_calc(flash.hit_vectorXfront, flash.weight_hit_vectorXfront, w_meanXfront)/pow(flash.weight_hit_vectorXfront.size(), 0.5);
         }
         std::vector<float> xmean_vec;
         std::vector<float> zmean_vec;
         std::vector<float> ymean_vec;
         std::vector<float> ymean_vec;
         std::vector<float> xweight_vec;
         std::vector<float> zweight_vec;
         if (w_errXR > 0){
            xmean_vec.push_back(w_meanXR);
            xweight_vec.push_back(w_errXR);
         }
         if (w_errXL > 0){
            xmean_vec.push_back(w_meanXL);
            xweight_vec.push_back(w_errXL);
         }
         if (w_errXback > 0){
            xmean_vec.push_back(w_meanXback);
            xweight_vec.push_back(w_errXback);
         }
         if (w_errXfront > 0){
            xmean_vec.push_back(w_meanXfront);
            xweight_vec.push_back(w_errXfront);
         }
         if (w_errZR > 0){
            zmean_vec.push_back(w_meanZR);
            zweight_vec.push_back(w_errZR);
         }
         if (w_errZL > 0){
            zmean_vec.push_back(w_meanZL);
            zweight_vec.push_back(w_errZL);
         }
         if (w_errZC > 0){
            zmean_vec.push_back(w_meanZC);
            zweight_vec.push_back(w_errZC);
         }
         if (w_errYC > 0){
            ymean_vec.push_back(w_meanYC);
            yweight_vec.push_back(w_errYC);
         }
         if (w_errYback > 0){
            ymean_vec.push_back(w_meanYback);
            yweight_vec.push_back(w_errYback);
         }
         if (w_errYfront > 0){
            ymean_vec.push_back(w_meanYfront);
            yweight_vec.push_back(w_errYfront);
         }
         float combined_x = combined_mean_calc(xmean_vec, xweight_vec);
         float combined_y = combined_mean_calc(ymean_vec, yweight_vec);
         float combined_z = combined_mean_calc(zmean_vec, zweight_vec);
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
            combined_yrrx = rms_calc(ymean_vec, yweight_vec, combined_y)/pow(ymean_vec.size(), 0.5);
         }
         else if(ymean_vec.size() == 1){
            combined_yrrx = yweight_vec.at(0);

         }
         if (zmean_vec.size() > 1){
            combined_errz = rms_calc(zmean_vec, zweight_vec, combined_z)/pow(zmean_vec.size(), 0.5);
         }
         else if(zmean_vec.size() == 1){
            combined_errz = zweight_vec.at(0);

         }

//          std::cout << "X: " << combined_x << " " << combined_errx << std::endl;
//          std::cout << "Z: " << combined_z << " " << combined_errz << std::endl;
//          std::cout << " Right wall X: " << " " << w_meanXR << " " << w_errXR << " " << flash.hit_vectorXR.size() << std::endl;
//          std::cout << " Left wall X: " << " " << w_meanXL << " " << w_errXL << " " << flash.hit_vectorXL.size() << std::endl;
//             std::cout << " Cathode Z: " << " " << w_meanZC << " " << w_errZC << " " << flash.hit_vectorZC.size() << std::endl;
//             std::cout << " Right wall Z: " << " " << w_meanZR << " " << w_errZR << " " << flash.hit_vectorZR.size() << std::endl;
//             std::cout << " Left wall Z: " << " " << w_meanZL << " " << w_errZL << " " << flash.hit_vectorZL.size() << std::endl;
//             std::cout << "Cathode Y: " << w_meanYC << " " << w_errYC << " " << flash.hit_vectorYC.size() << std::endl;
            if(combined_errx > 0){
               newres_x = abs(combined_x - TrueX);
               hist -> Fill(newres_x);
            }
            if(combined_erry > 0){
                  newres_y = abs(combined_y - TrueY);
              //    vector_newresy.push_back(newres_y);
                  hist1 -> Fill(newres_y);
                  RMS_vectorY.push_back(rms_calc(ymean_vec, yweight_vec, combined_y));
                  TrueY_vector.push_back(TrueY);
                  TrueX_vector.push_back(TrueX);

               }
               else{
                  hist2D -> Fill(TrueY, TrueZ);
                  hist2D2 -> Fill(TrueX, TrueY);
                  hist6 -> Fill(flash.sumC + flash.sumL +flash.sumR + flash.sumback + flash.sumfront);



            }
            if(combined_errz > 0){
                  newres_z = abs(combined_z - TrueZ);
                  hist2 -> Fill(newres_z);
                  RMS_vectorZ.push_back(rms_calc(zmean_vec, zweight_vec, combined_z));
                  TrueZ_vector.push_back(TrueZ);
               }
            else{
               hist2D1 -> Fill(TrueY, TrueZ);
               hist2D3 -> Fill(TrueX, TrueY);

            }

            oldres_x = abs(oldw_wmeanX - TrueX);
            oldres_y = abs(oldw_wmeanY - TrueY);
            oldres_z = abs(oldw_wmeanZ - TrueZ);
            hist3 -> Fill(oldres_x);
            hist4 -> Fill(oldres_y);
            hist5 -> Fill(oldres_z);
//             if (flash.min_dist_C != 1000000){
//             hist1 -> Fill(flash.min_dist_C);
//             }
//             res_C = pow(pow(w_meanYC - TrueY, 2) + pow(combined_z - TrueZ, 2), 0.5);
//             res_R_wall = pow(pow(w_meanXR - TrueX, 2) + pow(w_meanZR - TrueZ, 2), 0.5);
//             res_L_wall = pow(pow(w_meanXL - TrueX, 2) + pow(w_meanZL - TrueZ, 2), 0.5);
//             res_3D = pow(pow(w_meanYC - TrueY, 2) + pow(combined_z - TrueZ, 2), 0.5);
//          std::cout << w_meanXR << " " << w_meanXL << std::endl;
//          std::cout << w_meanZR << " " << w_meanZL << std::endl;
//          std::cout << w_meanYC << " " << w_meanZC << std::endl;
//             hist -> Fill(newres_x);
//             hist1 -> Fill(newres_y);
//             hist2 -> Fill(newres_z);

      }
//             for(int j=0; j<RecoXVector->size();j++){
//                oldres_x = abs(RecoXVector->at(j) - TrueX);
//                oldres_y = abs(YCenterVector->at(j) - TrueY);
//                oldres_z = abs(ZCenterVector->at(j) - TrueZ);
//                hist -> Fill(oldres_x);
//                hist1 -> Fill(oldres_y);
//                hist2 -> Fill(oldres_z);
//             }

   }
//    std::cout << w_meanY << std::endl;

   //Plot new 2D reco
   auto c = new TCanvas("histogram", "True yz, erry = 0");
//    TGraph *gr4 = new TGraph(vector_newresy.size(), &vector_newresy[0], &RMS_vectorY[0]);
//    gr4 -> SetTitle("RMS Y vs Resolution Y");
//    gr4 -> GetXaxis() -> SetTitle("Resolution Y (cm)");
//    gr4 -> GetYaxis() -> SetTitle("RMS Y");
//    gr4 -> Draw("AP");

      hist2 -> GetXaxis()->SetTitle("Resolution (cm)");
      hist2 -> SetTitle("OpSlicer resolution for 2D reconstruction");
      hist2 -> Draw();
      hist1 -> Draw("SAME");
      hist1 -> GetXaxis()->SetTitle("Resolution (cm)");
      hist1 -> SetLineColor(kRed);
      hist -> SetLineColor(kGreen);
      hist -> GetXaxis()->SetTitle("Resolution (cm)");
      hist -> Draw("SAME");
      // hist2D -> Draw("COLZ");

   auto legend1 = new TLegend();
   legend1 -> AddEntry(hist, Form("New mean X:  %5.2f (cm)", hist->GetMean())); //Form("New mean: " %5.2f, std::to_string(hist->GetMean()).c_str())
   legend1 -> AddEntry(hist1, Form("New mean Y:  %5.2f (cm)", hist1->GetMean()));
   legend1 -> AddEntry(hist2, Form("New mean Z:  %5.2f (cm)", hist2->GetMean()));
   legend1 -> Draw();

   //Plot old 3D reco

//    auto h = new TCanvas("histogram3", "True yz, errz = 0");
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
   TGraph *gr2 = new TGraph(RMS_vectorY.size(), &TrueX_vector[0], &RMS_vectorY[0]);
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
