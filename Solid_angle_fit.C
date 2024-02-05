#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/FitResult.h"
#include "Fit/DataOptions.h"
#include "Fit/FitConfig.h"
#include "TRandom3.h"
#include "TF1.h"
#include <limits>
#include "Math/WrappedMultiTF1.h"

//Can we do a complete 3D Fit?
//Cathode (Y,Z)
//Walls X (and some Z)

// void AddPD(int id, double xin, double yin);
//Gaussian clustering ... https://github.com/Ransaka/GMM-from-scratch

//class/structure with all arapucas
// 0s and then fill with hits..


class Fit3D
{
public:
  static  double PDx[], PDy[], PDz[];

  Fit3D(){
   
    hit=new double[N];
    std::cout<<"constructed ok"<<std::endl;
    Clear();
  };
  
  /*  PDhitlist(){
    //straight copy..
    }*/
  //destructor
  
  void Addhit(int i, double h){
    if(i<N){
      hit[i]=hit[i]+h;
//       std::cout << hit[i] << std::endl;
    }else{
      std::cout<<"Problem bin doesn't exist! "<<std::endl;
    }
  }
  
  void Clear(){ //empty hit array

    for(int i=0;i<N;i++)
      {
	hit[i]=0.0;
      }
  }
  
  float WRMS(int Np, double vals[], double weights[]){
   float value = 0;
   float sum_weights = 0;
   float mean = TMath::Mean(N, vals, weights);
   for(int k=0; k<Np;k++){
      value += pow((vals[k] - mean), 2)*weights[k];
      sum_weights += weights[k];
   }
   if (sum_weights > 0){
      value = pow(value/sum_weights, 0.5);
   }
   else{
      value = 0;
   }
   return value;

}

  //x[0] -- could be detector number [NOT AN INT!]
 static Double_t funcVol(Double_t *x, Double_t *par)
  {//par 0 1 2 = x y z
    //calculate distance from par[0..2] to x[0..2]
    //par[3] Nph

    int index=(int)x[0];
    //PDx[index] PDy[index] PDz[index]
    
    double res[1];
    double distance = TMath::Sqrt( pow(par[0]-PDx[index],2) +
				   pow(par[1]-PDy[index],2) + 
				   pow(par[2]-PDz[index],2));
    double theta=1.0;
    
    double area=60*60; 
    //orientation of Arapucas!
    //angle 60x60 square -- orientation ..cathode/wall
    //    std::cout<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;

    if(index>=56){
       theta=TMath::ASin(abs(par[0]-Fit3D::PDx[index])/distance);
         
    }else{//wall
	 theta = TMath::ASin(abs(par[1]-Fit3D::PDy[index])/distance);
  }
    

     return par[3]*area*TMath::Sin(theta)/(TMath::Pi()*4*pow(distance,2));

  }
  /*---------------------------------------*/
  //combined fit? x,y,z constrained by cathode and walls..
  // -- take x,y,z to be real space point?
  // and then solid angle to each Arapuca? [assuming no diffusion?]
  //
  //x is drift, y side, z- beam direction
   /*---------------------------------------*/
 
  
  //miniut?
  
 void Fit(){

    ROOT::Fit::BinData mydata(N, 1, ROOT::Fit::BinData::ErrorType::kNoError);
    double sum=0.0;
    double amp=0.0;

    double sumwall=0.0;
    
    for(int i=0;i<N;i++) //index i 0 - 167 inc
      {
	mydata.Add(i, hit[i]);
	sum+=hit[i];
	if(hit[i]>amp){amp=hit[i];}
	if(N<56){sumwall+=hit[i];}
      }
    
    double Xrms = 60; //WRMS(56, PDx, hit);
    double Yrms = 60;//WRMS(N, PDy, hit);
    double Zrms = 60;//WRMS(N, PDz, hit);

    double Xmean = TMath::Mean(56,PDx, hit);
    if(sumwall<1){//no hits on wall...
      Xmean = 150;//?
    }
    
    TF1 *ffunc = new TF1("ffunc",funcVol, 0, N,4);

    double initialParams[4] = {Xmean,TMath::Mean(N,PDy, hit),TMath::Mean(N,PDz,hit), sum};
    double step[4]={Xrms , Yrms, Zrms , sum};
    ROOT::Math::WrappedMultiTF1 fitFunction( *ffunc, 1 );
    fitter.SetFunction( fitFunction, false);

    std::cout<<"initial params : ";
    for(int i=0;i<4;i++)
      {
	std::cout<<initialParams[i]<<" ";
      }
    std::cout<<std::endl;
    /* initial params */

    fitter.Config().SetMinimizer("Minuit");
  
        
    fitter.Config().SetParamsSettings(4,initialParams, step);

    
    fitter.Config().ParSettings(0).SetLimits(rxmin, rxmax);
    fitter.Config().ParSettings(1).SetLimits(rymin, rymax);
    fitter.Config().ParSettings(2).SetLimits(rzmin, rzmax);
    fitter.Config().ParSettings(3).SetLimits(0,sum*10000); //need to work on this normalisation...
 


  
    
    fitter.Config().SetUpdateAfterFit();
    
    fitter.LikelihoodFit(mydata, true);
    

    
    r=fitter.Result();
    r.Print(std::cout);
    

    }
  
  ROOT::Fit::FitResult FitResult()
  {/* or merge with Fit function */
    return r;
  }
private:
  static int N;
  double  *hit;
  static double rxmin, rxmax, rymin, rymax, rzmin, rzmax;

  
  
  ROOT::Fit::Fitter fitter;

  ROOT::Fit::FitResult r;
};

// initialize the static data members  
 
int Fit3D::N = 168;   
double Fit3D::rxmin = -327;
double Fit3D::rxmax = 327;
double Fit3D::rymin = -736.9;
double Fit3D::rymax = 736.9;
double Fit3D::rzmin = 0;
double Fit3D::rzmax = 2000;

double Fit3D::PDx[168]={  285.07, 285.07, 285.07, 285.07, 285.07, 285.07, 285.07, 285.07, 285.07, 285.07, 285.07, 285.07, 285.07, 285.07, 210.07, 210.07, 210.07, 210.07, 210.07, 210.07, 210.07, 210.07, 210.07, 210.07, 210.07, 210.07, 210.07, 210.07, 135.07, 135.07, 135.07, 135.07, 135.07, 135.07, 135.07, 135.07, 135.07, 135.07, 135.07, 135.07, 135.07, 135.07, 60.07, 60.07, 60.07, 60.07, 60.07, 60.07, 60.07, 60.07, 60.07, 60.07, 60.07, 60.07, 60.07, 60.07, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26, -327.26};

double Fit3D::PDy[168]={ 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 732.062, -732.062, 538.102, 202.701, -132.701, -468.102, 373.402, 38.0006, -297.401, -632.802, 628.802, 293.401, -42.0006, -377.402, 464.102, 128.701, -206.701, -542.102, 538.102, 202.701, -132.701, -468.102, 373.402, 38.0006, -297.401, -632.802, 628.802, 293.401, -42.0006, -377.402, 464.102, 128.701, -206.701, -542.102, 538.102, 202.701, -132.701, -468.102, 373.402, 38.0006, -297.401, -632.802, 628.802, 293.401, -42.0006, -377.402, 464.102, 128.701, -206.701, -542.102, 538.102, 202.701, -132.701, -468.102, 373.402, 38.0006, -297.401, -632.802, 628.802, 293.401, -42.0006, -377.402, 464.102, 128.701, -206.701, -542.102, 538.102, 202.701, -132.701, -468.102, 373.402, 38.0006, -297.401, -632.802, 628.802, 293.401, -42.0006, -377.402, 464.102, 128.701, -206.701, -542.102, 538.102, 202.701, -132.701, -468.102, 373.402, 38.0006, -297.401, -632.802, 628.802, 293.401, -42.0006, -377.402, 464.102, 128.701, -206.701, -542.102, 538.102, 202.701, -132.701, -468.102, 373.402, 38.0006, -297.401, -632.802, 628.802, 293.401, -42.0006, -377.402, 464.102, 128.701, -206.701, -542.102}; 

double Fit3D::PDz[168]={1939.46, 1939.46, 1641.62, 1641.62, 1343.78, 1343.78, 1045.94, 1045.94, 748.1, 748.1, 450.26, 450.26, 152.42, 152.42, 1939.46, 1939.46, 1641.62, 1641.62, 1343.78, 1343.78, 1045.94, 1045.94, 748.1, 748.1, 450.26, 450.26, 152.42, 152.42, 1939.46, 1939.46, 1641.62, 1641.62, 1343.78, 1343.78, 1045.94, 1045.94, 748.1, 748.1, 450.26, 450.26, 152.42, 152.42, 1939.46, 1939.46, 1641.62, 1641.62, 1343.78, 1343.78, 1045.94, 1045.94, 748.1, 748.1, 450.26, 450.26, 152.42, 152.42, 2044.46, 2044.46, 2044.46, 2044.46, 1973.46, 1973.46, 1973.46, 1973.46, 1898.46, 1898.46, 1898.46, 1898.46, 1827.46, 1827.46, 1827.46, 1827.46, 1746.62, 1746.62, 1746.62, 1746.62, 1675.62, 1675.62, 1675.62, 1675.62, 1600.62, 1600.62, 1600.62, 1600.62, 1529.62, 1529.62, 1529.62, 1529.62, 1448.78, 1448.78, 1448.78, 1448.78, 1377.78, 1377.78, 1377.78, 1377.78, 1302.78, 1302.78, 1302.78, 1302.78, 1231.78, 1231.78, 1231.78, 1231.78, 1150.94, 1150.94, 1150.94, 1150.94, 1079.94, 1079.94, 1079.94, 1079.94, 1004.94, 1004.94, 1004.94, 1004.94, 933.94, 933.94, 933.94, 933.94, 853.1, 853.1, 853.1, 853.1, 782.1, 782.1, 782.1, 782.1, 707.1, 707.1, 707.1, 707.1, 636.1, 636.1, 636.1, 636.1, 555.26, 555.26, 555.26, 555.26, 484.26, 484.26, 484.26, 484.26, 409.26, 409.26, 409.26, 409.26, 338.26, 338.26, 338.26, 338.26, 257.42, 257.42, 257.42, 257.42, 186.42, 186.42, 186.42, 186.42, 111.42, 111.42, 111.42, 111.42, 40.42, 40.42, 40.42, 40.42};


// could make test now...

void testSA(int index, double xoffset, double yoffset, double zoffset)
{
  double x[1]={(double) index};
  double theta=0.0;
  double phi=0.0;
  double p[4]={Fit3D::PDx[index]+xoffset, Fit3D::PDy[index]+yoffset, Fit3D::PDz[index]+zoffset};
  double distance = TMath::Sqrt( pow(p[0]-Fit3D::PDx[index],2) +
				   pow(p[1]-Fit3D::PDy[index],2) + 
				   pow(p[2]-Fit3D::PDz[index],2));
  double area=60*60; 
 
       if(index>=56){
       theta=TMath::ASin(abs(p[0]-Fit3D::PDx[index])/distance);
         
  }else{//wall
	 theta = TMath::ASin(abs(p[1]-Fit3D::PDy[index])/distance);
    std::cout<<"wall"<<std::endl;
  }
       std::cout<<index<<" angle "<<TMath::Sin(theta)<<" "<<std::endl;
    //hmm not sure -- check this... 
  std::cout<<"SAF: "<<area*TMath::Sin(theta)/(TMath::Pi()*4*pow(distance,2))<<std::endl;
  }

/*
double xmin = 0;
  double ymin=0;
  double ymax = 5;
  double xmax = 5;
 
 TH2D htemp("htemp","",5,xmin,xmax,5, ymin, ymax);

void test(double meanx, double meany, double sigma)
{


  int Npds=50; //5x5 array
  double xpos[Npds], ypos[Npds];
 
  
  int binx, biny,binz;
  
  for(int i=0;i<Npds;i++)
    {

      htemp.GetBinXYZ(i,binx,biny,binz);
   
      xpos[i]=htemp.GetXaxis()->GetBinCenter(binx);
      ypos[i]=htemp.GetYaxis()->GetBinCenter(biny);
      
    }

 
 
  Fit3D myversion(Npds,xpos,ypos,xmin,xmax,ymin,ymax);
  TRandom3 *grand= new TRandom3();

  double dumx,dumy;
  int ig;
  for(int i=0;i<100;i++)
    {
      dumx = grand->Gaus(meanx,sigma);
      dumy = grand->Gaus(meany,sigma);
      ig =  htemp.FindBin(dumx,dumy);
      //   std::cout<<dumx<<" "<<dumy<<" "<<ig<<std::endl;
      if(ig<Npds){
	myversion.Addhit(ig, 1.0);
	htemp.Fill(dumx,dumy,1.0);
      }else{
	std::cout<<"problem -> "<<dumx<<" "<<dumy<<" "<<ig<<std::endl;
      }
    }

  std::cout<<htemp.GetEntries()<<std::endl;
  TCanvas *c1 = new TCanvas("c1","",600,600);
  htemp.Draw("colZ");
  htemp.SetDirectory(0);
  
  myversion.Fit();

}
*/
