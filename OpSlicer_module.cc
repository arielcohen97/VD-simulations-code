// Dan Pershey
//
// OpHit clusterer, using y,z,t information to isolate OpHit's from a common
// origin.  Algorithm uses DBScan - but slightly modified, by calculating a
// cluster centroid, to allow for better clustering of mutiple-density
// clusters
//
// Using OpFlashFinder_module.cc as template
//
    
    
    #ifndef OpSlicer_H
    #define OpSlicer_H 1
    
    // LArSoft includes
    #include "larcore/Geometry/Geometry.h"
    #include "lardataobj/RecoBase/OpFlash.h"
   #include "lardataobj/RecoBase/OpHit.h"
   #include "lardata/Utilities/AssociationUtil.h"
   #include "lardata/DetectorInfoServices/DetectorClocksService.h"
   #include "lardataalg/DetectorInfo/DetectorClocks.h"
    
    
    // Framework includes
    #include "art/Framework/Core/EDProducer.h"
   #include "art/Framework/Core/ModuleMacros.h"
    #include "art/Framework/Principal/Event.h"
    #include "fhiclcpp/ParameterSet.h"
    #include "art/Framework/Principal/Handle.h"
    #include "canvas/Persistency/Common/Ptr.h"
    #include "canvas/Persistency/Common/PtrVector.h"
    #include "art/Framework/Services/Registry/ServiceHandle.h"
    
    // ROOT includes
    
    // C++ Includes
    #include <vector>
    #include <string>
    #include <memory>
    #include <limits>
    
    namespace opdet {

    enum Plane {
        kPlaneXZR,
        kPlaneXZL,
        kPlaneYZ
    };

    std::map<Plane, std::vector<int>> plane_map = {
        {kPlaneXZR, {0, 2}},
        {kPlaneXZL, {0, 2}},
        {kPlaneYZ, {1, 2}}
    };
 
    bool sortOpHitByTime(const art::Ptr<recob::OpHit> &left,
                        const art::Ptr<recob::OpHit> &right)
    {
    return left->PeakTime() < right->PeakTime();
    }

    class OpSlicer : public art::EDProducer{
    public:
    
    // Standard constructor and destructor for an ART module.
    explicit OpSlicer(const fhicl::ParameterSet&);
    virtual ~OpSlicer();

    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& pset);

    // The producer routine, called once per event.
    void produce(art::Event&);

    private:

    double YZDist(art::Ptr<recob::OpHit>, art::Ptr<recob::OpHit>, const Plane &plane ) const;
    double Dist(art::Ptr<recob::OpHit>, art::Ptr<recob::OpHit>, const Plane &plane) const;
    int YZCentroid(std::vector< art::Ptr<recob::OpHit> >,
                    std::vector<int>, const Plane& plane) const;

    void GetHitYZ(std::vector< art::Ptr<recob::OpHit> >,
                    std::vector<int>,
                    std::vector<double>&,
                    std::vector<double>&,
                    const Plane &plane) const;
    void ClusterHits(std::vector<art::Ptr<recob::OpHit>>,
                        std::vector<recob::OpFlash>&,
                        std::vector< std::vector<int> >&,
                        detinfo::DetectorClocksData const&,
                        const Plane& plane)const;


    // The parameters we'll read from the .fcl file.
    std::string fOpHitModuleLabel; // Input tag for OpHit collection

    double fTScale;
    double fRScale;
    double fR0;
    double fR0_RW;
    double fR0_LW;
    double fBreakTime;
    int fMinN;

    double fTrigCoinc;
        
    std::map< Plane, float>  dbstep;


    art::ServiceHandle<geo::Geometry> geo;
    };
    
    }
    
    namespace opdet {
      DEFINE_ART_MODULE(OpSlicer)
    }
   
   #endif
   
   namespace opdet {
   
     //--------------------------------------------------------------------------
     // Constructor
     OpSlicer::OpSlicer(const fhicl::ParameterSet & pset) : EDProducer{pset}
     {
   
       reconfigure(pset);
   
       produces< std::vector< recob::OpFlash > >();
       produces< art::Assns< recob::OpFlash, recob::OpHit > >();
   
     }
   
     //--------------------------------------------------------------------------
     void OpSlicer::reconfigure(fhicl::ParameterSet const& pset)
     {
   
       // Indicate that the Input Module comes from .fcl
       fOpHitModuleLabel = pset.get< std::string >("OpHitModuleLabel");
   
       fTScale    = pset.get<double>("TScale");
       fRScale    = pset.get<double>("RScale");
       fR0        = pset.get<double>("R0");
       fR0_RW     = pset.get<double>("R0_RW");
       fR0_LW     = pset.get<double>("R0_LW");
       fBreakTime = pset.get<double>("BreakTime");
       fMinN      = pset.get<int>   ("MinN");
   
       fTrigCoinc      = pset.get< double >("TrigCoinc");
   
     }
   
     //--------------------------------------------------------------------------
     // Destructor
     OpSlicer::~OpSlicer()
     {
     }
     //--------------------------------------------------------------------------
     void OpSlicer::beginJob()
     {
     }
   
     //--------------------------------------------------------------------------
     void OpSlicer::endJob()
     {
     }
   
     //--------------------------------------------------------------------------
     void OpSlicer::produce(art::Event& evt)
     {
   
       // These are the storage pointers we will put in the event
       std::unique_ptr< std::vector< recob::OpFlash > >
                         flashPtr(new std::vector< recob::OpFlash >);
       std::unique_ptr< art::Assns< recob::OpFlash, recob::OpHit > >
                         assnPtr(new art::Assns< recob::OpFlash, recob::OpHit >);
   
       // This will keep track of what flashes will assoc to what ophits
       // at the end of processing
       
      std::vector< std::vector< int > > assoc_LW;
      std::vector< std::vector< int > > assoc_RW;
      std::vector< std::vector< int > > assoc_C;
      std::map<Plane, std::vector<std::vector< int >> > assocList;
      assocList[kPlaneXZL] = assoc_LW;
      assocList[kPlaneXZR] = assoc_RW;
      assocList[kPlaneYZ] = assoc_C;
       //std::vector< std::vector< int > > assocList;
   
       auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
   
       // Get OpHits from the event
       auto opHitHandle = evt.getHandle< std::vector< recob::OpHit > >(fOpHitModuleLabel);
   
       std::vector< art::Ptr<recob::OpHit> > ohits;
       for (int i = 0; i < int(opHitHandle->size()); i++){
         art::Ptr<recob::OpHit> opHitPtr(opHitHandle,i);
         ohits.push_back(opHitPtr);
       }

        std::vector< art::Ptr<recob::OpHit> > ophits_C;
        std::vector< art::Ptr<recob::OpHit> > ophits_RW;
        std::vector< art::Ptr<recob::OpHit> > ophits_LW;

        for (int i = 0; i < int(ohits.size()); i++){
            int channel = ohits[i]->OpChannel();
            auto pos_xyz = geo->OpDetGeoFromOpChannel(channel).GetCenter();
            double cx = pos_xyz.X();
            double cy = pos_xyz.Y();
            if (cx < -300){
                ophits_C.push_back(ohits[i]);
            }
            if (cy < -700){
                ophits_LW.push_back(ohits[i]);
            }
            if (cy > 700){
                ophits_RW.push_back(ohits[i]);
            }
        }

        std::map<Plane, std::vector< art::Ptr<recob::OpHit> >> ophits_by_plane;
        ophits_by_plane[kPlaneXZL] = ophits_LW;
        ophits_by_plane[kPlaneXZR] = ophits_RW;
        ophits_by_plane[kPlaneYZ] = ophits_C;

        std::cout << ophits_LW.size() << " " << ophits_RW.size() << " " << ophits_C.size() << std::endl; 

        dbstep = {
          {kPlaneXZR, fR0_RW},
          {kPlaneXZL, fR0_LW},
          {kPlaneYZ, fR0},
};
   
        for(auto const& [plane, hits] : ophits_by_plane){
            ClusterHits(hits, *flashPtr, assocList[plane], clockData, plane);
        }
   
       // Run the clustering
    //    ClusterHits(ohits, *flashPtr, assocList, clockData);
   
       int assn_idx = 0; 
       
      for(auto const& [plane, hits] : ophits_by_plane){
            
    
      // Make the associations which we noted we need
       for (size_t i = 0; i != assocList[plane].size(); ++i)
       {
         art::PtrVector< recob::OpHit > opHitPtrVector;
         for (int const& hitIndex : assocList[plane].at(i))
         {
           art::Ptr< recob::OpHit > opHitPtr(hits[hitIndex]);
           opHitPtrVector.push_back(opHitPtr);
         }
   
         util::CreateAssn(*this, evt, *flashPtr, opHitPtrVector,
                          *(assnPtr.get()), assn_idx++);
       }
}
   
       // Store results into the event
       evt.put(std::move(flashPtr));
       evt.put(std::move(assnPtr));
     }
   
     double OpSlicer::YZDist(art::Ptr<recob::OpHit> a,
                             art::Ptr<recob::OpHit> b,
                             const Plane &plane) const
     {
       // First need to ask the geometry for the y-z location of both OpHits
       std::vector<int> axes = plane_map[plane];
       int channela = a->OpChannel();
       auto xyza_point = geo->OpDetGeoFromOpChannel(channela).GetCenter();
       double xyza[3];
       xyza_point.GetCoordinates(xyza);
       double ay = xyza[axes[0]];
       double az = xyza[axes[1]];
   
       int channelb = b->OpChannel();
       auto xyzb_point = geo->OpDetGeoFromOpChannel(channelb).GetCenter();
       double xyzb[3];
       xyzb_point.GetCoordinates(xyzb);
       double by = xyzb[axes[0]];
       double bz = xyzb[axes[1]];
   
       double r2 = pow((ay-by),2) + pow((az-bz),2);
   
       double r = sqrt(r2);
       return r;
     }
     double OpSlicer::Dist(art::Ptr<recob::OpHit> a,
                           art::Ptr<recob::OpHit> b,
                           const Plane &plane) const
     {
       double r = YZDist(a,b, plane);
       return sqrt(pow((a->PeakTime()-b->PeakTime())/fTScale,2)+pow(r/fRScale,2));
     }
   
     int OpSlicer::YZCentroid(std::vector< art::Ptr<recob::OpHit> > ohits,
                              std::vector<int> curN,
                              const Plane & plane) const
     {
       double maxDens = 0;
       int maxIdx = -1;
       for (int i = 0; i < int(curN.size()); i++){
         double dens = 0;
         for (int j = 0; j < int(curN.size()); j++){
           double r = YZDist(ohits[i],ohits[j], plane);
           r = sqrt(r)/100;
           dens += ohits[curN[j]]->PE() * exp(-r*r);
         }
         if (dens > maxDens){
           maxDens = dens;  maxIdx = curN[i];
         }
       }
       return maxIdx;
     }
   
   
   
     void OpSlicer::GetHitYZ(std::vector< art::Ptr<recob::OpHit> > ohits,
                             std::vector<int> curN,
                             std::vector<double> &ys,
                             std::vector<double> &zs,
                             const Plane &plane) const
     {
        std::vector<int> axes = plane_map[plane];
       for (int cur : curN){
         art::Ptr<recob::OpHit> oh = ohits[cur];
         int channel = oh->OpChannel();
         auto xyz_point = geo->OpDetGeoFromOpChannel(channel).GetCenter();
         double xyz[3];
         xyz_point.GetCoordinates(xyz);
         ys.push_back(xyz[axes[0]]);
         zs.push_back(xyz[axes[1]]);
       }
     }
     
     void OpSlicer::ClusterHits(std::vector< art::Ptr<recob::OpHit> > ohits,
                                std::vector< recob::OpFlash>& oflashes,
                                std::vector< std::vector<int> >& assoc,
                                detinfo::DetectorClocksData const &ts,
                                const Plane &plane)const
     {
       float min_dist = dbstep.at(plane);

       std::sort(ohits.begin(),ohits.end(),sortOpHitByTime);
   
       std::vector<bool> isClust;
       for (int i = 0; i < int(ohits.size()); i++) isClust.push_back(false);
   
       std::vector<int> neigh;
   
       for (int i = 0; i < int(ohits.size()); i++){
   
         if (isClust[i]) continue; // Don't base clusts off of clustered hits!
         neigh.erase(neigh.begin(),neigh.end()); // Start from scratch every time
   
         // check nearby hits in time for coincidence
         for (int j = i-1; j > 0; j--){
           if ( YZDist(ohits[i],ohits[j], plane) < min_dist && !isClust[j] ){
             neigh.push_back(j);
           }
           if (abs(ohits[i]->PeakTimeAbs()-ohits[j]->PeakTimeAbs())>2) break;
         }
         for (int j = i+1; j < int(ohits.size()); j++){
           if ( YZDist(ohits[i],ohits[j], plane) < min_dist && !isClust[j] ){
             neigh.push_back(j);
           }
           if (abs(ohits[i]->PeakTimeAbs()-ohits[j]->PeakTimeAbs())>2) break;
         }
         if (int(neigh.size())<fMinN) continue;
         neigh.erase(neigh.begin(),neigh.end());
   
   
         std::vector<int> cands;
         for (int j = i; j < int(ohits.size()); j++){
           if (isClust[j]) continue;
           if (ohits[j]->PeakTimeAbs()-ohits[i]->PeakTimeAbs()>fBreakTime) break;
           cands.push_back(j);
         }
         int centroidIdx = YZCentroid(ohits,cands, plane);
         if (centroidIdx < 0) // No centroid found
           continue;
         
         art::Ptr<recob::OpHit> centroid = ohits[centroidIdx];
         neigh.push_back(centroidIdx);
   
         std::vector<int> curN;
         curN.push_back(centroidIdx);
   
         // check nearby hits in time for coincidence
         for (int j = centroidIdx-1; j > 0; j--){
           if ( YZDist(ohits[centroidIdx],ohits[j], plane) < min_dist && !isClust[j] ){
             neigh.push_back(j);
             curN.push_back(j);
           }
           if (abs(ohits[centroidIdx]->PeakTimeAbs()-ohits[j]->PeakTimeAbs())>2) break;
         }
         for (int j = centroidIdx+1; j < int(ohits.size()); j++){
           if ( YZDist(ohits[centroidIdx],ohits[j], plane) < min_dist && !isClust[j] ){
             neigh.push_back(j);
             curN.push_back(j);
           }
           if (abs(ohits[centroidIdx]->PeakTimeAbs()-ohits[j]->PeakTimeAbs())>2) break;
         }
         double totPE = 0;
         for (int idx : curN) totPE += ohits[idx]->PE();
         if (int(curN.size())<fMinN) continue;
   
         // Loop through neighboring hits, chck if it's a core hit
         while (neigh.size() > 0){
           std::vector<int> curNeigh;  curNeigh.push_back(neigh[0]);
           for (int j = neigh[0]-1; j > 0; j--){
             if ( YZDist(ohits[neigh[0]],ohits[j], plane) < min_dist && !isClust[j] ){
               curNeigh.push_back(j);
             }
             if (abs(ohits[centroidIdx]->PeakTimeAbs()-ohits[j]->PeakTimeAbs())>2)
               break;
           }
           for (int j = neigh[0]+1; j < int(ohits.size()); j++){
             if ( YZDist(ohits[neigh[0]],ohits[j], plane) < min_dist && !isClust[j] ){
               curNeigh.push_back(j);
             }
             if (abs(ohits[centroidIdx]->PeakTimeAbs()-ohits[j]->PeakTimeAbs())>2)
               break;
           }
           // If this is a core point, add in all reachable hits to neighborhood
           if (int(curNeigh.size())>=fMinN){
             for (int cur : curNeigh){
               if (std::find(curN.begin(),curN.end(),cur)==curN.end()){
                 curN.push_back(cur);
                 if (YZDist(ohits[centroidIdx],ohits[cur], plane) < min_dist)
                   neigh.push_back(cur);
               }
             }
           }
           neigh.erase(neigh.begin());
         }
   
         if (int(curN.size())<fMinN) continue;

         std::vector<int> axes = plane_map[plane];
   
         // Time to make the OpFlash;
         // Y-Z coordinates come from the centroid
         int channelcentroid = centroid->OpChannel();
         auto xyzcentroid_point = geo->OpDetGeoFromOpChannel(channelcentroid).GetCenter();
         double xyzcentroid[3];
         xyzcentroid_point.GetCoordinates(xyzcentroid);
         double yCenter = xyzcentroid[axes[0]];
         double zCenter = xyzcentroid[axes[1]];
         double tCenter = centroid->PeakTimeAbs();
   
   
         // Now that we have centroid coordinates, include ana delayed light
         for (int j = i; j < int(ohits.size()); j++){
           if (std::find(curN.begin(),curN.end(),j)!=curN.end()) continue;
           double r = YZDist(ohits[j],centroid, plane);
           if ( r < fRScale*min_dist){
             curN.push_back(j);
           }
           if (abs(ohits[j]->PeakTimeAbs()-tCenter)>fBreakTime) break;
         }
   
         double finE = 0;
         for (int idx : curN) finE += ohits[idx]->PE();
   
         // Grab the y-z information from the geometry
         std::vector<double> ys;
         std::vector<double> zs;
         GetHitYZ(ohits,curN,ys,zs, plane);
   
         double minT = std::numeric_limits<double>::max(); double maxT = -std::numeric_limits<double>::max();
         double minY = 1e6; double maxY = -1e6;
         double minZ = 1e6; double maxZ = -1e6;
   
         std::vector<double> PEs (geo->MaxOpChannel() + 1,0.0);
         std::vector<double> PE2s (geo->MaxOpChannel() + 1,0.0);
         double fastToTotal = 0;
         for (int hIdx = 0; hIdx < int(ys.size()); hIdx++){
           int cIdx = curN[hIdx];
   
           minT = std::min(minT,ohits[cIdx]->PeakTimeAbs());
           maxT = std::max(maxT,ohits[cIdx]->PeakTimeAbs());
           minY = std::min(minY,ys[hIdx]);
           maxY = std::min(maxY,ys[hIdx]);
           minZ = std::min(minZ,zs[hIdx]);
           maxZ = std::min(maxZ,zs[hIdx]);
           PEs[ohits[cIdx]->OpChannel()] += ohits[cIdx]->PE();
           PE2s[ohits[cIdx]->OpChannel()] += ohits[cIdx]->PE();
           fastToTotal += ohits[hIdx]->FastToTotal();
         }
         double yWidth = maxY-minY;
         double zWidth = maxZ-minZ;
   
         double tot1 = 0;
         double tot2 = 0;
         for (double PE : PEs) tot1 += PE;
         for (double PE : PE2s) tot2 += PE;
   
         // From OpFlashAlg
         int Frame = ts.OpticalClock().Frame(tCenter - 18.1);
         if (Frame == 0) Frame = 1;
   
         int BeamFrame = ts.OpticalClock().Frame(ts.TriggerTime());
         bool InBeamFrame = false;
         if (!(ts.TriggerTime() < 0)) InBeamFrame = (Frame == BeamFrame);
   
         double tWidth = (maxT-minT)/2;
   
         int OnBeamTime = 0;
         if (InBeamFrame && (std::abs(tCenter) < fTrigCoinc)) OnBeamTime = 1;

         double xcenter = -999;
         double xwidth = -999;
         double ycenter = -999;
         double ywidth = -999;
         double zcenter = -999;
         double zwidth = -999;

         switch (plane){
            case kPlaneXZR:
                xcenter = yCenter;
                xwidth = yWidth;
                ycenter = zCenter;
                ywidth = zWidth;
                break;
            case kPlaneYZ:
                ycenter = yCenter;
                ywidth = yWidth;
                zcenter = zCenter;
                zwidth = zWidth;
                break;
            case kPlaneXZL:
                xcenter = yCenter;
                xwidth = yWidth;
                zcenter = zCenter;
                zwidth = zWidth;
                break;
         }
   
         oflashes.emplace_back(tCenter,tWidth,tCenter,Frame,
                               PEs,InBeamFrame,OnBeamTime,fastToTotal,
                               xcenter, xwidth, ycenter,ywidth,zcenter,zwidth);
         assoc.emplace_back(curN);
   
   
         // And finally, indicate that current hits have been clustered
         for (int cur : curN) isClust[cur] = true;
   
       }
     
    }//I need to move this bracket to encompass the ClusterHits function
   
   } // namespace opdet
