#ifndef MICHEL_H
#define MICHEL_H


#include "CVUniverse.h"


#include "Cluster.h"
// A michel object is defined by a cluster of hits AND by the vertex to which
// it is matched.
//
// The michel tool allows one cluster to be matched to multiple vertices, so
// there can be multiple michel objects associated with a single cluster of
// hits.
//
// Specifically, michel object includes the info:
// * which vertex and track it's associated with
// * the distance from the endpoint
// * whether the fit distance is acceptable or not

class Michel;

// c++11 way to do a typedef
//using MichelMap = std::map<int, Michel>;

//// old school typedef that still works
//typedef std::map<int, Michel> MichelMap;

class Michel
{  
 public:
  //constructors
  Michel(const CVUniverse& univ, int ci);

  // Functions
  Michel(){};  // fill in most basic stuff in "generic info"
  void Clear(Michel* &current_michel);
  void DoMoreInitializationAndProcessing(); // fill in more complicated stuff in "generic info"
  void DoMatching();                          // fill in "best matching stuff"
  void DoesMichelMatchVtx(const CVUniverse& univ, Michel* &current_michel); //GEts info for Vtx Match
  void DoesMichelMatchClus(const CVUniverse& univ, Michel* &current_michel); // Gets info for ClusterMatch
  void GetBestMatch(Michel* &current_michel); // get enum for best match out of all four saved matches
  //std::vector<Michel*> CreateMichels(CVUniverse& univ);
  // generic info
  std::vector<double> up_location; // upstream location 0 X 1 U 2 V 3 Z
  std::vector<double> down_location; //downstream location
  double energy = -999.;
  double time = -999.;
  int is_fitted = -1;
  std::vector<double> up_to_vertex_dist2D;
  std::vector<double> down_to_vertex_dist2D;
  double up_to_vertex_dist3D = 9999.;
  double down_to_vertex_dist3D = 9999.;
  std::vector<Cluster> cluster_to_up_match;
  std::vector<Cluster> cluster_to_down_match;
  double up_to_cluster_dist3D = 9999.;
  double down_to_cluster_dist3D = 9999.;
  
  int nclusters = 0;
  
  // best matching stuff
  //enum *best_cluster_match; // just tells you which of the four matches are the best match
  //enum BestClusterMatch {kUpVtxMatch, kDownVtxMatch, kUpClusMatch, kDownClusMatch, kNClusterMatches};
  
  //the following is in place until i figure out how to use the enum.
  int BestMatch = 0; // 0 = null match, 1= kUpVtxMatch, 2 = kDownVtxMatch, 3 = kUpClusMatch, 4=  kDownClusMatch, 
  
  double Best3Ddist = 9999.;
                
};

void Michel::Clear(Michel* &current_michel){
  up_location.clear(); // upstream location 0 X 1 U 2 V 3 Z
  down_location.clear(); //downstream location
  energy = -999.;
  time = -999.;
  is_fitted = -1;
  up_to_vertex_dist2D.clear();
  down_to_vertex_dist2D.clear();
  up_to_vertex_dist3D = 9999.;
  down_to_vertex_dist3D = 9999.;
  cluster_to_up_match.clear();
  cluster_to_down_match.clear();
  up_to_cluster_dist3D = 9999.;
  down_to_cluster_dist3D = 9999.;
  
  // best matching stuff
  //enum *best_cluster_match; // just tells you which of the four matches are the best match
  //enum BestClusterMatch {kUpVtxMatch, kDownVtxMatch, kUpClusMatch, kDownClusMatch, kNClusterMatches};
  
  //the following is in place until i figure out how to use the enum.
  BestMatch = 0; // 0 = null match, 1= kUpVtxMatch, 2 = kDownVtxMatch, 3 = kUpClusMatch, 4=  kDownClusMatch, 
  
  Best3Ddist = 9999.;
}

void Michel::GetBestMatch(Michel* &current_michel){
     int upvtxmatch = 0;
     int downvtxmatch = 0;
     int upclusmatch = 0;
     int downclusmatch = 0;
     
     //if (current_michel->up_to_vertex_dist3D != NULL && )
     
     //std::cout << "GET BEST MATCH FOR THIS MICHEL " << std::endl;
     if (current_michel->up_to_vertex_dist3D < current_michel->down_to_vertex_dist3D){

        upvtxmatch = 1;
     }
     else if (current_michel->up_to_vertex_dist3D >= current_michel->down_to_vertex_dist3D){

        downvtxmatch = 1;
     }
     if (current_michel->up_to_cluster_dist3D < current_michel->down_to_cluster_dist3D){

        upclusmatch = 1;
     }
     else if (current_michel->up_to_cluster_dist3D >= current_michel->down_to_cluster_dist3D){

        downclusmatch = 1;
     }
     
     if (upvtxmatch == 1 && (upclusmatch == 1 || downclusmatch == 1 )){
        if (current_michel->up_to_vertex_dist3D < current_michel->up_to_cluster_dist3D) {
           //std::cout << "UPVTX IS BEST " << std::endl;

           current_michel->BestMatch = 1;
           current_michel->Best3Ddist = current_michel->up_to_vertex_dist3D;
        }
        else if (current_michel->up_to_vertex_dist3D >= current_michel->up_to_cluster_dist3D)
        {
           //std::cout << "UPCLUS IS BEST " << std::endl;

           current_michel->BestMatch = 3;
           current_michel->Best3Ddist = current_michel->up_to_cluster_dist3D;
        }
        else if (current_michel->up_to_vertex_dist3D < current_michel->down_to_cluster_dist3D) {
           //std::cout << "UPVTX IS BEST " << std::endl;

           current_michel->BestMatch = 1;
           current_michel->Best3Ddist = current_michel->up_to_vertex_dist3D;

        }
        else if (current_michel->up_to_vertex_dist3D >= current_michel->down_to_cluster_dist3D)
        {
           //std::cout << "DOWNCLUSISBEST" << std::endl;

           current_michel->BestMatch = 4;
           current_michel->Best3Ddist = current_michel->down_to_cluster_dist3D;

        }
     }
     else if (downvtxmatch == 1 && (upclusmatch == 1 || downclusmatch == 1 )){
       if (current_michel->down_to_vertex_dist3D < current_michel->up_to_cluster_dist3D) {
           //std::cout << "DOWNVTX IS BEST " << std::endl;

           current_michel->BestMatch = 2;
           
           //std::cout << "Setting Best Dist" << std::endl;

           current_michel->Best3Ddist = current_michel->down_to_vertex_dist3D;

       }
       else if (current_michel->down_to_vertex_dist3D >= current_michel->up_to_cluster_dist3D)
       {
           //std::cout << "UPCLUS IS BEST " << std::endl;

           current_michel->BestMatch = 3;
           current_michel->Best3Ddist = current_michel->up_to_cluster_dist3D;

       }
       else if (current_michel->down_to_vertex_dist3D < current_michel->down_to_cluster_dist3D)
       {        
           //std::cout << "DOWNVTX IS BEST " << std::endl;

           current_michel->BestMatch = 2;
           
           //std::cout << "Setting Best Dist" << std::endl;

           current_michel->Best3Ddist = current_michel->down_to_vertex_dist3D;

       }
       else if (current_michel->down_to_vertex_dist3D >= current_michel->down_to_cluster_dist3D)
       {
       
           //std::cout << "DOWNCLUSISBEST" << std::endl;

           current_michel->BestMatch = 4;
           current_michel->Best3Ddist = current_michel->down_to_cluster_dist3D;

       }
     
     }
     else{
      current_michel->BestMatch = 0;
     }
     
}

Michel::Michel(const CVUniverse& univ, int ci) 
{ 
  energy = univ.GetVecElem("FittedMichel_michel_energy", ci);
  time = univ.GetVecElem("FittedMichel_michel_time", ci)/pow(10,3);
  is_fitted = univ.GetVecElem("FittedMichel_michel_fitPass", ci);
  up_location.push_back(univ.GetVecElem("FittedMichel_michel_x1", ci));
  up_location.push_back(univ.GetVecElem("FittedMichel_michel_u1", ci));
  up_location.push_back(univ.GetVecElem("FittedMichel_michel_v1", ci));
  up_location.push_back(univ.GetVecElem("FittedMichel_michel_z1", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_x2", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_u1", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_u2", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_v2", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_z2", ci));
  
  nclusters = univ.GetInt("FittedMichel_cluster_view_sz");

  //double up_to_vertex_dist3d;
  //double down_to_vertex_dist3d;
  //std::vector<Cluster *> cluster_to_up_match;
  //std::vector<Cluster *> cluster_to_down_match;
  //double up_to_cluster_dist3d;
  //double down_to_cluster_dist3d;
}
/*
std::vector<Cluster*> CreateClusters(const CVUniverse& univ)
{
  vector<Cluster*> return_clusters;
  int nclusters = univ.GetInt("FittedMichel_cluster_view_sz");
  for (int i = 0; i < nclusters; ++i)
  {
    Cluster *c = new Cluster(univ, i);
    if (c->energy < 2. || c->energy > 1000.) continue;
    return_clusters.push_back(c);
    
  }
  return return_clusters;
}
*/


void Michel::DoesMichelMatchVtx(const CVUniverse& univ, Michel* &current_michel){

  //std::cout << "GETTING VTX MATCH FOR MICHEL " << std::endl;
  

  double vtx_x = univ.GetVertex().X(); //mm
  double vtx_y = univ.GetVertex().Y(); //mm
  double vtx_z = univ.GetVertex().Z(); //mm
  double vtx_t = univ.GetVertex().T()/pow(10, 3); //mus
  double vtx_u = (0.5 * (vtx_x - sqrt(3.) * vtx_y));
  double vtx_v = (0.5 * (vtx_x + sqrt(3.) * vtx_y));
  
  //std::cout << "VTX POSITION is (x, u , v, y, z) (" << vtx_x << " , " << vtx_u << " , " << vtx_v << " , " << vtx_y << " , " << vtx_z << std::endl;

  double zdiff1 = vtx_z - current_michel->up_location[3];
  double zdiff2 = vtx_z - current_michel->down_location[3];
  double xdiff = 9999.;
  double udiff = 9999.;
  double vdiff = 9999.;
  double XZdist = 9999.; 
  double UZdist = 9999.; 
  double VZdist = 9999.; 

  double michely1 = 9999.;
  double michely2 = 9999.;
  double michelx1 = 9999.;
  double michelx2 = 9999.;
  
  double timediff = current_michel->time - vtx_t;
  
  //std::cout << "DONE INITIALIZING VARIABLES" << std::endl;

  if (timediff < 0.400 || (zdiff1 > 1000. && zdiff2 > 1000.)){
  current_michel->up_to_vertex_dist3D = 9999.;
  current_michel->down_to_vertex_dist3D = 9999.;
  }
  else{
  
  if (zdiff1 <= zdiff2){
  
    //if (zdiff1 > 500.) break;
  
     //std::cout << "GETTING MICHEL UPSTREAM VARS" << std::endl;
  
    xdiff = abs(vtx_x - current_michel->up_location[0]);
    udiff = abs(vtx_u - current_michel->up_location[1]);
    vdiff = abs(vtx_v - current_michel->up_location[2]);
    
    XZdist = sqrt(xdiff*xdiff + zdiff1*zdiff1);
    UZdist = sqrt(udiff*udiff + zdiff1*zdiff1);
    VZdist = sqrt(vdiff*vdiff + zdiff1*zdiff1);

    michelx1 = current_michel->up_location[0];
    
    //std::cout << "COMBINATORIX TO GET BEST  X AND Y calculations" << std::endl;

    if (XZdist <= UZdist && XZdist < VZdist && UZdist <= VZdist){ 
        michely1 = (1. / sqrt(3)) * (current_michel->up_location[0] - 2*current_michel->up_location[1]);
        //std::cout << "XZdist <= UZdist && XZdist < VZdist && UZdist <= VZdist" << std::endl;

    }
    else if (XZdist <= UZdist && XZdist < VZdist && UZdist > VZdist){
        //std::cout << "XZdist <= UZdist && XZdist < VZdist && UZdist > VZdist" << std::endl;
        michely1 = (1. / sqrt(3)) * (2*current_michel->up_location[2] - current_michel->up_location[0]);
        
        }
    else if (UZdist <= XZdist && UZdist < VZdist && VZdist <= XZdist) {
        //std::cout << "UZdist <= XZdist && UZdist < VZdist && VZdist <= XZdist" << std::endl;
       michely1 = (1. / sqrt(3)) * (current_michel->up_location[2] - current_michel->up_location[1]);
       michelx1 = current_michel->up_location[1] + current_michel->up_location[2];
    }
    else if (UZdist <= XZdist && UZdist < VZdist && VZdist > XZdist) 
    {   //std::cout << "UZdist <= XZdist && UZdist < VZdist && VZdist > XZdist" << std::endl;
    
        michely1 = (1. / sqrt(3)) * (current_michel->up_location[0] - 2*current_michel->up_location[1]);}
    else if (VZdist <= XZdist && VZdist < UZdist && XZdist <= UZdist)
    { 
        //std::cout << "VZdist <= XZdist && VZdist < UZdist && XZdist <= UZdist" << std::endl;
    
       michely1 = (1. / sqrt(3)) * (2 * current_michel->up_location[2] - current_michel->up_location[0]);}
    else if (VZdist <= XZdist && VZdist < UZdist && XZdist > UZdist){
     
      //std::cout << "VZdist <= XZdist && VZdist < UZdist && XZdist > UZdist" << std::endl;
    
      michely1 = (1. / sqrt(3)) * (current_michel->up_location[2] - current_michel->up_location[1]);
      michelx1 = current_michel->up_location[1] + current_michel->up_location[2];
    }
    
    //std::cout << "SETTING 2D DISTANCE FOR MICHEL-VTX" << XZdist << "  " << UZdist << "  " << VZdist << "    " << std::endl;

    current_michel->up_to_vertex_dist2D.push_back(XZdist);
    current_michel->up_to_vertex_dist2D.push_back(UZdist);
    current_michel->up_to_vertex_dist2D.push_back(VZdist);

  }
  else if (zdiff2 < zdiff1 && zdiff2 < 1000.){
    //if (zdiff2 > 500.) break;

    //std::cout << "SETTING MICHEL VARS FOR DOWN STREAM" << std::endl;
    xdiff = abs(vtx_x - current_michel->down_location[0]);
    udiff = abs(vtx_u - current_michel->down_location[1]);
    vdiff = abs(vtx_v - current_michel->down_location[2]);
    XZdist = sqrt(xdiff*xdiff + zdiff1*zdiff2);
    UZdist = sqrt(udiff*udiff + zdiff1*zdiff2);
    VZdist = sqrt(vdiff*vdiff + zdiff1*zdiff2);
    michelx2 = current_michel->down_location[0];
    
    if (XZdist < UZdist && XZdist < VZdist && UZdist <= VZdist){
        //std::cout << "XZdist <= UZdist && XZdist < VZdist && UZdist <= VZdist" << std::endl;
    
        michely2 = (1. / sqrt(3)) * (current_michel->up_location[0] - 2*current_michel->up_location[1]);}
    else if (XZdist < UZdist && XZdist < VZdist && UZdist > VZdist){
        //std::cout << "XZdist <= UZdist && XZdist < VZdist && UZdist > VZdist" << std::endl;
        michely2 = (1. / sqrt(3)) * (2*current_michel->up_location[2] - current_michel->up_location[0]);}
    else if (UZdist < XZdist && UZdist < VZdist && VZdist <= XZdist) {
       //std::cout << "UZdist <= XZdist && UZdist < VZdist && VZdist <= XZdist" << std::endl;
    
       michely2 = (1. / sqrt(3)) * (current_michel->up_location[2] - current_michel->up_location[1]);
       michelx2 = current_michel->up_location[1] + current_michel->up_location[2];
    }
    else if (UZdist < XZdist && UZdist < VZdist && VZdist > XZdist){
        //std::cout << "UZdist <= XZdist && UZdist < VZdist && VZdist > XZdist" << std::endl;
        michely2 = (1. / sqrt(3)) * (current_michel->up_location[0] - 2*current_michel->up_location[1]);
     }
    else if (VZdist < XZdist && VZdist < UZdist && XZdist <= UZdist)
    { 
        //std::cout << "VZdist <= XZdist && VZdist < UZdist && XZdist <= UZdist" << std::endl;
        michely2 = (1. / sqrt(3)) * (2 * current_michel->up_location[2] - current_michel->up_location[0]);}
    else if (VZdist < XZdist && VZdist < UZdist && XZdist > UZdist){
      //std::cout << "VZdist <= XZdist && VZdist < UZdist && XZdist > UZdist" << std::endl;
    
      michely2 = (1. / sqrt(3)) * (current_michel->up_location[2] - current_michel->up_location[1]);
      michelx2 = current_michel->up_location[1] + current_michel->up_location[2];
    }
    
    //std::cout << "SETTING 2D DISTANCE FOR MICHEL-VTX" << XZdist << "  " << UZdist << "  " << VZdist << "    " << std::endl;
    current_michel->down_to_vertex_dist2D.push_back(XZdist);
    current_michel->down_to_vertex_dist2D.push_back(UZdist);
    current_michel->down_to_vertex_dist2D.push_back(VZdist);
    
  }
 
  
  //std::cout << "CALCULATING THE 3D DISTANCES" << std::endl;


  double xdiff1 = abs(vtx_x - michelx1);
  double xdiff2 = abs(vtx_x - michelx2);
  double ydiff1 = abs(vtx_y - michely1);
  double ydiff2 = abs(vtx_y - michely2);

  double dist1 = sqrt(zdiff1*zdiff1 + xdiff1*xdiff1 + ydiff1*ydiff1);
  double dist2 = sqrt(zdiff2*zdiff2 + xdiff2*xdiff2 + ydiff2*ydiff2);


  current_michel->up_to_vertex_dist3D = dist1;
  current_michel->down_to_vertex_dist3D = dist2;
  
  //std::cout << "Vtx Dist 1 " << dist1 << " Vtx Dist 2 " << dist2 << std::endl;
  
  }

  //if (dist1 < dist2 && dist1 < 10.2) {
  //  match = true;
  //  dist3D = dist1;
  //}
  //else if (dist1 >= dist2 && dist2 < 10.2)
  //{
  //   match = true;
  //   dist3D = dist2;
  //}
  
  ////std::cout << "END OF LOOP TO FIND VTX MATCH" << std::endl;


}




void Michel::DoesMichelMatchClus(const CVUniverse& univ, Michel* &current_michel){

  //This is where the function for Cluster Matching goes
  
  ////std::cout << "STARTING SEARCH FOR CLUSTER MATCH " << std::endl;


  //int nclusters = univ.GetInt("FittedMichel_cluster_view_sz");
  //CVUniverse* univ;
  int nclusters = current_michel->nclusters;
  

  //std::map<int, Cluster> match_clusters;
  
 
  
 
  double closestdistance1x = 9999.;
  double closestdistance1u = 9999.;
  double closestdistance1v = 9999.;
  double closestdistance1z = 9999.;
 
  double closestdistance2x = 9999.;
  double closestdistance2u = 9999.;
  double closestdistance2v = 9999.;
  double closestdistance2z = 9999.;


  double michelx1 = current_michel->up_location[0];
  double michelx2 = current_michel->down_location[0];
  double michelu1 = current_michel->up_location[1];
  double michelu2 = current_michel->down_location[1];
  double michelv1 = current_michel->up_location[2];
  double michelv2 = current_michel->down_location[2];
  double michelz1 = current_michel->up_location[3];
  double michelz2 = current_michel->down_location[3];
  double michely1 = 9999.;
  double michely2 = 9999.;
  
  double micheltime = current_michel->time;
  
  //std::cout << "Michel position 1 is " << michelx1 << " , " << michelu1 << " , " << michelv1 << " , " << michelz1 << std::endl;
  
  //std::cout << "Michel position 2 is " << michelx2 << " , " << michelu2 << " , " << michelv2 << " , " << michelz2 << std::endl;

  std::vector<Cluster> endpoint1_clus;
  std::vector<Cluster> endpoint2_clus;

// Get the closest distance for each view

  ////std::cout << "STARTING LOOP OVER CLUSTERS " << std::endl;

  for (int i = 0; i < nclusters; i++){

    Cluster current_cluster = Cluster(univ, i);

    double energy = current_cluster.energy;
    double time = current_cluster.time;  
    double pos = current_cluster.pos;   
    double zpos = current_cluster.zpos;  
    int view = current_cluster.view;
    double timediff = micheltime - time;


    if (energy < 2.) continue;
    if (timediff < .400) continue; //only get clusters that occur BEFORE the michel time
    
    //std::cout << "printing cluster info " << "energy " << energy << " time " << time << " pos " << pos << " zpos " << zpos << std::endl;

    double zdiff1 = abs(zpos - michelz1);
    double zdiff2 = abs(zpos - michelz2);
    
    //if (zdiff1 > 1000. && zdiff2 > 1000. ) continue;
    
    //if (zdiff1 > 500. && zdiff2 > 500.) continue;
    
    ////std::cout << "GETTING 2D Distances for each cluster view" << std::endl;


    if (view == 1)
    {

      double xdiff1 = abs(pos - michelx1);

      double x2Ddistance1 = sqrt(xdiff1 * xdiff1 + zdiff1 * zdiff1);

      if (x2Ddistance1 <= closestdistance1x ) closestdistance1x = x2Ddistance1;

      double xdiff2 = abs(pos - michelx2);

      double x2Ddistance2 = sqrt(xdiff2 * xdiff2 + zdiff2 * zdiff2);

      if (x2Ddistance2 <= closestdistance2x) closestdistance2x = x2Ddistance2;
    }
    else if (view == 2)
    {
      double udiff1 = abs(pos - michelu1);

      double u2Ddistance1 = sqrt(udiff1 * udiff1 + zdiff1 * zdiff1);

      if (u2Ddistance1 <= closestdistance1u ) closestdistance1u = u2Ddistance1;

      double udiff2 = abs(pos - michelu2);

      double u2Ddistance2 = sqrt(udiff2 * udiff2 + zdiff2 * zdiff2);

      if (u2Ddistance2 <= closestdistance2u) closestdistance2u = u2Ddistance2;

    }
    else if (view == 3)
    {
      double vdiff1 = abs(pos - michelv1);

      double v2Ddistance1 = sqrt(vdiff1 * vdiff1 + zdiff1 * zdiff1);

      if (v2Ddistance1 <= closestdistance1v ) closestdistance1v = v2Ddistance1;

      double vdiff2 = abs(pos - michelv2);

      double v2Ddistance2 = sqrt(vdiff2 * vdiff2 + zdiff2 * zdiff2);

      if (v2Ddistance2 <= closestdistance2v) closestdistance2v = v2Ddistance2;
    }
  }
  
  double closesum1 = closestdistance1v + closestdistance1u + closestdistance1x;
  double closesum2 = closestdistance2v + closestdistance2u + closestdistance2x;

//Now store the closest X, u, v clusters for each Michel Endpoint based on the above closest distance

  //std::cout << "Looping over clusters again to match to closest distance in each view" << std::endl;
  
  //if (closestdistance1v < 100. || closestdistance1x < 100. || closestdistance1u < 100.){
  //std::cout << "Closest distance x1 " << closestdistance1x << " u1 " << closestdistance1u << " v1 " << closestdistance1v << std::endl;
  //} 
  //if (closestdistance2v < 100. || closestdistance2x < 100. || closestdistance2u < 100.){
  //std::cout << "Closest distance x2 " << closestdistance2x << " u2 " << closestdistance2u << " v2 " << closestdistance2v << std::endl; 
  //}
  for (int i = 0; i < nclusters; i++){
    //if (closesum1 > 1000. && closesum2 > 1000.) continue;
    Cluster current_cluster = Cluster(univ, i);

    double energy = current_cluster.energy;
    double time = current_cluster.time;  
    double pos = current_cluster.pos;   
    double zpos = current_cluster.zpos;  
    int view = current_cluster.view;
    double timediff = micheltime - time;

    if (energy < 2.) continue;
    if (timediff < 0.) continue; 
   

    double zdiff1 = abs(zpos - michelz1);
    double zdiff2 = abs(zpos - michelz2);
    
    if (zdiff1 > 1000. && zdiff2 > 1000.) continue;

    
    //std::cout << "Pushing back clusters ONLY equal to closest" << std::endl;


    if (view == 1)
    {

      double xdiff1 = abs(pos - michelx1);

      double x2Ddistance1 = sqrt(xdiff1 * xdiff1 + zdiff1 * zdiff1);

      if (x2Ddistance1 == closestdistance1x ){
        endpoint1_clus.push_back(current_cluster);
        current_michel->cluster_to_up_match.push_back(current_cluster);
        
      }

      double xdiff2 = abs(pos - michelx2);

      double x2Ddistance2 = sqrt(xdiff2 * xdiff2 + zdiff2 * zdiff2);

      if (x2Ddistance2 == closestdistance2x ){
        endpoint2_clus.push_back(current_cluster);
        current_michel->cluster_to_down_match.push_back(current_cluster);
      }

    }
    else if (view == 2)
    {
      double udiff1 = abs(pos - michelu1);

      double u2Ddistance1 = sqrt(udiff1 * udiff1 + zdiff1 * zdiff1);

      if (u2Ddistance1 == closestdistance1u ){
        endpoint1_clus.push_back(current_cluster);
        current_michel->cluster_to_up_match.push_back(current_cluster);

      }
      double udiff2 = abs(pos - michelu2);

      double u2Ddistance2 = sqrt(udiff2 * udiff2 + zdiff2 * zdiff2);

      if (u2Ddistance2 == closestdistance2u)
      {
        endpoint2_clus.push_back(current_cluster);
        current_michel->cluster_to_down_match.push_back(current_cluster);
      }
    }
    else if (view == 3)
    {
      double vdiff1 = abs(pos - michelv1);

      double v2Ddistance1 = sqrt(vdiff1 * vdiff1 + zdiff1 * zdiff1);

      if (v2Ddistance1 == closestdistance1v ){
        endpoint1_clus.push_back(current_cluster);
        current_michel->cluster_to_up_match.push_back(current_cluster);
      }

      double vdiff2 = abs(pos - michelv2);

      double v2Ddistance2 = sqrt(vdiff2 * vdiff2 + zdiff2 * zdiff2);

       if (v2Ddistance2 == closestdistance2v ){
        endpoint1_clus.push_back(current_cluster);
        current_michel->cluster_to_down_match.push_back(current_cluster);
      }
    }

    
  }
   
 //std::cout << "Pushing back clusters ONLY equal to closest" << std::endl;


 std::vector<double> clusx1;
 std::vector<double> clusx2;

 std::vector<double> clusu1;
 std::vector<double> clusu2;

 std::vector<double> clusv1;
 std::vector<double> clusv2;

 std::vector<double> matchclus1;
 std::vector<double> matchclus2;
 
 //std::cout << "LOOPING OVER ENDPOINT1 CLUSTERS" << std::endl;

 double XZdist1 = 9999.;
 double UZdist1 = 9999.;
 double VZdist1 = 9999.;

  // get 3D point from the enpoint1 clusters

  for (unsigned int i = 0; i < endpoint1_clus.size(); i++){
  
    //std::cout << "Looking at " << i << "endpoint1 cluster" << std::endl;

    Cluster iclus = endpoint1_clus[i];
    double energy = iclus.energy;
    double time = iclus.time;  
    double pos = iclus.pos;   
    double zpos = iclus.zpos;  
    int view = iclus.view;
    
    //std::cout << "FILLING IN INFO FOR FINAL END POINT1 CLUS INFORMATION" << std::endl;

 
    if (view == 1){
     // XZdist = sqrt(pow((pos - michelx1), 2) + pow((zpos - michelz1), 2));
      clusx1.push_back(pos);
      clusx1.push_back(zpos);
      }
    else if (view == 2 ) {
      // UZdist = sqrt(pow((pos - michelu1), 2) + pow((zpos - michelz1), 2));
       clusu1.push_back(pos);
       clusu1.push_back(zpos);
    }
    else if (view == 3){
      //VZdist = sqrt(pow((pos - michelv1), 2) + pow((zpos - michelz1), 2));
      clusv1.push_back(pos);
      clusv1.push_back(zpos);
      }

  }


if (!clusx1.empty()){
    double xdif = michelx1 - clusx1[0];
    double zdif = michelz1 - clusx1[1];
    XZdist1 = sqrt(xdif*xdif + zdif*zdif);}
if (!clusu1.empty()){
    double udif = michelx1 - clusu1[0];
    double zdif = michelz1 - clusu1[1];
    UZdist1 = sqrt(udif*udif + zdif*zdif);}
if (!clusv1.empty()){
    double vdif = michelx1 - clusv1[0];
    double zdif = michelz1 - clusv1[1];
    VZdist1 = sqrt(vdif*vdif + zdif*zdif);}

//std::cout << "GET 3D Information for Clusters - ENDPOINT1" << std::endl;


  if (XZdist1 < UZdist1 && XZdist1 < VZdist1 && UZdist1 <= VZdist1){
    if (!clusu1.empty() && !clusx1.empty()){
    michely1 = (1. / sqrt(3.)) * (current_michel->up_location[0] - 2*current_michel->up_location[1]);
    double yclus = (1. / sqrt(3.)) * (clusx1[0] - 2 * clusu1[0]);
    matchclus1.push_back(clusx1[0]);
    matchclus1.push_back(yclus); // y point of match 3D point
    matchclus1.push_back(clusx1[1]); // seting the cluster 3D point z to be of the closest view
    }
  }
  else if (XZdist1 < UZdist1 && XZdist1 < VZdist1 && UZdist1 > VZdist1){
  
    if (!clusv1.empty() && !clusx1.empty()){
  
    michely1 = (1. / sqrt(3.)) * (2*current_michel->up_location[2] - current_michel->up_location[0]);
    double yclus =  (1. / sqrt(3.)) * (2 * clusv1[0] - clusx1[0]);
    matchclus1.push_back(clusx1[0]);
    matchclus1.push_back(yclus);
    matchclus1.push_back(clusx1[1]); // seting the cluster 3D point z to be of the closest view
    }
  }
  else if (UZdist1 < XZdist1 && UZdist1 < VZdist1 && VZdist1 <= XZdist1) {
    if (!clusv1.empty() && !clusu1.empty()){

     michely1 = (1. / sqrt(3.)) * (current_michel->up_location[2] - current_michel->up_location[1]);
     michelx1 = current_michel->up_location[1] + current_michel->up_location[2];
     double yclus = (1. / sqrt(3.)) * (clusv1[0] - clusu1[0]);
     double xclus = clusu1[0] + clusv1[0];
     matchclus1.push_back(xclus);
     matchclus1.push_back(yclus);
     matchclus1.push_back(clusu1[1]);  // seting the cluster 3D point z to be of the closest view
    }
  }
  else if (UZdist1 < XZdist1 && UZdist1 < VZdist1 && VZdist1 > XZdist1) { 
    if (!clusu1.empty() && !clusx1.empty()){

    michely1 = (1. / sqrt(3.)) * (current_michel->up_location[0] - 2*current_michel->up_location[1]);
    double yclus = (1. / sqrt(3.)) * (clusx1[0] - 2*clusu1[0]);
    matchclus1.push_back(clusx1[0]);
    matchclus1.push_back(yclus);
    matchclus1.push_back(clusu1[1]);  // seting the cluster 3D point z to be of the closest view
    }
  }
  else if (VZdist1 < XZdist1 && VZdist1 < UZdist1 && XZdist1 <= UZdist1) {
    if (!clusv1.empty() && !clusx1.empty()){

    michely1 = (1. / sqrt(3.)) * (2 * current_michel->up_location[2] - current_michel->up_location[0]);
    double yclus = ((1. / sqrt(3.)) * (2 * clusv1[0] - clusx1[0]));
    matchclus1.push_back(clusx1[0]);
    matchclus1.push_back(yclus);
    matchclus1.push_back(clusv1[1]);  // seting the cluster 3D point z to be of the closest view
    }
  }
  else if (VZdist1 < XZdist1 && VZdist1 < UZdist1 && XZdist1 > UZdist1){
    if (!clusu1.empty() && !clusv1.empty()){

    michely1 = (1. / sqrt(3.)) * (current_michel->up_location[2] - current_michel->up_location[1]);
    michelx1 = current_michel->up_location[1] + current_michel->up_location[2];
    double xclus = (1. / sqrt(3.)) * (clusv1[0] - clusu1[0]);
    double yclus = clusu1[0] + clusv1[0];
    matchclus1.push_back(xclus);
    matchclus1.push_back(yclus);
    matchclus1.push_back(clusv1[1]);  // seting the cluster 3D point z to be of the closest view
    }
  }

  // get 3D point from the enpoint2 clusters
  

  double XZdist2 = 9999.;
  double UZdist2 = 9999.;
  double VZdist2 = 9999.;

  for (unsigned int i = 0; i < endpoint2_clus.size(); i++)
  {
    Cluster iclus = endpoint2_clus[i];
    double energy = iclus.energy;
    double time = iclus.time;
    double pos = iclus.pos;
    double zpos = iclus.zpos;
    int view = iclus.view;

    if (view == 1)
    {
      //XZdist2 = sqrt(pow((pos - michelx2), 2) + pow((zpos - michelz2), 2));
      clusx2.push_back(pos);
      clusx2.push_back(zpos);
    }
    else if (view == 2)
    {
      //UZdist2 = sqrt(pow((pos - michelu2), 2) + pow((zpos - michelz2), 2));
      clusu2.push_back(pos);
      clusu2.push_back(zpos);
    }
    else if (view == 3)
    {
      //VZdist2 = sqrt(pow((pos - michelv2), 2) + pow((zpos - michelz2), 2));
      clusv2.push_back(pos);
      clusv2.push_back(zpos);
    }
    
  }
  
     
  if (!clusx2.empty()){
      double xdif = michelx2 - clusx2[0];
      double zdif = michelz2 - clusx2[1];
      XZdist2 = sqrt(xdif*xdif + zdif*zdif);}
  if (!clusu2.empty()){
      double udif = michelx2 - clusu2[0];
      double zdif = michelz2 - clusu2[1];
      UZdist2 = sqrt(udif*udif + zdif*zdif);}
  if (!clusv2.empty()){
    double vdif = michelx2 - clusv2[0];
    double zdif = michelz2 - clusv2[1];
    VZdist2 = sqrt(vdif*vdif + zdif*zdif);}

  
  //std::cout << "GET 3D Information for Clusters - ENDPOINT2" << std::endl;

  
  if (XZdist2 < UZdist2 && XZdist2 < VZdist2 && UZdist2 <= VZdist2){
    //std::cout << "GET 3D Information CASE1" << std::endl;    
    if (!clusu2.empty() && !clusx2.empty()){
    
    michely2 = (1. / sqrt(3.)) * (current_michel->down_location[0]- 2*current_michel->down_location[1]);
    double yclus = (1. / sqrt(3.)) * (clusx2[0] - 2 * clusu2[0]);
    matchclus2.push_back(clusx2[0]);
    matchclus2.push_back(yclus);
    matchclus2.push_back(clusx2[1]); // seting the cluster 3D point z to be of the closest view
    }
  }
  else if (XZdist2 < UZdist2 && XZdist2 < VZdist2 && UZdist2 > VZdist2){
    //std::cout << "GET 3D Information CASE2" << std::endl;
    if (!clusv2.empty() && !clusx2.empty()){

    michely2 = (1. / sqrt(3.)) * (2*current_michel->down_location[2] - current_michel->down_location[0]);
    double yclus = (1. / sqrt(3.)) * (2 * clusv2[0] - clusx2[0]);
    matchclus2.push_back(clusx2[0]);
    matchclus2.push_back(yclus);
    matchclus2.push_back(clusx2[1]); // seting the cluster 3D point z to be of the closest view
    }
  }
  else if (UZdist2 < XZdist2 && UZdist2 < VZdist2 && VZdist2 <= XZdist2) {
     //std::cout << "GET 3D Information CASE3" << std::endl;
     if (!clusu2.empty() && !clusv2.empty()){

     michely2 = (1. / sqrt(3.)) * (current_michel->down_location[2] - current_michel->down_location[1]);
     michelx2 = current_michel->down_location[1] + current_michel->down_location[2];
     double xclus = clusu2[0] + clusv2[0];
     double yclus = (1. / sqrt(3.)) * (clusv2[0] - clusu2[0]);
     matchclus2.push_back(xclus);
     matchclus2.push_back(yclus);
     matchclus2.push_back(clusu2[1]);  // seting the cluster 3D point z to be of the closest view
    }
  }
  else if (UZdist2 < XZdist2 && UZdist2 < VZdist2 && VZdist2 > XZdist2) { 
    //std::cout << "GET 3D Information CASE4" << std::endl;
    if (!clusu2.empty() && !clusx2.empty()){
    michely2 = (1. / sqrt(3.)) * (current_michel->down_location[0] - 2*current_michel->up_location[1]);
    double yclus =(1. / sqrt(3.)) * (clusx2[0] - 2*clusu2[0]);
    matchclus2.push_back(clusx2[0]);
    matchclus2.push_back(yclus);
    matchclus2.push_back(clusu2[1]);  // seting the cluster 3D point z to be of the closest view
  
    }
  }
  else if (VZdist2 < XZdist2 && VZdist2 < UZdist2 && XZdist2 <= UZdist2) {
    //std::cout << "GET 3D Information CASE5" << std::endl;
    if (!clusv2.empty() && !clusx2.empty()){
    michely2 = (1. / sqrt(3.)) * (2 * current_michel->down_location[2] - current_michel->down_location[0]);
    double yclus = ((1. / sqrt(3.)) * (2 * clusv2[0] - clusx2[0]));
    matchclus2.push_back(clusx2[0]);
    matchclus2.push_back(yclus);
    matchclus2.push_back(clusv2[1]);  // seting the cluster 3D point z to be of the closest view
  
    }
  }
  else if (VZdist2 < XZdist2 && VZdist2 < UZdist2 && XZdist2 > UZdist2){
    //std::cout << "GET 3D Information CASE6" << std::endl;
    if (!clusu2.empty() && !clusv2.empty()){
    michely2 = (1. / sqrt(3.)) * (current_michel->down_location[2] - current_michel->down_location[1]);
    michelx2 = current_michel->down_location[1] + current_michel->down_location[2];
    double xclus = (1. / sqrt(3.)) * (clusv2[0] - clusu2[0]);
    double yclus = clusu2[0] + clusv2[0];
    matchclus2.push_back(xclus);
    matchclus2.push_back(yclus);
    matchclus2.push_back(clusv2[1]);  // seting the cluster 3D point z to be of the closest view
  
    }
  }
  
  //std::cout << "GETTING 3D DISTANCES FOR THE CUSTERS " << std::endl;
  double clusx1diff = 9999.;
  double clusy1diff = 9999.;
  double clusz1diff = 9999.;
  if (!matchclus1.empty()){
    clusx1diff = michelx1 - matchclus1[0];
    clusy1diff = michely1 - matchclus1[1];
    clusz1diff = michelz1 - matchclus1[2];
  }
  double clusx2diff = 9999.;
  double clusy2diff = 9999.;
  double clusz2diff = 9999.;
  
  if (!matchclus2.empty()){
  clusx2diff = michelx2 - matchclus2[0];
  clusy2diff = michely2 - matchclus2[1];
  clusz2diff = michelz2 - matchclus2[2];
  }
  double dist1 = sqrt(pow(clusx1diff, 2) + pow(clusy1diff, 2) + pow(clusz1diff, 2));
  double dist2 = sqrt(pow(clusx2diff, 2) + pow(clusy2diff, 2) + pow(clusz2diff, 2));


 if (dist1 <= dist2) current_michel->up_to_cluster_dist3D = dist1;
 else if (dist1 > dist2) current_michel->down_to_cluster_dist3D = dist2;
 
 //if (dist1 < 100. || dist2 < 100.) //std::cout << "Clus Dist 1 " << dist1 << " Clus Dist 2 " << dist2 << endl;



}






#endif // MICHEL_H
