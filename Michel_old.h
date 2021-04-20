#ifndef MICHEL_H
#define MICHEL_H

#ifndef __CINT__

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
using MichelMap = std::map<int, Michel>;

//// old school typedef that still works
//typedef std::map<int, Michel> MichelMap;

class Michel
{
 public:
  //constructors
  Michel(const CVUniverse& univ, int ci);

  // Functions
  Michel(CVUniverse univ){};                  // fill in most basic stuff in "generic info"
  void DoMoreInitializationAndProcessing(); // fill in more complicated stuff in "generic info"
  void DoMatching();                          // fill in "best matching stuff"
  void DoesMichelMatchVtx(const CVUniverse& univ, Michel current_michel); //GEts info for Vtx Match
  void DoesMichelMatchClus(const CVUniverse& univ, Michel current_michel); // Gets info for ClusterMatch

  // generic info
  std::vector<double> up_location; // upstream location 0 X 1 U 2 V 3 Z
  std::vector<double> down_location; //downstream location
  double energy;
  double time;
  int is_fitted;
  std::vector<double> up_to_vertex_dist2d;
  std::vector<double> down_to_vertex_dist2d;
  double up_to_vertex_dist3d;
  double down_to_vertex_dist3d;
  std::vector<Cluster *> cluster_to_up_match;
  std::vector<Cluster *> cluster_to_down_match;
  double up_to_cluster_dist3d;
  double down_to_cluster_dist3d;
  
  
  // best matching stuff
  enum BestClusterMatch {kUpMatch, kDownMatch, kNClusterMatches}; // just tells you which of the four matches are the best match
};

Michel::Michel(const CVUniverse& univ, int ci) 
  : 
    energy(univ.GetVecElem("FittedMichel_michel_energy", ci)),
    time(univ.GetVecElem("FittedMichel_michel_time", ci)),
    is_fitted(univ.GetVecElem("FittedMichel_michel_fitPass", ci)) 
{   up_location.push_back(univ.GetVecElem("FittedMichel_michel_x1",ci));
    up_location.push_back(univ.GetVecElem("FittedMichel_michel_u1", ci));
    up_location.push_back(univ.GetVecElem("FittedMichel_michel_v1", ci));
    up_location.push_back(univ.GetVecElem("FittedMichel_michel_z1", ci));
    down_location.push_back(univ.GetVecElem("FittedMichel_michel_x2",ci));
    down_location.push_back(univ.GetVecElem("FittedMichel_michel_u1",ci));
    down_location.push_back(univ.GetVecElem("FittedMichel_michel_u2",ci));
    down_location.push_back(univ.GetVecElem("FittedMichel_michel_v2", ci));
    down_location.push_back(univ.GetVecElem("FittedMichel_michel_z2", ci));
                 
}

vector<Cluster *> CreateClusters(CVUniverse &univ)
{
  std::vector<Cluster *> return_clusters;
  int nclusters = univ.GetInt("cluster_sz");
  for (int i = 0; i < ncluster; ++i)
  {
    Cluster *c = new Cluster(univ, i);
    return_cluster.push_back(c);
  }
  return return_clusters;
}


//==============================================================================
// Collect event's quality michels and require at most one-to-one
// cluster-to-vertex matching.
//
// An interaction vertex can have at most one michel cluster match (this is a
// restriction built into the michel tool).
// 
// But a michel cluster can be matched to more than one vertex.
//
// This function selects michel clusters that are well-matched to track
// endpoints.
//==============================================================================


void DoesMichelMatchVtx(const CVUniverse& univ, Michel current_michel){

  double vtx_x = univ.GetVertex().X() / 10.; //cm
  double vtx_y = univ.GetVertex().Y() / 10.; //cm
  double vtx_z = univ.GetVertex().Z() / 10.; //cm
  double vtx_t = univ.GetVertex().T(); //cm
  double vtx_u = (0.5 * (vtx_x - sqrt(3) * vtx_y))/10.;
  double vtx_v = (0.5 * (vtx_x + sqrt(3) * vtx_y))/10.;

  double zdiff1 = abs(vtx_z - current_michel.endpointz1);
  double zdiff2 = abs(vtx_z - current_michel.endpointz2);
  double xdiff;
  double udiff;
  double vdiff;
  double XZdist; 
  double UZdist; 
  double VZdist; 

  double michely1;
  double michely2;
  double michelx1;
  double michelx2;

  if (zdiff1 <= zdiff2){
    xdiff = abs(vtx_x - current_michel.endpointx1);
    udiff = abs(vtx_u - current_michel.endpointu1);
    vdiff = abs(vtx_v - current_michel.endpointv1);
    
    XZdist = sqrt(xdiff*xdiff + zdiff1*zdiff1);
    UZdist = sqrt(udiff*udiff + zdiff1*zdiff1);
    VZdist = sqrt(vdiff*vdiff + zdiff1*zdiff1);

    michelx1 = current_michel.endpointx1;

    if (XZdist < UZdist && XZdist < VZdist && UZdist <= VZdist) michely1 = (1. / sqrt(3)) * (current_michel.endpointx1 - 2*current_michel.endpointu1);
    else if (XZdist < UZdist && XZdist < VZdist && UZdist > VZdist) michely1 = (1. / sqrt(3)) * (2*curren_michel.endpointv1 - current_michel.endpointx1);
    else if (UZdist < XZdist && UZdist < VZdist && VZdist <= XZdist) {
       michely1 = (1. / sqrt(3)) * (current_michel.endpointv1 - current_michel.endpointu1);
       michelx1 = current_michel.endpointu1 + current_michel.endpointv1;
    }
    else if (UZdist < XZdist && UZdist < VZdist && VZdist > XZdist) michely1 = (1. / sqrt(3)) * (current_michel.endpointx1 - 2*current_michel.endpointu1);
    else if (VZdist < XZdist && VZdist < UZdist && XZdist <= UZdist) michely1 = (1. / sqrt(3)) * (2 * curren_michel.endpointv1 - current_michel.endpointx1);
    else if (VZdist < XZdist && VZdist < UZdist && XZdist > UZdist){
      michely1 = (1. / sqrt(3)) * (current_michel.endpointv1 - current_michel.endpointu1);
      michelx1 = current_michel.endpointu1 + current_michel.endpointv1;
    }

    up_to_vertex_dist2D.push_back(XZdist,UZdist,VZdist);

  }
  else if (zdiff2 < zdiff1){
    xdiff = abs(vtx_x - current_michel.endpointx2);
    udiff = abs(vtx_u - current_michel.endpointu2);
    vdiff = abs(vtx_v - current_michel.endpointv2);
    XZdist = sqrt(xdiff*xdiff + zdiff1*zdiff2);
    UZdist = sqrt(udiff*udiff + zdiff1*zdiff2);
    VZdist = sqrt(vdiff*vdiff + zdiff1*zdiff2);
    michelx2 = current_michel.endpointx2;
    if (XZdist < UZdist && XZdist < VZdist && UZdist <= VZdist) michely2 = (1. / sqrt(3)) * (current_michel.endpointx1 - 2*current_michel.endpointu1);
    else if (XZdist < UZdist && XZdist < VZdist && UZdist > VZdist) michely2 = (1. / sqrt(3)) * (2*curren_michel.endpointv1 - current_michel.endpointx1);
    else if (UZdist < XZdist && UZdist < VZdist && VZdist <= XZdist) {
       michely2 = (1. / sqrt(3)) * (current_michel.endpointv1 - current_michel.endpointu1);
       michelx2 = current_michel.endpointu1 + current_michel.endpointv1;
    }
    else if (UZdist < XZdist && UZdist < VZdist && VZdist > XZdist) michely2 = (1. / sqrt(3)) * (current_michel.endpointx1 - 2*current_michel.endpointu1);
    else if (VZdist < XZdist && VZdist < UZdist && XZdist <= UZdist) michely2 = (1. / sqrt(3)) * (2 * curren_michel.endpointv1 - current_michel.endpointx1);
    else if (VZdist < XZdist && VZdist < UZdist && XZdist > UZdist){
      michely2 = (1. / sqrt(3)) * (current_michel.endpointv1 - current_michel.endpointu1);
      michelx2 = current_michel.endpointu1 + current_michel.endpointv1;
    }

  }

  double xdiff1 = abs(vtx_x - michelx1);
  double xdiff2 = abs(vtx_x - michelx2);
  double ydiff1 = abs(vtx_y - michely1);
  double ydiff2 = abs(vtx_y - michely2);

  double dist1 = sqrt(zdiff1*zdiff1 + xdiff1*xdiff1 + ydiff1*ydiff1);
  double dist2 = sqrt(zdiff2*zdiff2 + xdiff2*xdiff2 + ydiff2*ydiff2);

  up_to_vertex_dist3d = dist1;
  down_to_vertex_dist3d = dist2;

  vtxdist3D1 = dist1;
  vtxdist3D2 = dist2;

  //if (dist1 < dist2 && dist1 < 10.2) {
  //  match = true;
  //  dist3D = dist1;
  //}
  //else if (dist1 >= dist2 && dist2 < 10.2)
  //{
  //   match = true;
  //   dist3D = dist2;
  //}

}


void DoesMichelMatchClus(const CVUniverse& univ, Michel current_michel){

  //This is where the function for Cluster Matching goes

  int nclusters = univ.Getint("FittedMichel_cluster_view_sz");

  std::map<int, Cluster> match_clusters;
 
  double closestdistance1x = 9999.;
  double closestdistance1u = 9999.;
  double closestdistance1v = 9999.;
  double closestdistance1z = 9999.;
 
  double closestdistance2x = 9999.;
  double closestdistance2u = 9999.;
  double closestdistance2v = 9999.;
  double closestdistance2z = 9999.;


  double michelx1 = current_michel.endpointx1;
  double michelx2 = current_michel.endpointx2;
  double michelu1 = current_michel.endpointu1;
  double michelu2 = current_michel.endpointu2;
  double michelv1 = current_michel.endpointv1;
  double michelv2 = current_michel.endpointv2;
  double michelz1 = current_michel.endpointz1;
  double michelz2 = current_michel.endpointz2;
  double michely1;
  double michely2;

  std::vector<Cluster> endpoint1_clus;
  std::vector<Cluster> endpoint2_clus;

// Get the closest distance for each view
  for (int i = 0; i < nclusters; i++){

    Cluster current_cluster = Cluster(univ, i);

    double energy = current_cluster.energy;
    double time = current_cluster.time;  
    double pos = current_cluster.pos;   
    double zpos = current_cluster.zpos;  
    int view = current_cluster.view;

    if (energy < 2.) continue;

    double zdiff1 = abs(zpos - michelz1);
    double zdiff2 = abs(zpos - michelz2);

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

//Now store the closest X, u, v clusters for each Michel Endpoint based on the above closest distance
  for (int i = 0; i < nclusters; i++){

    Cluster current_cluster = Cluster(univ, i);

    double energy = current_cluster.energy;
    double time = current_cluster.time;  
    double pos = current_cluster.pos;   
    double zpos = current_cluster.zpos;  
    int view = current_cluster.view;

    if (energy < 2.) continue;

    double zdiff1 = abs(zpos - michelz1);
    double zdiff2 = abs(zpos - michelz2);

    if (view == 1)
    {

      double xdiff1 = abs(pos - michelx1);

      double x2Ddistance1 = sqrt(xdiff1 * xdiff1 + zdiff1 * zdiff1);

      if (x2Ddistance1 == closestdistance1x ){
        endpoint1_clus.push_back(current_cluster)
          cluster_to_up_match.push_back(current_cluster);
        
      }

      double xdiff2 = abs(pos - michelx2);

      double x2Ddistance2 = sqrt(xdiff2 * xdiff2 + zdiff2 * zdiff2);

      if (x2Ddistance2 == closestdistance2x ){
        endpoint2_clus.push_back(current_cluster);
        cluster_to_down_match.push_back(current_cluster);
      }

    }
    else if (view == 2)
    {
      double udiff1 = abs(pos - michelu1);

      double u2Ddistance1 = sqrt(udiff1 * udiff1 + zdiff1 * zdiff1);

      if (u2Ddistance1 == closestdistance1u ){
        endpoint1_clus.push_back(current_cluster);
        cluster_to_up_match.push_back(current_cluster);

      }
      double udiff2 = abs(pos - michelu2);

      double u2Ddistance2 = sqrt(udiff2 * udiff2 + zdiff2 * zdiff2);

      if (u2Ddistance2 == closestdistance2u)
      {
        endpoint2_clus.push_back(current_cluster);
        cluster_to_down_match.push_back(current_cluster);
      }
    }
    else if (view == 3)
    {
      double vdiff1 = abs(pos - michelv1);

      double v2Ddistance1 = sqrt(vdiff1 * vdiff1 + zdiff1 * zdiff1);

      if (v2Ddistance1 == closestdistance1v ){
        endpoint1_clus.push_back(current_cluster);
        cluster_to_up_match.push_back(current_cluster);
      }

      double vdiff2 = abs(pos - michelv2);

      double v2Ddistance2 = sqrt(vdiff2 * vdiff2 + zdiff2 * zdiff2);

       if (v2Ddistance2 == closestdistance2v ){
        endpoint1_clus.push_back(current_cluster);
        cluster_to_down_match.push_back(current_cluster);
      }
    }

    
  }

 std::vector<double> clusx1;
 std::vector<double> clusx2;

 std::vector<double> clusu1;
 std::vector<double> clusu2;

 std::vector<double> clusv1;
 std::vector<double> clusv2;

 std::vector<double> matchclus1;
 std::vector<double> matchclus2;


  // get 3D point from the enpoint1 clusters

  for (int i = 0; i < endpoint1_clus.size(); i++){
    Cluster iclusx = endpoint1_clus[i];
    double energy = iclus.energy;
    double time = iclus.time;  
    double pos = iclus.pos;   
    double zpos = iclus.zpos;  
    int view = iclus.view;
 
    if (view == 1){
      XZdist = sqrt(pow((pos - michelx1), 2) + pow((zpos - michelz1), 2));
      clusx1[0] = pos;
      clusx1[1] = zpos;
      }
    else if (view == 2 ) {
       UZdist = sqrt(pow((pos - michelu1), 2) + pow((zpos - michelz1), 2));
       clusu1[0] = pos;
       clusu1[1] = zpos;
    }
    else if (view == 3){
      VZdist = sqrt(pow((pos - michelv1), 2) + pow((zpos - michelz1), 2));
      clusv1[0] = pos;
      clusv1[1] = zpos;
      }

  }


double XZdist1 = sqrt(clusx1[0]*clusx1[0] + clusx1[1]*clusx1[1]);
double UZdist1 = sqrt(clusu1[0]*clusu1[0] + clusu1[1]*clusu1[1]);
double VZdist1 = sqrt(clusv1[0]*clusv1[0] + clusv1[1]*clusv1[1]);

  if (XZdist1 < UZdist1 && XZdist1 < VZdist1 && UZdist1 <= VZdist1){
    michely1 = (1. / sqrt(3)) * (current_michel.endpointx1 - 2*current_michel.endpointu1);
    matchclus1[0] = clusx1[0];
    matchclus1[1] = (1. / sqrt(3)) * (clusx1[0] - 2 * clusu1[0]);
    matchclus1[2] = clusx1[1]; // seting the cluster 3D point z to be of the closest view
    }
  else if (XZdist1 < UZdist1 && XZdist1 < VZdist1 && UZdist1 > VZdist1){
    michely1 = (1. / sqrt(3)) * (2*current_michel.endpointv1 - current_michel.endpointx1);
    matchclus1[0] = clusx1[0];
    matchclus1[1] = (1. / sqrt(3)) * (2 * clusv1[0] - clusx1[0]);
    matchclus1[2] = clusx1[1]; // seting the cluster 3D point z to be of the closest view

  }
  else if (UZdist1 < XZdist1 && UZdist1 < VZdist1 && VZdist1 <= XZdist1) {
     michely1 = (1. / sqrt(3)) * (current_michel.endpointv1 - current_michel.endpointu1);
     michelx1 = current_michel.endpointu1 + current_michel.endpointv1;
     matchclus1[0] = clusu1[0] + clusv1[0];
     matchclus1[1] = (1. / sqrt(3)) * (clusv1[0] - clusu1[0]);
     matchclus1[2] = clusu1[1];  // seting the cluster 3D point z to be of the closest view
  }
  else if (UZdist1 < XZdist1 && UZdist1 < VZdist1 && VZdist1 > XZdist1) { 
    michely1 = (1. / sqrt(3)) * (current_michel.endpointx1 - 2*current_michel.endpointu1);
    matchclus1[0] = clusx1[0];
    matchclus1[1] = (1. / sqrt(3)) * (clusx1[0] - 2*clusu1[0]);
    matchclus1[2] = clusu1[1];  // seting the cluster 3D point z to be of the closest view
  }
  else if (VZdist1 < XZdist1 && VZdist1 < UZdist1 && XZdist1 <= UZdist1) {
    michely1 = (1. / sqrt(3)) * (2 * curren_michel.endpointv1 - current_michel.endpointx1);
    matchclus1[0] = clusx1[0];
    matchclus1[1] = ((1. / sqrt(3)) * (2 * clusv1[0] - clusx1[0]);
    matchclus1[2] = clusv1[1];  // seting the cluster 3D point z to be of the closest view
  }
  else if (VZdist1 < XZdist1 && VZdist1 < UZdist1 && XZdist1 > UZdist1){
    michely1 = (1. / sqrt(3)) * (current_michel.endpointv1 - current_michel.endpointu1);
    michelx1 = current_michel.endpointu1 + current_michel.endpointv1;
    matchclus1[1] = (1. / sqrt(3)) * (clusv1[0] - clusu1[0]);
    matchclus1[0] = clusu1[0] + clusv1[0];
    matchclus1[2] = clusv1[1];  // seting the cluster 3D point z to be of the closest view
  }

  // get 3D point from the enpoint2 clusters

  for (int i = 0; i < endpoint2_clus.size(); i++)
  {
    Cluster iclusx = endpoint2_clus[i];
    double energy = iclus.energy;
    double time = iclus.time;
    double pos = iclus.pos;
    double zpos = iclus.zpos;
    int view = iclus.view;

    if (view == 1)
    {
      XZdist2 = sqrt(pow((pos - michelx2), 2) + pow((zpos - michelz2), 2));
      clusx2[0] = pos;
      clusx2[1] = zpos;
    }
    else if (view == 2)
    {
      UZdist2 = sqrt(pow((pos - michelu2), 2) + pow((zpos - michelz2), 2));
      clusu2[0] = pos;
      clusu2[1] = zpos;
    }
    else if (view == 3)
    {
      VZdist2 = sqrt(pow((pos - michelv2), 2) + pow((zpos - michelz2), 2));
      clusv2[0] = pos;
      clusv2[1] = zpos;
    }
  }


  double XZdist2 = sqrt(clusx2[0]*clusx2[0] + clusx2[1]*clusx2[1]);
  double UZdist2 = sqrt(clusu2[0]*clusu2[0] + clusu2[1]*clusu2[1]);
  double VZdist2 = sqrt(clusv2[0]*clusv2[0] + clusv2[1]*clusv2[1]);

  if (XZdist2 < UZdist2 && XZdist2 < VZdist2 && UZdist2 <= VZdist2){
    michely2 = (1. / sqrt(3)) * (current_michel.endpointx2- 2*current_michel.endpointu2);
    matchclus2[0] = clusx2[0];
    matchclus2[1] = (1. / sqrt(3)) * (clusx2[0] - 2 * clusu2[0]);
    matchclus2[2] = clusx2[1]; // seting the cluster 3D point z to be of the closest view
    }
  else if (XZdist2 < UZdist2 && XZdist2 < VZdist2 && UZdist2 > VZdist2){
    michely2 = (1. / sqrt(3)) * (2*current_michel.endpointv2 - current_michel.endpointx2);
    matchclus2[0] = clusx2[0];
    matchclus2[1] = (1. / sqrt(3)) * (2 * clusv2[0] - clusx2[0]);
    matchclus2[2] = clusx2[1]; // seting the cluster 3D point z to be of the closest view

  }
  else if (UZdist2 < XZdist2 && UZdist2 < VZdist2 && VZdist2 <= XZdist2) {
     michely2 = (1. / sqrt(3)) * (current_michel.endpointv2 - current_michel.endpointu2);
     michelx2 = current_michel.endpointu2 + current_michel.endpointv2;
     matchclus2[0] = clusu2[0] + clusv2[0];
     matchclus2[1] = (1. / sqrt(3)) * (clusv2[0] - clusu2[0]);
     matchclus2[2] = clusu2[1];  // seting the cluster 3D point z to be of the closest view
  }
  else if (UZdist2 < XZdist2 && UZdist2 < VZdist2 && VZdist2 > XZdist2) { 
    michely2 = (1. / sqrt(3)) * (current_michel.endpointx2 - 2*current_michel.endpointu12;
    matchclus2[0] = clusx1[0];
    matchclus2[1] = (1. / sqrt(3)) * (clusx2[0] - 2*clusu2[0]);
    matchclus2[2] = clusu2[1];  // seting the cluster 3D point z to be of the closest view
  }
  else if (VZdist2 < XZdist2 && VZdist2 < UZdist2 && XZdist2 <= UZdist2) {
    michely2 = (1. / sqrt(3)) * (2 * curren_michel.endpointv2 - current_michel.endpointx2);
    matchclus2[0] = clusx2[0];
    matchclus2[1] = ((1. / sqrt(3)) * (2 * clusv2[0] - clusx2[0]);
    matchclus2[2] = clusv2[1];  // seting the cluster 3D point z to be of the closest view
  }
  else if (VZdist2 < XZdist2 && VZdist2 < UZdist2 && XZdist2 > UZdist2){
    michely2 = (1. / sqrt(3)) * (current_michel.endpointv2 - current_michel.endpointu2);
    michelx2 = current_michel.endpointu2 + current_michel.endpointv2;
    matchclus2[1] = (1. / sqrt(3)) * (clusv2[0] - clusu2[0]);
    matchclus2[0] = clusu2[0] + clusv2[0];
    matchclus2[2] = clusv2[1];  // seting the cluster 3D point z to be of the closest view
  }

  double clusx1diff = michelx1 - matchclus1[0];
  double clusx2diff = michelx2 - matchclus2[0];
  double clusy1diff = michely1 - matchclus1[1];
  double clusy2diff = michely2 - matchclus2[1];
  double clusz1diff = michelz1 - matchclus1[2];
  double clusz2diff = michelz2 - matchclus2[2];
  double dist1 = sqrt(pow(clusx1diff, 2) + pow(clusy1diff, 2) + pow(clusz1diff, 2));
  double dist2 = sqrt(pow(clusx2diff, 2) + pow(clusy2diff, 2) + pow(clusz2diff, 2));

  if (dist1 < 10.2  && dist2 < 10.2 ){
    match = true;
    if (dist1 <= dist2) up_to_cluster_dist3d = dist1;
    else if (dist1 > dist2) down_to_cluster_dist3d = dist2;
  }
  else{
    match = false;
    if (dist1 <= dist2) up_to_cluster_dist3d = dist1;
    else if (dist1 > dist2) down_to_cluster_dist3d = dist2;
  } 

}

#endif // __CINT__

#endif // MICHEL_H
