#ifndef STDAFX_H
#define STDAFX_H

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <map>
#include <stdint.h>
#include "c_utils.h"
#include <iostream>
using namespace Eigen;




//hash bucket
typedef uint KeyType;

typedef struct hashnode_i
{
  KeyType key;
  void *data;
  struct hashnode_i *next;
} hashnode_i ;


typedef struct HSHTBL_i
{
  size_t size;
  struct hashnode_i **nodes;
  size_t (*hashfunc)(uint);
} hashtable_int;

typedef struct THash
{
  int id;
  int i, ppfInd;
  float angle;
} THash;

    // pose.alpha = alpha;
    // pose.modelIndex = refIndMax;
    // pose.numVotes = maxVotes;
    // pose.pose = rawPose;


typedef struct Pose_3DM2E
{
    double alpha;//
    double residual;
    double angle;
    size_t modelIndex;//
    size_t numVotes;//
    Eigen::Matrix4f pose;//
    Eigen::Vector3f t;
    Eigen::Vector4f q;
    void computeqandt(){
      
      Eigen::Vector4f q_temp;
      Eigen::Matrix3f R =pose.block(0,0,3,3);
      dcmToQuat(R, q_temp);
      t=pose.block(0, 3, 3, 1);
      q=q_temp;
     
    }
    void computeangle(){
      Eigen::Matrix3f R =pose.block(0,0,3,3);
      Vec3f axis;
      double  angle1;
      dcmToAA(R, axis, angle1);
      angle=angle1;
      
    }

    void computepose(){
      Eigen::Matrix3f R ;
      // std::cout<<"ok11"<<std::endl;
      quatToDCM(q, R);
      // std::cout<<"ok"<<std::endl;
      rtToPose(R,t,pose);
      // std::cout<<"ok22"<<std::endl;
    }
    Pose_3DM2E& operator = (const  Pose_3DM2E& a){
      std::cout<<"init";
      Pose_3DM2E b;
      b.alpha=a.alpha;
      b.residual=a.residual;
      b.angle=a.angle;
      b.modelIndex=a.modelIndex;//
      b.numVotes=a.numVotes;//
     std::cout<<"a.pose:"<<a.pose<<std::endl;
      b.t=a.t;
      b.q=a.q;
      b.computepose();
      std::cout<<"b.pose:"<<b.pose<<std::endl;
    
      return b;

    }



} Pose_3DM2E;



typedef struct ClusterM2E
{
    std::vector<Pose_3DM2E> poses;
    size_t accu_votes;
    ClusterM2E(){
      accu_votes=0;
      // std::cout<<"调用ClusterM2E()"<<std::endl;
    }

   
} ClusterM2E;



#endif // STDAFX_H
