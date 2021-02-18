#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/point_cloud.h>
#include <pcl/console/parse.h>
#include <pcl/common/transforms.h>	//	pcl::transformPointCloudWithNormals 用到这个头文件
#include <pcl/visualization/pcl_visualizer.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <thread>
#include <pcl/features/ppf.h>
#include <pcl/io/pcd_io.h>
#include <pcl/features/normal_3d.h>
#include <pcl/registration/ppf_registration.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include </home/lcy/Documents/pcl_demo/PPFfeaturelcy.hpp>
#include</home/lcy/Documents/pcl_demo/PPFlcy.hpp>
#include <pcl/visualization/cloud_viewer.h>       
#include <pcl/filters/voxel_grid.h>
#include <string>
#include "hash_murmur.hpp"
#include "hash_murmur64.hpp"
#include "t_hash_int.hpp"
#include "stdafx.h"
#include <pcl/features/normal_3d.h>
#include <pcl/surface/gp3.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/octree/octree.h>
#include <pcl/octree/octree_pointcloud_changedetector.h>
#include<pcl/common/geometry.h>
#include <pcl/surface/convex_hull.h>
#include <pcl/io/pcd_io.h>
#include <pcl/common/impl/io.hpp>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <fstream>
#include <math.h>
#include <time.h>//



using namespace pcl;
using namespace std;
using namespace Eigen;



#define PI 3.1415926
#define EPS     1e-8

hashtable_int* hash_tableM2E = NULL;//最大的hash表
THash* hash_nodesM2E = NULL;


static const size_t PPF_LENGTH = 5;
typedef Eigen::Matrix<float, 5, 1> Vector5f;


//hashPPFM2E  angle_step_radiansM2E, distanceStep
const double angle_step_relativeM2E = 15;//30
const double angle_step_radiansM2E = (360.0/angle_step_relativeM2E)*M_PI/180.0;
const double distanceStep = 0.008;

//CalculateNormals
const float search_radius = 0.1f;
// const float leaf = 0.05;//0.02


const double sampling_step_relativeM2E = 0.05;
const double distance_step_relativeM2E = 0.05;
const double scene_sample_stepM2E = (int)(1/0.04);
const double angle_stepM2E = angle_step_radiansM2E;

const int numAngles = (int) (floor (2 * M_PI / angle_stepM2E));
const double relativeSceneSampleStep = 1.0/2.0; // 重要
const double relativeSceneDistance = 0.03;





//use_weighted_avg
const bool use_weighted_avg=true;
int model_Index0[30];
int model_Index1[30];
int model_Index2[30];
int temp0=0;
int temp1=0;
int temp2=0;

//求 (nx,ny,nz) 、(mx,my,mz) 角度
double  get_angle(float nx,float ny,float nz,float mx,float my,float mz){
    Eigen::Vector3f v1{nx,ny,nz};
    Eigen::Vector3f v2{mx,my,mz};
    v1.normalized();
    v2.normalized();
    double angle= (v1.dot(v2))/(v1.norm()*v2.norm());
    angle=acos(angle);
    angle=180*angle/PI;
    return angle;
}
double get_volume(PointCloud<PointXYZ>::Ptr& cloud){
  if(cloud->points.size()==1) return 0;
   pcl::PointXYZ min_p, max_p;
	pcl::getMinMax3D(*cloud, min_p, max_p);
  double volume=(max_p.x-min_p.x)*(max_p.y-min_p.y)*(max_p.z-min_p.z);
  volume=volume*1000000;
  return volume;
}

//这里有一个角度的判断;懒得改了   
//计算映射在yoz平面的角度的计算
double compute_a(pcl::PointNormal point1,pcl::PointNormal point2){
  Eigen::Affine3f transform_2 = Eigen::Affine3f::Identity();
  transform_2.translation() << -point1.x,-point1.y,-point1.z;
  // std::cout << transform_2.matrix() << std::endl;
  pcl::PointCloud<pcl::PointNormal>::Ptr transformed_cloud (new pcl::PointCloud<pcl::PointNormal> ());
  pcl::PointCloud<pcl::PointNormal>::Ptr source_cloud (new pcl::PointCloud<pcl::PointNormal> ());
  source_cloud->push_back(point1);
  source_cloud->push_back(point2);
  pcl::transformPointCloudWithNormals (*source_cloud, *transformed_cloud, transform_2);
  Eigen::Vector3f vecbefore;
	vecbefore << point1.normal_x,point1.normal_y,point1.normal_z;
	Eigen::Vector3f vecafter;
	vecafter << 1,0,0;

  Eigen::Matrix4f transform2;
  double lalal=get_angle(point1.normal_x,point1.normal_y,point1.normal_z,1,0,0);
  	//【1】求两个向量间的旋转角angle(点积)
  if(lalal>90){
    lalal=180-lalal;
    double angle=lalal*PI/180;
    //【2】求旋转轴（叉积）
    Eigen::Vector3f axis1 = vecbefore.cross(vecafter);
    // std::cout << "求旋转轴： " << axis1 << std::endl;
    Eigen::Vector3f axis2 = vecafter.cross(vecbefore);
    // std::cout << "求旋转轴： " << axis2 << std::endl;
    //【3】求旋转矩阵
    transform_2 = Eigen::Affine3f::Identity();
    // Define a translation of 2.5 meters on the x axis.
    transform_2.translation() << 0, 0, 0;
    // The same rotation matrix as before; theta radians arround Z axis
    transform_2.rotate(Eigen::AngleAxisf(angle, axis2.normalized()));
    // Print the transformation
  //	std::cout << transform_2.matrix() << std::endl;
    transform2=transform_2.matrix();
  }
  else{
    double angle=lalal*PI/180;
    //【2】求旋转轴（叉积）
    Eigen::Vector3f axis1 = vecbefore.cross(vecafter);
    // std::cout << "求旋转轴： " << axis1 << std::endl;
    Eigen::Vector3f axis2 = vecafter.cross(vecbefore);
    // std::cout << "求旋转轴： " << axis2 << std::endl;
    //【3】求旋转矩阵
    transform_2 = Eigen::Affine3f::Identity();
    // Define a translation of 2.5 meters on the x axis.
    transform_2.translation() << 0, 0, 0;
    // The same rotation matrix as before; theta radians arround Z axis
    transform_2.rotate(Eigen::AngleAxisf(angle, axis1.normalized()));
    // Print the transformation
    // std::cout << transform_2.matrix() << std::endl;
    transform2=transform_2.matrix();
    }
    // cout<<endl;
   // cout<<transform2<<endl;
    pcl::PointCloud<pcl::PointNormal>::Ptr transformed_tmp(new pcl::PointCloud<pcl::PointNormal>);

    pcl::transformPointCloudWithNormals(*transformed_cloud, *transformed_tmp, transform2);//.inverse()

    // cout<<"transformed_tmp->clouds[0]"<<transformed_tmp->points[0]<<endl;
    // cout<<"transformed_tmp->clouds[1]"<<transformed_tmp->points[1]<<endl;

    // double a=get_angle(transformed_tmp->points[1].x,transformed_tmp->points[1].y,transformed_tmp->points[1].z,0,1,0);
    double a=get_angle(0,transformed_tmp->points[1].y,transformed_tmp->points[1].z,0,1,0);
    a=a*PI/180;
    return a;
}


void compute_transform_x(pcl::PointNormal point1,Eigen::Matrix4f &Tmg){
  Eigen::Affine3f transform_2 = Eigen::Affine3f::Identity();
  transform_2.translation() << -point1.x,-point1.y,-point1.z;
  // std::cout << transform_2.matrix() << std::endl;
  Eigen::Matrix4f tmg=transform_2.matrix();
  pcl::PointCloud<pcl::PointNormal>::Ptr transformed_cloud (new pcl::PointCloud<pcl::PointNormal> ());
  pcl::PointCloud<pcl::PointNormal>::Ptr source_cloud (new pcl::PointCloud<pcl::PointNormal> ());
  source_cloud->push_back(point1);
  pcl::transformPointCloudWithNormals (*source_cloud, *transformed_cloud, transform_2);
  Eigen::Vector3f vecbefore;
	vecbefore << point1.normal_x,point1.normal_y,point1.normal_z;
	Eigen::Vector3f vecafter;
	vecafter << 1,0,0;

  Eigen::Matrix4f transform2;
  double lalal=get_angle(point1.normal_x,point1.normal_y,point1.normal_z,1,0,0);
  	//【1】求两个向量间的旋转角angle(点积)
  if(lalal>90){
    lalal=180-lalal;
    double angle=lalal*PI/180;
    //【2】求旋转轴（叉积）
    Eigen::Vector3f axis1 = vecbefore.cross(vecafter);
    // std::cout << "求旋转轴： " << axis1 << std::endl;
    Eigen::Vector3f axis2 = vecafter.cross(vecbefore);
    // std::cout << "求旋转轴： " << axis2 << std::endl;
    //【3】求旋转矩阵
    transform_2 = Eigen::Affine3f::Identity();
    // Define a translation of 2.5 meters on the x axis.
    transform_2.translation() <<  0,0,0;
    // The same rotation matrix as before; theta radians arround Z axis
    transform_2.rotate(Eigen::AngleAxisf(angle, axis2.normalized()));
    // Print the transformation
  //	std::cout << transform_2.matrix() << std::endl;
    transform2=transform_2.matrix();
  }
  else{
    double angle=lalal*PI/180;
    //【2】求旋转轴（叉积）
    Eigen::Vector3f axis1 = vecbefore.cross(vecafter);
    // std::cout << "求旋转轴： " << axis1 << std::endl;
    Eigen::Vector3f axis2 = vecafter.cross(vecbefore);
    // std::cout << "求旋转轴： " << axis2 << std::endl;
    //【3】求旋转矩阵
    transform_2 = Eigen::Affine3f::Identity();
    // Define a translation of 2.5 meters on the x axis.
    transform_2.translation() << 0,0,0;
    // The same rotation matrix as before; theta radians arround Z axis
    transform_2.rotate(Eigen::AngleAxisf(angle, axis1.normalized()));
    // Print the transformation
    // std::cout << transform_2.matrix() << std::endl;
    transform2=transform_2.matrix();
    }
    Tmg=transform2*tmg;
    

}

void computeUnitX_Rotation(double angle, Eigen::Matrix4f& Talpha)
{
  const double sx = sin(angle);
  const double cx = cos(angle);

  //Mat(Rx.eye()).copyTo(Rx);
  Talpha(0, 0) = 1, Talpha(1, 0) = 0, Talpha(2, 0) = 0; Talpha(3, 0) = 0;
  Talpha(0, 1) = 0, Talpha(1, 1) = 1, Talpha(2, 1) = 0; Talpha(3, 1) = 0;
  Talpha(0, 2) = 0, Talpha(1, 2) = 0, Talpha(2, 2) = 1; Talpha(3, 2) = 0;
  Talpha(0, 3) = 0, Talpha(1, 3) = 0, Talpha(2, 3) = 0; Talpha(3, 3) = 1;

  Talpha(1, 1) = cx;
  Talpha(1, 2) = -sx;
  Talpha(2, 1) = sx;
  Talpha(2, 2) = cx;
}


//ppf p1 (distance)  p2 (reference points   lianxian)  p3 ( points   lianxian) p4 (nor1 nor2) p5(alpha)
Vector5f computePPFfeature(pcl::PointNormal point1,pcl::PointNormal point2){
  Vector5f ppf;
  ppf << 1.0 , 2.0 , 3.0, 4.0, 5.0;
  ppf[0]=sqrt(pow((point1.x-point2.x),2)+pow((point1.y-point2.y),2)+pow((point1.z-point2.z),2));
  double angle1=get_angle(point1.normal_x,point1.normal_y,point1.normal_z,(point1.x-point2.x),(point1.y-point2.y),(point1.z-point2.z));
  double angle2=get_angle((point1.x-point2.x),(point1.y-point2.y),(point1.z-point2.z),point2.normal_x,point2.normal_y,point2.normal_z);
  double angle3=get_angle(point1.normal_x,point1.normal_y,point1.normal_z,point2.normal_x,point2.normal_y,point2.normal_z);
  //转换为弧度
  ppf[1]=angle1*PI/180;
  ppf[2]=angle2*PI/180;
  ppf[3]=angle3*PI/180;
  ppf[4]=compute_a(point1,point2);
  return ppf;

}


void  computePPFfeature1(pcl::PointNormal point1,pcl::PointNormal point2,  Vector5f &ppf){
  ppf[0]=sqrt(pow((point1.x-point2.x),2)+pow((point1.y-point2.y),2)+pow((point1.z-point2.z),2));
  double angle1=get_angle(point1.normal_x,point1.normal_y,point1.normal_z,(point1.x-point2.x),(point1.y-point2.y),(point1.z-point2.z));
  double angle2=get_angle((point1.x-point2.x),(point1.y-point2.y),(point1.z-point2.z),point2.normal_x,point2.normal_y,point2.normal_z);
  double angle3=get_angle(point1.normal_x,point1.normal_y,point1.normal_z,point2.normal_x,point2.normal_y,point2.normal_z);
  ppf[1]=angle1*PI/180;;
  ppf[2]=angle2*PI/180;;
  ppf[3]=angle3*PI/180;;
  ppf[4]=compute_a(point1,point2);
}
float distance_point(pcl::PointXYZ point1,pcl::PointXYZ point2){
  float distance =sqrt(pow((point1.x-point2.x),2)+pow((point1.y-point2.y),2)+pow((point1.z-point2.z),2));
  return distance;

}


void FPS(PointCloud<PointXYZ>::Ptr& cloud_in,PointCloud<PointXYZ>::Ptr& cloud_out,int number,int index_first=0){
  std::vector<int> index_point;
  index_point.push_back(index_first);
  int count=0;
  cloud_out->push_back(cloud_in->points[ index_first]);
  // cout<<cloud_in->points[ index_first].x<<"   "<<cloud_in->points[ index_first].y<<cloud_in->points[ index_first].z<<endl;
  while(count<=number){

    float max_distance=-1;
    int index_max=0;
    int index_min=0;      
    for(int i=0;i<cloud_in->points.size();i++){
        float min_distance=999999;
        float distance=0;
        for(int j=0;j<index_point.size();j++){

          vector<int>::iterator ret;
          ret=std::find(index_point.begin(),index_point.end(),i);
          if(ret==index_point.end()){
            distance=distance_point(cloud_in->points[i],cloud_in->points[index_point[j]]);
            if(distance <min_distance){
              min_distance=distance;
              index_min=i;
            }
        }
      }
      if(min_distance > max_distance&&min_distance!=999999){
          index_max=index_min;
          max_distance=min_distance;
        }
    }  
    index_point.push_back(index_max);
    cloud_out->push_back(cloud_in->points[ index_max ]); 
    count++;
    if(count%5 == 0)
    {
        std::cout << "sampleing_trained: " << std::floor(100*count/number) << "%" << std::endl;
    }
  }
}

//leaf  search_radius
PointCloud<PointNormal>::Ptr FPS_CalculateNormals (PointCloud<PointXYZ>::Ptr& cloud,int number)
{
  pcl::PointCloud<pcl::PointXYZ>::Ptr out (new pcl::PointCloud<pcl::PointXYZ> ());
  // pcl::VoxelGrid<pcl::PointXYZ> down_filter;
  // down_filter.setLeafSize(leaf, leaf, leaf);
  // down_filter.setInputCloud(cloud);
  // down_filter.filter(*out);
  // cout<<cloud->points.size()<<endl;
 
  FPS(cloud,out,number-2);

  PointCloud<Normal>::Ptr cloud_normals (new PointCloud<Normal> ());
  NormalEstimation<PointXYZ, Normal> normal_estimation_filter;
  normal_estimation_filter.setInputCloud (out);
  search::KdTree<PointXYZ>::Ptr search_tree (new search::KdTree<PointXYZ>);
  normal_estimation_filter.setSearchMethod (search_tree);
  normal_estimation_filter.setRadiusSearch (search_radius);
  normal_estimation_filter.compute (*cloud_normals);

  PointCloud<PointNormal>::Ptr cloud_with_normals (new PointCloud<PointNormal> ());
  concatenateFields (*out, *cloud_normals, *cloud_with_normals);
  return cloud_with_normals;
}

PointCloud<PointNormal>::Ptr Voxel_CalculateNormals (PointCloud<PointXYZ>::Ptr& cloud,float leaf)
{
  pcl::PointCloud<pcl::PointXYZ>::Ptr out (new pcl::PointCloud<pcl::PointXYZ> ());
  pcl::VoxelGrid<pcl::PointXYZ> down_filter;
  down_filter.setLeafSize(leaf, leaf, leaf);
  down_filter.setInputCloud(cloud);
  down_filter.filter(*out);
  cout<<cloud->points.size()<<endl;
 
  // FPS(cloud,out,29);

  PointCloud<Normal>::Ptr cloud_normals (new PointCloud<Normal> ());
  NormalEstimation<PointXYZ, Normal> normal_estimation_filter;
  normal_estimation_filter.setInputCloud (out);
  search::KdTree<PointXYZ>::Ptr search_tree (new search::KdTree<PointXYZ>);
  normal_estimation_filter.setSearchMethod (search_tree);
  normal_estimation_filter.setRadiusSearch (search_radius);
  normal_estimation_filter.compute (*cloud_normals);

  PointCloud<PointNormal>::Ptr cloud_with_normals (new PointCloud<PointNormal> ());
  concatenateFields (*out, *cloud_normals, *cloud_with_normals);
  return cloud_with_normals;
}


static KeyType hashPPFM2E(Vector5f& f, const double AngleStep, const double DistanceStep)
{
    KeyType *key = new KeyType[4];

    key[0] = (int)(f[0] / DistanceStep);
    key[1] = (int)(f[1] / AngleStep);
    key[2] = (int)(f[2] / AngleStep);
    key[3] = (int)(f[3] / AngleStep);
    //  cout<<" "<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<f[3]<<"  ";
    // cout<<key[0]<<" "<<key[1]<<" "<<key[2]<<" "<<key[3]<<"  ";
    KeyType hashKey = 0;
    murmurHash(key, 4*sizeof(int), 42, &hashKey);
    return hashKey;
}
void inverse_Matrix44(Eigen::Matrix4f &Tsg,Eigen::Matrix4f &TsgInv){

    Eigen::Matrix3f R_Tsg;
    Eigen::Vector3f t_Tsg;
    R_Tsg=Tsg.block(0, 0, 3, 3) ;//Index startRow, Index startCol, Index blockRows, Index blockCols
    t_Tsg=Tsg.block(0, 3, 3, 1) ;//Index startRow, Index startCol, Index blockRows, Index blockCols
    R_Tsg = R_Tsg.transpose().eval();
    t_Tsg=t_Tsg*(-1);

    Eigen::Matrix4f TsgInv_R;
    Eigen::Matrix4f TsgInv_t;
    TsgInv_R.setIdentity();
    TsgInv_t.setIdentity();
    TsgInv_R.block(0, 0, 3, 3) = R_Tsg;
    TsgInv_t.block(0, 3, 3, 1) = t_Tsg;

    TsgInv=TsgInv_R*TsgInv_t;

}
void SelectSort(std::vector<Pose_3DM2E> &a,std::vector<Pose_3DM2E> &b,int len){
  Pose_3DM2E temp;
  bool flag[50000] = {false};
  
  
  for(int i=0;i<len;i++){
    int max=-1;
    int maxvaule=0;
   
    for(int j=0;j<len;j++){
  
      if(a[j].numVotes > maxvaule && flag[j]==false){
       max=j;
       maxvaule=a[j].numVotes;
     }

    }
    flag[max]=true;
    if(max!=i){
      
      b.push_back(a[max]);
  
    }else{
       
      b.push_back(a[i]);

    }
    // for(int i=0;i<b.size();i++){
    //   cout<<b[i].numVotes<<"  ";
    // }
  }

}
static bool pose_3DCompareM2E(const Pose_3DM2E& a, const Pose_3DM2E& b)
{
    // cout<<"len:"<<sizeof(a)<<endl;
    // long a_number=a.numVotes;
    // long b_number=b.numVotes;
    // cout<<" a_number: "<<a_number<<"  "<<" b_number: "<<b_number<<endl;
    // bool result=( a.numVotes > b.numVotes);
    // cout<<"result: "<<result<<endl;
 
    return ( a.numVotes > b.numVotes);
}
// static bool pose_3DCompareM2E(Pose_3DM2E  a, Pose_3DM2E b)
// {
//   cout<<"len:"<<sizeof(a)<<endl;
//     long a_number=a.numVotes;
//     long b_number=b.numVotes;
//     cout<<" a_number: "<<a_number<<"  "<<" b_number: "<<b_number<<endl;
//     bool result=( a.numVotes > b.numVotes);
//     cout<<"result: "<<result<<endl;
//     return ( a.numVotes > b.numVotes);
// }


static int Pose_3DClustersM2E(const ClusterM2E& a, const ClusterM2E& b)
{
  return ( a.accu_votes > b.accu_votes);
}


bool matchPoseM2E(const Pose_3DM2E& sourcePose, const Pose_3DM2E& targetPose,double rotation_threshold,double position_threshold)
{
    Eigen::Vector3f dv = targetPose.t - sourcePose.t;
    double dNorm = dv.norm();
    const double phi = fabs ( sourcePose.angle - targetPose.angle );

    return (phi < rotation_threshold && dNorm < position_threshold);
}
void matrix2angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle)
{
	double ax, ay, az;
	if (result_trans(2, 0) == 1 || result_trans(2, 0) == -1)
	{
		az = 0;
		double dlta;
		dlta = atan2(result_trans(0, 1), result_trans(0, 2));
		if (result_trans(2, 0) == -1)
		{
			ay = M_PI / 2;
			ax = az + dlta;
		}
		else
		{
			ay = -M_PI / 2;
			ax = -az + dlta;
		}
	}
	else
	{
		ay = -asin(result_trans(2, 0));
		ax = atan2(result_trans(2, 1) / cos(ay), result_trans(2, 2) / cos(ay));
		az = atan2(result_trans(1, 0) / cos(ay), result_trans(0, 0) / cos(ay));
	}
	result_angle << ax, ay, az;

	cout << "x轴旋转角度：" << ax << endl;
	cout << "y轴旋转角度：" << ay << endl;
	cout << "z轴旋转角度：" << az << endl;
}




void  clusterPosesM2E(std::vector<Pose_3DM2E>& poseList1,std::vector<Pose_3DM2E>& poseList, int numPoses, std::vector<Pose_3DM2E> &finalPoses,double rotation_threshold,double position_threshold){
    std::vector<ClusterM2E> poseClusters;
    // std::sort(poseList.begin(), poseList.end(), pose_3DCompareM2E);
  

    SelectSort(poseList1,poseList,poseList1.size());
    cout<<"begin:-----"<<poseList.size()<<endl;
   
    // poseList=poseList1;
    for(int i=0;i<poseList.size();i++){
      cout<<poseList[i].numVotes<<"  ";
    }
    for (int i=0; i<numPoses; i++)
    {
        Pose_3DM2E pose = poseList[i];
        bool assigned = false;

        // search all clusters
        for (size_t j=0; j<poseClusters.size() && !assigned; j++)
        {
            //const Pose_3D poseCenter = poseClusters[j]->poseList[0];
            const Pose_3DM2E poseCenter = poseClusters[j].poses[0];
            if (matchPoseM2E(pose, poseCenter,rotation_threshold,position_threshold))
            {
                //poseClusters[j]->addPose(pose);
                poseClusters[j].poses.push_back(pose);
                poseClusters[j].accu_votes+=pose.numVotes;
                cout<<j<<"   accu_votes:  "<<poseClusters[j].accu_votes<<endl;
                assigned = true;
                break;
            }
        }

        if (!assigned)
        {
            ClusterM2E poseCluster;
            poseCluster.poses.push_back(pose);
            poseCluster.accu_votes+=pose.numVotes;
            poseClusters.push_back(poseCluster);
        }
    }
    cout<<"okkk"<<endl;
    std::sort(poseClusters.begin(), poseClusters.end(), Pose_3DClustersM2E);
    for(int i=0;i<poseClusters.size();i++){
      cout<<poseClusters[i].accu_votes<<"   ";
    }
    cout<<"okkk"<<endl;
    cout<<" Pose clusters:  "<<poseClusters.size()<<endl;
    if (use_weighted_avg)
    {

        // uses weighting by the number of votes
        for (int i=0; i<(poseClusters.size()); i++)
        {
            // We could only average the quaternions. So I will make use of them here
            Eigen::Vector4f qAvg(0,0,0,0);
            Eigen::Vector3f tAvg(0,0,0);

            // Perform the final averaging
            ClusterM2E curCluster = poseClusters[i];
        
            std::vector<Pose_3DM2E> curPoses = curCluster.poses;
            int curSize = (int)curPoses.size();
            size_t numTotalVotes = 0;

            for (int j=0; j<curSize; j++)
                numTotalVotes += curPoses[j].numVotes;

            double wSum=0;

            for (int j=0; j<curSize; j++)
            {
                const double w = (double)curPoses[j].numVotes / (double)numTotalVotes;
                if(i==0){
                 model_Index0[temp0]=curPoses[j].modelIndex;
                 temp0++;
                }
                if(i==1){
                 model_Index1[temp1]=curPoses[j].modelIndex;
                 temp1++;
                }
                if(i==2){
                 model_Index2[temp2]=curPoses[j].modelIndex;
                 temp2++;
                }
                curPoses[j].q.normalize();
                // curPoses[j].t.normalize();

                qAvg += w * curPoses[j].q;
                tAvg += w * curPoses[j].t;
                wSum += w;
            }
            

            tAvg *= 1.0 / wSum;
            qAvg *= 1.0 / wSum;
            qAvg.normalize();
            
          // cout<<tAvg<<endl;
            Pose_3DM2E pose_tmp;
            pose_tmp.q=qAvg;
            pose_tmp.t=tAvg;
            pose_tmp.numVotes=curCluster.accu_votes;
            pose_tmp.computepose();
            // cout<<pose_tmp.q<<"  "<<pose_tmp.q<<"  pose: "<<pose_tmp.pose<<endl;
            // cout<<"i:"<<i<<endl;
            finalPoses.push_back(pose_tmp);

            

            
        }
    }
    else
    {

        for (int i=0; i<static_cast<int>(poseClusters.size()); i++)
        {
            // We could only average the quaternions. So I will make use of them here
            Eigen::Vector4f qAvg(0,0,0,0);
            Eigen::Vector3f tAvg(0,0,0);

            // Perform the final averaging
            ClusterM2E curCluster = poseClusters[i];
            std::vector<Pose_3DM2E> curPoses = curCluster.poses;
            const int curSize = (int)curPoses.size();

            for (int j=0; j<curSize; j++)
            {   

                curPoses[j].q.normalize();
                // curPoses[j].t.normalize();
                qAvg += curPoses[j].q;
                tAvg += curPoses[j].t;
            }

            tAvg *= 1.0 / curSize;
            qAvg *= 1.0 / curSize;
            qAvg.normalize();

            //curPoses[0]->updatePoseQuat(qAvg, tAvg);
           
            Pose_3DM2E pose_tmp;
            pose_tmp.q=qAvg;
            pose_tmp.t=tAvg;
            pose_tmp.computepose();
            pose_tmp.numVotes=curCluster.accu_votes;
            finalPoses.push_back(pose_tmp);
        }
    }



   




}

//  const double rotation_threshold=0.21;//12
// const double position_threshold=0.01;//1cm


int main (int argc, char** argv)
{
  pcl::PointCloud<pcl::PointXYZ>::Ptr source_cloud (new pcl::PointCloud<pcl::PointXYZ> ());
  pcl::PointCloud<pcl::PointXYZ>::Ptr model_cloud (new pcl::PointCloud<pcl::PointXYZ> ());
  PointCloud<PointNormal>::Ptr source_with_normals (new PointCloud<PointNormal> ());
  PointCloud<PointNormal>::Ptr model_with_normals (new PointCloud<PointNormal> ());
  double rotation_threshold=atof(argv[2]);//12
  double position_threshold=atof(argv[3]);//1cm
  string scene_path="../data/"+string(argv[1])+".ply";
  pcl::io::loadPLYFile ("../data/bunny.ply", *model_cloud);
  pcl::io::loadPLYFile (scene_path, *source_cloud);
  ofstream oFile;
  string file_ods=string(argv[1])+"_"+string(argv[4])+".ods";
  oFile.open(file_ods,ios::out|ios::trunc);
  oFile<<"rx"<<","<<"ry"<<","<<"rz"<<","<<"tx"<<","<<"ty"<<","<<"tz"<<","<<"bili"<<","<<"score"<<","<<"angle_score"<<","<<"distance_score"<<","<<"sample_points"<<","<<"sample_method"<<","<<"PoseList_number"
    <<","<<"PoseList_max"<<","<<"PoseList_avg"<<","<<"clusterPoses_num"<<","<<"clusterPoses_max_votes"<<","<<"clusterPoses_avg_votes"<<","<<"angle_threshold"<<","
    <<"distance_threshold"<<","<<"key_number"<<","<<"key_volumn"<<","<<"time"<<endl;
  bool is_FPS=false;
  int cloud_number[5]={25,100,250,300,400};
  float leaf[5]={0.04,0.027,0.017,0.012,0.007};
  int sum_number=0;
  double k1[6]={60,50,90,0.05,0.07,0.09};
  double k2[6]={180,0,0,0.12,0.01,0.06};
  double k3[6]={28,23,32,0.2,0.4,0.6};
  double k4[6]={10,12,16,-0.3,0.01,0.02};
  double *k;
  if(string(argv[1])=="bunny_tran1"){
    k=k1;
  }else if(string(argv[1])=="bunny_tran2"){
    k=k2;
  }else if(string(argv[1])=="bunny_tran3"){
    k=k3;
  }else{
    k=k4;
  }
  

  while(sum_number<10){
    clock_t start = clock();
    if(sum_number<5){
      is_FPS=false;

    }else{
      is_FPS=true;

    }
  if(is_FPS){
    model_with_normals=FPS_CalculateNormals(model_cloud,cloud_number[sum_number-5]);
    source_with_normals=FPS_CalculateNormals(source_cloud,cloud_number[sum_number-5]);

  }else{
    model_with_normals=Voxel_CalculateNormals(model_cloud,leaf[sum_number]);
    source_with_normals=Voxel_CalculateNormals(source_cloud,leaf[sum_number]);

  }

  
  string model_with_normals_name="../data/"+string(argv[1])+"/"+string(argv[4])+"/"+to_string(sum_number)+"_bunny_with_normals"+".ply";
  string source_with_normals_name="../data/"+string(argv[1])+"/"+string(argv[4])+"/"+to_string(sum_number)+"bunny_scene_with_normals"+".ply";//
  pcl::io::savePLYFile (model_with_normals_name, *model_with_normals);
  pcl::io::savePLYFile (source_with_normals_name, *source_with_normals);
  cout<<"model_with_normals->points.size()"<<model_with_normals->points.size()<<endl;

  PCL_INFO ("Training models ...\n");
  int model_points=model_with_normals->points.size();
  int size=model_points*model_points;
  hashtable_int* hashTable = hashtableCreate(size, NULL);
  hash_nodesM2E = (THash*)calloc(size, sizeof(THash));
  int count=0;
   for(int i=0;i<model_with_normals->points.size();i++){
    for(int j=0;j<model_with_normals->points.size();j++){
      if(j!=i){
        Vector5f ppf_temp=computePPFfeature(model_with_normals->points[i],model_with_normals->points[j]);
        KeyType hashValue = hashPPFM2E(ppf_temp, angle_step_radiansM2E, distanceStep);
        uint ppfInd = i*model_points+j;
        THash* hashNode = &hash_nodesM2E[i*model_points+j];
        hashNode->id = hashValue;
        hashNode->i = i;
        hashNode->ppfInd = ppfInd;
        hashNode->angle = (float)ppf_temp[4];
        hashtableInsertHashed(hashTable, hashValue, (void*)hashNode);
        count++;

      }
    }
     if(i%10 == 0)
        {
            std::cout << "trained: " << std::floor(100*i/model_points) << "%" << std::endl;
        }
  }
  cout<<"count(node insert numbers) :  "<<count<<endl;
  hash_tableM2E=hashTable;
  Vector5f ppf_temp;
  uint* accumulator = (uint*)calloc(numAngles*model_points, sizeof(uint));
  int scene_points=source_with_normals->points.size();

  int scene_sample_stepM2E = (int)(1.0/relativeSceneSampleStep);//10个点来寻找
  std::vector<Pose_3DM2E> poseList;
  int sceneSamplingStep = scene_sample_stepM2E;
  poseList.reserve((scene_points/sceneSamplingStep)+4);//多少个pose 相关

  for(int i=0;i<scene_points;i+=sceneSamplingStep){
    uint refIndMax = 0, alphaIndMax = 0;
    uint maxVotes = 0;
    for(int j=0;j<scene_points;j++){
      if(j!=i){
        
        computePPFfeature1(source_with_normals->points[i],source_with_normals->points[j],ppf_temp);
        KeyType hashValue = hashPPFM2E(ppf_temp, angle_step_radiansM2E, distanceStep);
        hashnode_i* node = hashtableGetBucketHashed(hash_tableM2E, (hashValue));
       
        while (node)
        {
            THash* tData = (THash*) node->data;
            int corrI = (int)tData->i;
            int ppfInd = (int)tData->ppfInd;
            double alpha_model = (float)tData->angle;
            double alpha = alpha_model - (float)ppf_temp[4];
            int alpha_index = (int)(numAngles*(alpha + 2*M_PI) / (4*M_PI));
            uint accIndex = corrI * numAngles + alpha_index;
            accumulator[accIndex]++;
            node = node->next;
            
        }

      }
     
    }
   
    for (uint k = 0; k < model_points; k++)
    {
        for (int j = 0; j < numAngles; j++)
        {

          const uint accInd = k*numAngles + j;
          const uint accVal = accumulator[ accInd ];
          if (accVal > maxVotes)
          {
              maxVotes = accVal;
              refIndMax = k;
              alphaIndMax = j;
          }
          accumulator[accInd ] = 0;

        }
    }
    cout<<" maxVotes: "<<maxVotes<<" "<<" refIndMax:(model) "<<refIndMax<<" "<<" alphaIndMax: "<<alphaIndMax<<"scene_ i: "<<i<<endl;


    // invert Tsg : Luckily rotation is orthogonal:
    // Inverse = Transpose.
    // We are not required to invert
    
    //local coordinates(mr ,alpha)
    pcl::PointNormal point1;
    point1=model_with_normals->points[refIndMax];//model
    pcl::PointNormal point2;
    point2=source_with_normals ->points[i];//scene

    int alpha_index = alphaIndMax;
    double alpha = (alpha_index*(4*M_PI))/numAngles-2*M_PI;

    Eigen::Matrix4f Tmg;
    Eigen::Matrix4f Tsg;
    Eigen::Matrix4f TsgInv;
    Eigen::Matrix4f Talpha;

    compute_transform_x(point1,Tmg);
    compute_transform_x(point2,Tsg);
    inverse_Matrix44(Tsg,TsgInv);
    computeUnitX_Rotation(alpha,Talpha);


    Eigen::Matrix4f rawPose = TsgInv * (Talpha * Tmg);
    
    // string plyname="transformed_cloud"+to_string(i)+".ply";
    // PointCloud<PointNormal>::Ptr transformed_cloud (new PointCloud<PointNormal> ());
    // pcl::transformPointCloudWithNormals (*model_with_normals, *transformed_cloud, rawPose);
     // string plyname="transformed_cloud"+to_string(i)+".ply";
    // pcl::io::savePLYFile (plyname, *transformed_cloud);

    


  //这里和opencv的实现就不太一样了，根本没有用四元数之类的;
  //看论文，看两个的实现方法，看opencv pcl的实现方法选择哪一个才是真正正确的;
    Pose_3DM2E pose;
    pose.alpha = alpha;
    pose.modelIndex = refIndMax;
    pose.numVotes = maxVotes;
    pose.pose = rawPose;
    pose.computeqandt();
    pose.computeangle();
    // cout<<"  q;t;angle:  "<<pose.q<<"  "<<pose.t<<"  "<<pose.angle<<endl;
    poseList.push_back(pose);
 
    


  
  }

 
  int numPosesAdded = scene_points/sceneSamplingStep;
  std::cout << "numPoses = " << numPosesAdded << std::endl;
  std::vector<Pose_3DM2E> results;
  double avg_votes=0;
  double max_votes=0;
 
  
  cout<<"poseList:  "<<poseList.size()<<endl;
  for(int i=0;i<poseList.size();i++){
    cout<<poseList[i].numVotes<<"   ";

    // cout<<poseList[i]<<endl;
    
    if(poseList[i].numVotes> max_votes){
      max_votes=poseList[i].numVotes;
    }
    avg_votes+=poseList[i].numVotes;
  }
  avg_votes/=poseList.size();
  cout<<"avg_votes  "<<avg_votes<<endl;
  cout<<"max_votes  "<<max_votes<<endl;
  cout<<endl;
  cout<<endl;

  cout<<"okk"<<endl;
   std::vector<Pose_3DM2E> poseList_result;
  clusterPosesM2E(poseList, poseList_result,numPosesAdded, results,rotation_threshold,position_threshold);
  // sort(results.begin(),results.end(),pose_3DCompareM2E);


  cout<<"okk"<<endl;
  cout<<"final_pose-count:"<<results.size()<<endl;
  double clusterPoses_avg_votes=0;
  double clusterPoses_max_votes=0;
  for(int i=0;i<results.size();i++){
    cout<<" i: "<<results[i].numVotes<<"    ";
     if(results[i].numVotes> max_votes){
      clusterPoses_max_votes=results[i].numVotes;
    }
    clusterPoses_avg_votes+=results[i].numVotes;
  }
   clusterPoses_avg_votes/=results.size();
   cout<<"avg_votes_result  "<<clusterPoses_avg_votes<<endl;
  cout<<"max_votes _result "<<clusterPoses_max_votes<<endl;
  cout<<endl;
  cout<<endl;
  pcl::PointXYZ min_p, max_p;
	pcl::getMinMax3D(*model_cloud, min_p, max_p);
  double volume=(max_p.x-min_p.x)*(max_p.y-min_p.y)*(max_p.z-min_p.z);
  double bili=(pow(volume*3/4*M_PI,1.0 / 3))*(2*M_PI/360.0);
  
  double distance_score=0;
  double angle_score=0;

 
  clock_t sac_time = clock();
  double time=(double)(sac_time - start) / (double)CLOCKS_PER_SEC;
	// cout << "sac time: " << (double)(sac_time - start) / (double)CLOCKS_PER_SEC << " s" << endl;
    
   

  cout << "Visualization_PointCloud...\n";

  // for(int i=0;i<poseList_result.size();i++){
  //   pcl::PointCloud<pcl::PointXYZ>::Ptr transform_model2scene(new pcl::PointCloud<pcl::PointXYZ> ());
  //   Eigen::Matrix4f tmp=poseList_result[i].pose;
  //   pcl::transformPointCloud(*model_cloud, *transform_model2scene, tmp);
  //   string plyname="transformed"+to_string(i)+".ply";
  //   pcl::io::savePLYFile (plyname, *transform_model2scene);

  // }
  pcl::PointCloud<pcl::PointXYZ>::Ptr transform_model2scene(new pcl::PointCloud<pcl::PointXYZ> ());
  Eigen::Matrix4f tmp=results[0].pose; 
  // Eigen::Matrix4f tmp=poseList_result[16].pose;
  

  // cout<<"angles0  "<<results[0].angle<<endl;
  // cout<<"q0    "<<results[0].q<<endl;
  // cout<<"t0    "<<results[0].t<<endl;
  Eigen::Vector3f result_angle;
  cout<<"tmp:"<<tmp<<endl;
  pcl::transformPointCloud(*model_cloud, *transform_model2scene, tmp);
  // pcl::io::savePLYFile ("match_save3.ply", *transform_model2scene);
  // pcl::io::savePLYFile ("match_model3.ply", *model_cloud);
   matrix2angle(tmp, result_angle);
  cout<<"6D POSE"<<endl;
  cout<<result_angle[0]*180/M_PI<<"   "<<result_angle[1]*180/M_PI<<"   "<<result_angle[2]*180/M_PI<<endl;
  cout<<tmp(0,3)<<"   "<<tmp(1,3)<<"   "<<tmp(2,3)<<endl;
  angle_score=0;
  angle_score+=abs(result_angle[0]*180/M_PI-k[0]);
  angle_score+=abs(result_angle[1]*180/M_PI-k[1]);
  angle_score+=abs(result_angle[2]*180/M_PI-k[2]);
  distance_score=0;
  distance_score+=abs(tmp(0,3)-k[3]);
  distance_score+=abs(tmp(1,3)-k[4]);
  distance_score+=abs(tmp(2,3)-k[5]);

    pcl::PointCloud<pcl::PointXYZ>::Ptr points_im(new pcl::PointCloud<pcl::PointXYZ> ());
    for(int i=0;i<temp0;i++){
      points_im->push_back(transform_model2scene->points[model_Index0[i]]);

    }
  oFile<<to_string(result_angle[0]*180/M_PI)<<","<<to_string(result_angle[1]*180/M_PI)<<","<<to_string(result_angle[2]*180/M_PI)<<","<<
to_string(tmp(0,3))<<","<<to_string(tmp(1,3))<<","<<to_string(tmp(2,3))<<","<<bili<<","<<to_string((angle_score*bili)+(distance_score))<<","<<
angle_score<<","<<distance_score<<","<<to_string(model_with_normals->points.size())<<","<<((is_FPS ==true) ? "FPS" : "Voxel")<<","<<
to_string(poseList_result.size())<<","<<to_string(max_votes)<<","<<to_string(avg_votes)<<","<<to_string(results.size())<<","<<to_string(clusterPoses_max_votes)
<<","<<to_string(clusterPoses_avg_votes)<<","<<to_string(rotation_threshold)<<","<<to_string(position_threshold)<<","<<to_string(points_im->points.size())<<
","<<to_string(get_volume(points_im))<<","<<time<<endl;


// oFile<<"rx"<<","<<"ry"<<","<<"rz"<<","<<"tx"<<","<<"ty"<<","<<"tz"<<","<<"bili"<<","<<"score"<<","<<"angle_score"<<","<<"distance_score"<<","<<"sample_points"<<","<<"sample_method"<<","<<"PoseList_number"
//     <<","<<"PoseList_max"<<","<<"PoseList_avg"<<","<<"clusterPoses_num"<<","<<"clusterPoses_max_votes"<<","<<"clusterPoses_avg_votes"<<","<<"angle_threshold"<<","
//     <<"distance_threshold"<<endl;
  
    

  //=================1=================================================
  pcl::PointCloud<pcl::PointXYZ>::Ptr transform_model2scene1(new pcl::PointCloud<pcl::PointXYZ> ());
   Eigen::Matrix4f  tmp1=results[1].pose;
  // Eigen::Matrix4f tmp1=poseList_result[32].pose;
    cout<<"tmp1:"<<tmp1<<endl;
  pcl::transformPointCloud(*model_cloud, *transform_model2scene1, tmp1);
  matrix2angle(tmp1, result_angle);
  cout<<"6D POSE"<<endl;
  cout<<result_angle[0]*180/M_PI<<"   "<<result_angle[1]*180/M_PI<<"   "<<result_angle[2]*180/M_PI<<endl;
  cout<<tmp1(0,3)<<"   "<<tmp1(1,3)<<"   "<<tmp1(2,3)<<endl;
  angle_score=0;
  angle_score+=abs(result_angle[0]*180/M_PI-k[0]);
  angle_score+=abs(result_angle[1]*180/M_PI-k[1]);
  angle_score+=abs(result_angle[2]*180/M_PI-k[2]);
  distance_score=0;
  distance_score+=abs(tmp1(0,3)-k[3]);
  distance_score+=abs(tmp1(1,3)-k[4]);
  distance_score+=abs(tmp1(2,3)-k[5]);
     pcl::PointCloud<pcl::PointXYZ>::Ptr points_im1(new pcl::PointCloud<pcl::PointXYZ> ());
    for(int i=0;i<temp1;i++){
      points_im1->push_back(transform_model2scene1->points[model_Index1[i]]);

    }
  
   oFile<<to_string(result_angle[0]*180/M_PI)<<","<<to_string(result_angle[1]*180/M_PI)<<","<<to_string(result_angle[2]*180/M_PI)<<","<<
to_string(tmp1(0,3))<<","<<to_string(tmp1(1,3))<<","<<to_string(tmp1(2,3))<<","<<bili<<","<<to_string((angle_score*bili)+(distance_score))<<","<<
angle_score<<","<<distance_score<<","<<to_string(model_with_normals->points.size())<<","<<((is_FPS ==true) ? "FPS" : "Voxel")<<","<<
to_string(poseList_result.size())<<","<<to_string(max_votes)<<","<<to_string(avg_votes)<<","<<to_string(results.size())<<","<<to_string(clusterPoses_max_votes)
<<","<<to_string(clusterPoses_avg_votes)<<","<<to_string(rotation_threshold)<<","<<to_string(position_threshold)<<","<<to_string(points_im1->points.size())<<
","<<to_string(get_volume(points_im1))<<","<<time<<endl;

   
//==================2==============================================
    pcl::PointCloud<pcl::PointXYZ>::Ptr transform_model2scene2(new pcl::PointCloud<pcl::PointXYZ> ());
   Eigen::Matrix4f  tmp2=results[2].pose;
  // Eigen::Matrix4f tmp2=poseList_result[3].pose;
     cout<<"tmp2:"<<tmp2<<endl;
  pcl::transformPointCloud(*model_cloud, *transform_model2scene2, tmp2);
 matrix2angle(tmp2, result_angle);
  cout<<"6D POSE"<<endl;
  cout<<result_angle[0]*180/M_PI<<"   "<<result_angle[1]*180/M_PI<<"   "<<result_angle[2]*180/M_PI<<endl;
  cout<<tmp2(0,3)<<"   "<<tmp2(1,3)<<"   "<<tmp2(2,3)<<endl;
  angle_score=0;
  angle_score+=abs(result_angle[0]*180/M_PI-k[0]);
  angle_score+=abs(result_angle[1]*180/M_PI-k[1]);
  angle_score+=abs(result_angle[2]*180/M_PI-k[2]);
  distance_score=0;
  distance_score+=abs(tmp2(0,3)-k[3]);
  distance_score+=abs(tmp2(1,3)-k[4]);
  distance_score+=abs(tmp2(2,3)-k[5]);

       pcl::PointCloud<pcl::PointXYZ>::Ptr points_im2(new pcl::PointCloud<pcl::PointXYZ> ());
    for(int i=0;i<temp2;i++){
      points_im2->push_back(transform_model2scene2->points[model_Index2[i]]);

    }
    oFile<<to_string(result_angle[0]*180/M_PI)<<","<<to_string(result_angle[1]*180/M_PI)<<","<<to_string(result_angle[2]*180/M_PI)<<","<<
to_string(tmp2(0,3))<<","<<to_string(tmp2(1,3))<<","<<to_string(tmp2(2,3))<<","<<bili<<","<<to_string((angle_score*bili)+(distance_score))<<","<<
angle_score<<","<<distance_score<<","<<to_string(model_with_normals->points.size())<<","<<((is_FPS ==true) ? "FPS" : "Voxel")<<","<<
to_string(poseList_result.size())<<","<<to_string(max_votes)<<","<<to_string(avg_votes)<<","<<to_string(results.size())<<","<<to_string(clusterPoses_max_votes)
<<","<<to_string(clusterPoses_avg_votes)<<","<<to_string(rotation_threshold)<<","<<to_string(position_threshold)<<","<<to_string(points_im2->points.size())<<
","<<to_string(get_volume(points_im2))<<","<<time<<endl;
  
  string path_transform_model2scene="../data/"+string(argv[1])+"/"+string(argv[4])+"/"+to_string(sum_number)+"_transform_model2scene"+".ply";
  string path_transform_model2scene1="../data/"+string(argv[1])+"/"+string(argv[4])+"/"+to_string(sum_number)+"_transform_model2scene1"+".ply";
  string path_transform_model2scene2="../data/"+string(argv[1])+"/"+string(argv[4])+"/"+to_string(sum_number)+"_transform_model2scene2"+".ply";
  string path_points_im="../data/"+string(argv[1])+"/"+string(argv[4])+"/"+to_string(sum_number)+"_points_im"+".ply";
  string path_points_im1="../data/"+string(argv[1])+"/"+string(argv[4])+"/"+to_string(sum_number)+"_points_im1"+".ply";
  string path_points_im2="../data/"+string(argv[1])+"/"+string(argv[4])+"/"+to_string(sum_number)+"_points_im2"+".ply";
  pcl::io::savePLYFile (path_transform_model2scene, *transform_model2scene);
  pcl::io::savePLYFile (path_transform_model2scene1, *transform_model2scene1);
  pcl::io::savePLYFile (path_transform_model2scene2, *transform_model2scene2);
  pcl::io::savePLYFile (path_points_im, *points_im);
  pcl::io::savePLYFile (path_points_im1, *points_im1);
  pcl::io::savePLYFile (path_points_im2, *points_im2);
  cout<<path_transform_model2scene<<endl;
  cout<<path_transform_model2scene1<<endl;
  cout<<path_transform_model2scene2<<endl;
  cout<<path_points_im<<endl;
  cout<<path_points_im1<<endl;
  cout<<path_points_im2<<endl;


 

	// boost::shared_ptr<pcl::visualization::PCLVisualizer > viewer(new pcl::visualization::PCLVisualizer("Ransac"));
	// viewer->setBackgroundColor(0, 0, 0);
	// //创建窗口
	// int vp;
	// viewer->createViewPort(0.0, 0.0, 1.0, 1.0, vp);
	// //设置点云颜色
	// pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> source_color(source_cloud, 255, 255, 255);
	// pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> target_color(transform_model2scene, 255, 0, 0);
  // pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> target_color1(transform_model2scene1, 0, 255, 0);
  // pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> target_color2(transform_model2scene2, 0, 0, 255);
  // pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> point_im_color(points_im, 158, 157, 131);
  //  pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> point_im_color1(points_im1, 244, 208, 0);
  //   pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> point_im_color2(points_im2, 222 ,125, 44);
 
	// viewer->addPointCloud<pcl::PointXYZ>(source_cloud, source_color, "source", vp);
	// viewer->addPointCloud<pcl::PointXYZ>(transform_model2scene, target_color, "target", vp);
  // viewer->addPointCloud<pcl::PointXYZ>(transform_model2scene1, target_color1, "target1", vp);
  // viewer->addPointCloud<pcl::PointXYZ>(transform_model2scene2, target_color2, "target2", vp);
  // viewer->addPointCloud<pcl::PointXYZ>(points_im, point_im_color, "point_im", vp);
  //  viewer->addPointCloud<pcl::PointXYZ>(points_im1, point_im_color1, "point_im1", vp);
  //     viewer->addPointCloud<pcl::PointXYZ>(points_im2, point_im_color2, "point_im2", vp);
 
	// viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "source");
	// viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "target");
  // viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "target1");
  // viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "target2");
  //  viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "point_im");
  //   viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "point_im1");
  //   viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "point_im2");
	// viewer->spin();
  oFile<<endl;

  sum_number++;
  }
   
 oFile.close();







  



  
          
 
  return 0;
}

