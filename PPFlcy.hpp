#include <pcl/features/feature.h>
#include <pcl/features/boost.h>
using namespace pcl; 
using namespace std; 

template <typename PointInT, typename PointNT, typename PointOutT>
  class PPFEstimationlcy : public FeatureFromNormals<PointInT, PointNT, PointOutT>
  {
    public:
      typedef boost::shared_ptr<PPFEstimationlcy<PointInT, PointNT, PointOutT> > Ptr;
      typedef boost::shared_ptr<const PPFEstimationlcy<PointInT, PointNT, PointOutT> > ConstPtr;
      using PCLBase<PointInT>::indices_;
      using Feature<PointInT, PointOutT>::input_;
      using Feature<PointInT, PointOutT>::feature_name_;
      using Feature<PointInT, PointOutT>::getClassName;
      using FeatureFromNormals<PointInT, PointNT, PointOutT>::normals_;

      typedef pcl::PointCloud<PointOutT> PointCloudOut;

      /** \brief Empty Constructor. */
      PPFEstimationlcy ();


    private:
      /** \brief The method called for actually doing the computations
        * \param[out] output the resulting point cloud (which should be of type pcl::PPFSignature);
        * its size is the size of the input cloud, squared (i.e., one point for each pair in
        * the input cloud);
        */
      void
      computeFeature (PointCloudOut &output);
  };

template <typename PointInT, typename PointNT, typename PointOutT>
PPFEstimationlcy<PointInT, PointNT, PointOutT>::PPFEstimationlcy ()
    : FeatureFromNormals <PointInT, PointNT, PointOutT> ()
{
  feature_name_ = "PPFEstimationlcy";
  // Slight hack in order to pass the check for the presence of a search method in Feature::initCompute ()
  Feature<PointInT, PointOutT>::tree_.reset (new pcl::search::KdTree <PointInT> ());
  Feature<PointInT, PointOutT>::search_radius_ = 0.05f;//radius m
}


//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointInT, typename PointNT, typename PointOutT> void
PPFEstimationlcy<PointInT, PointNT, PointOutT>::computeFeature (PointCloudOut &output)
{
  // Initialize output container - overwrite the sizes done by Feature::initCompute ()
  output.points.resize (indices_->size () * input_->points.size ());
  output.height = 1;
  output.width = static_cast<uint32_t> (output.points.size ());
  output.is_dense = true;

  // Compute point pair features for every pair of points in the cloud
  for (size_t index_i = 0; index_i < indices_->size (); ++index_i)
  {
    size_t i = (*indices_)[index_i];
    for (size_t j = 0 ; j < input_->points.size (); ++j)
    {
      PointOutT p;
      if (i != j)
      {
        if (//pcl::computePPFPairFeature
            pcl::computePairFeatures (input_->points[i].getVector4fMap (),
                                      normals_->points[i].getNormalVector4fMap (),
                                      input_->points[j].getVector4fMap (),
                                      normals_->points[j].getNormalVector4fMap (),
                                      p.f1, p.f2, p.f3, p.f4))
        {
          // Calculate alpha_m angle
          Eigen::Vector3f model_reference_point = input_->points[i].getVector3fMap (),
                          model_reference_normal = normals_->points[i].getNormalVector3fMap (),
                          model_point = input_->points[j].getVector3fMap ();
          float rotation_angle = acosf (model_reference_normal.dot (Eigen::Vector3f::UnitX ()));
          bool parallel_to_x = (model_reference_normal.y() == 0.0f && model_reference_normal.z() == 0.0f);
          Eigen::Vector3f rotation_axis = (parallel_to_x)?(Eigen::Vector3f::UnitY ()):(model_reference_normal.cross (Eigen::Vector3f::UnitX ()). normalized());
          Eigen::AngleAxisf rotation_mg (rotation_angle, rotation_axis);
          Eigen::Affine3f transform_mg (Eigen::Translation3f ( rotation_mg * ((-1) * model_reference_point)) * rotation_mg);

          Eigen::Vector3f model_point_transformed = transform_mg * model_point;
          float angle = atan2f ( -model_point_transformed(2), model_point_transformed(1));
          if (sin (angle) * model_point_transformed(2) < 0.0f)
            angle *= (-1);
          p.alpha_m = -angle;
        }
        else
        {
          PCL_ERROR ("[pcl::%s::computeFeature] Computing pair feature vector between points %u and %u went wrong.\n", getClassName ().c_str (), i, j);
          p.f1 = p.f2 = p.f3 = p.f4 = p.alpha_m = std::numeric_limits<float>::quiet_NaN ();
          output.is_dense = false;
        }
      }
      // Do not calculate the feature for identity pairs (i, i) as they are not used
      // in the following computations
      else
      {
        p.f1 = p.f2 = p.f3 = p.f4 = p.alpha_m = std::numeric_limits<float>::quiet_NaN ();
        output.is_dense = false;
      }

      output.points[index_i*input_->points.size () + j] = p;
    }
  }
}

