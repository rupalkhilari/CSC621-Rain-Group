/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include <vector>
#include <unordered_map>

#include "itkConnectedThresholdImageFilter.h"
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkRGBPixel.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkBinomialBlurImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

typedef   float           InternalPixelType;
const     unsigned int    Dimension = 3;
typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;
typedef float                           OutputPixelType;
typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

class RegionGrowingNoThreshold {
 public:
  std::vector<std::vector<int>> GetCentroids(char filename[], int seed_x, int seed_y, int seed_z);
  void WriteImageWithCentroids(char filename[], char outfilename[], std::vector<std::vector<int>> centroids);
 protected:
  void ComputeStats(InternalImageType::Pointer& image,
                    std::unordered_map<long long, InternalImageType::IndexType>& region_points,
                     double* mgv, double* upper_dev, double* lower_dev);
  void ComputeThreshold(InternalImageType::Pointer& image,
                        std::unordered_map<long long, InternalImageType::IndexType>& region_points,
                        double* upper_threshold, double* lower_threshold);
  void ComputeFinalThreshold(InternalImageType::Pointer& image,
                             std::unordered_map<long long, InternalImageType::IndexType>& region_points,
                             double* upper_threshold, double* lower_threshold, double* upper_dev, double* lower_dev);
  long long ComputePointHash(InternalImageType::IndexType& point);
  bool CheckBounds(InternalImageType::IndexType& point,
                   InternalImageType::SizeType size);
  bool GrowRegions(InternalImageType::Pointer& image,
                   InternalImageType::IndexType& seed,
                   double* upper_threshold,
                   double* lower_threshold);
  std::vector<std::vector<int>> ComputeCentroids(InternalImageType::Pointer& image,
    double* max_image_x_dist_ratio, double* max_image_y_dist_ratio, double* max_image_z_dist_ratio, bool final);

  void ComputePerSliceAvgAndDev(InternalImageType::Pointer& image, double* avg, double* dev);

};
