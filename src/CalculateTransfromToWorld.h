#pragma once
#include <opencv2/opencv.hpp>

#include <fstream>
#include <iostream>
#include <vector>

#include <Eigen/SVD>
#include <Eigen/Dense>

Eigen::Matrix4d CalTransformToWorldWhenKnowTranslationDir(std::vector<Eigen::Vector3d>& pts_fixed, std::vector<Eigen::Vector3d>& pts_aligned);

Eigen::Matrix4d PoseEstimation3D3D_SVD(std::vector<Eigen::Vector3d> ptsFixed, std::vector<Eigen::Vector3d> ptsAligned);
Eigen::Matrix4d PoseEstimation3D3D_DLT(std::vector<Eigen::Vector3d> ptsFixed, std::vector<Eigen::Vector3d> ptsAligned);
Eigen::Matrix4d PoseEstimation3D3D_CrossProduct(std::vector<Eigen::Vector3d> ptsFixed, std::vector<Eigen::Vector3d> ptsAligned, double s);
Eigen::Matrix4d PoseEstimation3D3D_CrossProduct_point(std::vector<Eigen::Vector3d> ptsFixed, std::vector<Eigen::Vector3d> ptsScaledAligned, double s);

