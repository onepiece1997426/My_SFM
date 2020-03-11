#pragma once

#include <opencv2/opencv.hpp>

#include <fstream>
#include <iostream>
#include <vector>

#include <Eigen/SVD>
#include <Eigen/Dense>

#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include "CalMatchPoints.h"

#include <omp.h>

using namespace std;

#define COUTENDL(x) std::cout << x << std::endl

#define IMG_WIDTH 4096
#define IMG_HEIGHT 3000

#define USE_KERNEL_BA true
#define LINE_SPLIT " "
#define POINT_SPLIT ","


//global
static omp_lock_t omp_lock;  //全局互斥锁

//一个点对应了一个Tracks
/*
	该结构描述了某张图片上的某个特征点对应的一连串的track，这些track会在别的图片上出现

*/
struct PointTrack
{
	int lineIdx;
	Eigen::Vector2d point;
	int trackID;   //这个point对应的trackID tracks里面的trackID和这个值相同
	bool hasBeenAddInStructure;
	int featIdx;

	PointTrack() {}
	PointTrack(Eigen::Vector2d _p, int lidx = -1, int fidx = -1) 
	{
		point = _p;
		lineIdx = lidx;
		trackID = -1;   //没有被查找的track
		hasBeenAddInStructure = false;
		featIdx = fidx;
	}
};

//每张图片的信息，包含了位姿，包含的Track
class View
{
public:
	//informations
	std::string path;
    std::string name;
	int viewID;
	std::vector<int> neighborsID;
	//me.point -> neighbor.line
    std::map<uint, Eigen::Matrix3d> F_Me2Neighbors; //元素顺序和neigiborID要一一对应

	//poses
	Eigen::Matrix3d rotation;
	Eigen::Vector3d center;
	Eigen::Matrix3d K;
	Eigen::Matrix3d K_Inv;
	std::vector<double> distortionParams;
	Eigen::Matrix3d t_;
	Eigen::Vector3d t;

	std::vector<PointTrack> pointTracks;  
	
	double* rotation_array_aa; //轴角
	double* translatrion_array;
	double scale;
	View() {}
	View(std::string path_, int id,Eigen::Matrix3d r, Eigen::Vector3d c, Eigen::Matrix3d _K, std::vector<double> _distort)
	{
		this->path = path_;
		this->viewID = id;
        size_t sz = path.find_last_of('/') + 1;
        this->name = path.substr(sz);
		this->rotation = r;
		this->center = c;
		this->K = _K;
		this->distortionParams = _distort;

		this->K_Inv = this->K.inverse();
		t = -r * c;
		t_ << 0, -t(2), t(1),
			t(2), 0, -t(0),
			-t(1), t(0), 0;

		this->rotation_array_aa = new double[3];
		Eigen::AngleAxisd aa(this->rotation);
		Eigen::Vector3d v = aa.angle() * aa.axis();
		rotation_array_aa[0] = v.x();
		rotation_array_aa[1] = v.y();
		rotation_array_aa[2] = v.z();

		this->translatrion_array = new double[3];
		translatrion_array[0] = t.x();
		translatrion_array[1] = t.y();
		translatrion_array[2] = t.z();
		this->scale = 1.0f;   //和world的尺度
	}

	void updateParam()
	{
		t = -1 * rotation * center;

		Eigen::AngleAxisd aa(this->rotation);
		Eigen::Vector3d v = aa.angle() * aa.axis();
		rotation_array_aa[0] = v.x();
		rotation_array_aa[1] = v.y();
		rotation_array_aa[2] = v.z();

		translatrion_array[0] = t.x();
		translatrion_array[1] = t.y();
		translatrion_array[2] = t.z();		
	}

	void updateParam(Eigen::Vector3d t_)
	{
		t = t_;
		center = -1 * rotation.transpose() * t_;

		Eigen::AngleAxisd aa(this->rotation);
		Eigen::Vector3d v = aa.angle() * aa.axis();
		rotation_array_aa[0] = v.x();
		rotation_array_aa[1] = v.y();
		rotation_array_aa[2] = v.z();

		translatrion_array[0] = t.x();
		translatrion_array[1] = t.y();
		translatrion_array[2] = t.z();
	}

};


typedef std::map< uint,std::vector<std::vector<Eigen::Vector2d>> > VID_Lines;
class Reconstruction
{
public:

	struct Observation
	{
		//这个3D点能被哪些像素看到
		int viewID;
		Eigen::Vector2d pixel;
		Eigen::Vector2d pixel_norm;
		uint pixelID; 

		Observation(){}
		Observation(int id, Eigen::Vector2d p, uint p_id = 0)
		{
			viewID = id;
			pixel = p;
			pixel_norm = Eigen::Vector2d(p.x() / IMG_WIDTH, p.y() / IMG_HEIGHT);
			pixelID = p_id;
		}
	};

	struct Structure
	{
		Eigen::Vector3d position;
		Eigen::Vector3d colors; //rgb
		double* positions_array;
		std::vector<Observation> obs;
		
		int lineIdx;  //这个结构属于那一条直线（简易的mesh要用）
		bool validStructure;

		Structure() {}
		Structure(Eigen::Vector3d _p)
		{
			position = _p;
			colors = Eigen::Vector3d::Zero();
			positions_array = new double[3];
			positions_array[0] = _p.x();
			positions_array[1] = _p.y();
			positions_array[2] = _p.z();
			validStructure = true;
		}

		
		//line idx 排序 为了简便。。
		static bool sortmethod_lineIdx(Structure s1, Structure s2)
		{
			if (s1.lineIdx < s2.lineIdx)
			{
				return true;
			}
			else
				return false;
		}

		static bool sortmethod_z(Structure s1, Structure s2)
		{           

			if (s1.position.z() < s2.position.z())
			{
				return true;
			}
			else
				return false;
		}

	};


	std::vector<View> Views;
	std::vector<Structure> vec_structures;  //新 2020.01.10

	uint structureCntBeforeExpand;  //mesh用

	int trackCnt;
	
	std::map<uint, std::vector<Matches>> map_viewID_Matches;  //ID-vec_Matches,而且这后面的match的图片ID只会比前面的ViewIDX小

    std::map<uint,Eigen::Vector3d> debug_lineColor;

    Reconstruction()
    {
        trackCnt = 0;
    }
	Reconstruction(std::vector<View> vs)
	{
		Views = vs;
		trackCnt = 0;
	}

    bool Init(std::vector<string> vec_imageFile, std::vector<Eigen::Matrix<double,3,4>> vec_PMatrix, std::vector<Eigen::Matrix3d> vec_KMatrix,std::vector<std::vector<double>> vec_distort,std::vector<string> vec_featureFilePath, std::vector<string> vec_matchesFile, std::vector<string> vec_linesPointsFile);
    bool Init(std::vector<string> vec_imageFile, std::vector<Eigen::Matrix<double,3,4>> vec_PMatrix, std::vector<Eigen::Matrix3d> vec_KMatrix,std::vector<std::vector<double>> vec_distort,std::vector<std::vector<Eigen::Vector2d>> vec_vec_featurePoints, std::vector<Matches> vec_matches, VID_Lines mp_vec_vec_linesPoints);


    bool AddView(View view);
    bool UpdateNeighbor(uint hostID,uint neighborID);
    bool UpdatePointTrackFromMatch(std::map<uint, std::vector<Eigen::Vector3d> > map_VID_vecFeatures, std::vector<Matches> vec_matches);


	void findTrackInOneView(uint ViewID);
    void findLineTrackInOneView_WithKnownMatch(uint ViewID);   //直线的匹配关系已知,就是最简单的对应关系,pointtrack的序号是对应的
    bool findTrack(uint ViewID, std::map<uint,Eigen::Matrix3d> F_neighbor2this, PointTrack& fromPointTrack); //一个点在所有neighbor中的track
	void findTracksInViews();
    void FindCorelatedPointByF(double& minDst, int& minDstIdx, Eigen::Matrix3d F, uint neighborId);
	void GetStructure();
	void SavePLYFile_PointCloud(std::string filePath);
    void SavePLYFile_Mesh(std::string filePath);  //8 divide
    void SavePLYFile_Mesh_6(std::string filePath);  //6 divide
    void SavePLYFile_Mesh_Rectangle(std::string filePath);
    void SavePLYFile_Mesh_UpOneMM(std::string filePath);
	void SaveTXTFile_PointCloud(std::string filePath);

	void CalculateStructure_Init_DLT(std::map<uint, Eigen::Matrix<double, 3, 4>>& map_RtMatrix, bool fixCenter = false);  
	void CalculateStructure_Ceres(std::map<uint, Eigen::Matrix<double, 3, 4>>& map_RtMatrix, std::map<uint, Eigen::Matrix<double, 3, 4>>& map_PMatrix);
	void CalculateStructure_Ceres_AdjustRt();
	void CalculateStructure_Ceres_Adjust_t(std::map<uint, Eigen::Matrix<double, 3, 4>>& map_RtMatrix);
	void CalculateStructure_Ceres_AdjustR_FixC();
	void CalculateStructure_FirstXThenRt(std::map<uint, Eigen::Matrix<double, 3, 4>>& map_RtMatrix, std::map<uint, Eigen::Matrix<double, 3, 4>>& map_PMatrix);
    void CalculateStructure_FirstXThenR_FixCenter(std::map<uint, Eigen::Matrix<double, 3, 4>>& map_RtMatrix, std::map<uint, Eigen::Matrix<double, 3, 4>>& map_PMatrix);

    void LerpVoltaageLines(uint cntBetween2Point);  //对电线进行差值
    void ExpandStructure();  //8 division
    void ExpandStructure_4(); //6_division
    void ExpandStructure_6(); //6_division
    void ExpandStructure_Rectangle();
    void ExpandStructure_UpOneMM(); //向上扩1mm
	void ApplyTransform(Eigen::Matrix4d T);
	bool isInSide(uint viewID, Eigen::Vector3d X, Eigen::Vector2d x, Eigen::Matrix<double,3,4> T); //判断投影是否在图片中
	void CalTranslationByCeres(std::map<uint, Eigen::Matrix<double, 3, 4>>& map_RtMatrix);
	void CheckValidStructure();  //检查这个 structure 是否合理
	void EraseInvalidStructure(); //删除不合理的structure


};

Eigen::Matrix3d calFundamentalMatrix(View from, View to);
inline Eigen::Matrix3d toAntisymmetricMatrix(Eigen::Vector3d t);

void Test_reconstruct_with_KnownPoses(std::string workingFolder, std::vector<string> vec_imageFile, std::vector<string> vec_featureFilePath, std::vector<Eigen::Matrix<double,3,4>> vec_PMatrix, std::vector<Eigen::Matrix3d> vec_KMatrix,std::vector<std::vector<double>> vec_distort, std::vector<string> vec_matchesFile = std::vector<string>());

//BA
struct Ceres_Triangulate
{
	//2d 点
	const Eigen::Vector2d x;	
	const Eigen::Matrix<double, 3, 4> P;  //投影矩阵

	Ceres_Triangulate(Eigen::Vector2d x_,  Eigen::Matrix<double, 3, 4> P_):x(x_),P(P_) {}

	template<typename T>
	bool operator()(const T* const ceres_X, T* residual) const
	{

		T PX0 = T(P(0, 0)) * ceres_X[0] + T(P(0, 1)) * ceres_X[1] + T(P(0, 2)) * ceres_X[2] + T(P(0, 3));
		T PX1 = T(P(1, 0)) * ceres_X[0] + T(P(1, 1)) * ceres_X[1] + T(P(1, 2)) * ceres_X[2] + T(P(1, 3));
		T PX2 = T(P(2, 0)) * ceres_X[0] + T(P(2, 1)) * ceres_X[1] + T(P(2, 2)) * ceres_X[2] + T(P(2, 3));

		PX0 = PX0 / PX2;
		PX1 = PX1 / PX2;	

		residual[0] = T(x.x()) - PX0;
		residual[1] = T(x.y()) - PX1;

		return true;
	}

};


struct Ceres_Triangulate_AdjustRt
{
	const Eigen::Vector2d x;
	const Eigen::Matrix<double, 3, 3> K;  //投影矩阵

	Ceres_Triangulate_AdjustRt(Eigen::Vector2d x_, Eigen::Matrix<double, 3, 3> K_ ) :x(x_), K(K_){}

	template<typename T>
	bool operator()(const T* const ceres_X, const T* const ceres_angleAxis, const T* const ceres_t, T* residual) const
	{	
		T PX[3];
		PX[0] = ceres_X[0];
		PX[1] = ceres_X[1];
		PX[2] = ceres_X[2];

		T PX_r[3];
		ceres::AngleAxisRotatePoint(ceres_angleAxis, PX, PX_r);

		T PX0 = T(K(0, 0)) * (PX_r[0] + ceres_t[0]) + T(K(0, 2)) * (PX_r[2] + ceres_t[2]);
		T PX1 = T(K(1, 1)) * (PX_r[1] + ceres_t[1]) + T(K(1, 2)) * (PX_r[2] + ceres_t[2]);
		T PX2 = (PX_r[2] + ceres_t[2]);

		PX0 = PX0 / PX2;
		PX1 = PX1 / PX2;

		residual[0] = T(x.x()) - PX0;
		residual[1] = T(x.y()) - PX1;
		return true;
	}

};

struct Ceres_Triangulate_withScale_R
{
	const Eigen::Vector2d x;
	const Eigen::Matrix<double, 3, 3> K;  //投影矩阵
	const double scaleToWorld;
	const Eigen::Vector3d c;

	Ceres_Triangulate_withScale_R(Eigen::Vector2d x_, Eigen::Matrix<double, 3, 3> K_, double scaleToWorld_ , Eigen::Vector3d c_) :x(x_), K(K_), scaleToWorld(scaleToWorld_),c(c_) {}

	template<typename T>
	bool operator()(const T* const ceres_X, const T* const ceres_angleAxis, T* residual) const
	{
		T PX[3];
		PX[0] = ceres_X[0] * T(scaleToWorld);
		PX[1] = ceres_X[1] * T(scaleToWorld);
		PX[2] = ceres_X[2] * T(scaleToWorld);

		T PX_C[3];
		PX_C[0] = PX[0] - T(c.x());
		PX_C[1] = PX[1] - T(c.y());
		PX_C[2] = PX[2] - T(c.z());

		T PX_r[3];
		ceres::AngleAxisRotatePoint(ceres_angleAxis, PX_C, PX_r);
		

		T PX0 = T(K(0, 0)) * (PX_r[0]) + T(K(0, 2)) * (PX_r[2]);
		T PX1 = T(K(1, 1)) * (PX_r[1]) + T(K(1, 2)) * (PX_r[2]);
		T PX2 = PX_r[2] ;

		PX0 = PX0 / PX2;
		PX1 = PX1 / PX2;

		residual[0] = T(x.x()) - PX0;
		residual[1] = T(x.y()) - PX1;
		return true;
	}
};

struct Ceres_AdjustT
{
	const Eigen::Vector2d x;
	const Eigen::Matrix<double, 3, 3> K;  //投影矩阵
	double* angel_axis_aa;

	Ceres_AdjustT(Eigen::Vector2d x_, Eigen::Matrix<double, 3, 3> K_, double* angle_axis_aa_) :x(x_), K(K_)  
	{
		angel_axis_aa = angle_axis_aa_;
	}

	template<typename T>
	bool operator()(const T* const ceres_X, const T* const ceres_t, T* residual) const
	{
		T PX[3];
		PX[0] = ceres_X[0];
		PX[1] = ceres_X[1];
		PX[2] = ceres_X[2];

		T PX_r[3];
		T ceres_angel_axis_aa[3];
		ceres_angel_axis_aa[0] = T(angel_axis_aa[0]);
		ceres_angel_axis_aa[1] = T(angel_axis_aa[1]);
		ceres_angel_axis_aa[2] = T(angel_axis_aa[2]);

		ceres::AngleAxisRotatePoint(ceres_angel_axis_aa, PX, PX_r);

		T PX0 = T(K(0, 0)) * (PX_r[0] + ceres_t[0]) + T(K(0, 2)) * (PX_r[2] + ceres_t[2]);
		T PX1 = T(K(1, 1)) * (PX_r[1] + ceres_t[1]) + T(K(1, 2)) * (PX_r[2] + ceres_t[2]);
		T PX2 = (PX_r[2] + ceres_t[2]);

		PX0 = PX0 / PX2;
		PX1 = PX1 / PX2;

		residual[0] = T(x.x()) - PX0;
		residual[1] = T(x.y()) - PX1;
		return true;
	}

};


struct Ceres_globalTranslation
{
	const Eigen::Matrix3d R_relative;  //i: I|0   j: R_relative|t_relative
	double* r_relative_array;

	Ceres_globalTranslation(Eigen::Matrix3d R_relative_) : R_relative(R_relative_)
	{
		r_relative_array = new double[3];
		Eigen::AngleAxisd aa(R_relative);
		Eigen::Vector3d tmp = aa.angle() * aa.axis();
		r_relative_array[0] = tmp.x();
		r_relative_array[1] = tmp.y();
		r_relative_array[2] = tmp.z();
	}

	template<typename T> 
	bool operator()(const T* const ti, const T* const tj, const T* const lambda, T* residual) const
	{
		//Eigen::Vector3d t_relative = map_RtMatrix.at(neighborID).block(0, 3, 3, 1) - R_relative * map_RtMatrix.at(hostID).block(0, 3, 3, 1);
		T t_relative[3];	
		T t_r[3]; 		

		T r[3];
		r[0] = T(r_relative_array[0]);
		r[1] = T(r_relative_array[1]);
		r[2] = T(r_relative_array[2]);
		ceres::AngleAxisRotatePoint(r, ti, t_r);
		t_relative[0] = tj[0] - t_r[0];
		t_relative[1] = tj[1] - t_r[1];
		t_relative[2] = tj[2] - t_r[2];

		
		/*residual[0] = tj[0] - t_r[0] - lambda[0] * t_relative[0];
		residual[1] = tj[1] - t_r[1] - lambda[0] * t_relative[1];
		residual[2] = tj[2] - t_r[2] - lambda[0] * t_relative[2];*/

		residual[0] = ti[0] - tj[0] - lambda[0] * t_relative[0];
		residual[1] = ti[1] - tj[1] - lambda[0] * t_relative[1];
		residual[2] = ti[2] - tj[2] - lambda[0] * t_relative[2];

		return true;
	}

};


struct Ceres_globalTranslation_KnownXInit
{
	const Eigen::Matrix3d R_relative;
	const Eigen::Vector3d tao_relative;  //相对平移方向
	const Eigen::Vector3d X_init;
	const Eigen::Vector2d x;
	const Eigen::Matrix3d K;

	Ceres_globalTranslation_KnownXInit(Eigen::Vector3d X_init_, Eigen::Vector2d x_, Eigen::Matrix3d K_, Eigen::Matrix3d R_relative_, Eigen::Vector3d tao_relative_) : X_init(X_init_), x(x_), K(K_), R_relative(R_relative_), tao_relative(tao_relative_)
	{

	}

	template<typename T>
	bool operator()(const T* const lambda, T* residual) const
	{

		Eigen::Vector3d RX = R_relative * X_init;
		
		T ceres_X[3];
		ceres_X[0] = T(RX(0));
		ceres_X[1] = T(RX(1));
		ceres_X[2] = T(RX(2));

		T XC0 = ceres_X[0] + lambda[0] * T(tao_relative(0));
		T XC1 = ceres_X[1] + lambda[0] * T(tao_relative(1));
		T XC2 = ceres_X[2] + lambda[0] * T(tao_relative(2));
		
		T PX0 = T(K(0, 0)) * XC0 + T(K(0, 2)) * XC2;
		T PX1 = T(K(1, 1)) * XC1 + T(K(1, 2)) * XC2;
		T PX2 = XC2;
		PX0 /= PX2;
		PX1 /= PX2;

		residual[0] = T(x.x()) - PX0;
		residual[1] = T(x.y()) - PX1;

		return true;
	}


};
