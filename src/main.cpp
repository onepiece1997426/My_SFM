#pragma once
#include "Reconstruction.h"
#include "CalculateTransfromToWorld.h"


/*
	知道直线顺序的情况下进行3D重建，简化了在对极线上找点的情况
*/
int main()
{
	std::string filePath0 = "E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\001_l_Undistort.png";
	std::string filePath1 = "E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\001_m_Undistort.png";
	std::string filePath2 = "E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\001_r_Undistort.png";

	/*
	//find points 用二值化的曲线 相当于找到了直线上的点 顺序 我可以通过某种方法根据横坐标的增加找直线
	cv::Mat mat = cv::imread("E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\mask_l_center.png");
	cv::cvtColor(mat, mat,cv::COLOR_BGR2GRAY);
	std::ofstream ofs("E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\point_l.txt", std::ios::out);
	int cnt = 0;
	while (true)
	{
		for (uint i = 0; i < mat.rows; i++)
		{
			for (uint j = 0; j < mat.cols; j++)
			{
				if (mat.ptr<uchar>(i)[j] > 200)
				{
					cnt++;
					ofs << j << " " << i << std::endl;  //每行只找最小的x，找完再把她置零，这样可以按x顺序找到不同的直线
					mat.ptr<uchar>(i)[j] = 0;
					break;
				}
			}
		}
		ofs << std::endl;
		if (cnt == 0)
		{
			break;
		}
		cnt = 0;
	}
	ofs.close();
	ofs.flush();
	return 0;
	*/
	std::vector<PointTrack> pointsTracks;
	std::vector<PointTrack> pointsTracks2;
	std::vector<PointTrack> pointsTracks3;
	//init line points
	std::ifstream ifs("E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\point_l.txt", std::ios::in);
	if (!ifs)
	{		
		COUTENDL("Err in file");
		return -1;
	}
	std::string str;
	uint lineIdx = 0;
	while (std::getline(ifs, str))
	{
		if (str == "")
		{
			lineIdx++;
			continue;
		}
			
		size_t pos = str.find_first_of(' ', 0);
		if (pos == str.npos)
			continue;
		double x = std::atof(str.substr(0, pos).c_str());
		double y = std::atof(str.substr(pos+1).c_str());
		Eigen::Vector2d X(x, y);
		PointTrack pt(X,lineIdx);
		pointsTracks.emplace_back(pt);
	}
	ifs.close();

	lineIdx = 0;
	std::ifstream ifs2("E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\point_m.txt", std::ios::in);
	if (!ifs2)
	{
		COUTENDL("Err in file");
		return -1;
	}	
	while (std::getline(ifs2, str))
	{
		if (str == "")
		{
			lineIdx++;
			continue;
		}
		size_t pos = str.find_first_of(' ', 0);
		if (pos == str.npos)
			continue;
		double x = std::atof(str.substr(0, pos).c_str());
		double y = std::atof(str.substr(pos + 1).c_str());
		Eigen::Vector2d X(x, y);
		PointTrack pt(X,lineIdx);
		pointsTracks2.emplace_back(pt);
	}
	ifs.close();	


	lineIdx = 0;
	std::ifstream ifs3("E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\point_r.txt", std::ios::in);	
	if (!ifs3)
	{
		COUTENDL("Err in file");
		return -1;
	}
	while (std::getline(ifs3, str))
	{
		if (str == "")
		{
			lineIdx++;
			continue;
		}
		size_t pos = str.find_first_of(' ', 0);
		if (pos == str.npos)
			continue;
		double x = std::atof(str.substr(0, pos).c_str());
		double y = std::atof(str.substr(pos + 1).c_str());
		Eigen::Vector2d X(x, y);
		PointTrack pt(X, lineIdx);
		pointsTracks3.emplace_back(pt);
	}
	ifs.close();
	//----------------finish init points

	//R: 绝对旋转，到world的I的旋转矩阵 这里还没写

	Eigen::Matrix3d R0;
	R0 = Eigen::Matrix3d::Identity();//cv_compute F
	//R0 << 0.9999985164015987, -0.0016852721912964716, 0.0003564438848244992,
	//	0.0016852537762661407, 0.9999985786090674, 0.00005195723735700068,
	//	-0.0003565309402656449, -0.000051356461870418339, 0.9999999351240991;  // open mvg
	Eigen::Vector3d c0;
	//c0 << 3761.29, 2045.79,1479.62;
	c0 << 0, 0, 0; //cv_compute F
	//c0 << -0.015210322515032325, -0.010315869704961822, -0.0369290944359687; //open mvg
	Eigen::Matrix3d K0;
	K0 << 3.7612906067774788e+003, 0., 2.0457995013785016e+003, 0.,
		3.7159736130516289e+003, 1.4796241830921265e+003, 0., 0., 1.;
	std::vector<double> distort0{ -7.5267404892772435e-002, 1.2710728604878160e-001,
	   -2.5774289715705485e-003, -3.6861402083794338e-005,
	   2.1573505064134832e-001 };
	View view0(filePath0, 0, R0, c0,K0,distort0);
	view0.pointTracks = pointsTracks;
	view0.scale = 1.0f;

	Eigen::Matrix3d R1 ;
	//R1 = Eigen::Matrix3d::Identity();
	R1 << 0.9997877399536809, 0.001533996409847495, 0.02054560520708353,
	-0.0009714381574838053, 0.9996250910635111, -0.02736299735359312,
	-0.02057987721578888, 0.02733723049764072, 0.9994144007780266;   //cv compute f
	//R1 << 0.9998191962875127, -0.0006365931544684841, 0.01900445958578646,
	//	0.0011306431573686433, 0.9996613784675856, -0.025997116115776925,
	//	-0.018981474680403587, 0.026013903002859389, 0.9994813557388227;  //open mvg
	Eigen::Vector3d c1;
	c1 << 0.999264977438602, -0.02242492054618143, -0.03109063851270102; //cv compute f
	//c1 << 0.012314168371418967, -0.010701346358933224, -0.03601467825925361; //open mvg
	/*Eigen::Vector3d t10;
	t10 << -0.9989705126323998, 0.02035518830274818, 0.04054110506782228;
	c1 = -R1.transpose() * t10;*/
	Eigen::Matrix3d K1;
	K1 << 3.7516261809333914e+003, 0., 2.0321668349416393e+003, 0.,
		3.7075460080696898e+003, 1.4757973862050546e+003, 0., 0., 1.;
	std::vector<double> distort1{ -8.9706636166420994e-002, 2.2503707166823739e-001,
	   -2.7655590432857850e-003, 1.1201754344822839e-004,
	   -8.2725767124464909e-003 };
	View view1(filePath1, 1, R1, c1, K1, distort1);
	view1.pointTracks = pointsTracks2;
	view1.scale = 1.0f;

	Eigen::Matrix3d R2;
	R2 << 0.9993297494017308, -0.0104452284852808, 0.03508488510119413,
		0.01157459428116065, 0.9994162839037561, -0.03214218777919282,
		-0.03472867299365523, 0.03252673776905168, 0.9988673238234408; //cv compute f
	//R2 << 0.99930614288737, -0.012159130892335009, 0.03520494742945117,
	//	0.013310805482134959, 0.9993776830569874, -0.03266602309838801,
	//	-0.03478584834360263, 0.03311196375316205, 0.9988461055695344;  //open mvg
	Eigen::Vector3d c2;
	c2 << 0.998090453592277, 0.00187607723956934, 0.06174080321921585;  //cv compute F
	//c2 << 0.029342534293331978, -0.009267891856573398, -0.030697318565627349; //open mvg
	/*Eigen::Vector3d t20;
	t20 << -0.9996867425646876,-0.02464782888913964,-0.004347559252713634;
	c2 = -R2.transpose() * t20;*/

	Eigen::Matrix3d K2;
	K2 << 3.7604405383652875e+003, 0., 2.0415438075955501e+003, 0.,
		3.7160779331549415e+003, 1.4845370511494227e+003, 0., 0., 1.;
	std::vector<double> distort2{ -8.9249034584200818e-002, 2.4976814960048890e-001,
	   -2.3008498513017160e-003, 1.4399452868398138e-003,
	   -4.7341695340444627e-002 };
	View view2(filePath2, 2, R2, c2, K2, distort2);
	view2.pointTracks = pointsTracks3;
	view2.scale = 1.0f;

	//neighborsID 为 viewvector中的index
	//ADD neighbors informations
	Eigen::Matrix3d F1, F2, F3;
	view0.neighborsID.emplace_back(1);	
	F1 << -1.148289396524095e-09, -1.443613059989445e-07, 0.0004640568681676394,
	8.849411915859721e-08, 7.943834668178488e-08, 0.01018649933100411,
	-0.0003520947904074312, -0.01029941539296455, 1;   //l->m cv compute F
	//Eigen::Matrix3d R_relative = R1 * R0.transpose();
	//Eigen::Vector3d t_relative = view1.t - R_relative * view0.t;
	//F1 = view1.K_Inv.transpose() * toAntisymmetricMatrix(t_relative) * R_relative * view0.K_Inv; // 1->m openmvg F	
	//F1 /= F1(2, 2);
	view0.F_Me2Neighbors.emplace_back(F1);	
	view0.neighborsID.emplace_back(2);
	F2 << 1.313892622952379e-09, 4.98403599615894e-08, -0.0001619510423509496,
	-1.15110500440471e-07, 5.883297886126666e-08, 0.007141016183111804,
	0.0001665841826363756, -0.007202233298315958, 1; //cv compute F
	/*R_relative = R2 * R0.transpose();
	t_relative = view2.t - R_relative * view0.t;
	F2 = view2.K_Inv.transpose() * toAntisymmetricMatrix(t_relative) * R_relative * view0.K_Inv;
	F2 /= F2(2, 2);*/
	view0.F_Me2Neighbors.emplace_back(F2);   //l->r

	view1.neighborsID.emplace_back(0);
	view1.F_Me2Neighbors.emplace_back(F1.transpose()); //m->l 
	view1.neighborsID.emplace_back(2);
	F3 << -5.197542569214979e-09, -5.8243541334072e-07, 0.0005972575469230523,
	4.536056796174939e-07, 1.627600770017867e-08, 0.02502001577511281,
	-0.0006681380558518357, -0.02495766868377103, 1;
	view1.F_Me2Neighbors.emplace_back(F3); //m->r

	view2.neighborsID.emplace_back(0);
	view2.F_Me2Neighbors.emplace_back(F2.transpose()); //r->l
	view2.neighborsID.emplace_back(1);
	view2.F_Me2Neighbors.emplace_back(F3.transpose()); //r->m

	//通过已知三个世界相机的c和分解出的c求解一个到世界系的变换 
	//Eigen::Vector3d cw0, cw1, cw2;
	//cw0 << 945.115, 1024.913, 52.990;
	//cw1 << 949.300, 1021.623, 52.989;
	//cw2 << 951.961, 1019.501, 52.985;
	//std::vector<Eigen::Vector3d> pts_wordC{ cw0,cw1,cw2 };
	//std::vector<Eigen::Vector3d> pts_needAlignedC{ c0,c1,c2 };
	//Eigen::Matrix4d T = CalTransformToWorldWhenKnowTranslationDir(pts_wordC, pts_needAlignedC);
	//view0.center = pts_needAlignedC[0];
	//view0.updateParam();
	//view1.center = pts_needAlignedC[1];
	//view1.updateParam();
	//view2.center = pts_needAlignedC[2]; //按照local三个点的坐标进行缩放之后的点，不是local2world，没有*T呢还
	//view2.updateParam();

	//更新参数： 哪个好？ 先变换(P->PT^-1)，再直接做BA，还是先做BA（更新C后），在去做旋转平移？
	//Eigen::Matrix4d T_1 = T.inverse();
	//Eigen::Matrix3d sR1 = view0.rotation * T_1.block(0, 0, 3, 3);
	//Eigen::Vector3d t1 = view0.rotation * T_1.block(0, 3, 3, 1) + view0.t;
	////Eigen::Vector3d c_test = -sR1.inverse() * t1;
	//double scale1 = std::pow(sR1.determinant(), 0.3333333);
	//view0.rotation = sR1 / scale1;
	//view0.t = t1; 
	//view0.scale = scale1;
	//view0.updateParam();
	//Eigen::Matrix3d sR2 = view1.rotation * T_1.block(0, 0, 3, 3);
	//Eigen::Vector3d t2 = view1.rotation * T_1.block(0, 3, 3, 1) + view1.t;
	//double scale2 = std::pow(sR2.determinant(), 0.3333333);
	//view1.rotation = sR2 / scale2;
	//view1.t = t2;
	//view1.scale = scale2;
	//view1.updateParam();
	//Eigen::Matrix3d sR3 = view2.rotation * T_1.block(0, 0, 3, 3);
	//Eigen::Vector3d t3 = view2.rotation * T_1.block(0, 3, 3, 1) + view2.t;
	//double scale3 = std::pow(sR3.determinant(), 0.3333333);
	//view2.rotation = sR3 / scale3;
	//view2.t = t3;
	//view2.scale = scale3;
	//view2.updateParam();

	
	//std::vector<View> views{ view0,view1 };

	/*Reconstruction reconstruction(views);
	reconstruction.findTrackInOneView(0);*/
	//reconstruction.findTracksInViews();

	/*test all kind structure.*/	
	std::string file_l_Points = "E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\data-total\\data-total\\001_l_Undistort.txt";
	std::string file_m_Points = "E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\data-total\\data-total\\001_m_Undistort.txt";
	std::string file_r_Points = "E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\data-total\\data-total\\001_r_Undistort.txt";
	std::vector<Eigen::Vector2d> points_l = getPoint(file_l_Points, ",");
	uint linePointSize1 = view0.pointTracks.size();
	uint linePointSize2 = view1.pointTracks.size();
	uint linePointSize3 = view2.pointTracks.size();
	for (uint i = 0; i < points_l.size(); i++)
	{
		PointTrack pt(points_l[i], -1, i + linePointSize1);
		view0.pointTracks.emplace_back(pt);
	}
	points_l = getPoint(file_m_Points, ",");
	for (uint i = 0; i < points_l.size(); i++)
	{
		PointTrack pt(points_l[i], -1, i + linePointSize2);
		view1.pointTracks.emplace_back(pt);
	}
	points_l = getPoint(file_r_Points, ",");
	for (uint i = 0; i < points_l.size(); i++)
	{
		PointTrack pt(points_l[i], -1, i + linePointSize3);
		view2.pointTracks.emplace_back(pt);
	}
	std::string matchFile_lm = "E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\data-total\\data-total\\001_l_Undistort.txt_001_m_Undistort.txt";
	std::string matchFile_lr = "E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\data-total\\data-total\\001_l_Undistort.txt_001_r_Undistort.txt";
	std::string matchFile_mr = "E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\data-total\\data-total\\001_m_Undistort.txt_001_r_Undistort.txt";
	Matches matches_lm = getMatchesFromFile(matchFile_lm, 0, 1, linePointSize1, linePointSize2);
	Matches matches_lr = getMatchesFromFile(matchFile_lr, 0, 2, linePointSize1, linePointSize3);
	Matches matches_mr = getMatchesFromFile(matchFile_mr, 1, 2, linePointSize2, linePointSize3);
	std::vector<Matches> view0Matches{ matches_lm ,matches_lr };
	std::vector<Matches> view1Matches{ matches_mr };
	//
	
	std::vector<View> views{ view0,view1,view2 };
	Reconstruction reconstruction(views);
	reconstruction.map_viewID_Matches[0] = view0Matches;
	reconstruction.map_viewID_Matches[1] = view1Matches;
	reconstruction.findTrackInOneView(0);
	//reconstruction.findTracksInViews();

	COUTENDL("begin sTRUCTURE....");
	reconstruction.GetStructure();
	reconstruction.SavePLYFile_PointCloud("E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\test_mvg_orig.ply");
	reconstruction.SaveTXTFile_PointCloud("E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\test_mvg_orig.txt");
	//reconstruction.ExpandStructure();
	//reconstruction.SavePLYFile_PointCloud("E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\test_mvg_expand.ply");
	//reconstruction.ApplyTransform(T);
	//reconstruction.SavePLYFile_PointCloud("E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\test_mvg_orig_transfrom.ply");
	//reconstruction.SavePLYFile_Mesh("E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\test\\test1_Mesh.ply");

	return 0;
}


