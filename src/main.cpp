#pragma once
#include "Reconstruction.h"
#include "CalculateTransfromToWorld.h"


/*
	知道直线顺序的情况下进行3D重建，简化了在对极线上找点的情况
*/
int main()
{
    std::string filePath0 = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/part1/huzhou/undistort/001_l_Undistort.png";
    std::string filePath1 = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/part1/huzhou/undistort/001_m_Undistort.png";
    std::string filePath2 = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/part1/huzhou/undistort/001_r_Undistort.png";
//    std::string filePath0 = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/part3/wafangdian_chooseName/undistort/001_l_Undistort.png";
//    std::string filePath1 = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/part3/wafangdian_chooseName/undistort/001_m_Undistort.png";
//    std::string filePath2 = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/part3/wafangdian_chooseName/undistort/001_r_Undistort.png";


    std::string pointDataFolder = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/data-total/data-total";
    std::string featureFile0 = pointDataFolder + "/001_l_Undistort.txt";
    std::string featureFile1 = pointDataFolder + "/001_m_Undistort.txt";
    std::string featureFile2 = pointDataFolder + "/001_r_Undistort.txt";
    std::string matchFile_lm = pointDataFolder + "/001_l_Undistort.txt_001_m_Undistort.txt";
    std::string matchFile_lr = pointDataFolder + "/001_l_Undistort.txt_001_r_Undistort.txt";
    std::string matchFile_mr = pointDataFolder + "/001_m_Undistort.txt_001_r_Undistort.txt";

    std::string lineDataFolder = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/final";
    std::string lineDataFile0 = lineDataFolder + "/point_l.txt" ;
    std::string lineDataFile1 = lineDataFolder + "/point_m.txt" ;
    std::string lineDataFile2 = lineDataFolder + "/point_r.txt" ;
//    std::string lineDataFile0 = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/jiajie/huzhou/huzhou_data/huzhou_data/match_points_l.txt" ;
//    std::string lineDataFile1 = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/jiajie/huzhou/huzhou_data/huzhou_data/match_points_m.txt" ;

//    std::string lineDataFile0 = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/part3/wafangdian_chooseName/dalian_lines/jaijie2/dalian_data/point_l.txt" ;
//    std::string lineDataFile1 = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/part3/wafangdian_chooseName/dalian_lines/jaijie2/dalian_data/point_m.txt" ;
//    std::string lineDataFile2 = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/part3/wafangdian_chooseName/dalian_lines/jaijie2/dalian_data/point_r.txt" ;

     /*
	//find points 用二值化的曲线 相当于找到了直线上的点 顺序 我可以通过某种方法根据横坐标的增加找直线
    cv::Mat mat = cv::imread("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/data-total/data-total\\mask_l_center.png");
	cv::cvtColor(mat, mat,cv::COLOR_BGR2GRAY);
    std::ofstream ofs("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/data-total/data-total\\point_l.txt", std::ios::out);
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

	//R: 绝对旋转，到world的I的旋转矩阵 这里还没写

	Eigen::Matrix3d R0;
//    R0 = Eigen::Matrix3d::Identity();//cv_compute F
//	R0 << 0.9999985164015987, -0.0016852721912964716, 0.0003564438848244992,
//		0.0016852537762661407, 0.9999985786090674, 0.00005195723735700068,
//		-0.0003565309402656449, -0.000051356461870418339, 0.9999999351240991;  // open mvg
    R0 << 0.775071,  -0.631781, -0.0108743,
            0.0352285,  0.0603884,  -0.997553,
             0.630892,   0.772791,   0.069062;    //world_pnp huzhou
//    R0 << -0.75596,    -0.654554,   0.00914511  ,
//            -0.0120758, -2.38882e-05,    -0.999927 ,
//              0.654507,    -0.756015,   -0.0078862 ; //dalian
    Eigen::Vector3d c0;
//    c0 << 0, 0, 0; //cv_compute F
	//c0 << -0.015210322515032325, -0.010315869704961822, -0.0369290944359687; //open mvg
    c0 <<  943.42, 1019.13, 54.7952; //world pnp huzhou
//    c0 << 974.676, 1021.41, 16.9945;

	Eigen::Matrix3d K0;
    K0 << 3.7612906067774788e+003, 0., 2.0457995013785016e+003, 0.,
        3.7159736130516289e+003, 1.4796241830921265e+003, 0., 0., 1.;  //huzhou
//    K0 << 3.7343664904677757e+003, 0., 2.0458971603779855e+003, 0.,
//            3.6840622549040932e+003, 1.5309521054590580e+003, 0., 0., 1.;  //dalian

    std::vector<double> distort0{ -7.5267404892772435e-002, 1.2710728604878160e-001,
       -2.5774289715705485e-003, -3.6861402083794338e-005,
       2.1573505064134832e-001 };  //huzhou
//    std::vector<double> distort0{ -9.5087746874849985e-002, 1.9459308485232937e-001,7.1079938650945612e-002 ,
//                                  8.4600547566187389e-005, -1.1331658818119056e-004}; //dalian

	Eigen::Matrix3d R1 ;
	//R1 = Eigen::Matrix3d::Identity();
//    R1 << 0.9997844543318415, 0.001318258126926052, 0.02071972663652263,
//    -0.0007616346938534441, 0.9996392000591341, -0.02684938766019824,
//    -0.02074764538386343, 0.02682781952834184, 0.9994247361909652;   //cv compute f
	//R1 << 0.9998191962875127, -0.0006365931544684841, 0.01900445958578646,
	//	0.0011306431573686433, 0.9996613784675856, -0.025997116115776925,
	//	-0.018981474680403587, 0.026013903002859389, 0.9994813557388227;  //open mvg
    R1 << 0.785548, -0.618716, -0.010268,
          0.0190466, 0.0407612, -0.998987,
           0.618508,  0.784557, 0.0438043;  //world pnp huzhou
//    R1 << -0.740306,   -0.672267, -0.00202373 ,
//            0.00930091, -0.00723215,   -0.999931,
//              0.672205 ,  -0.740274 ,  0.0116067 ;   //dalian
	Eigen::Vector3d c1;
//    c1 << 0.9988942870564429, -0.02090994236007995, -0.04210674051144664; //cv compute f
	//c1 << 0.012314168371418967, -0.010701346358933224, -0.03601467825925361; //open mvg
    c1 << 947.823, 1016.18, 54.5612; //world pnp huzhou
//    c1 << 970.431, 1017.26,  15.493;   //dalian
	Eigen::Matrix3d K1;
    K1 << 3.7516261809333914e+003, 0., 2.0321668349416393e+003, 0.,
        3.7075460080696898e+003, 1.4757973862050546e+003, 0., 0., 1.; //HUZHOU
//    K1 <<3.6807406030598072e+003, 0., 2.0575809009145160e+003, 0.,
//            3.6316348410018245e+003, 1.5323338632609398e+003, 0., 0., 1.;  //dalian

    std::vector<double> distort1{ -8.9706636166420994e-002, 2.2503707166823739e-001, -8.2725767124464909e-003 ,
                                  -2.7655590432857850e-003, 1.1201754344822839e-004 }; //huzhou
//    std::vector<double> distort1{ -1.0203338374479372e-001, 2.1315285352229257e-001,1.6984551476590334e-002,
//                                  -3.0681132919279703e-004, 1.1225988223251450e-005};  //dalian


	Eigen::Matrix3d R2;
    R2 << 0.9993811778822265, -0.01078803104811036, 0.03347954122801314,
    0.01167160140057348, 0.9995857163480499, -0.02630911237008151,
    -0.03318184768031093, 0.02668359156973558, 0.9990930641964541; //cv compute f
	//R2 << 0.99930614288737, -0.012159130892335009, 0.03520494742945117,
	//	0.013310805482134959, 0.9993776830569874, -0.03266602309838801,
	//	-0.03478584834360263, 0.03311196375316205, 0.9988461055695344;  //open mvg
    R2 <<    0.79421,   -0.60764, 0.00190083,
             0.0225506,  0.0263483,  -0.999398,
              0.607224,   0.793775,  0.0346288; //world pnp huzhou
//    R2 <<  -0.762675,   -0.646221,  -0.0269323,
//            0.0197389,   0.0183655,   -0.999636,
//              0.64648,    -0.76293, -0.00125122;   //dalian

	Eigen::Vector3d c2;
//    c2 << 0.9995034302815931, -0.00769056588222664, -0.03055729130239231; //cv compute F
    c2 << 951.094, 1014.74, 54.7585;  //huzhou
	//c2 << 0.029342534293331978, -0.009267891856573398, -0.030697318565627349; //open mvg
    //c2 << 966.922, 1014.19, 17.0275; //dalian


	Eigen::Matrix3d K2;
    K2 << 3.7604405383652875e+003, 0., 2.0415438075955501e+003, 0.,
        3.7160779331549415e+003, 1.4845370511494227e+003, 0., 0., 1.; //huzhou
//    K2 << 3.6902009246169255e+003, 0., 2.0795207000036444e+003, 0.,
//                3.6414589667424284e+003, 1.5276047177099535e+003, 0., 0., 1. ; //dalian

    std::vector<double> distort2{ -8.9249034584200818e-002, 2.4976814960048890e-001,
       -2.3008498513017160e-003, 1.4399452868398138e-003,
       -4.7341695340444627e-002 };  //huzhou
//    std::vector<double> distort2{ -1.0623327573034792e-001, 2.2930636719207204e-001,1.6003332119788841e-003
//                                  -1.0970447984199379e-003, 3.9033060894630466e-004}; //dalian

    //init reconstruction
    std::vector<string> vec_imageFile{filePath0,filePath1,filePath2};
    std::vector<string> vec_featureFilePath{ featureFile0, featureFile1, featureFile2};
//    std::vector<string> vec_imageFile{filePath0,filePath1};
//    std::vector<string> vec_featureFilePath{ featureFile0, featureFile1};
    std::vector<Eigen::Matrix3d> vec_R{R0,R1,R2};
    std::vector<Eigen::Vector3d> vec_c{c0,c1,c2};
//    std::vector<Eigen::Matrix3d> vec_R{R0,R1};
//    std::vector<Eigen::Vector3d> vec_c{c0,c1};
    std::vector<Eigen::Matrix<double,3,4>> vec_P(vec_R.size());
    for(int i =0;i<vec_P.size();i++)
    {
        vec_P[i].block(0,0,3,3) = vec_R[i];
        vec_P[i].block(0,3,3,1) = -vec_R[i] * vec_c[i];
    }
    std::vector<Eigen::Matrix3d> vec_KMatrix{K0,K1,K2};
    std::vector<std::vector<double>> vec_distort{distort0,distort1,distort2};
    //std::vector<std::vector<double>> vec_distort{distort0,distort1};
    std::vector<string> vec_matchesFile{matchFile_lm,matchFile_lr,matchFile_mr};
    //std::vector<string> vec_matchesFile{matchFile_lm,matchFile_lr};
    std::vector<string> vec_linePointsFile{lineDataFile0, lineDataFile1, lineDataFile2};
//    std::vector<string> vec_linePointsFile{lineDataFile0, lineDataFile1};

//    std::vector<string> vec_linePointsFile{lineDataFile0, lineDataFile1};
//    Test_reconstruct_with_KnownPoses("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/jiajie",vec_imageFile,vec_featureFilePath,vec_P,vec_KMatrix,vec_distort);
//    return 0;
    Reconstruction reconstruction;
    if(!reconstruction.Init(vec_imageFile,vec_P,vec_KMatrix,vec_distort,vec_featureFilePath,vec_matchesFile,vec_linePointsFile))
            return -1;
//    if(!reconstruction.Init(vec_imageFile,vec_P,vec_KMatrix,vec_distort,std::vector<string>(),std::vector<string>(),vec_linePointsFile))
//        return -1;

    //-----------------------

	//neighborsID 为 viewvector中的index
	//ADD neighbors informations
//	Eigen::Matrix3d F1, F2, F3;
//	view0.neighborsID.emplace_back(1);
////    F1 << -1.036742949829026e-09, -1.728890788495411e-07, 0.0004844618069332896,
////    1.169605847868869e-07, 7.469060661146156e-08, 0.009948625752603912,
////    -0.0003722534903969066, -0.01006154179399343, 1;   //l->m cv compute F
//	//Eigen::Matrix3d R_relative = R1 * R0.transpose();
//	//Eigen::Vector3d t_relative = view1.t - R_relative * view0.t;
//	//F1 = view1.K_Inv.transpose() * toAntisymmetricMatrix(t_relative) * R_relative * view0.K_Inv; // 1->m openmvg F
//	//F1 /= F1(2, 2);
//    F1 <<  1.56085e-09,  1.90842e-07, -0.000659478,
//           -2.36214e-07,  6.59792e-08 ,   0.0100108,
//             0.00070144,   -0.0101018 ,           1;   //WORLD pnp  huzhou
////    F1 << -1.31688e-08, -9.39254e-08 ,  0.00240304,
////          1.3876e-07,  4.92075e-08,  -0.00937679,
////         -0.00249231,   0.00899161,            1;   //dalian
//	view0.F_Me2Neighbors.emplace_back(F1);
//	view0.neighborsID.emplace_back(2);
////    F2 << -1.182109417917045e-09, -1.341513336584367e-07, 0.0001809402313467379,
////    6.417165844798604e-08, 3.756743735613564e-08, 0.007672951659434601,
////    -0.0001444273635979876, -0.00764985144987862, 1; //cv compute F
//	/*R_relative = R2 * R0.transpose();
//	t_relative = view2.t - R_relative * view0.t;
//	F2 = view2.K_Inv.transpose() * toAntisymmetricMatrix(t_relative) * R_relative * view0.K_Inv;
//	F2 /= F2(2, 2);*/
//    F2 <<  3.09045e-09,  2.19503e-07, -0.000424799,
//           -2.71411e-07,  6.02334e-08,   0.00661553,
//            0.000390681,  -0.00671668,            1;   //world pnp huzhou
////    F2 <<   3.2867e-09, -1.13409e-07, -0.000105977,
////            8.79383e-08,  2.02944e-08,  -0.00929817,
////            -0.00019418,   0.00918697 ,           1; //dalian
//	view0.F_Me2Neighbors.emplace_back(F2);   //l->r

//	view1.neighborsID.emplace_back(0);
//	view1.F_Me2Neighbors.emplace_back(F1.transpose()); //m->l
//	view1.neighborsID.emplace_back(2);
////	F3 << 1.884761703884511e-08, 1.343253350740144e-06, -0.00268961709437078,
////		-1.394848637065609e-06, 8.755146818652295e-08, 0.01535013030723664,
////		0.002509511883667218, -0.01565728246774467, 1;
//    F3 <<  1.01735e-08,  1.15511e-06, -0.000953044,
//           -1.22126e-06,  5.79408e-08 ,   0.0200092,
//            0.000774415,   -0.020152,            1;   //world pnp huzhou
////    F3 << 2.92439e-08, 4.99717e-08,  0.00309925,
////          4.66116e-08, 3.12012e-08,  0.00899255,
////          -0.00315183, -0.00939096,           1;  //dalian
//	view1.F_Me2Neighbors.emplace_back(F3); //m->r

//	view2.neighborsID.emplace_back(0);
//	view2.F_Me2Neighbors.emplace_back(F2.transpose()); //r->l
//	view2.neighborsID.emplace_back(1);
//	view2.F_Me2Neighbors.emplace_back(F3.transpose()); //r->m

	//通过已知三个世界相机的c和分解出的c求解一个到世界系的变换 
//    Eigen::Vector3d cw0, cw1, cw2;
//    cw0 << 945.115, 1024.913, 52.990;
//    cw1 << 949.300, 1021.623, 52.989;
//    cw2 << 951.961, 1019.501, 52.985;
//    std::vector<Eigen::Vector3d> pts_wordC{ cw0,cw1,cw2 };
//    std::vector<Eigen::Vector3d> pts_needAlignedC{ c0,c1,c2 };
//    Eigen::Matrix4d T = CalTransformToWorldWhenKnowTranslationDir(pts_wordC, pts_needAlignedC);
//    view0.center = pts_needAlignedC[0];
//    view0.updateParam();
//    view1.center = pts_needAlignedC[1];
//    view1.updateParam();
//    view2.center = pts_needAlignedC[2]; //按照local三个点的坐标进行缩放之后的点，不是local2world，没有*T呢还
//    view2.updateParam();

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

	/*test all kind structure.*/	
//    std::string file_l_Points = pointDataFile + "/001_l_Undistort.txt";
//    std::string file_m_Points = pointDataFile + "/001_m_Undistort.txt";
//    std::string file_r_Points = pointDataFile + "/001_r_Undistort.txt";
//    uint linePointSize1 = view0.pointTracks.size();
//    uint linePointSize2 = view1.pointTracks.size();
//    uint linePointSize3 = view2.pointTracks.size();
//    std::vector<Eigen::Vector2d> points_l = getPoint(file_l_Points, ",");
//    for (uint i = 0; i < points_l.size(); i++)
//    {
//        PointTrack pt(points_l[i], -1, i + linePointSize1);
//        view0.pointTracks.emplace_back(pt);
//    }
//    std::vector<Eigen::Vector2d> points_m = getPoint(file_m_Points, ",");
//    for (uint i = 0; i < points_m.size(); i++)
//    {
//        PointTrack pt(points_m[i], -1, i + linePointSize2);
//        view1.pointTracks.emplace_back(pt);
//    }
//    std::vector<Eigen::Vector2d> points_r = getPoint(file_r_Points, ",");
//    for (uint i = 0; i < points_r.size(); i++)
//    {
//        PointTrack pt(points_r[i], -1, i + linePointSize3);
//        view2.pointTracks.emplace_back(pt);
//    }
//    std::string matchFile_lm = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/data-total/data-total/001_l_Undistort.txt_001_m_Undistort.txt";
//    std::string matchFile_lr = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/data-total/data-total/001_l_Undistort.txt_001_r_Undistort.txt";
//    std::string matchFile_mr = "/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/data-total/data-total/001_m_Undistort.txt_001_r_Undistort.txt";
//    Matches matches_lm = getMatchesFromFile(matchFile_lm, 0, 1, linePointSize1, linePointSize2);
//    Matches matches_lr = getMatchesFromFile(matchFile_lr, 0, 2, linePointSize1, linePointSize3);
//    Matches matches_mr = getMatchesFromFile(matchFile_mr, 1, 2, linePointSize2, linePointSize3);
//    //geometric filter:
//    GeometricFilter(matches_lm, F1, points_l, points_m, linePointSize1, linePointSize2);
//    GeometricFilter(matches_lr, F2, points_l, points_r, linePointSize1, linePointSize3);
//    GeometricFilter(matches_mr, F3, points_m, points_r, linePointSize2, linePointSize3);
//    //----------------
//    std::vector<Matches> view0Matches{ matches_lm ,matches_lr };
//    //std::vector<Matches> view0Matches{ matches_lr };
//    std::vector<Matches> view1Matches{ matches_mr };
	
//	std::vector<View> views{ view0,view1,view2 };
//	Reconstruction reconstruction(views);
//    reconstruction.map_viewID_Matches[0] = view0Matches;
//    reconstruction.map_viewID_Matches[1] = view1Matches;
    //reconstruction.findTrackInOneView(0);
    reconstruction.findTracksInViews();   //20200222 dalian
    //reconstruction.findLineTrackInOneView_WithKnownMatch(0); //20200226 huzhou

	COUTENDL("begin sTRUCTURE....");
	reconstruction.GetStructure();
    reconstruction.SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/test_orig_fixRT.ply");
	//reconstruction.ExpandStructure();
    //reconstruction.SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/data-total/data-total/test_mvg_expand.ply");
	//reconstruction.ApplyTransform(T);
    //reconstruction.SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/data-total/data-total/test_mvg_orig_transfrom.ply");
    //reconstruction.SavePLYFile_Mesh("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/data-total/data-total/test1_Mesh.ply");

    //reconstruction.ExpandStructure_4();
    //reconstruction.ExpandStructure_6();
    //reconstruction.SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/jiajie/dalian/test_orig_fixRT_expand4.ply");

    //reconstruction.ExpandStructure_UpOneMM();
    //reconstruction.SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/jiajie/dalian/test_orig_fixRT_expand_after4.ply");
    //reconstruction.SavePLYFile_Mesh_UpOneMM("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/jiajie/dalian/test_orig_fixRT_expand_rectangle.ply");



    return 0;
}


