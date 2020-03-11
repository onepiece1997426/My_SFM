#pragma once
#include "Reconstruction.h"

#define MAX_POINTLINE_DST  0.005
#define MAX_DEPTH 500
#define MIN_DEPTH 5
#define MAX_REPROJ_ERR 20


//from 的某个点 映射到 to的对极线
Eigen::Matrix3d calFundamentalMatrix(View from, View to)
{
	//R*from+t=to 
	Eigen::Matrix3d Rtf = to.rotation * from.rotation.transpose();
	Eigen::Vector3d t_tf = to.t - Rtf * from.t;	
	Eigen::Matrix3d E_tf = toAntisymmetricMatrix(t_tf) * Rtf;
	Eigen::Matrix3d F = to.K_Inv.transpose() * E_tf * from.K_Inv;
	return F ;
}

inline Eigen::Matrix3d toAntisymmetricMatrix(Eigen::Vector3d t)
{
	Eigen::Matrix3d _t;
	_t << 0, -t(2), t(1),
		t(2), 0, -t(0),
		-t(1), t(0), 0;
	return _t;
}


//vec_featureFilePath with vec_matchesFile : featurePoints with known matchesFile
//vec_linesPointsFile linesPointsFile,fine Correlate by F
bool Reconstruction::Init(std::vector<string> vec_imageFile, std::vector<Eigen::Matrix<double,3,4>> vec_PMatrix, std::vector<Eigen::Matrix3d> vec_KMatrix,std::vector<std::vector<double>> vec_distort,std::vector<string> vec_featureFilePath, std::vector<string> vec_matchesFile, std::vector<string> vec_linesPointsFile)
{
    //init views
    std::vector<View> views(vec_PMatrix.size());
//    for(int i=0;i<views.size();i++)
//    {
//        Eigen::Matrix3d R = vec_PMatrix[i].block(0,0,3,3);
//        Eigen::Vector3d t = vec_PMatrix[i].block(0,3,3,1);
//        views[i] = View(vec_imageFile[i],i,R, -R.transpose() * t,vec_KMatrix[i], vec_distort[i]);
//    }

    //init lines
    srand(time(0));
    //std::vector<uint> vec_linePointCnt(views.size()); //default value 0
    VID_Lines vid_lines;
    if(vec_linesPointsFile.size() > 0)
    {
        for(int i=0;i<views.size();i++)
        {
            std::vector<std::vector<Eigen::Vector2d>> vec_vec_linePoints;
            std::vector<Eigen::Vector2d> vec_linePoints;
            std::ifstream ifs(vec_linesPointsFile[i], std::ios::in);
            if (!ifs)
            {
                COUTENDL("Err in line file");
                return false;
            }
            std::vector<PointTrack> pointsTracks;
            std::string str;
            uint lineIdx = 0;
            while (std::getline(ifs, str))
            {
                if (str == "" || str == "\r")
                {
                    Eigen::Vector3d color;
                    color << std::ceil(random()%255+1), std::ceil(random()%255 +1), std::ceil(random()%255 + 1);
                    if(!debug_lineColor.count(lineIdx))
                        debug_lineColor[lineIdx] = (color);
                    lineIdx++;
                    vec_vec_linePoints.emplace_back(vec_linePoints);
                    vec_linePoints.clear();
                    continue;
                }

                size_t pos = str.find_first_of(LINE_SPLIT, 0);  //use LINE_SPLIT to split string
                if (pos == str.npos)
                    continue;
                double x = std::atof(str.substr(0, pos).c_str());
                double y = std::atof(str.substr(pos+1).c_str());
                Eigen::Vector2d X(x, y);
                vec_linePoints.emplace_back(X);
//                PointTrack pt(X,lineIdx);
//                pointsTracks.emplace_back(pt);
            }
            ifs.close();
//            views[i].pointTracks = pointsTracks;
//            vec_linePointCnt[i] = pointsTracks.size();
            vid_lines[i] = vec_vec_linePoints;
        }
    }

    //init featurePOints
    std::vector<std::vector<Eigen::Vector2d>> vec_vec_featurePoints(vec_featureFilePath.size());
    for(uint i=0;i< vec_featureFilePath.size();i++)
    {
        vec_vec_featurePoints[i] = getPoint(vec_featureFilePath[i],POINT_SPLIT);
//        for(uint k=0;k<vec_featurePoints[i].size();k++)
//        {
//            PointTrack pt(vec_featurePoints[i][k], -1, k + vec_linePointCnt[i]);
//            views[i].pointTracks.emplace_back(pt);
//        }
    }


    //relative correspondences
    int matchesIdx = 0;
    std::vector<Matches> vec_I_Matches;
    for(int i=0;i< views.size();i++)
    {
        for(int j= i + 1;j < views.size();j++)
        {
            //init neighbor
//            views[i].neighborsID.emplace_back(j);

            //init F
//            Eigen::Matrix3d R_relative = views[j].rotation * (views[i].rotation.transpose());
//            Eigen::Vector3d t_relative = views[j].t - R_relative * views[i].t;
//            Eigen::Matrix3d F = views[j].K_Inv.transpose() * toAntisymmetricMatrix(t_relative) * R_relative * views[i].K_Inv;
//            F /= F(2,2);
//            views[i].F_Me2Neighbors.emplace_back(F);


            if(vec_featureFilePath.size() == 0)
                continue;
            Matches matches_ij(i,j);
            if(vec_matchesFile.size() > 0)
            {
                matches_ij = getMatchesFromFile(vec_matchesFile[matchesIdx],i,j);
                matchesIdx++;
            }

            //GeometricFilter(matches_ij, F, vec_featurePoints[i], vec_featurePoints[j]);
            vec_I_Matches.emplace_back(matches_ij);
        }
        //->map_viewID_Matches[i] = vec_I_Matches;
    }


//    this->Views = views;

    return Init(vec_imageFile,vec_PMatrix,vec_KMatrix,vec_distort,vec_vec_featurePoints,vec_I_Matches,vid_lines);
}

//same as before, but input  from vector<vector>
bool Reconstruction::Init(std::vector<string> vec_imageFile, std::vector<Eigen::Matrix<double,3,4>> vec_PMatrix, std::vector<Eigen::Matrix3d> vec_KMatrix,std::vector<std::vector<double>> vec_distort,std::vector<std::vector<Eigen::Vector2d>> vec_vec_featurePoints, std::vector<Matches> vec_matches, VID_Lines mp_vec_vec_linesPoints)
{
    //init views
    std::vector<View> views;
    if(this->Views.empty())
        views =  std::vector<View>(vec_PMatrix.size());
    else
        views = this->Views;
    for(int i=0;i<views.size();i++)
    {
        Eigen::Matrix3d R = vec_PMatrix[i].block(0,0,3,3);
        Eigen::Vector3d t = vec_PMatrix[i].block(0,3,3,1);
        views[i] = View(vec_imageFile[i],i,R, -R.transpose() * t,vec_KMatrix[i], vec_distort[i]);
    }

    //init lines
    srand(time(0));
    std::vector<uint> vec_linePointCnt(views.size()); //default value 0
    if(mp_vec_vec_linesPoints.size() > 0)
    {
        for(VID_Lines::iterator iter = mp_vec_vec_linesPoints.begin();iter != mp_vec_vec_linesPoints.end();iter++)
        {
            std::vector<PointTrack> pointsTracks;
            uint viewID = iter->first;
            std::vector<std::vector<Eigen::Vector2d>> vec_vec_linePoints = iter->second;  //这张图片中所有的线
            for(uint id=0;id<vec_vec_linePoints.size();id++)  //不同的lineid
            {
                for(uint j=0;j<vec_vec_linePoints[id].size();j++)  //这个ID上所有的点
                {
                    Eigen::Vector2d X(vec_vec_linePoints[id][j]);
                    PointTrack pt(X,id);
                    pointsTracks.emplace_back(pt);
                }
            }
            views[viewID].pointTracks = pointsTracks;
            vec_linePointCnt[viewID] = pointsTracks.size();
        }
    }

    //inti featurePOints
    for(uint i=0;i< vec_vec_featurePoints.size();i++)
    {
        for(uint k=0;k<vec_vec_featurePoints[i].size();k++)
        {
            PointTrack pt(vec_vec_featurePoints[i][k], -1, k + vec_linePointCnt[i]);
            views[i].pointTracks.emplace_back(pt);
        }
    }


    //relative correspondences
    int matchesIdx = 0;
    for(int i=0;i< views.size();i++)
    {
        //std::vector<Matches> vec_I_Matches;
        for(int j= i + 1;j < views.size();j++)
        {
            //init neighbor
            views[i].neighborsID.emplace_back(j);

            //init F
            Eigen::Matrix3d R_relative = views[j].rotation * (views[i].rotation.transpose());
            Eigen::Vector3d t_relative = views[j].t - R_relative * views[i].t;
            Eigen::Matrix3d F = views[j].K_Inv.transpose() * toAntisymmetricMatrix(t_relative) * R_relative * views[i].K_Inv;
            F /= F(2,2);
            views[i].F_Me2Neighbors[views[j].viewID] = (F);
        }
    }

    //matches
    for(uint i=0;i<vec_matches.size();i++)
    {
        Matches& matches_ij = vec_matches[i];
        uint hostID = matches_ij.viewIdx1;
        uint neighborID = matches_ij.viewIdx2;
        Eigen::Matrix3d F = views[hostID].F_Me2Neighbors[neighborID];
        if(vec_matches.size() == 0)
            continue;
        GeometricFilter(matches_ij, F, vec_vec_featurePoints[hostID], vec_vec_featurePoints[neighborID]);
        Matches match(hostID,neighborID);
        uint offsetHost = vec_linePointCnt[hostID];
        uint offsetNeighbor = vec_linePointCnt[neighborID];
        for(std::map<uint,uint>::iterator iter = matches_ij.map_featIdx.begin();iter != matches_ij.map_featIdx.end();iter++)
        {
            match.map_featIdx[iter->first + offsetHost] = iter->second + offsetNeighbor;
        }
        if(this->map_viewID_Matches.count(hostID))
        {
            this->map_viewID_Matches.at(hostID).emplace_back(match);
        }
        else
        {
            std::vector<Matches> vec_m{match};
            this->map_viewID_Matches[hostID] = vec_m;
        }
    }


    this->Views = views;

    return true;
}

//add neighbor to host(neighbor must be in views)
bool Reconstruction::UpdateNeighbor(uint hostID,uint neighborID)
{
    this->Views[hostID].neighborsID.emplace_back(neighborID);
    return true;
}


bool Reconstruction::UpdatePointTrackFromMatch(std::map<uint, std::vector<Eigen::Vector3d> > map_VID_vecFeatures, std::vector<Matches> vec_matches)
{

    for(std::map<uint, std::vector<Eigen::Vector3d> >::iterator iter = map_VID_vecFeatures.begin();iter != map_VID_vecFeatures.end();iter++)
    {
        for(uint i=0;i<iter->second.size();i++)
        {
            Eigen::Vector2d feat = Eigen::Vector2d(iter->second[i].x(),iter->second[i].y());
            PointTrack pt(feat,-1,Views[iter->first].pointTracks.size());
            Views[iter->first].pointTracks.emplace_back(pt);
        }
    }

    for(uint i=0;i<vec_matches.size();i++)
    {
        uint hostID = vec_matches[i].viewIdx1;
        uint neighborID = vec_matches[i].viewIdx2;

        if ( std::find(Views[hostID].neighborsID.begin(), Views[hostID].neighborsID.end(), neighborID) == Views[hostID].neighborsID.end() )
        {
            UpdateNeighbor(hostID,neighborID);
        }

        if(this->map_viewID_Matches.count(hostID))
            this->map_viewID_Matches[hostID].emplace_back(vec_matches[i]);
        else
        {
            vector<Matches> vec_matches_{vec_matches[i]};
            this->map_viewID_Matches[hostID] = vec_matches_;
        }
    }

    return true;
}




//所有图的track
void Reconstruction::findTracksInViews()
{
    omp_init_lock(&omp_lock);
	for (uint i = 0; i < Views.size() - 1; i++)
	{
        std::cout << "begin find corr lines in " << i << std::endl;
		findTrackInOneView(i);
	}
    omp_destroy_lock(&omp_lock); //销毁互斥器
}

//一张图的track
void Reconstruction::findTrackInOneView(uint ViewID)
{
	View& view = this->Views[ViewID];
	/*std::vector<Eigen::Matrix3d> F_neighbor2this;
	for (uint i = 0; i < view.neighborsID.size(); i++)
	{
		View neighbor = Views[view.neighborsID[i]];
		F_neighbor2this.push_back(calFundamentalMatrix(neighbor, view));
	}*/

	//test 30
    //uint cnt = 0;
    uint interval = 1;

    #pragma omp for
    for (uint pointID = 0; pointID < view.pointTracks.size(); pointID += interval )
	{			
		if (pointID > view.pointTracks.size())
			break;
        if(view.pointTracks[pointID].lineIdx < 0)
        {
            interval = 1;
        }
        else
        {
            interval = 1;
        }
		
		//每个point对应一个Tracks
		if (findTrack(ViewID, view.F_Me2Neighbors, view.pointTracks[pointID]))
		{
			//COUTENDL(view.pointTracks[pointID].lineIdx << "th line, " << pointID << " point has no point in neighbor");
            omp_set_lock(&omp_lock);
            trackCnt++;  //structures idx
            omp_unset_lock(&omp_lock);
		}
		
//		cnt++;
//		if (cnt > 100)
//		{
//			cnt = 0;
//			printf("\r %.3f", (double)pointID / view.pointTracks.size());
//		}
	}
	printf("\n");
	std::cout << " finish find corres lines" << std::endl;
	
	/*test*/
	/*
    cv::Mat test = cv::imread("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/out_corres.jpg");
	test.setTo(0);
	for (uint i = 0; i < view.pointTracks.size(); i++)
	{
		uint col = view.pointTracks[i].point.x();
		uint row = view.pointTracks[i].point.y();

		test.ptr<uchar>(row)[3 * col] = 0;
		test.ptr<uchar>(row)[3 * col + 1] = 0;
		test.ptr<uchar>(row)[3 * col + 2] = 255;

		if (view.pointTracks[i].tracks.size() > 0)
		{
			uint other_col = view.pointTracks[i].tracks[0].position.x();
			uint other_row = view.pointTracks[i].tracks[0].position.y();
			test.ptr<uchar>(other_row)[3 * other_col] = 255;
			test.ptr<uchar>(other_row)[3 * other_col + 1] = 255;
			test.ptr<uchar>(other_row)[3 * other_col + 2] = 255;
		}
		
	}
    cv::imwrite("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/out_corres.jpg", test);
	*/

}


void Reconstruction::findLineTrackInOneView_WithKnownMatch(uint ViewID)
{
    View& view = this->Views[ViewID];
    for (uint pointID = 0; pointID < view.pointTracks.size(); pointID++)
    {
        if (pointID > view.pointTracks.size())
            break;
        //每个point对应一个Tracks
        PointTrack& fromPointTrack = view.pointTracks[pointID];
        int lineId = fromPointTrack.lineIdx;
        Structure structure(Eigen::Vector3d::Ones());  //直接 一个track对应一个structure  至少得两个点
        structure.obs.emplace_back(ViewID, fromPointTrack.point, fromPointTrack.featIdx);
        structure.lineIdx = lineId;
        fromPointTrack.trackID = trackCnt;
        for (uint i = 0; i < view.neighborsID.size(); i++)
        {
            uint neighborId = view.neighborsID[i];
            Views[neighborId].pointTracks[pointID].trackID = trackCnt;
            structure.obs.emplace_back(neighborId, Views[neighborId].pointTracks[pointID].point, Views[neighborId].pointTracks[pointID].featIdx);
        }

        vec_structures.emplace_back(structure);
        trackCnt++;  //structures idx
    }
    printf("\n");
    std::cout << " finish find corres lines by known matches" << std::endl;

}


void testEippolalline(Eigen::Vector3d epipolarLine);


//一组track
bool Reconstruction::findTrack(uint ViewID, std::map<uint, Eigen::Matrix3d> F_neighbor2this, PointTrack& fromPointTrack)
{
	//if (fromPointTrack.trackID > 0)
	//{
	//	return false; //这个点已经有了对应的track了
	//}

	View view = Views[ViewID];
	Eigen::Vector3d x(fromPointTrack.point.x(), fromPointTrack.point.y(), 1);

    //旧版本
    {
	//for (uint i = 0; i < view.neighborsID.size(); i++)
	//{
	//	Eigen::Matrix3d F = F_neighbor2this[i];
	//	Eigen::Vector3d epipolarLine = F * x;
	//	uint neighborId = view.neighborsID[i];
	//	uint lineId = fromPointTrack.lineIdx;
	//	Track selfTrackToNeighbor(ViewID, fromPointTrack.point, trackCnt);
	//	
	//	/*test*/
	//	//x << 1633, 945, 1;
	//	F << -3.255013967077802e-13, 1.937330888047119e-09, -3.227454484534114e-06,
	//		-9.006384206573653e-10, -1.403213647902775e-09, -0.0001867188722359071,
	//		1.378189959106581e-06, 0.0001886452384087879, -0.01839238892919742;
	//	epipolarLine = F * x;
	//	/*testEippolalline(epipolarLine);
	//	Eigen::Vector3d x2(1454, 846, 1);
	//	COUTENDL(x2.transpose() * epipolarLine);*/
	//	double minDst = 9999;
	//	int minDstIdx = -1;
	//	for (uint j = 0; j < Views[neighborId].pointTracks.size(); j++)  //遍历所有的线点
	//	{			
	//		//属于不同图片下的同一条线才进行判定
	//		if (Views[neighborId].pointTracks[j].lineIdx == lineId)
	//		{
	//			Eigen::Vector3d x(Views[neighborId].pointTracks[j].point.x(), Views[neighborId].pointTracks[j].point.y(), 1.0f);
	//			double dst = std::abs(x.transpose() * epipolarLine);
	//			if (dst < MAX_POINTLINE_DST && dst < minDst)
	//			{
	//				minDst = dst;
	//				minDstIdx = j;   //neighbor对应的点的index
	//			}
	//		}
	//	}
	//	if (minDstIdx > 0)
	//	{
	//		Track t(neighborId, Views[neighborId].pointTracks[minDstIdx].point, trackCnt);
	//		fromPointTrack.tracks.emplace_back(t);	
	//		fromPointTrack.trackID = trackCnt; //这句得加上 不然出BUG emplaceBack最后一个参数trackID没赋值进去。。
	//		Views[neighborId].pointTracks[minDstIdx].trackID = trackCnt;
	//		Views[neighborId].pointTracks[minDstIdx].tracks.emplace_back(selfTrackToNeighbor);
	//	}
	//}	
	//if (fromPointTrack.tracks.size() > 0)
	//	return true;
	//else
	//	return false;
    }


    int lineId = fromPointTrack.lineIdx;
	Structure structure(Eigen::Vector3d::Ones());  //直接 一个track对应一个structure  至少得两个点
	structure.obs.emplace_back(ViewID, fromPointTrack.point, fromPointTrack.featIdx);
	structure.lineIdx = lineId;
	if (lineId >= 0)
	{
		//对线来说，前面找到的后面也能找到
		if (fromPointTrack.trackID >= 0)
		{
			return false; //这个点已经有了对应的track了
		}
		//通过对极约束找线
		for (uint i = 0; i < view.neighborsID.size(); i++)
		{
			if (view.neighborsID[i] < ViewID)   //只会在id比自己大的neighbor中找（小的neighbor都已经先找过了）
				continue;
            uint neighborId = view.neighborsID[i];
            Eigen::Matrix3d F = F_neighbor2this[neighborId];
			Eigen::Vector3d epipolarLine = F * x;

			/*test*/
			//x << 1633, 945, 1;
			///*F << -7.386202811949957e-10, -1.085529734493583e-07, 0.0003779755959086861,
			//	7.390186733841628e-08, 6.78692511334697e-08, 0.009940129693771615,
			//	-0.00031341315781519, -0.01004092592294054, 1;*/  //l->m
			//	/*F << 1.701859298628866e-09, 7.739500238132806e-09, -0.0001969476749336674,
			//		-7.413446245958539e-08, 3.523620291614568e-08, 0.007524715309168872,
			//		0.0001968570437233252, -0.007502674980432, 1;*///l->r
			//epipolarLine = F * x;
			//testEippolalline(epipolarLine);
			//Eigen::Vector3d x2(1454, 846, 1);   //m
			////Eigen::Vector3d x2(1378, 813, 1); //r
			//COUTENDL(x2.transpose() * epipolarLine);

 			double minDst = 9999;
			int minDstIdx = -1;
            //在一个Neighbor中通过对级约束找对应点
			for (uint j = 0; j < Views[neighborId].pointTracks.size(); j++)  //遍历所有的线点
			{
				//属于不同图片下的同一条线才进行判定
                if (Views[neighborId].pointTracks[j].lineIdx == lineId)
				{
                    Eigen::Vector3d x_(Views[neighborId].pointTracks[j].point.x(), Views[neighborId].pointTracks[j].point.y(), 1.0f);
                    double dst = std::abs(x_.transpose() * epipolarLine);
					if (dst < MAX_POINTLINE_DST && dst < minDst)
					{
						minDst = dst;
						minDstIdx = j;   //neighbor对应的点的index
					}
				}
			}

            if (minDstIdx > 0)  //找到合理的对应点
			{
				fromPointTrack.trackID = trackCnt; //这句得加上 。。
				Views[neighborId].pointTracks[minDstIdx].trackID = trackCnt;
				structure.obs.emplace_back(neighborId, Views[neighborId].pointTracks[minDstIdx].point);
			}
		}
		if (structure.obs.size() > 1)   //至少有一个邻居
		{
            omp_set_lock(&omp_lock);
			vec_structures.emplace_back(structure);
            omp_unset_lock(&omp_lock);
			return true;
		}
		else
			return false;

	}
	else
	{
		bool newTrack = (fromPointTrack.trackID < 0);

		//通过已经计算的匹配点找对应点		
        if(!this->map_viewID_Matches.count(ViewID))
            return false;
        std::vector<Matches> vec_matches = this->map_viewID_Matches.at(ViewID);  //和其neighbor对应的matches
		for (uint i = 0; i < view.neighborsID.size(); i++)
		{
			uint neighborId = view.neighborsID[i];
			if (neighborId < ViewID)
				continue;
			
			//找到对应neighbor的matches
			for (uint j = 0; j < vec_matches.size(); j++)
			{
				if (vec_matches[j].viewIdx2 != neighborId || vec_matches[j].viewIdx1 != ViewID)
					continue;
				if (!vec_matches[j].map_featIdx.count(fromPointTrack.featIdx))
				{
					//std::cout << "err feat idx" << std::endl;  //这个特征点没有匹配
					continue;
				}
					
				uint neihborPointIdx = vec_matches[j].map_featIdx.at(fromPointTrack.featIdx);
				PointTrack& pt = Views[neighborId].pointTracks[neihborPointIdx];
				if (pt.trackID > 0)
					continue;
				Eigen::Vector2d p = pt.point;
				structure.obs.emplace_back(neighborId, p,pt.featIdx);
				if (newTrack)
				{
					fromPointTrack.trackID = trackCnt; //这句得加上 。。
				}
				else
				{
					vec_structures[fromPointTrack.trackID].obs.emplace_back(neighborId, p, pt.featIdx);
				}
				pt.trackID = fromPointTrack.trackID;
				break;
			}

		}		
		if ((structure.obs.size() > 1) && newTrack)   //至少有一个邻居
		{
            omp_set_lock(&omp_lock);
			vec_structures.emplace_back(structure);
            omp_unset_lock(&omp_lock);
			return true;
		}
		else 
		{
			return false;
		}
	}
}

void Reconstruction::FindCorelatedPointByF(double& minDst, int& minDstIdx, Eigen::Matrix3d F, uint neighborId)
{


}

//求3D点
void Reconstruction::GetStructure()
{	
	//利用解析解做的，但是一般精确解不合理。。。在图片多了的时候，后面可以改为用Ceres做非线性优化
 	uint invalidDepth = 0;
	std::map<uint, Eigen::Matrix<double, 3, 4>> map_PMatrix;
	std::map<uint, Eigen::Matrix<double, 3, 4>> map_RtMatrix;

	for (uint viewIdx = 0; viewIdx < Views.size(); viewIdx++)
	{		
		Eigen::Matrix<double, 3, 4> P;
		P.block(0, 0, 3, 3) = Views[viewIdx].rotation;
		P.block(0, 3, 3, 1) = Views[viewIdx].t;
		map_RtMatrix[viewIdx] = P;
		map_PMatrix[viewIdx] = Views[viewIdx].K * P;
	}


    //cal lines when known pose
    CalculateStructure_Init_DLT(map_RtMatrix,true);
    SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/test_orig_fixRT_DLT.ply");
    CalculateStructure_Ceres(map_RtMatrix, map_PMatrix);

    this->CalculateStructure_FirstXThenRt(map_RtMatrix, map_PMatrix);
    //this->CalculateStructure_FirstXThenR_FixCenter(map_RtMatrix,map_PMatrix);  //fix center

	this->EraseInvalidStructure();
	COUTENDL("total structure: " << vec_structures.size());
}

void Reconstruction::SavePLYFile_PointCloud(std::string filePath)
{
	ofstream ofs(filePath);
	if (!ofs)
	{
		COUTENDL("err in create ply!");
		return;
	}
	else
	{
		ofs << "ply " << endl << "format ascii 1.0" << endl;
		//ofs << "element vertex " << this->structures.size() + this->Views.size() << endl;//old
		ofs << "element vertex " << this->vec_structures.size() + this->Views.size() << endl;
		ofs << "property float x" << endl << "property float y" << endl << "property float z" << endl;
		ofs << "property uchar blue                   { start of vertex color }"
			<< endl << "property uchar green"
			<< endl << "property uchar red" << endl;
		ofs << "end_header" << endl;
		
		/*old*/
		/*for (std::map<uint, Structure>::iterator iter = this->structures.begin(); iter != this->structures.end(); iter++)
		{
			ofs << iter->second.position.x() << " " << iter->second.position.y() << " " << iter->second.position.z()
				<< " " << iter->second.colors.x() << " " << iter->second.colors.y() << " " << iter->second.colors.z() << endl;

		}*/
		
		for (uint i = 0; i < vec_structures.size(); i++)
		{
			if(vec_structures[i].validStructure)
				ofs << vec_structures[i].position.x() << " " << vec_structures[i].position.y() << " " << vec_structures[i].position.z()
					<< " " << vec_structures[i].colors.x() << " " << vec_structures[i].colors.y() << " " << vec_structures[i].colors.z() << endl;
			else
				ofs << vec_structures[i].position.x() << " " << vec_structures[i].position.y() << " " << vec_structures[i].position.z()
				<< " " << 0 << " " << 0 << " " << 255 << endl;
		}
	

		for (uint i = 0; i < this->Views.size(); i++)
		{
			ofs << this->Views[i].center.x() << " " << this->Views[i].center.y() << " " << this->Views[i].center.z()
				<< " 0 255 0" << endl;
		}

	}
	ofs.close();
	ofs.flush();
	COUTENDL("FINISH SAVE PLY");
}

void Reconstruction::SavePLYFile_Mesh(std::string filePath)
{
	ofstream ofs(filePath);
	if (!ofs)
	{
		COUTENDL("err in create ply!");
		return;
	}
	else
	{
		ofs << "ply " << endl << "format ascii 1.0" << endl;
		ofs << "element vertex " << this->vec_structures.size() << endl;
		ofs << "property float x" << endl << "property float y" << endl << "property float z" << endl;
		ofs << "property uchar blue                   { start of vertex color }"
			<< endl << "property uchar green"
			<< endl << "property uchar red" << endl;
		ofs << "element face " << (this->structureCntBeforeExpand - 4) * 8 << endl;
		ofs << "property list uchar int vertex_index" << endl;
		ofs << "end_header" << endl;

		/*old*/
		/*for (std::map<uint, Structure>::iterator iter = this->structures.begin();iter != this->structures.end();iter++)
		{
			ofs << iter->second.position.x() << " " << iter->second.position.y() << " " << iter->second.position.z()
				<< " " << iter->second.colors.x() << " " << iter->second.colors.y() << " " << iter->second.colors.z() << endl;
		}*/
			
		for (uint i = 0; i < vec_structures.size(); i++)
		{
			ofs << vec_structures[i].position.x() << " " << vec_structures[i].position.y() << " " << vec_structures[i].position.z()
				<< " " << vec_structures[i].colors.x() << " " << vec_structures[i].colors.y() << " " << vec_structures[i].colors.z() << endl;

		}

		uint preMeshLineIdx = 0;
		uint cnt = 0;
		for (uint i = 0; i < structureCntBeforeExpand; i++)
		{
			if (this->vec_structures[i + 1].lineIdx != preMeshLineIdx)
			{
				preMeshLineIdx = this->vec_structures[i + 1].lineIdx;
				continue;
			}
			
			//一个面
			for (int k = 0; k < 7; k++)
			{			
				ofs << "4 "
					<< structureCntBeforeExpand + 8 * i + k + 1 << " " << structureCntBeforeExpand + 8 * (i + 1) + k + 1 << " " << structureCntBeforeExpand + 8 * (i + 1) + k
					<< " " << structureCntBeforeExpand + 8 * i + k << endl;
			}
			ofs << "4 "
				<< structureCntBeforeExpand + 8 * i << " " << structureCntBeforeExpand + 8 * (i + 1) << " " << structureCntBeforeExpand + 8 * (i + 1) + 7
				<< " " << structureCntBeforeExpand + 8 * i + 7 << endl;
			cnt += 8;
		}
		/*
		while (cnt < vec_structures.size())
		{
			if (this->vec_structures[cnt + 1].lineIdx != preMeshLineIdx)
			{
				preMeshLineIdx = this->vec_structures[cnt + 1].lineIdx;
				continue;
			}

			for (int k = 0; k < 7; k++)
			{
				ofs << "4 "
					<<  cnt + k + 1 << " " <<  (cnt + 1) + k + 1 << " " <<  (cnt + 1) + k
					<< " " << cnt + k << endl;
			}
			ofs << "4 "
				<< cnt << " " << (cnt + 1) << " " << (cnt + 1) + 7
				<< " " << cnt + 7 << endl;

			cnt += 8;
		}
		COUTENDL(cnt);
		*/
	}
	ofs.close();
	ofs.flush();
	COUTENDL("FINISH SAVE PLY");

}


void Reconstruction::SavePLYFile_Mesh_Rectangle(std::string filePath)
{
    ofstream ofs(filePath);
    if (!ofs)
    {
        COUTENDL("err in create ply!");
        return;
    }
    else
    {
        ofs << "ply " << endl << "format ascii 1.0" << endl;
        ofs << "element vertex " << this->vec_structures.size() << endl;
        ofs << "property float x" << endl << "property float y" << endl << "property float z" << endl;
        ofs << "property uchar blue                   { start of vertex color }"
            << endl << "property uchar green"
            << endl << "property uchar red" << endl;
        ofs << "element face " << (this->structureCntBeforeExpand - 14) * 4 << endl;
        ofs << "property list uchar int vertex_index" << endl;
        ofs << "end_header" << endl;

        /*old*/
        /*for (std::map<uint, Structure>::iterator iter = this->structures.begin();iter != this->structures.end();iter++)
        {
            ofs << iter->second.position.x() << " " << iter->second.position.y() << " " << iter->second.position.z()
                << " " << iter->second.colors.x() << " " << iter->second.colors.y() << " " << iter->second.colors.z() << endl;
        }*/

        for (uint i = 0; i < vec_structures.size(); i++)
        {
            ofs << vec_structures[i].position.x() << " " << vec_structures[i].position.y() << " " << vec_structures[i].position.z()
                << " " << vec_structures[i].colors.x() << " " << vec_structures[i].colors.y() << " " << vec_structures[i].colors.z() << endl;

        }

        uint preMeshLineIdx = 0;
        uint cnt = 0;
        for (uint i = 0; i < structureCntBeforeExpand; i++)
        {
            if (this->vec_structures[i + 1].lineIdx != preMeshLineIdx)
            {
                preMeshLineIdx = this->vec_structures[i + 1].lineIdx;
                continue;
                //break;
            }

            //一个面
            for (int k = 0; k < 3; k++)
            {
                ofs << "4 "
                    << structureCntBeforeExpand + 4 * i + k + 1 << " " << structureCntBeforeExpand + 4 * (i + 1) + k + 1 << " " << structureCntBeforeExpand + 4 * (i + 1) + k
                    << " " << structureCntBeforeExpand + 4 * i + k << endl;
            }
            ofs << "4 "
                << structureCntBeforeExpand + 4 * i << " " << structureCntBeforeExpand + 4 * (i + 1) << " " << structureCntBeforeExpand + 4 * (i + 1) + 3
                << " " << structureCntBeforeExpand + 4 * i + 3 << endl;
            cnt += 4;
        }
        /*
        while (cnt < vec_structures.size())
        {
            if (this->vec_structures[cnt + 1].lineIdx != preMeshLineIdx)
            {
                preMeshLineIdx = this->vec_structures[cnt + 1].lineIdx;
                continue;
            }

            for (int k = 0; k < 7; k++)
            {
                ofs << "4 "
                    <<  cnt + k + 1 << " " <<  (cnt + 1) + k + 1 << " " <<  (cnt + 1) + k
                    << " " << cnt + k << endl;
            }
            ofs << "4 "
                << cnt << " " << (cnt + 1) << " " << (cnt + 1) + 7
                << " " << cnt + 7 << endl;

            cnt += 8;
        }
        COUTENDL(cnt);
        */
    }
    ofs.close();
    ofs.flush();
    COUTENDL("FINISH SAVE PLY");

}

void Reconstruction::SavePLYFile_Mesh_UpOneMM(std::string filePath)
{
    ofstream ofs(filePath);
    if (!ofs)
    {
        COUTENDL("err in create ply!");
        return;
    }
    else
    {
        uint lineCnt = this->vec_structures[structureCntBeforeExpand-1].lineIdx+1;
        ofs << "ply " << endl << "format ascii 1.0" << endl;
        ofs << "element vertex " << this->vec_structures.size() << endl;
        ofs << "property float x" << endl << "property float y" << endl << "property float z" << endl;
        ofs << "property uchar blue                   { start of vertex color }"
            << endl << "property uchar green"
            << endl << "property uchar red" << endl;
        ofs << "element face " << (this->vec_structures.size() / 2 ) - 1 - lineCnt - 1 << endl;
        ofs << "property list uchar int vertex_index" << endl;
        ofs << "end_header" << endl;

        /*old*/
        /*for (std::map<uint, Structure>::iterator iter = this->structures.begin();iter != this->structures.end();iter++)
        {
            ofs << iter->second.position.x() << " " << iter->second.position.y() << " " << iter->second.position.z()
                << " " << iter->second.colors.x() << " " << iter->second.colors.y() << " " << iter->second.colors.z() << endl;
        }*/

        for (uint i = 0; i < vec_structures.size(); i++)
        {
            ofs << vec_structures[i].position.x() << " " << vec_structures[i].position.y() << " " << vec_structures[i].position.z()
                << " " << vec_structures[i].colors.x() << " " << vec_structures[i].colors.y() << " " << vec_structures[i].colors.z() << endl;

        }

        uint preMeshLineIdx = 0;
        for (uint i = 0; i < structureCntBeforeExpand; i++)
        {
            if (this->vec_structures[i + 1].lineIdx != preMeshLineIdx)
            {
                preMeshLineIdx = this->vec_structures[i + 1].lineIdx;
                continue;
                //break;
            }
            //一个面
            ofs << "4 "
                <<  i << " " << structureCntBeforeExpand + i  << " " << structureCntBeforeExpand + i + 1
                << " " << i + 1 << endl;
        }
    }
    ofs.close();
    ofs.flush();
    COUTENDL("FINISH SAVE PLY");

}


void Reconstruction::SaveTXTFile_PointCloud(std::string filePath)
{
	ofstream ofs(filePath);
	if (!ofs)
	{
		COUTENDL("err in create txt!");
		return;
	}
	else
	{
		for (uint i = 0; i < this->vec_structures.size(); i++)
		{
			ofs << this->vec_structures[i].position.x() << ' ' << this->vec_structures[i].position.y() << ' ' << this->vec_structures[i].position.z() << std::endl;
		}
	}
	ofs.close();
	COUTENDL("finish save");
}


//void Reconstruction::CalculateStructure_DLT(std::map<uint, Eigen::Matrix<double, 3, 4>> map_PMatrix)
//{
//	for (uint viewIdx = 0; viewIdx < Views.size(); viewIdx++)
//	{
//		cv::Mat img = cv::imread(Views[viewIdx].path);
//
//		//求在K1坐标系下的深度
//		/*slam14
//		//std::map<uint, Eigen::Matrix3d> map_R_relative;
//		//std::map<uint, Eigen::Vector3d> map_t_relative;
//		//min s x2^ A x1 + x2^ b s为深度 （slam14 154页）
//		//std::map<uint, Eigen::Matrix3d> map_Amartix;
//		//std::map<uint, Eigen::Vector3d> map_bVector;
//
//		//for (uint i = 0; i < Views[viewIdx].neighborsID.size(); i++)
//		//{
//		//	uint neighborID = Views[viewIdx].neighborsID[i];
//		//	Eigen::Matrix3d r_re = Views[viewIdx].rotation * Views[neighborID].rotation.transpose();
//		//	Eigen::Vector3d t_re = Views[viewIdx].t - r_re * Views[neighborID].t;
//		//	//map_R_relative[neighborID] = r_re;
//		//	//map_t_relative[neighborID] = t_re;
//		//	map_Amartix[neighborID] = Views[neighborID].K * r_re * Views[viewIdx].K_Inv;
//		//	map_bVector[neighborID] = Views[neighborID].K * t_re;
//
//		//}
//
//
//		//sx = KX
//		for (uint i = 0; i < Views[viewIdx].pointTracks.size(); i++)
//		{
//			PointTrack pointTrack = Views[viewIdx].pointTracks[i];
//			Eigen::Vector3d x(pointTrack.point.x(), pointTrack.point.y(), 1);
//			if (pointTrack.tracks.size() == 0)
//				continue;
//			for (uint k = 0; k < pointTrack.tracks.size(); k++)
//			{
//				Track track = pointTrack.tracks[k];
//				Eigen::Vector3d x2(track.position.x(), track.position.y(), 1);
//				Eigen::Matrix3d x2_ = toAntisymmetricMatrix(x2);
//				Eigen::Vector3d vec1 = x2_ * map_Amartix.at(track.viewIdx) * x;
//				Eigen::Vector3d vec2 = x2_ * map_bVector.at(track.viewIdx);
//				double s = -(vec2.dot(vec1)) / vec1.dot(vec1);
//
//				if (s <= 0)
//				{
//					invalidDepth++;
//					continue;
//				}
//
//				//反投影
//				  //没有图片宽高归一化的结果
//				Eigen::Vector3d Xrelative = s * Views[viewIdx].K_Inv * x;
//				Eigen::Vector3d X_world = Views[viewIdx].rotation.transpose() * (Xrelative - Views[viewIdx].t);
//
//
//				Structure structure(X_world);
//				structure.colors << img.ptr<uchar>((int)x.y())[3 * (int)x.x()], img.ptr<uchar>((int)x.y())[3 * (int)x.x() + 1], img.ptr<uchar>((int)x.y())[3 * (int)x.x() + 2];
//				structures.push_back(structure);
//			}
//
//			//深度值容易为负值 不好解 我的方程或者揭发可能有问题
//
//		}
//		*/
//		//教材 线形三角化 p218
//		
//		for (uint i = 0; i < Views[viewIdx].pointTracks.size(); i++)
//		{
//			PointTrack pointTrack = Views[viewIdx].pointTracks[i];
//			if (pointTrack.tracks.size() == 0)
//				continue;
//
//			Structure structure(Eigen::Vector3d::Ones());
//			Eigen::Vector3d x(pointTrack.point.x(), pointTrack.point.y(), 1);
//			if (pointTrack.tracks.size() == 0)
//				continue;
//			Eigen::Matrix<double, Eigen::Dynamic, 4> A;
//			A.resize(2 * pointTrack.tracks.size() + 2, 4);
//			Eigen::Matrix<double, 3, 4> P = map_PMatrix[viewIdx];
//			A.block(0, 0, 1, 4) = (x.x() * P.block(2, 0, 1, 4) - P.block(0, 0, 1, 4));
//			A.block(1, 0, 1, 4) = (x.y() * P.block(2, 0, 1, 4) - P.block(1, 0, 1, 4));
//			structure.obs.emplace_back(Observation(viewIdx, pointTrack.point));
//			for (uint k = 0; k < pointTrack.tracks.size(); k++)
//			{
//				Eigen::Matrix<double, 3, 4> P_neighbor = map_PMatrix[pointTrack.tracks[k].viewIdx];
//				Eigen::Vector2d x2 = pointTrack.tracks[k].position;
//				A.block(2 + 2 * k, 0, 1, 4) = (x2.x() * P_neighbor.block(2, 0, 1, 4) - P_neighbor.block(0, 0, 1, 4));
//				A.block(2 + 2 * k + 1, 0, 1, 4) = (x2.y() * P_neighbor.block(2, 0, 1, 4) - P_neighbor.block(0, 0, 1, 4));
//				structure.obs.emplace_back(Observation(pointTrack.tracks[k].viewIdx, x2));
//			}
//
//			//一个点解一个SVD...
//			Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinV);
//			Eigen::VectorXd vn = svd.matrixV().block(0, 3, 4, 1);
//			vn /= vn(3);
//
//			Eigen::Vector3d X_world(vn.x(), vn.y(), vn.z());
//			structure.position = X_world;
//			structure.colors << img.ptr<uchar>((int)x.y())[3 * (int)x.x()], img.ptr<uchar>((int)x.y())[3 * (int)x.x() + 1], img.ptr<uchar>((int)x.y())[3 * (int)x.x() + 2];
//			structure.positions_array[0] = vn.x();
//			structure.positions_array[1] = vn.y();
//			structure.positions_array[2] = vn.z();
//			structures[pointTrack.trackID] = (structure);			
//		}
//		
//	}
//
//}


//FixCenter:
//如果已知相机外参数或者是已经有一个光心的初值了，就可以不用在这里修改t的初值了（fixcenter置true）
void Reconstruction::CalculateStructure_Init_DLT(std::map<uint, Eigen::Matrix<double, 3, 4>>& map_RtMatrix, bool fixCenter)
{
	std::map<std::pair<uint, uint>, double> m_viewPair_lambda;
	std::vector<bool> vec_hasInitTranslation;

	for (uint i = 0; i < Views.size(); i++)
	{
		for (uint j = i + 1; j < Views.size(); j++)
		{
			m_viewPair_lambda[std::pair<uint, uint>(i, j)] = 1.0f;
		}
		vec_hasInitTranslation.emplace_back(false);
	}
	vec_hasInitTranslation[0] = vec_hasInitTranslation[1] = true;
	//CalTranslationByCeres(map_RtMatrix);

	ceres::Problem problem;
	ceres::LossFunction* lossFunc = USE_KERNEL_BA ? new ceres::HuberLoss(2.0f) : NULL;
	//给点云一个初值
	for (uint i = 0; i < vec_structures.size(); i++)
	{
		std::vector<Observation> obs = vec_structures[i].obs;
		Eigen::Matrix<double, 4, 4> A;
		//A.resize(2 * obs.size(), 4);	
		uint hostID = obs[0].viewID;
		uint neighborID = obs[1].viewID;

		Eigen::Vector3d X_init;
		Eigen::Matrix3d R_relative = map_RtMatrix.at(neighborID).block(0, 0, 3, 3) * (map_RtMatrix.at(hostID).block(0, 0, 3, 3).transpose());
		Eigen::Vector3d t_relative = map_RtMatrix.at(neighborID).block(0, 3, 3, 1) - R_relative * map_RtMatrix.at(hostID).block(0, 3, 3, 1);
		Eigen::Matrix<double, 3, 4> P1, P2;
		P1 = map_RtMatrix.at(hostID);
		P2 = map_RtMatrix.at(neighborID);
		/*P1.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity();
		P1.block(0, 3, 3, 1) = Eigen::Vector3d::Zero();*/
		P1 = this->Views[hostID].K * P1;
		/*P2.block(0, 0, 3, 3) = R_relative;
		P2.block(0, 3, 3, 1) = t_relative;*/
		P2 = this->Views[neighborID].K * P2;
		A.block(0, 0, 1, 4) =
			(obs[0].pixel.x()) * P1.row(2) - P1.row(0);
		A.block(1, 0, 1, 4) =
			(obs[0].pixel.y()) * P1.row(2) - P1.row(1);
		A.block(2, 0, 1, 4) =
			(obs[1].pixel.x()) * P2.row(2) - P2.row(0);
		A.block(3, 0, 1, 4) = 
			(obs[1].pixel.y()) * P2.row(2) - P2.row(1);
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinV);
		Eigen::Vector4d X_ = svd.matrixV().col(3);
		X_init << X_.x() / X_.w(), X_.y() / X_.w(), X_.z() / X_.w();
		if (X_init.z() < MIN_DEPTH || X_init.z() > MAX_DEPTH)
		{
			vec_structures[i].validStructure = false;
			continue;
		}

		Eigen::Vector3d X_init_local = map_RtMatrix.at(hostID).block(0, 0, 3, 3) * X_init + map_RtMatrix.at(hostID).block(0, 3, 3, 1);
		
		if(neighborID != 1 && !fixCenter)
		{ 
			ceres::CostFunction* cost_function =
			new ceres::AutoDiffCostFunction<Ceres_globalTranslation_KnownXInit, 2, 1>(new Ceres_globalTranslation_KnownXInit(X_init_local, obs[1].pixel, Views[neighborID].K, R_relative, t_relative));
			problem.AddResidualBlock(cost_function, lossFunc, &(m_viewPair_lambda.at(std::pair<uint, uint>(hostID, neighborID))));
			
		}
		
		//ceres::CostFunction* cost_function =
		//	new ceres::AutoDiffCostFunction<Ceres_globalTranslation_KnownXInit, 2, 1>(new Ceres_globalTranslation_KnownXInit(X_init_local, obs[1].pixel, Views[neighborID].K, R_relative, t_relative));
		//problem.AddResidualBlock(cost_function, NULL, &(m_viewPair_lambda.at(std::pair<uint, uint>(hostID, neighborID))));
		//if (hostID == 0 && neighborID == 1)
		//	problem.SetParameterBlockConstant(&(m_viewPair_lambda.at(std::pair<uint, uint>(hostID, neighborID))));  //通过这个函数可以固定某个地址中的变量不被更新



		for (uint j = 2; j < obs.size(); j++)
		{		
			neighborID = obs[j].viewID;
			R_relative = map_RtMatrix.at(neighborID).block(0, 0, 3, 3) * (map_RtMatrix.at(hostID).block(0, 0, 3, 3).transpose());
			t_relative = map_RtMatrix.at(neighborID).block(0, 3, 3, 1) - R_relative * map_RtMatrix.at(hostID).block(0, 3, 3, 1);
			Eigen::Vector3d x;
			x << obs[j].pixel.x(), obs[j].pixel.y(), 1;			
			

			if (neighborID != 1 && !fixCenter )
			{
				ceres::CostFunction* cost_function =
					new ceres::AutoDiffCostFunction<Ceres_globalTranslation_KnownXInit, 2, 1>(new Ceres_globalTranslation_KnownXInit(X_init_local, obs[j].pixel, Views[neighborID].K, R_relative, t_relative));
				problem.AddResidualBlock(cost_function, lossFunc, &(m_viewPair_lambda.at(std::pair<uint, uint>(hostID, neighborID))));
			}
			//ceres::CostFunction* cost_function =
			//	new ceres::AutoDiffCostFunction<Ceres_globalTranslation_KnownXInit, 2, 1>(new Ceres_globalTranslation_KnownXInit(X_init_local, obs[j].pixel, Views[neighborID].K, R_relative,t_relative));
			//problem.AddResidualBlock(cost_function, NULL, &(m_viewPair_lambda.at(std::pair<uint, uint>(hostID, neighborID))));
			//if (hostID == 0 && neighborID == 1)
			//	problem.SetParameterBlockConstant(&(m_viewPair_lambda.at(std::pair<uint, uint>(hostID, neighborID))));  //通过这个函数可以固定某个地址中的变量不被更新

			/*A.block(0, 0, 1, 4) =
				(obs[0].pixel.x()) * P1.row(2) - P1.row(0);
			A.block(1, 0, 1, 4) =
				(obs[0].pixel.y()) * P1.row(2) - P1.row(1);
			A.block(2, 0, 1, 4) =
				(obs[j].pixel.x()) * P2.row(2) - P2.row(0);
			A.block(3, 0, 1, 4) =
				(obs[j].pixel.y()) * P2.row(2) - P2.row(1);
			Eigen::JacobiSVD<Eigen::MatrixXd> svd_(A, Eigen::ComputeThinV);
			Eigen::Vector4d X_ = svd_.matrixV().col(3);*/
		}


		vec_structures[i].position = X_init;
		vec_structures[i].positions_array[0] = X_init.x();
		vec_structures[i].positions_array[1] = X_init.y();
		vec_structures[i].positions_array[2] = X_init.z();
	}
	
	if (!fixCenter)
	{
		ceres::Solver::Options options;
		options.linear_solver_type = ceres::DENSE_QR;
		options.max_num_iterations = 20;
		options.minimizer_progress_to_stdout = true;
		ceres::Solver::Summary summary;
		ceres::Solve(options, &problem, &summary);
		std::cout << "init translation: " << std::endl;
		std::cout << summary.BriefReport() << "\n";
	

		//别忘了更新t
		for (std::map<std::pair<uint, uint>, double>::iterator iter = m_viewPair_lambda.begin(); iter != m_viewPair_lambda.end(); iter++)
		{
			uint hostID = iter->first.first;
			uint neighborID = iter->first.second;
			if (neighborID == 1)
				continue;
			if (!vec_hasInitTranslation[neighborID])
			{
				Eigen::Matrix3d R_relative = map_RtMatrix.at(neighborID).block(0, 0, 3, 3) * (map_RtMatrix.at(hostID).block(0, 0, 3, 3).transpose());
				Eigen::Vector3d t_relative = map_RtMatrix.at(neighborID).block(0, 3, 3, 1) - R_relative * map_RtMatrix.at(hostID).block(0, 3, 3, 1);
			
				map_RtMatrix.at(neighborID).block(0, 3, 3, 1) = iter->second * t_relative + R_relative * map_RtMatrix[hostID].block(0, 3, 3, 1);
				this->Views[neighborID].updateParam(map_RtMatrix[neighborID].block(0, 3, 3, 1));
				vec_hasInitTranslation[neighborID] = true;
			}
		
		}
	}

	for (uint i = 0; i < vec_structures.size(); i++)
	{
		for (uint j = 0; j < vec_structures[i].obs.size(); j++)
		{
			bool valid = true;
			if (!isInSide(vec_structures[i].obs[j].viewID, vec_structures[i].position, vec_structures[i].obs[j].pixel, map_RtMatrix.at(vec_structures[i].obs[j].viewID)))
			{
				vec_structures[i].validStructure = false;
				valid = false;
				break;
			}
			if (!valid)
				break;
		}

	}


	COUTENDL("Finish Init Points");
}

void Reconstruction::CalculateStructure_Ceres(std::map<uint, Eigen::Matrix<double, 3, 4>>& map_RtMatrix, std::map<uint, Eigen::Matrix<double, 3, 4>>& map_PMatrix)
{
	ceres::Problem problem;
	ceres::LossFunction* lossFunc = USE_KERNEL_BA ? new ceres::HuberLoss(2.0f) : NULL;
	//旧版本
	//int structureID;
	//for (uint viewIdx = 0; viewIdx < Views.size(); viewIdx++)
	//{
	//	cv::Mat img = cv::imread(Views[viewIdx].path);
	//	for (uint i = 0; i < Views[viewIdx].pointTracks.size(); i++)
	//	{
	//		PointTrack& pointTrack = Views[viewIdx].pointTracks[i];
	//		structureID = pointTrack.trackID;
	//		if (structureID < 0)
	//			continue;
	//		Eigen::Vector2d x = pointTrack.point;
	//		Structure structure(Eigen::Vector3d::Ones());
	//		structure.colors << img.ptr<uchar>((int)x.y())[3 * (int)x.x()], img.ptr<uchar>((int)x.y())[3 * (int)x.x() + 1], img.ptr<uchar>((int)x.y())[3 * (int)x.x() + 2];
	//		structure.obs.emplace_back(Observation(viewIdx, x));
	//		this->structures[structureID] = (structure);
	//		pointTrack.hasBeenAddInStructure = true;
	//		ceres::CostFunction* cost_function =
	//			new ceres::AutoDiffCostFunction<Ceres_Triangulate, 2, 3>(new Ceres_Triangulate(x, map_PMatrix[viewIdx]));
	//		problem.AddResidualBlock(cost_function, NULL, this->structures.at(structureID).positions_array);
	//		for (uint k = 0; k < pointTrack.tracks.size(); k++)
	//		{
	//			Track& t = pointTrack.tracks[k];
	//			this->structures.at(structureID).obs.emplace_back(t.viewIdx, t.position);
	//			ceres::CostFunction* cost_function_track =
	//				new ceres::AutoDiffCostFunction<Ceres_Triangulate, 2, 3>(new Ceres_Triangulate(t.position, map_PMatrix[t.viewIdx]));
	//			problem.AddResidualBlock(cost_function_track, NULL, this->structures.at(structureID).positions_array);
	//			
	//			t.hasBeenReconstruct = true; //
	//			
	//		}
	//	}
	//}
	cv::Mat img = cv::imread(Views[0].path);
	map_PMatrix.clear();
	for (uint i = 0; i < Views.size(); i++)
	{
		map_PMatrix[Views[i].viewID] = Views[i].K * map_RtMatrix[Views[i].viewID];
	}


	for (uint structureIdx = 0; structureIdx < vec_structures.size(); structureIdx++)
	{
		if (!vec_structures[structureIdx].validStructure)
			continue;
		for (uint i = 0; i < vec_structures[structureIdx].obs.size(); i++)
		{
			Observation ob = vec_structures[structureIdx].obs[i];
            //if (ob.viewID == 0)
            //{
                //vec_structures[structureIdx].colors << img.ptr<uchar>((int)ob.pixel.y())[3 * (int)ob.pixel.x()], img.ptr<uchar>((int)ob.pixel.y())[3 * (int)ob.pixel.x() + 1], img.ptr<uchar>((int)ob.pixel.y())[3 * (int)ob.pixel.x() + 2];
                vec_structures[structureIdx].colors << 203,196,186;
                //vec_structures[structureIdx].colors = debug_lineColor[vec_structures[structureIdx].lineIdx];
            //}
			ceres::CostFunction* cost_function =
				new ceres::AutoDiffCostFunction<Ceres_Triangulate, 2, 3>(new Ceres_Triangulate(ob.pixel, map_PMatrix[ob.viewID]));
			problem.AddResidualBlock(cost_function, lossFunc, vec_structures[structureIdx].positions_array);
		}

	}

	ceres::Solver::Options options;
	options.max_num_iterations = 20;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	options.minimizer_progress_to_stdout = true;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.FullReport() << "\n";

	this->CheckValidStructure();

	structureCntBeforeExpand = this->vec_structures.size();

	COUTENDL("finish ceres BA");
}


//void Reconstruction::CalculateStructure_DLTandCeres(std::map<uint, Eigen::Matrix<double, 3, 4>> map_PMatrix)
//{
//	CalculateStructure_DLT(map_PMatrix);
//
//	//SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/init_ply.ply");
//
//	ceres::Problem problem;
//	for (std::map<uint, Structure>::iterator iter = this->structures.begin(); iter != this->structures.end(); iter++)
//	{
//		for (uint j = 0; j < iter->second.obs.size(); j++)
//		{
//			Observation ob = iter->second.obs[j];
//			ceres::CostFunction* cost_function =
//				new ceres::AutoDiffCostFunction<Ceres_Triangulate, 2, 3>(new Ceres_Triangulate(ob.pixel, map_PMatrix[ob.viewID]));
//			problem.AddResidualBlock(cost_function, NULL, iter->second.positions_array);
//		}
//	}
//
//	ceres::Solver::Options options;
//	options.max_num_iterations = 100;
//	options.linear_solver_type = ceres::DENSE_SCHUR;
//	options.minimizer_progress_to_stdout = true;
//	ceres::Solver::Summary summary;
//	ceres::Solve(options, &problem, &summary);
//	std::cout << summary.FullReport() << "\n";
//
//	for (std::map<uint, Structure>::iterator iter = this->structures.begin();iter != this->structures.end();iter++)
//	{
//		iter->second.position << iter->second.positions_array[0], iter->second.positions_array[1], iter->second.positions_array[2];
//	}
//
//	COUTENDL("finish ceres BA");
//}


void Reconstruction::CalculateStructure_Ceres_AdjustRt()
{
	ceres::Problem problem;
	ceres::LossFunction* lossFunc = USE_KERNEL_BA ? new ceres::HuberLoss(2.0f) : NULL;
	/*old*/
	//int structureID = 0;
	//for (uint viewIdx = 0; viewIdx < Views.size()-1; viewIdx++)
	//{
	//	cv::Mat img = cv::imread(Views[viewIdx].path);
	//	for (uint i = 0; i < Views[viewIdx].pointTracks.size(); i++)
	//	{
	//		PointTrack pointTrack = Views[viewIdx].pointTracks[i];
	//		structureID = pointTrack.trackID;
	//		if (structureID < 0)
	//			continue;
	//		Eigen::Vector2d x = pointTrack.point;
	//		Structure structure(Eigen::Vector3d::Ones());
	//		structure.colors << img.ptr<uchar>((int)x.y())[3 * (int)x.x()], img.ptr<uchar>((int)x.y())[3 * (int)x.x() + 1], img.ptr<uchar>((int)x.y())[3 * (int)x.x() + 2];
	//		structure.obs.emplace_back(Observation(viewIdx, x));
	//		structure.lineIdx = pointTrack.lineIdx;
	//		this->structures[structureID] = (structure);
	//		ceres::CostFunction * cost_function =
	//			new ceres::AutoDiffCostFunction<Ceres_Triangulate_AdjustRt, 2, 3, 3, 3>(new Ceres_Triangulate_AdjustRt(x, Views[viewIdx].K));
	//		problem.AddResidualBlock(cost_function, NULL, this->structures.at(structureID).positions_array, Views[viewIdx].rotation_array_aa, Views[viewIdx].translatrion_array);
	//		
	//		for (uint k = 0; k < pointTrack.tracks.size(); k++)
	//		{
	//			Track& track = pointTrack.tracks[k];
	//			this->structures.at(structureID).obs.emplace_back(track.viewIdx, track.position);
	//			ceres::CostFunction* cost_function_track =
	//				new ceres::AutoDiffCostFunction<Ceres_Triangulate_AdjustRt, 2, 3, 3, 3>(new Ceres_Triangulate_AdjustRt(track.position, Views[track.viewIdx].K));
	//			problem.AddResidualBlock(cost_function, NULL, this->structures.at(structureID).positions_array, Views[track.viewIdx].rotation_array_aa, Views[track.viewIdx].translatrion_array);
	//			track.hasBeenReconstruct = true; //
	//		}
	//	}
	//}
	cv::Mat img = cv::imread(Views[0].path);
	for (uint structureIdx = 0; structureIdx < vec_structures.size(); structureIdx++)
	{
		if (!vec_structures[structureIdx].validStructure)
			continue;

		for (uint i = 0; i < vec_structures[structureIdx].obs.size(); i++)
		{
			Observation ob = vec_structures[structureIdx].obs[i];
			if (ob.viewID == 0)
			{
				vec_structures[structureIdx].colors << img.ptr<uchar>((int)ob.pixel.y())[3 * (int)ob.pixel.x()], img.ptr<uchar>((int)ob.pixel.y())[3 * (int)ob.pixel.x() + 1], img.ptr<uchar>((int)ob.pixel.y())[3 * (int)ob.pixel.x() + 2];
			}
			ceres::CostFunction* cost_function_track =
				new ceres::AutoDiffCostFunction<Ceres_Triangulate_AdjustRt, 2, 3, 3, 3>(new Ceres_Triangulate_AdjustRt(ob.pixel, Views[ob.viewID].K));
			problem.AddResidualBlock(cost_function_track, lossFunc, vec_structures[structureIdx].positions_array, Views[ob.viewID].rotation_array_aa, Views[ob.viewID].translatrion_array);

		}
	}

	ceres::Solver::Options options;
	options.max_num_iterations = 20;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	options.minimizer_progress_to_stdout = true;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.FullReport() << "\n";

	for (uint structureIdx = 0; structureIdx < vec_structures.size(); structureIdx++)
	{
		if (!vec_structures[structureIdx].validStructure)
			continue;
		vec_structures[structureIdx].position << vec_structures[structureIdx].positions_array[0], vec_structures[structureIdx].positions_array[1], vec_structures[structureIdx].positions_array[2];
	}

    std::ofstream ofs("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/adjust_rotation.txt");
	for (uint i = 0; i < Views.size(); i++)
	{
		Eigen::AngleAxis<double> angelAxis;	
		angelAxis.angle() = std::sqrt(Views[i].rotation_array_aa[0] * Views[i].rotation_array_aa[0]  + Views[i].rotation_array_aa[1] * Views[i].rotation_array_aa[1] + Views[i].rotation_array_aa[2] * Views[i].rotation_array_aa[2]);
		Eigen::Vector3d axis(Views[i].rotation_array_aa[0], Views[i].rotation_array_aa[1], Views[i].rotation_array_aa[2]);
		axis.normalize();
		angelAxis.axis() = axis;
		Views[i].rotation = angelAxis.toRotationMatrix();

		Views[i].t = Eigen::Vector3d(Views[i].translatrion_array[0], Views[i].translatrion_array[1], Views[i].translatrion_array[2]);
		Views[i].center = -Views[i].rotation.transpose() * Views[i].t;

		ofs << Views[i].viewID << std::endl;
		ofs << Views[i].rotation << std::endl;
		ofs << "t: " << Views[i].t.transpose() << endl;
		ofs << "c: " << Views[i].center.transpose() << std::endl << std::endl;

	}
	ofs.close();

	this->CheckValidStructure();

	COUTENDL("finish ceres BA_rt");

}

void Reconstruction::CalculateStructure_Ceres_Adjust_t(std::map<uint, Eigen::Matrix<double, 3, 4>>& map_RtMatrix )
{
	ceres::Problem problem;
	ceres::LossFunction* lossFunc = USE_KERNEL_BA ? new ceres::HuberLoss(2.0f) : NULL;

	cv::Mat img = cv::imread(Views[0].path);
	for (uint structureIdx = 0; structureIdx < vec_structures.size(); structureIdx++)
	{
		if (!vec_structures[structureIdx].validStructure)
			continue;

		for (uint i = 0; i < vec_structures[structureIdx].obs.size(); i++)
		{
			Observation ob = vec_structures[structureIdx].obs[i];		
			if (ob.viewID <= 1)
				continue;
			ceres::CostFunction* cost_function_track =
				new ceres::AutoDiffCostFunction<Ceres_Triangulate_AdjustRt, 2, 3, 3, 3>(new Ceres_Triangulate_AdjustRt(ob.pixel, Views[ob.viewID].K));
			problem.AddResidualBlock(cost_function_track, lossFunc, vec_structures[structureIdx].positions_array, Views[ob.viewID].rotation_array_aa, Views[ob.viewID].translatrion_array);
			/*if (ob.viewID <= 1)
			{
				problem.SetParameterBlockConstant(Views[ob.viewID].translatrion_array);
			}*/
			problem.SetParameterBlockConstant(vec_structures[structureIdx].positions_array);
			problem.SetParameterBlockConstant(Views[ob.viewID].rotation_array_aa);
			
		}
	}

	ceres::Solver::Options options;
	options.max_num_iterations = 20;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	
	options.minimizer_progress_to_stdout = true;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.FullReport() << "\n";

	for (uint structureIdx = 0; structureIdx < vec_structures.size(); structureIdx++)
	{
		if (!vec_structures[structureIdx].validStructure)
			continue;
		vec_structures[structureIdx].position << vec_structures[structureIdx].positions_array[0], vec_structures[structureIdx].positions_array[1], vec_structures[structureIdx].positions_array[2];
	}

    std::ofstream ofs("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/adjust_translation.txt");
	for (uint i = 0; i < Views.size(); i++)
	{
		Eigen::AngleAxis<double> angelAxis;
		angelAxis.angle() = std::sqrt(Views[i].rotation_array_aa[0] * Views[i].rotation_array_aa[0] + Views[i].rotation_array_aa[1] * Views[i].rotation_array_aa[1] + Views[i].rotation_array_aa[2] * Views[i].rotation_array_aa[2]);
		Eigen::Vector3d axis(Views[i].rotation_array_aa[0], Views[i].rotation_array_aa[1], Views[i].rotation_array_aa[2]);
		if(axis.norm() != 0)
			axis.normalize();

		angelAxis.axis() = axis;
		Views[i].rotation = angelAxis.toRotationMatrix();

		Views[i].t = Eigen::Vector3d(Views[i].translatrion_array[0], Views[i].translatrion_array[1], Views[i].translatrion_array[2]);
		Views[i].center = -Views[i].rotation.transpose() * Views[i].t;

		ofs << Views[i].viewID << std::endl;
		ofs << Views[i].rotation << std::endl;
		ofs << "t: " << Views[i].t.transpose() << endl;
		ofs << "c: " << Views[i].center.transpose() << std::endl << std::endl;

		map_RtMatrix.at(Views[i].viewID).block(0, 3, 3, 1) = Views[i].t;

	}
	ofs.close();


	this->CheckValidStructure();
	COUTENDL("finish ceres BA_t");
}

void Reconstruction::CalculateStructure_Ceres_AdjustR_FixC()
{
	ceres::Problem problem;
	ceres::LossFunction* lossFunc = USE_KERNEL_BA ? new ceres::HuberLoss(2.0f) : NULL;
	cv::Mat img = cv::imread(Views[0].path);
	for (uint structureIdx = 0; structureIdx < vec_structures.size(); structureIdx++)
	{
		if (!vec_structures[structureIdx].validStructure)
			continue;
		for (uint i = 0; i < vec_structures[structureIdx].obs.size(); i++)
		{
			Observation ob = (vec_structures[structureIdx].obs[i]);
			ceres::CostFunction* cost_function_track =
				new ceres::AutoDiffCostFunction<Ceres_Triangulate_withScale_R, 2, 3, 3>(new Ceres_Triangulate_withScale_R(ob.pixel, Views[ob.viewID].K, 1, Views[ob.viewID].center));
			problem.AddResidualBlock(cost_function_track, lossFunc, vec_structures[structureIdx].positions_array, Views[ob.viewID].rotation_array_aa);
		}
	}

	ceres::Solver::Options options;
	options.max_num_iterations = 20;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	options.minimizer_progress_to_stdout = true;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.FullReport() << "\n";

	for (uint structureIdx = 0; structureIdx < vec_structures.size(); structureIdx++)
	{
		vec_structures[structureIdx].position << vec_structures[structureIdx].positions_array[0], vec_structures[structureIdx].positions_array[1], vec_structures[structureIdx].positions_array[2];
	}

    std::ofstream ofs("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/adjust_rotation.txt");
	for (uint i = 0; i < Views.size(); i++)
	{
		Eigen::AngleAxis<double> angelAxis;
		angelAxis.angle() = std::sqrt(Views[i].rotation_array_aa[0] * Views[i].rotation_array_aa[0] + Views[i].rotation_array_aa[1] * Views[i].rotation_array_aa[1] + Views[i].rotation_array_aa[2] * Views[i].rotation_array_aa[2]);
		Eigen::Vector3d axis(Views[i].rotation_array_aa[0], Views[i].rotation_array_aa[1], Views[i].rotation_array_aa[2]);
		axis.normalize();
		angelAxis.axis() = axis;
		Views[i].rotation = angelAxis.toRotationMatrix();
		Views[i].t = -Views[i].rotation * Views[i].center;
		//Views[i].t = Eigen::Vector3d(Views[i].translatrion_array[0], Views[i].translatrion_array[1], Views[i].translatrion_array[2]);
		//Views[i].center = -Views[i].rotation.transpose() * Views[i].t;

		ofs << Views[i].viewID << std::endl;
		ofs << Views[i].rotation << std::endl;
		ofs << "t: " << Views[i].t.transpose() << endl;
		ofs << "c: " << Views[i].center.transpose() << std::endl << std::endl;

	}
	ofs.close();

	this->CheckValidStructure();
	COUTENDL("finish ceres BA_fix c");

	structureCntBeforeExpand = this->vec_structures.size();
}

void Reconstruction::CalculateStructure_FirstXThenRt(std::map<uint, Eigen::Matrix<double, 3, 4>>& map_RtMatrix, std::map<uint, Eigen::Matrix<double, 3, 4>>& map_PMatrix)
{


    /*old*/
    {
	////first structure
	//ceres::Problem problem;
	//int structureID = 0;
	//for (uint viewIdx = 0; viewIdx < Views.size() - 1; viewIdx++)
	//{
	//	cv::Mat img = cv::imread(Views[viewIdx].path);
	//	for (uint i = 0; i < Views[viewIdx].pointTracks.size(); i++)
	//	{
	//		PointTrack pointTrack = Views[viewIdx].pointTracks[i];
	//		structureID = pointTrack.trackID;
	//		if (structureID < 0)
	//			continue;
	//		//debug
	//		/*if (viewIdx > 0)
	//		{
	//			COUTENDL(structureID);
	//		}*/
	//		Eigen::Vector2d x = pointTrack.point;
	//		//后面可以优化，现在相当于有重复添加点，或者是有重复操作
	//		/*if (structures.count(structureID))
	//		{
	//			structures.at(structureID).obs.emplace_back(Observation(viewIdx, x));
	//		}
	//		else
	//		{			
	//			Structure structure(Eigen::Vector3d::Ones());
	//			structure.colors << img.ptr<uchar>((int)x.y())[3 * (int)x.x()], img.ptr<uchar>((int)x.y())[3 * (int)x.x() + 1], img.ptr<uchar>((int)x.y())[3 * (int)x.x() + 2];
	//			structure.obs.emplace_back(Observation(viewIdx, x));
	//			this->structures[structureID] = (structure);
	//		}*/
	//		
	//		Structure structure(Eigen::Vector3d::Ones());
	//		structure.colors << img.ptr<uchar>((int)x.y())[3 * (int)x.x()], img.ptr<uchar>((int)x.y())[3 * (int)x.x() + 1], img.ptr<uchar>((int)x.y())[3 * (int)x.x() + 2];
	//		structure.obs.emplace_back(Observation(viewIdx, x));
	//		this->structures[structureID] = (structure);
	//		ceres::CostFunction * cost_function =
	//			new ceres::AutoDiffCostFunction<Ceres_Triangulate, 2, 3>(new Ceres_Triangulate(x, map_PMatrix[viewIdx]));
	//		problem.AddResidualBlock(cost_function, NULL, this->structures.at(structureID).positions_array);
	//		for (uint k = 0; k < pointTrack.tracks.size(); k++)
	//		{
	//			Track& t = pointTrack.tracks[k];
	//			this->structures.at(structureID).obs.emplace_back(t.viewIdx, t.position);
	//			ceres::CostFunction* cost_function_track =
	//				new ceres::AutoDiffCostFunction<Ceres_Triangulate, 2, 3>(new Ceres_Triangulate(t.position, map_PMatrix[t.viewIdx]));
	//			problem.AddResidualBlock(cost_function_track, NULL, this->structures.at(structureID).positions_array);
	//			t.hasBeenReconstruct = true; //
	//		}
	//	}
	//}
	//ceres::Solver::Options options;
	//options.max_num_iterations = 30;
	//options.linear_solver_type = ceres::DENSE_SCHUR;
	//options.minimizer_progress_to_stdout = true;
	//ceres::Solver::Summary summary;
	//ceres::Solve(options, &problem, &summary);
	//std::cout << summary.FullReport() << "\n";
	//for (std::map<uint, Structure>::iterator iter = this->structures.begin(); iter != this->structures.end(); iter++)
	//{
	//	iter->second.position << iter->second.positions_array[0], iter->second.positions_array[1], iter->second.positions_array[2];
	//}
	//COUTENDL("1. finish ceres BA_Structure, begin XRT BA");
	////------------------------------------
	//structureCntBeforeExpand = structures.size();
    //SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/test1_orig_temp.ply");
	//// second X RT
	//ceres::Problem problem_XRT;
	//for (std::map<uint, Structure>::iterator iter = this->structures.begin(); iter != this->structures.end(); iter++)
	//{
	//	std::vector<Observation> obs = iter->second.obs;
	//	for (uint i = 0; i < obs.size(); i++)
	//	{			
	//		ceres::CostFunction* cost_function_track =
	//			new ceres::AutoDiffCostFunction<Ceres_Triangulate_AdjustRt, 2, 3, 3, 3>(new Ceres_Triangulate_AdjustRt(obs[i].pixel, Views[obs[i].viewID].K));
	//		problem_XRT.AddResidualBlock(cost_function_track, NULL, iter->second.positions_array, Views[obs[i].viewID].rotation_array_aa, Views[obs[i].viewID].translatrion_array);
	//	}	
	//}
	//ceres::Solver::Options options_XRT;
	//options_XRT.max_num_iterations = 30;
	//options_XRT.linear_solver_type = ceres::DENSE_SCHUR;
	//options_XRT.minimizer_progress_to_stdout = true;
	//ceres::Solver::Summary summary_XRT;
	//ceres::Solve(options_XRT, &problem_XRT, &summary_XRT);
	//std::cout << summary_XRT.FullReport() << "\n";
	//for (std::map<uint, Structure>::iterator iter = this->structures.begin(); iter != this->structures.end(); iter++)
	//{
	//	iter->second.position << iter->second.positions_array[0], iter->second.positions_array[1], iter->second.positions_array[2];
	//}
    //std::ofstream ofs("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/adjust_rotation.txt");
	//for (uint i = 0; i < Views.size(); i++)
	//{
	//	Eigen::AngleAxis<double> angelAxis;
	//	angelAxis.angle() = std::sqrt(Views[i].rotation_array_aa[0] * Views[i].rotation_array_aa[0] + Views[i].rotation_array_aa[1] * Views[i].rotation_array_aa[1] + Views[i].rotation_array_aa[2] * Views[i].rotation_array_aa[2]);
	//	Eigen::Vector3d axis(Views[i].rotation_array_aa[0], Views[i].rotation_array_aa[1], Views[i].rotation_array_aa[2]);
	//	axis.normalize();
	//	angelAxis.axis() = axis;
	//	Views[i].rotation = angelAxis.toRotationMatrix();
	//	Views[i].t = Eigen::Vector3d(Views[i].translatrion_array[0], Views[i].translatrion_array[1], Views[i].translatrion_array[2]);
	//	Views[i].center = -Views[i].rotation * Views[i].t;
	//	ofs << Views[i].viewID << std::endl;
	//	ofs << Views[i].rotation << std::endl;
	//	ofs << "t: " << Views[i].t.transpose() << endl;
	//	ofs << "c: " << Views[i].center.transpose() << std::endl << std::endl;
	//}
	//ofs.close();
	//COUTENDL("2. finish ceres BA_XRT");
	//structureCntBeforeExpand = structures.size();
    //
    }

    COUTENDL("First begin cal Structure");
	CalculateStructure_Init_DLT(map_RtMatrix);
    SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/first_noRT_DLT.ply");
    CalculateStructure_Ceres_Adjust_t(map_RtMatrix);
    SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/first_t.ply");
    CalculateStructure_Ceres(map_RtMatrix, map_PMatrix);
    SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/first_noRT.ply");
	COUTENDL("Second begin cal Structure and Rt");
	CalculateStructure_Ceres_AdjustRt();
	COUTENDL("Finish Total BA");

}

void Reconstruction::CalculateStructure_FirstXThenR_FixCenter(std::map<uint, Eigen::Matrix<double, 3, 4>>& map_RtMatrix, std::map<uint, Eigen::Matrix<double, 3, 4>>& map_PMatrix)
{
	COUTENDL("First begin cal Structure");
	CalculateStructure_Init_DLT(map_RtMatrix,true);
    SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/first_noRT_DLT.ply");
	CalculateStructure_Ceres(map_RtMatrix, map_PMatrix);
    SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/first_noRT.ply");
	COUTENDL("Second begin cal Structure and R(fix c)");
	CalculateStructure_Ceres_AdjustR_FixC();
	COUTENDL("Finish Total BA");


}


static omp_lock_t lock;
void Reconstruction::LerpVoltaageLines(uint cntBetween2Point)
{
    //先对线按照lineIdx和z方向进行排序
    //先做变换到Views[0]坐标系下
    Eigen::Matrix3d R = Views[0].rotation;
    Eigen::Vector3d t = Views[0].t;
    for(uint i=0;i<vec_structures.size();i++)
    {
        vec_structures[i].position = R * vec_structures[i].position + t;
    }
    std::vector<uint> vec_linePointCnt;
    std::sort(vec_structures.begin(), vec_structures.end(), Structure::sortmethod_lineIdx);  //按照线的顺序排序
    uint preProcessLineIdx = 0;
    std::vector<Structure>::iterator lastItor = vec_structures.begin();
    uint linePointCnt = 0;
    for (std::vector<Structure>::iterator iter = vec_structures.begin(); ; iter++)
    {
        if (iter->lineIdx != preProcessLineIdx)
        {
            std::sort(lastItor, iter, Structure::sortmethod_z);  //按照线的顺序排序
            if(iter == vec_structures.end())
                break;
            lastItor = iter;
            preProcessLineIdx = iter->lineIdx;
            vec_linePointCnt.emplace_back(linePointCnt);
//            linePointCnt = 1;
        }
//        else
//        {
//            linePointCnt++;
//        }
    }

    //可以优化
//    uint newStructureSize=0;
//    std::vector<uint> vec_expandLinePointCnt(vec_linePointCnt.size());
//    for(uint i=0;i<vec_linePointCnt.size();i++)
//    {
//        vec_expandLinePointCnt[i] = vec_linePointCnt[i] + (vec_linePointCnt[i] - 1) * cntBetween2Point;
//        newStructureSize += vec_expandLinePointCnt[i];
//    }
//    std::vector<Structure> new_structure(newStructureSize);
//    uint currentLineIdx = 0;
//    uint totalProceedPointCnt = 0;
//    for(uint i=0;i<vec_expandLinePointCnt.size();i++)
//    {


//        totalProceedPointCnt += vec_expandLinePointCnt[i];
//    }

    omp_init_lock(&lock); // 初始化互斥锁

#pragma omp for
    for(uint i=0;i<vec_structures.size()-1;i++)
    {
        if(vec_structures[i+1].lineIdx != vec_structures[i].lineIdx)
            continue;
        for(uint j=1;j<cntBetween2Point;j++)
        {
            Eigen::Vector3d X_lerp = (((double)j) / (double)cntBetween2Point) * (vec_structures[i+1].position -  vec_structures[i].position) +  vec_structures[i].position;
            Structure structure(X_lerp);
            structure.obs = vec_structures[i].obs;
            structure.colors = vec_structures[i].colors;
            omp_set_lock(&lock); //获得互斥器
            vec_structures.emplace_back(structure);
            omp_unset_lock(&lock); //获得互斥器
        }

    }

    //再排序  (这方法效率比较低)
    preProcessLineIdx = 0;
    for (std::vector<Structure>::iterator iter = vec_structures.begin(); ; iter++)
    {
        if (iter->lineIdx != preProcessLineIdx)
        {
            std::sort(lastItor, iter, Structure::sortmethod_z);  //按照线的顺序排序
            if(iter == vec_structures.end())
                break;
            lastItor = iter;
            preProcessLineIdx = iter->lineIdx;
            vec_linePointCnt.emplace_back(linePointCnt);
//            linePointCnt = 1;
        }
//        else
//        {
//            linePointCnt++;
//        }
    }
    //变换回来
    for(uint i=0;i<vec_structures.size();i++)
    {
        vec_structures[i].position  = R.transpose() * (vec_structures[i].position - t);
    }

    omp_destroy_lock(&lock); //销毁互斥器

}

void Reconstruction::ExpandStructure()
{
	COUTENDL("Begin Expand");
	double r_size = 0.05f;
	uint origSize = vec_structures.size();
	uint num = structureCntBeforeExpand;
	
	std::sort(vec_structures.begin(), vec_structures.end(), Structure::sortmethod_lineIdx);  //按照线的顺序排序
	uint preProcessLineIdx = 0;
	std::vector<Structure>::iterator lastItor = vec_structures.begin();
    for (std::vector<Structure>::iterator iter = vec_structures.begin(); ; iter++)
	{
		if (iter->lineIdx != preProcessLineIdx)
		{
            std::sort(lastItor, iter, Structure::sortmethod_z);  //按照线的顺序排序
            if(iter == vec_structures.end())
                break;
			lastItor = iter;
			preProcessLineIdx = iter->lineIdx;
		}
	}

	for (uint i = 0; i < origSize; i++)
	{
		//在某个点的坐标系下生成八个点
		Eigen::Vector3d X_center = this->vec_structures[i].position;

		//深度都为正，Z方向都在相机Z方向
		for (uint j = 0; j < 8; j++)
		{
			Eigen::Vector3d direction(std::cos(j * 3.14f / 4), std::sin(j * 3.14f / 4), 0.0f);
			direction = Views[0].rotation * direction;
			Eigen::Vector3d X_j = X_center + r_size * direction;
			Structure structure(X_j);
			structure.colors = vec_structures[i].colors;
			structure.lineIdx = vec_structures[i].lineIdx;
			//obs?

			vec_structures.emplace_back(structure);
			num++;
		}

	}
	COUTENDL("Finish Expand");
}

void Reconstruction::ExpandStructure_4()
{
    COUTENDL("Begin Expand");
    //先做变换到Views[0]坐标系下
    Eigen::Matrix3d R = Views[0].rotation;
    Eigen::Vector3d t = Views[0].t;
    for(uint i=0;i<vec_structures.size();i++)
    {
        vec_structures[i].position = R * vec_structures[i].position + t;
    }

    std::vector<bool> vec_needExpandLine{true,false,true,false,true};
    double r_size = 0.225f; //400mm
    uint origSize = vec_structures.size();
    uint num = structureCntBeforeExpand;

    std::sort(vec_structures.begin(), vec_structures.end(), Structure::sortmethod_lineIdx);  //按照线的顺序排序
    uint preProcessLineIdx = 0;
    std::vector<Structure>::iterator lastItor = vec_structures.begin();
    for (std::vector<Structure>::iterator iter = vec_structures.begin(); ;iter++)
    {
        if (iter->lineIdx != preProcessLineIdx || iter == vec_structures.end())
        {
            std::sort(lastItor, iter, Structure::sortmethod_z);  //按照线的顺序排序
            if(iter == vec_structures.end())
                break;
            lastItor = iter;
            preProcessLineIdx = iter->lineIdx;
        }
    }

    uint lineIdxBeforeExpand = preProcessLineIdx;
    preProcessLineIdx = 0;
    uint baseLineIdx = lineIdxBeforeExpand + 1;
    for (uint i = 0; i < origSize; i++)
    {
        if(this->vec_structures[i].lineIdx < 0)
            continue;
        if(!vec_needExpandLine[this->vec_structures[i].lineIdx])
            continue;

        //在某个点的坐标系下生成6个点
        Eigen::Vector3d X_center = this->vec_structures[i].position;
        uint thisLineIdx = this->vec_structures[i].lineIdx;
        if(thisLineIdx != preProcessLineIdx)
        {
            baseLineIdx += 4;
            preProcessLineIdx = thisLineIdx;
        }
        //深度都为正，Z方向都在相机Z方向
        for (uint j = 0; j < 4; j++)
        {
            Eigen::Vector3d direction(std::cos(j * 3.14f / 2) + 3.14f/4, std::sin(j * 3.14f / 2) + 3.14f/4, 0.0f);
            Eigen::Vector3d X_j = X_center + r_size * direction;
            Structure structure(X_j);
            structure.colors = vec_structures[i].colors;
            structure.lineIdx = baseLineIdx + j;
            //obs?

            vec_structures.emplace_back(structure);
            num++;
        }

    }

    std::sort(vec_structures.begin(), vec_structures.end(), Structure::sortmethod_lineIdx);  //按照线的顺序排序
    preProcessLineIdx = 0;
    lastItor = vec_structures.begin();
    for (std::vector<Structure>::iterator iter = vec_structures.begin(); ;iter++)
    {
        if (iter->lineIdx != preProcessLineIdx || iter == vec_structures.end())
        {
            std::sort(lastItor, iter, Structure::sortmethod_z);  //按照线的顺序排序
            if(iter == vec_structures.end())
                break;
            lastItor = iter;
            preProcessLineIdx = iter->lineIdx;
        }
    }

    for(uint i=0;i<vec_structures.size();i++)
    {
        vec_structures[i].position  = R.transpose() * (vec_structures[i].position - t);
    }

    COUTENDL("Finish Expand 44444444444");
    this->structureCntBeforeExpand = this->vec_structures.size();
}


void Reconstruction::ExpandStructure_6()
{
    COUTENDL("Begin Expand");
    //先做变换到Views[0]坐标系下
    Eigen::Matrix3d R = Views[0].rotation;
    Eigen::Vector3d t = Views[0].t;
    for(uint i=0;i<vec_structures.size();i++)
    {
        vec_structures[i].position = R * vec_structures[i].position + t;
    }

    std::vector<bool> vec_needExpandLine{true,false,false,true};
    double r_size = 0.4f; //400mm
    uint origSize = vec_structures.size();
    uint num = structureCntBeforeExpand;

    std::sort(vec_structures.begin(), vec_structures.end(), Structure::sortmethod_lineIdx);  //按照线的顺序排序
    uint preProcessLineIdx = 0;
    std::vector<Structure>::iterator lastItor = vec_structures.begin();
    for (std::vector<Structure>::iterator iter = vec_structures.begin(); ;iter++)
    {
        if (iter->lineIdx != preProcessLineIdx || iter == vec_structures.end())
        {
            std::sort(lastItor, iter, Structure::sortmethod_z);  //按照线的顺序排序
            if(iter == vec_structures.end())
                break;
            lastItor = iter;
            preProcessLineIdx = iter->lineIdx;
        }
    }

    uint lineIdxBeforeExpand = preProcessLineIdx;
    preProcessLineIdx = 0;
    uint baseLineIdx = lineIdxBeforeExpand + 1;
    for (uint i = 0; i < origSize; i++)
    {
        if(this->vec_structures[i].lineIdx < 0)
            continue;
        if(!vec_needExpandLine[this->vec_structures[i].lineIdx])
            continue;

        //在某个点的坐标系下生成6个点
        Eigen::Vector3d X_center = this->vec_structures[i].position;
        uint thisLineIdx = this->vec_structures[i].lineIdx;
        if(thisLineIdx != preProcessLineIdx)
        {
            baseLineIdx += 6;
            preProcessLineIdx = thisLineIdx;
        }
        //深度都为正，Z方向都在相机Z方向
        for (uint j = 0; j < 6; j++)
        {
            Eigen::Vector3d direction(std::cos(j * 3.14f / 3), std::sin(j * 3.14f / 3), 0.0f);
            Eigen::Vector3d X_j = X_center + r_size * direction;
            Structure structure(X_j);
            structure.colors = vec_structures[i].colors;
            structure.lineIdx = baseLineIdx + j;
            //obs?

            vec_structures.emplace_back(structure);
            num++;
        }

    }

    std::sort(vec_structures.begin(), vec_structures.end(), Structure::sortmethod_lineIdx);  //按照线的顺序排序
//    preProcessLineIdx = 0;
//    Eigen::Vector3d color;
//    color << std::ceil(random()%255+1), std::ceil(random()%255 +1), std::ceil(random()%255 + 1);
//    for(uint i=0;i<vec_structures.size();i++)
//    {
//        if(vec_structures[i].lineIdx != preProcessLineIdx)
//        {
//            color << std::ceil(random()%255+1), std::ceil(random()%255 +1), std::ceil(random()%255 + 1);
//             preProcessLineIdx = vec_structures[i].lineIdx;
//        }
//        vec_structures[i].colors = color;
//    }
    preProcessLineIdx = 0;
    lastItor = vec_structures.begin();
    for (std::vector<Structure>::iterator iter = vec_structures.begin(); ;iter++)
    {
        if (iter->lineIdx != preProcessLineIdx || iter == vec_structures.end())
        {
            std::sort(lastItor, iter, Structure::sortmethod_z);  //按照线的顺序排序
            if(iter == vec_structures.end())
                break;
            lastItor = iter;
            preProcessLineIdx = iter->lineIdx;
        }
    }
    //SavePLYFile_PointCloud("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/jiajie/huzhou/test_orig_fixRT_expand_local.ply");

    for(uint i=0;i<vec_structures.size();i++)
    {
        vec_structures[i].position  = R.transpose() * (vec_structures[i].position - t);
    }

    COUTENDL("Finish Expand 66666");
    this->structureCntBeforeExpand = this->vec_structures.size();
}


void Reconstruction::ExpandStructure_Rectangle()
{
    COUTENDL("Begin Expand with known pos");
    double r_size = 0.01f;   //1mm
    uint origSize = vec_structures.size();
    uint num = structureCntBeforeExpand;

    //先做变换到Views[0]坐标系下
    Eigen::Matrix3d R = Views[0].rotation;
    Eigen::Vector3d t = Views[0].t;
    for(uint i=0;i<vec_structures.size();i++)
    {
        vec_structures[i].position = R * vec_structures[i].position + t;
    }

    //std::vector<Eigen::Vector3d> vec_colors;

    std::sort(vec_structures.begin(), vec_structures.end(), Structure::sortmethod_lineIdx);  //按照线的顺序排序
    uint preProcessLineIdx = 0;
    std::vector<Structure>::iterator lastItor = vec_structures.begin();
    for (std::vector<Structure>::iterator iter = vec_structures.begin(); ; iter++)
    {
        if (iter->lineIdx != preProcessLineIdx || iter == vec_structures.end())
        {
            std::sort(lastItor, iter, Structure::sortmethod_z);  //按照线的顺序排序
//            Eigen::Vector3d color;
//            color << std::ceil(random()%255+1), std::ceil(random()%255 +1), std::ceil(random()%255 + 1);
//            vec_colors.emplace_back(color);

            if(iter == vec_structures.end())
                break;
            lastItor = iter;
            preProcessLineIdx = iter->lineIdx;

        }
    }

    for (uint i = 0; i < origSize; i++)
    {
        //在某个点的坐标系下生成4个点
        Eigen::Vector3d X_center = this->vec_structures[i].position;

        //深度都为正，Z方向都在相机Z方向
        for (uint j = 0; j < 4; j++)
        {
            Eigen::Vector3d direction(std::cos((j * 3.14f / 2) + 3.14f / 4 ), std::sin( (j * 3.14f / 2) + 3.14f / 4), 0.0f);
            Eigen::Vector3d X_j = X_center + r_size * direction;
            Structure structure(X_j);
            structure.colors = vec_structures[i].colors;
            //structure.colors = vec_colors[vec_structures[i].lineIdx];
            structure.lineIdx = vec_structures[i].lineIdx;
            //obs?

            vec_structures.emplace_back(structure);
            num++;
        }

    }


    //变换回来
    for(uint i=0;i<vec_structures.size();i++)
    {
        vec_structures[i].position  = R.transpose() * (vec_structures[i].position - t);
    }

    COUTENDL("Finish Expand, size: " << vec_structures.size());

}

void Reconstruction::ExpandStructure_UpOneMM()
{
    COUTENDL("Begin Expand with known pos");
    double r_size = 0.01f;   //1mm
    uint origSize = vec_structures.size();
    uint num = structureCntBeforeExpand;

    //先做变换到Views[0]坐标系下
    Eigen::Matrix3d R = Views[0].rotation;
    Eigen::Vector3d t = Views[0].t;
    for(uint i=0;i<vec_structures.size();i++)
    {
        vec_structures[i].position = R * vec_structures[i].position + t;
    }

    //std::vector<Eigen::Vector3d> vec_colors;

    std::sort(vec_structures.begin(), vec_structures.end(), Structure::sortmethod_lineIdx);  //按照线的顺序排序
    uint preProcessLineIdx = 0;
    std::vector<Structure>::iterator lastItor = vec_structures.begin();
    for (std::vector<Structure>::iterator iter = vec_structures.begin(); ; iter++)
    {
        if (iter->lineIdx != preProcessLineIdx || iter == vec_structures.end())
        {
            std::sort(lastItor, iter, Structure::sortmethod_z);  //按照线的顺序排序
//            Eigen::Vector3d color;
//            color << std::ceil(random()%255+1), std::ceil(random()%255 +1), std::ceil(random()%255 + 1);
//            vec_colors.emplace_back(color);

            if(iter == vec_structures.end())
                break;
            lastItor = iter;
            preProcessLineIdx = iter->lineIdx;

        }
    }

    for (uint i = 0; i < origSize; i++)
    {
        //在某个点的坐标系下生成4个点
        Eigen::Vector3d X_center = this->vec_structures[i].position;

        //深度都为正，Z方向都在相机Z方向

        Eigen::Vector3d X_j = X_center + r_size * Eigen::Vector3d(0,1.0,0);
        Structure structure(X_j);
        structure.colors = vec_structures[i].colors;
        //structure.colors = vec_colors[vec_structures[i].lineIdx];
        structure.lineIdx = vec_structures[i].lineIdx;
        //obs?

        vec_structures.emplace_back(structure);
    }

    //变换回来
    for(uint i=0;i<vec_structures.size();i++)
    {
        vec_structures[i].position  = R.transpose() * (vec_structures[i].position - t);
    }

    COUTENDL("Finish Expand, size: " << vec_structures.size());

}

void Reconstruction::ApplyTransform(Eigen::Matrix4d T)
{
    ofstream ofs1("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/P_beforeTransform.txt");
	for (uint viewIdx = 0; viewIdx < Views.size(); viewIdx++)
	{
		Eigen::Matrix<double, 3, 4> P;
		P.block(0, 0, 3, 3) = Views[viewIdx].rotation;
		P.block(0, 3, 3, 1) = Views[viewIdx].t;
		ofs1 << "T" << viewIdx << ": " << std::endl << P << std::endl;
		P = Views[viewIdx].K * P;
		ofs1 << "P" << viewIdx << ": " << std::endl << P << std::endl;
	}
	ofs1.close();
	
	//apply transform
	Eigen::Matrix3d sR = T.block(0, 0, 3, 3);
	double s = std::pow(sR.determinant(), 0.33333333f);
	Eigen::Matrix3d R = sR / s;
	Eigen::Vector3d t = T.block(0, 3, 3, 1);
	for (uint i = 0; i < vec_structures.size(); i++)
	{
		vec_structures[i].position = sR * vec_structures[i].position + t;
		vec_structures[i].positions_array[0] = vec_structures[i].position(0);
		vec_structures[i].positions_array[1] = vec_structures[i].position(1);
		vec_structures[i].positions_array[2] = vec_structures[i].position(2);
	}
	for (uint i = 0; i < Views.size(); i++)
	{
		Views[i].center = sR * Views[i].center + t;
		Views[i].rotation = Views[i].rotation * R.transpose();
		Views[i].t = -Views[i].rotation * Views[i].center;
		Views[i].updateParam();
	}
	COUTENDL("Finish apply transform");

    ofstream ofs("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/P_afterTransform.txt");
	for (uint viewIdx = 0; viewIdx < Views.size(); viewIdx++)
	{
		Eigen::Matrix<double, 3, 4> P;
		P.block(0, 0, 3, 3) = Views[viewIdx].rotation;
		P.block(0, 3, 3, 1) = Views[viewIdx].t;
		ofs << "T" << viewIdx << ": " << std::endl << P << std::endl;
		P = Views[viewIdx].K * P;
		ofs << "P" << viewIdx << ": " << std::endl << P << std::endl;
	}
	ofs.close();
}

bool Reconstruction::isInSide(uint viewID, Eigen::Vector3d X, Eigen::Vector2d x, Eigen::Matrix<double, 3, 4> T)
{
	Eigen::Vector3d x_proj = Views[viewID].K * (T.block(0, 0, 3, 3) * X + T.block(0, 3, 3, 1));
	x_proj /= x_proj.z();
	
    double dst = (x_proj.block(0, 0, 2, 1) - x).norm();
    if (dst < MAX_REPROJ_ERR)
        return true;
    else
        return false;

//    if (x_proj.x() > 0 && x_proj.x() <= 2 * Views[viewID].K(0, 2)
//		&& x_proj.y() > 0 && x_proj.y() <= 2 * Views[viewID].K(1, 2)
//		)
//		return true;
//	else
//        return false;
}


//GLOBAL SFM paper functin (8)  不太一样 比那个简单些 
void Reconstruction::CalTranslationByCeres(std::map<uint, Eigen::Matrix<double, 3, 4>>& map_RtMatrix)
{
	ceres::Problem problem;
	std::map<std::pair<uint, uint>, double> m_lambda;

	for (uint i = 0; i < Views.size(); i++)
	{
		uint hostID = Views[i].viewID;
		for (uint j = i + 1; j < Views.size(); j++)
		{		
			uint neighborID = Views[j].viewID;
			m_lambda[std::pair<uint, uint>(hostID, neighborID)] = 1.0;
			Eigen::Matrix3d R_relative = map_RtMatrix.at(neighborID).block(0, 0, 3, 3) * (map_RtMatrix.at(hostID).block(0, 0, 3, 3).transpose());

			ceres::CostFunction* cost_function =
				new ceres::AutoDiffCostFunction<Ceres_globalTranslation, 3, 3, 3, 1>(new Ceres_globalTranslation(R_relative));
			problem.AddResidualBlock(cost_function, NULL, Views[i].translatrion_array,Views[j].translatrion_array, &(m_lambda.at(std::pair<uint, uint>(hostID, neighborID)) ));
			if (i == 0 && j == 1)
				problem.SetParameterBlockConstant(&(m_lambda.at(std::pair<uint, uint>(hostID, neighborID))));  //通过这个函数可以固定某个地址中的变量不被更新
		}
	}
	ceres::Solver::Options options;   
	options.linear_solver_type = ceres::DENSE_QR;   
	options.minimizer_progress_to_stdout = true;   
	ceres::Solver::Summary summary;   
	ceres::Solve(options, &problem, &summary);   
	std::cout << "init translation: " << std::endl;
	std::cout << summary.BriefReport() << "\n";
	COUTENDL("Finish init translation");
}

void Reconstruction::CheckValidStructure()
{
	for (uint structureIdx = 0; structureIdx < vec_structures.size(); structureIdx++)
	{
		if (!vec_structures[structureIdx].validStructure)
			continue;
		Eigen::Vector3d X;
		X << vec_structures[structureIdx].positions_array[0], vec_structures[structureIdx].positions_array[1], vec_structures[structureIdx].positions_array[2];
		for (uint i = 0; i < vec_structures[structureIdx].obs.size(); i++)
		{
			Eigen::Matrix<double, 3, 4> Rt;
			Rt.block(0, 0, 3, 3) = this->Views[vec_structures[structureIdx].obs[i].viewID].rotation;
			Rt.block(0, 3, 3, 1) = this->Views[vec_structures[structureIdx].obs[i].viewID].t;
			if (!isInSide(vec_structures[structureIdx].obs[i].viewID, X, vec_structures[structureIdx].obs[i].pixel, Rt))
			{
				vec_structures[structureIdx].validStructure = false;
				break;
			}
            Eigen::Vector3d X_cam = Rt.block(0, 0, 3, 3) * X + Rt.block(0, 3, 3, 1);
            if (X_cam.z() < MIN_DEPTH || X_cam.z() > MAX_DEPTH)
            {
                vec_structures[structureIdx].validStructure = false;
                continue;
            }
		}
		if (vec_structures[structureIdx].validStructure)
			vec_structures[structureIdx].position = X;
	}
}

void Reconstruction::EraseInvalidStructure()
{
	for (std::vector<Structure>::iterator iter = this->vec_structures.begin(); iter != this->vec_structures.end();)
	{
		if (!iter->validStructure)
		{
			iter = vec_structures.erase(iter);
			if (iter == vec_structures.end())
				break;
		}
		else
			iter++;
	}
	structureCntBeforeExpand = this->vec_structures.size();
}



std::string testImgPath = "E:\\ProjectsDuringPostgraduate\\HighVoltageLine\\images\\001_m_Undistort.png";
void testEippolalline(Eigen::Vector3d epipolarLine)
{
	/*test对极线*/
	cv::Mat out = cv::imread(testImgPath);
	cv::Point2i p1, p2;
	if (epipolarLine.x() == 0)
	{
		p1 = cv::Point2i(0, -epipolarLine.z() / epipolarLine.y());
		p2 = cv::Point2i(out.cols - 1, -epipolarLine.z() / epipolarLine.y());
		COUTENDL("k = inf");
	}
	else if (epipolarLine.y() == 0)
	{
		p1 = cv::Point2i(-epipolarLine.z()  / epipolarLine.x(), 0);
		p2 = cv::Point2i(-epipolarLine.z()  / epipolarLine.x(), out.rows - 1);
		COUTENDL("k = 0");
	}
	else
	{
		
		p1 = cv::Point2i(0, -epipolarLine.z()  / epipolarLine.y());
		p2 = cv::Point2i(out.cols - 1, -(epipolarLine.x() * (out.cols - 1) + epipolarLine.z()) / (epipolarLine.y()));
		
		COUTENDL("k = " << (-epipolarLine.x() / epipolarLine.y()));
	}
	COUTENDL(p1 << std::endl << p2);
	cv::line(out, p1, p2, cv::Scalar(0, 255, 0),1);
    cv::imwrite("/media/yangxuyuan/Planetarian/ProjectsDuringPostgraduate/HighVoltageLine/images/test/out_ep.jpg", out);
}


//20200224 测试工具,输入已知的匹配点,求出点云
void Test_reconstruct_with_KnownPoses(std::string workingFolder, std::vector<string> vec_imageFile, std::vector<string> vec_featureFilePath, std::vector<Eigen::Matrix<double,3,4>> vec_PMatrix, std::vector<Eigen::Matrix3d> vec_KMatrix, std::vector<std::vector<double>> vec_distort, std::vector<string> vec_matchesFile )
{
    std::vector<View> views(vec_featureFilePath.size());
    //init views
    std::vector<std::vector<Eigen::Vector2d>> vec_points(views.size());
    for(int i=0;i<views.size();i++)
    {
        Eigen::Matrix3d R = vec_PMatrix[i].block(0,0,3,3);
        Eigen::Vector3d t = vec_PMatrix[i].block(0,3,3,1);
        views[i] = View(vec_imageFile[i],i,R, -R.transpose() * t,vec_KMatrix[i], vec_distort[i]);
        std::vector<Eigen::Vector2d> points = getPoint(vec_featureFilePath[i], " ");
        vec_points[i] = points;
        for (uint k = 0; k < points.size(); k++)
        {
            PointTrack pt(points[k], -1, k);
            views[i].pointTracks.emplace_back(pt);
        }
    }
    //relative correspondences
    int matchesIdx = 0;
    std::vector<std::vector<Matches>> vec_vec_matches;
    for(int i=0;i< views.size() - 1 ;i++)
    {
        std::vector<Matches> vec_I_Matches;
        for(int j= i + 1;j < views.size();j++)
        {
            views[i].neighborsID.emplace_back(j);
            Eigen::Matrix3d R_relative = views[j].rotation * (views[i].rotation.transpose());
            Eigen::Vector3d t_relative = views[j].t - R_relative * views[i].t;
            Eigen::Matrix3d F = views[j].K_Inv.transpose() * toAntisymmetricMatrix(t_relative) * R_relative * views[i].K_Inv;
            F /= F(2,2);
            views[i].F_Me2Neighbors[j] = (F);

            Matches matches_ij(i,j);
            if(vec_matchesFile.size() > 0)
            {
                matches_ij = getMatchesFromFile(vec_matchesFile[matchesIdx],i,j);
                matchesIdx++;
            }
            else
            {
                if(views[i].pointTracks.size() != views[j].pointTracks.size())
                {
                    std::cout <<"not matching point size in " << i << " " << j << std::endl;
                    return;
                }
                for(int k=0;k<views[i].pointTracks.size();k++)
                {
                    matches_ij.map_featIdx[k] = k;
                }
            }
            GeometricFilter(matches_ij, F, vec_points[i], vec_points[j]);
            vec_I_Matches.emplace_back(matches_ij);
        }
        vec_vec_matches.emplace_back(vec_I_Matches);
    }

    Reconstruction reconstruction(views);
    //set matches
    for(int i=0;i<views.size()-1;i++)
    {
        reconstruction.map_viewID_Matches[i] = vec_vec_matches[i];
    }
    //reconstruction
    reconstruction.findTracksInViews();
    COUTENDL("begin sTRUCTURE....");
    reconstruction.GetStructure();
    reconstruction.SavePLYFile_PointCloud(workingFolder + "/test_orig_fixRT.ply");

}
