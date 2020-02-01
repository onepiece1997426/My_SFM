#pragma once
#include <opencv2/opencv.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

struct Matches
{
	uint viewIdx1, viewIdx2; //一般 idx1 < idx2
	//std::vector<uint> points1Idx, points2Idx;   //对应点在PointTracks中的序号
	std::map<uint, uint> map_featIdx;  //map_thisPointTrackIdx_anotherPointTrackIdx

	Matches() {}
	Matches(uint idx1, uint idx2)
	{
		viewIdx1 = idx1;
		viewIdx2 = idx2;
	}

};


std::vector<std::string> split(std::string str, std::string pattern);
std::vector<Eigen::Vector2d> getPoint(std::string file_path, std::string pattern);
Matches getMatchesFromFile(std::string filePath, uint viewID1, uint viewID2, uint offset1 = 0, uint offset2 = 0);
