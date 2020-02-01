#include "CalMatchPoints.h"


//�ָ��ַ���
std::vector<std::string> split(std::string str, std::string pattern)
{
	std::vector<std::string> ret;
	if (pattern.empty()) return ret;
	size_t start = 0, index = str.find_first_of(pattern, 0);
	while (index != str.npos)
	{
		if (start != index)
			ret.push_back(str.substr(start, index - start));
		start = index + 1;
		index = str.find_first_of(pattern, start);
	}
	if (!str.substr(start).empty())
		ret.push_back(str.substr(start));
	return ret;
}

std::vector<Eigen::Vector2d> getPoint(std::string file_path, std::string pattern)
{
	std::ifstream files;
	
	//������vector
	std::vector<Eigen::Vector2d> points;

	//����㣬��Ϊpush_back�Ķ���
	//�����Ҫע�⣬vector��ÿһ����Ҫ��Ϊһ��������pushback�������pushback���飬���pushback�ࡣ
	Eigen::Vector2d temp;

	//�������������������洢��ȡ������
	int temp1 = 0, temp2 = 0;

	files.open(file_path);

	std::string s;
	while (std::getline(files, s))
	{
		//ȥ��" ," д��vector��
		std::vector< std::string> result = split(s, pattern);
		temp1 = atof(result[0].c_str());
		temp2 = atof(result[1].c_str());

		temp.x() = temp1;
		temp.y() = temp2;	

		points.push_back(temp);
	}

	files.close();
	return points;
}

Matches getMatchesFromFile(std::string filePath, uint viewID1, uint viewID2,uint offset1, uint offset2)
{
	Matches matches(viewID1, viewID2);
	std::vector<Eigen::Vector2d> matchesPoints = getPoint(filePath, " ");
	/*matches.points1Idx.resize(matchesPoints.size() + offset1);
	matches.points2Idx.resize(matchesPoints.size() + offset2);*/
	for (uint i = 0; i < matchesPoints.size(); i++)
	{
		matches.map_featIdx[offset1 + matchesPoints[i].x()] = offset2 + matchesPoints[i].y();
	}

	return matches;
}
