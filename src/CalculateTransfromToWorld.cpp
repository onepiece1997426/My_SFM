#include "CalculateTransfromToWorld.h"


//�����������Ե��������Է���
Eigen::Matrix4d PoseEstimation3D3D_CrossProduct(std::vector<Eigen::Vector3d> ptsFixed, std::vector<Eigen::Vector3d> ptsAligned, double scale)
{
	if (ptsFixed.size() != ptsAligned.size() && ptsFixed.size() < 3)
		return Eigen::Matrix4d::Zero();
	
	//������  
	Eigen::Vector3d X0 = ptsFixed[1] - ptsFixed[0];
	X0.normalize();
	Eigen::Vector3d Y0 = ptsFixed[2] - ptsFixed[0];
	Y0 = Y0 - Y0.dot(X0) * X0;
	if (Y0.norm() == 0)
		return Eigen::Matrix4d::Zero();
	Y0.normalize();
	Eigen::Vector3d Z0 = X0.cross(Y0);
	
	Eigen::Vector3d X1 = ptsAligned[1] - ptsAligned[0];
	X1.normalize();
	Eigen::Vector3d Y1 = ptsAligned[2] - ptsAligned[0];
	Y1 = Y1 - Y1.dot(X1) * X1;
	if (Y1.norm() == 0)
		return Eigen::Matrix4d::Zero();
	Y1.normalize();
	Eigen::Vector3d Z1 = X1.cross(Y1);

	//����ϵ��������ı任��
	Eigen::Matrix3d R0, R1;
	R0.block(0, 0, 3, 1) = X0;
	R0.block(0, 1, 3, 1) = Y0;
	R0.block(0, 2, 3, 1) = Z0;
	R1.block(0, 0, 3, 1) = X1;
	R1.block(0, 1, 3, 1) = Y1;
	R1.block(0, 2, 3, 1) = Z1;
	Eigen::Matrix3d R = R0 * R1.transpose();

	/*Eigen::Vector3d t = Eigen::Vector3d::Zero();
	for (uint i = 0; i < ptsAligned.size(); i++)
	{
		t += ptsFixed[i] - R * scale * ptsAligned[i];
	}
	t /= ptsAligned.size();*/
	
	Eigen::Vector3d t = ptsFixed[0] - R * scale * ptsAligned[0];

	Eigen::Matrix4d T;

	//block(start,start,size,size)
	T.block(0, 0, 3, 3) = R;
	T.block(0, 3, 3, 1) = t;
	T.row(3) << 0, 0, 0, 1;
	std::cout << "Calculate T = " << std::endl << T << std::endl;

	return T;
}


Eigen::Matrix4d PoseEstimation3D3D_CrossProduct_point(std::vector<Eigen::Vector3d> ptsFixed, std::vector<Eigen::Vector3d> ptsScaledAligned,double s)
{
	if (ptsFixed.size() != ptsScaledAligned.size() && ptsFixed.size() < 3)
		return Eigen::Matrix4d::Zero();

	//�����ѡ����index ��Ϊԭ��, x�� ��һ��ͬƽ�����
	std::vector<uint> vec_initialIndexes;
	srand((uint)time(0));
	for (uint i = 0; i < 3; i++)
	{
		uint idx = rand() % (ptsFixed.size()); //ȡ��[a,b)�����������ʹ��(rand() % (b-a))+ a;
		bool repeatIdx = false;
		for (uint j = 0; j < i; j++)
		{
			if (vec_initialIndexes[j] == idx)
			{
				repeatIdx = true;
				break;
			}
		}

		if (!repeatIdx)
		{
			vec_initialIndexes.push_back(idx);
		}
		else
		{
			i--;
		}
	}
	//------------

	//���������໥��ֱ���� ����FΪFixed,AΪAligned �Ǻ� ������f1f2��Ӧ�ĵ㹹�ɵ�һ��x��,�����Ӧ��y,z��
	Eigen::Vector3d F1X = ptsFixed[vec_initialIndexes[1]] - ptsFixed[vec_initialIndexes[0]];
	Eigen::Vector3d F1F3 = ptsFixed[vec_initialIndexes[2]] - ptsFixed[vec_initialIndexes[0]];
	Eigen::Vector3d F1Y = F1X.cross(F1F3);
	F1Y.normalize();
	Eigen::Vector3d F1Z = F1X.cross(F1Y);
	F1Z.normalize();

	Eigen::Vector3d A1X = ptsScaledAligned[vec_initialIndexes[1]] - ptsScaledAligned[vec_initialIndexes[0]];
	Eigen::Vector3d A1A3 = ptsScaledAligned[vec_initialIndexes[2]] - ptsScaledAligned[vec_initialIndexes[0]];
	Eigen::Vector3d A1Y = A1X.cross(A1A3);
	A1Y.normalize();
	Eigen::Vector3d A1Z = A1X.cross(A1Y);
	A1Z.normalize();

	//���������  �ֱ���y,z����(ֻ��һ���Ǻ�,������worldϵ�µ�y,z��)�ĵ�λһ�� ȷ���������һ����������ƽ�����
	Eigen::Vector3d FY_point = F1Y * F1X.norm() + ptsFixed[vec_initialIndexes[0]];
	Eigen::Vector3d FZ_point = F1Z * F1X.norm() + ptsFixed[vec_initialIndexes[0]];
	ptsFixed.push_back(FY_point);
	ptsFixed.push_back(FZ_point);

	Eigen::Vector3d AY_point = A1Y * A1X.norm() + ptsScaledAligned[vec_initialIndexes[0]];
	Eigen::Vector3d AZ_point = A1Z * A1X.norm() + ptsScaledAligned[vec_initialIndexes[0]];
	ptsScaledAligned.push_back(AY_point);
	ptsScaledAligned.push_back(AZ_point);
	//------------------

	for (uint i = 0; i < ptsScaledAligned.size(); i++)
	{
		ptsScaledAligned[i] *= s;
	}

	//ֱ�Ӽ���SVD����
	//return PoseEstimation3D3D_SVD(ptsFixed, ptsScaledAligned);
	return PoseEstimation3D3D_DLT(ptsFixed, ptsScaledAligned);

}

Eigen::Matrix4d PoseEstimation3D3D_DLT(std::vector<Eigen::Vector3d> ptsFixed, std::vector<Eigen::Vector3d> ptsAligned)
{
	if (ptsFixed.size() != ptsAligned.size() && ptsFixed.size() < 4)
		return Eigen::Matrix4d::Zero();
	Eigen::VectorXd right_vec(ptsFixed.size() * 3);
	Eigen::MatrixXd left_Matrix(ptsAligned.size() * 3, 12);

	//R���ǶԳƾ���
	//��ʼ������
	for (uint i = 0; i < ptsAligned.size(); i++)
	{
		uint row_idx = i * 3;
		left_Matrix.row(row_idx) << ptsAligned[i].x(), ptsAligned[i].y(), ptsAligned[i].z(), 1, 0, 0, 0, 0, 0, 0, 0, 0;
		left_Matrix.row(row_idx + 1) << 0, 0, 0, 0, ptsAligned[i].x(), ptsAligned[i].y(), ptsAligned[i].z(), 1, 0, 0, 0, 0;
		left_Matrix.row(row_idx + 2) << 0, 0, 0, 0, 0, 0, 0, 0, ptsAligned[i].x(), ptsAligned[i].y(), ptsAligned[i].z(), 1;

		right_vec(row_idx) = ptsFixed[i].x();
		right_vec(row_idx + 1) = ptsFixed[i].y();
		right_vec(row_idx + 2) = ptsFixed[i].z();

	}

	//QR �ֽ�
	Eigen::VectorXd t = left_Matrix.colPivHouseholderQr().solve(right_vec);

	Eigen::Matrix4d T;
	T << t(0), t(1), t(2), t(3),
		t(4), t(5), t(6), t(7),
		t(8), t(9), t(10), t(11),
		0, 0, 0, 1;
	std::cout << "T: " << std::endl << T << std::endl;
	return T;

}

Eigen::Matrix4d PoseEstimation3D3D_SVD(std::vector<Eigen::Vector3d> ptsFixed, std::vector<Eigen::Vector3d> ptsAligned)
{
	//����ƥ������� �᲻����ƥ���������ڶ��� ��SVDЧ�����Щ
	if (ptsFixed.size() != ptsAligned.size())
	{
		std::cout << "invalid input" << std::endl;
		return Eigen::Matrix4d::Zero();
	}

	//------

	//��������
	Eigen::Vector3d p1, p2; //����
	uint N = ptsFixed.size();
	for (uint i = 0; i < N; i++)
	{
		p1 += ptsFixed[i];
		p2 += ptsAligned[i];
	}
	p1 /= N;
	p2 /= N;
	//---------

	//��ÿ���������ȥ����
	std::vector<Eigen::Vector3d> qFixed(N), qAligned(N);
	for (int i = 0; i < N; i++)
	{
		qFixed[i] = ptsFixed[i] - p1;
		qAligned[i] = ptsAligned[i] - p2;
	}
	//------

	//����q1q2^T  align��Fixed ����˻���˳����Ҫ��
	Eigen::Matrix3d W = Eigen::Matrix3d::Zero();
	for (uint i = 0; i < N; i++)
	{
		W += qFixed[i] * qAligned[i].transpose();
	}
	//std::cout << " W = " << std::endl << W <<std::endl;


	//SVD on W
	// U��V^T
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(W, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d U = svd.matrixU();
	Eigen::Matrix3d V = svd.matrixV();
	//    Eigen::Matrix3d S = U.inverse() * W * V.transpose().inverse(); //�Խ���
	//    std::cout << " U = " << std::endl << U <<std::endl;
	//    std::cout << " V = " << std::endl << V <<std::endl;
	//    std::cout << " S = " << std::endl << S <<std::endl;
	if (svd.rank() != 3)
	{
		//���ȵ�ʱ����Խ����R
		std::cout << "W rank " << svd.rank() << " is not 3" << std::endl;
		return Eigen::Matrix4d::Zero();
	}
	//---------
	//R t  ����û��scale
	Eigen::Matrix3d R = U * V.transpose();
	Eigen::Vector3d t = p1 - R * p2;

	Eigen::Matrix4d T;

	//block(start,start,size,size)
	T.block(0, 0, 3, 3) = R;
	T.block(0, 3, 3, 1) = t;
	T.row(3) << 0, 0, 0, 1;
	std::cout << "SVD_Calculate T = " << std::endl << T << std::endl;

	return T;
}

//ʹ�ò����ߵ�����
Eigen::Matrix4d PoseEstimation3D3D_CrossProduct_axis(std::vector<Eigen::Vector3d> & ptsFixed, std::vector<Eigen::Vector3d> & ptsScaledAligned)
{
	return Eigen::Matrix4d();
}

void CalScaleAndRescale(std::vector<Eigen::Vector3d> & ptsFixed, std::vector<Eigen::Vector3d> & ptsAligned)
{
	if (ptsFixed.size() != ptsAligned.size())
	{
		std::cout << "not the same scale" << std::endl;
		return;
	}

	double aveScale = 0;
	for (uint i = 1; i < ptsFixed.size(); i++)
	{
		Eigen::Vector3d v = ptsFixed[i] - ptsFixed[i - 1];
		Eigen::Vector3d u = ptsAligned[i] - ptsAligned[i - 1];

		double scale = v.norm() / u.norm();

		aveScale += scale;
	}
	aveScale /= (ptsFixed.size() - 1);

	for (uint i = 0; i < ptsAligned.size(); i++)
	{
		ptsAligned[i] *= aveScale;
	}
	std::cout << "scale: " << aveScale << std::endl;
}

void CalTranslationScale_DirectionBetweenPoints(std::vector<Eigen::Vector3d> & ptsFixed, std::vector<Eigen::Vector3d> & ptsAligned)
{
	if (ptsFixed.size() != ptsAligned.size())
	{
		std::cout << "not the same scale" << std::endl;
		return;
	}

	//��0,1֮��ĳ���Ϊ��λ1����������Ĳ��ֳ��ȵı�ֵ	
	std::map<uint, double> map_scaleLen;
	Eigen::Vector3d v0 = ptsFixed[1] - ptsFixed[0];
	Eigen::Vector3d u0 = ptsAligned[1] - ptsAligned[0];
	for (uint i = 2; i < ptsFixed.size(); i++)
	{
		Eigen::Vector3d v = ptsFixed[i] - ptsFixed[0];
		Eigen::Vector3d u = ptsAligned[i] - ptsAligned[0];

		double scale = (v.norm() * u0.norm()) / (u.norm() * v0.norm());

		map_scaleLen[i] = scale;
		ptsAligned[i] = ptsAligned[0] + scale * (ptsAligned[i] - ptsAligned[0]);
		std::cout << "c_i: " << ptsAligned[i].transpose() << std::endl;
	}
	std::cout << "finish first scale" << std::endl;
}

double CalTotalScale(std::vector<Eigen::Vector3d> ptsFixed, std::vector<Eigen::Vector3d> ptsAligned)
{
	double totalpartScale = 0;
	for (uint i = 1; i < ptsFixed.size(); i++)
	{
		totalpartScale += (ptsFixed[i] - ptsFixed[0]).norm() / (ptsAligned[i] - ptsAligned[0]).norm();
	}
	totalpartScale /= (ptsFixed.size() - 1);
	return totalpartScale;
}

//����E�õ���ʼ��R�ͷ���Ӻ�ͨ�������������任�������ı任��world�ϣ����󲿷�scale�����ڽ���֪���򣩣���������scale��
Eigen::Matrix4d CalTransformToWorldWhenKnowTranslationDir(std::vector<Eigen::Vector3d>& pts_fixed, std::vector<Eigen::Vector3d>& pts_aligned)
{
	/*Eigen::Vector3d c0,c1,c2;
	c0 << 945.115, 1024.913, 52.990;
	c1 << 949.300, 1021.623, 52.989;
	c2 << 951.961, 1019.501, 52.985;

	std::vector<Eigen::Vector3d> pts_fixed{ c0,c1,c2 };

	Eigen::Vector3d c0_, c1_, c2_;
	c0_ << 0, 0, 0;
	c1_ << 0.992218754067413, -0.0284178038593165, -0.121220346892435;
	c2_ << 0.994231543445891, -0.0116150178276532, -0.106624243856945;
	std::vector<Eigen::Vector3d> pts_aligned{ c0_,c1_,c2_ };*/

	//CalScaleAndRescale(pts_fixed, pts_aligned);
	CalTranslationScale_DirectionBetweenPoints(pts_fixed, pts_aligned);
	double s = CalTotalScale(pts_fixed, pts_aligned);

	//Eigen::Matrix4d T = PoseEstimation3D3D_CrossProduct(pts_fixed, pts_aligned,s);
	Eigen::Matrix4d T = PoseEstimation3D3D_CrossProduct_point(pts_fixed, pts_aligned, s);
	std::cout << "T: " << std::endl << T << std::endl;

	Eigen::Matrix3d R = T.block(0, 0, 3, 3);
	Eigen::Vector3d t = T.block(0, 3, 3, 1);
	std::cout << "s: " << s << std::endl;
	std::cout << "R: " << std::endl << R << std::endl << std::endl;
	T.block(0, 0, 3, 3) *= s;
	std::cout << "scaled_T: " << std::endl << T << std::endl;
	std::vector<Eigen::Vector3d> err_vec;
	double err_d = 0;
	for (uint i = 0; i < pts_aligned.size(); i++)
	{
		Eigen::Vector3d v = s * R * pts_aligned[i] + t;
		std::cout << "ci: " << v.transpose() << std::endl;
		err_vec.emplace_back(pts_fixed[i] - v);
		std::cout << "ei: " << err_vec[i].transpose() << std::endl << std::endl;
		err_d += err_vec[i].norm();
	}
	err_d /= err_vec.size();
	std::cout << "ave err: " << err_d << std::endl;

	return T;
}
