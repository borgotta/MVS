#include "StdAfx.h"
#include "Camera.h"
#include "MatLib.h"

#include <ios>
#include <fstream>
#include <opencv2\calib3d\calib3d.hpp>

using namespace MVS;

Camera::Camera(void)
{
}


Camera::~Camera(void)
{
}

void Camera::init(const std::string cname, const int maxLevel) {
	m_cname = cname;
	m_maxLevel = maxLevel;

	// initialize camera
	m_intrinsics.resize(6);
	m_extrinsics.resize(6);

	ifstream ifstr;
	ifstr.open(cname.c_str());

	string header;
	ifstr >> header;
	if (header != "CONTOUR")
	{
		cerr << "Unrecognizable txt format" << endl;
		exit (1);
	}
	//m_projection = Mat::eye(3,4,CV_32F);
	m_intrinsics = Mat::eye(3,3,CV_32F);
	m_extrinsics = Mat::eye(3,4,CV_32F);
	Mat rmat = Mat::eye(3,3,CV_32F);

	//instinsics 
	for (int i=0; i<3; i++) 
		for (int j=0; j<3; j++) 
			ifstr >> m_intrinsics.at<float>(i,j);
	//extrinsics
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			ifstr >> m_extrinsics.at<float>(i,j);
			rmat.at<float>(i,j) = m_extrinsics.at<float>(i,j);
		}
	}
	for (int i=0; i<3; i++) {
		ifstr >> m_extrinsics.at<float>(i,3);
		tvec[i] = m_extrinsics.at<float>(i,3);
	}
	//todo create rvec rodrigues
	Rodrigues(rmat,rvec);
	//setting projection matrix
	m_projection = m_intrinsics * m_extrinsics;

	//setting otical center 
	m_center = getOpticalCenter();

	//TODO fill projection matrix
	//for (int i=0; i<3; i++) {
	//	for (int j=0; j<4; j++) {
	//		ifstr >> temp.at<double>(i,j);
	//	}}
	//m_projection = temp.clone();
	//temp.release();

	//for (int i = 0; i < 6; ++i)
	//	ifstr >> m_intrinsics[i];
	//for (int i = 0; i < 6; ++i)
	//	ifstr >> m_extrinsics[i];

	ifstr.close();

	//----------------------------------------------------------------------
	//?? what is it i don't know
	//m_projection.resize(maxLevel*3);


	updateCamera();
}

void Camera::updateCamera(void) {

	//----------------------------------------------------------------------
	m_projection.row(2).copyTo(m_oaxis);
	m_oaxis[3] = 0.0;
	const float ftmp = norm(m_oaxis);
	m_oaxis[3] = m_projection.at<float>(2,3);
	m_oaxis /= ftmp;

	m_center = getOpticalCenter();

	m_zaxis = Vec3f(m_oaxis[0], m_oaxis[1], m_oaxis[2]);
	m_xaxis = Vec3f(m_projection.at<float>(0,0),
		m_projection.at<float>(0,1),
		m_projection.at<float>(0,2));
	m_yaxis = cross(m_zaxis, m_xaxis);
	unitize(m_yaxis);
	m_xaxis = cross(m_yaxis, m_zaxis);

	Vec4f xaxis;
	m_projection.row(0).copyTo(xaxis);  
	xaxis[3] = 0.0f;    
	Vec4f yaxis;
	m_projection.row(1).copyTo(yaxis);
	yaxis[3] = 0.0f;
	float ftmp2 = (norm(xaxis) + norm(yaxis)) / 2.0f;
	if (ftmp2 == 0.0f)
		ftmp2 = 1.0f;
	m_ipscale = ftmp2;
}
Vec4f Camera::getOpticalCenter(void) const {
  // orthographic case
  Vec4f ans;
  if (m_projection.at<float>(2,0) == 0.0 && m_projection.at<float>(2,1) == 0.0 &&
	  m_projection.at<float>(2,2) == 0.0) {
	Vec3f vtmp[2];
	for (int i = 0; i < 2; ++i)
	  for (int y = 0; y < 3; ++y)
	vtmp[i][y] = m_projection.at<float>(i,y);
	
	Vec3f vtmp2 = cross(vtmp[0], vtmp[1]);
	unitize(vtmp2);
	for (int y = 0; y < 3; ++y)
	  ans[y] = vtmp2[y];
	ans[3] = 0.0;
  }
  else {
	Mat A = Mat::eye(3,3,CV_32F);
	Mat b = Mat::eye(3,1,CV_32F);; //TODO maybe it is supposed to be int  & maybe it is needed to be 1x3
	for (int y = 0; y < 3; ++y) {
	  for (int x = 0; x < 3; ++x)
	A.at<float>(y,x) = m_projection.at<float>(y,x);
	  b.at<float>(y,0) = - m_projection.at<float>(y,3);
	}
	Mat iA = Mat::eye(3,3,CV_32F);
	invert(iA, A);
	b = iA * b;
	
	for (int y = 0; y < 3; ++y)
	  ans[y] = b.at<float>(y,0);
	ans[3] = 1.0;
  }
  return ans;
}
inline Vec3f Camera::project(const Vec4f& coord) const {
  Vec3f vtmp;    
  for (int i = 0; i < 3; ++i) {
		Vec4f tmp;
	  m_projection.row(i).copyTo(tmp);
	  vtmp[i] = mult(tmp, coord); //
  }

  if (vtmp[2] <= 0.0) {
	vtmp[0] = -0xffff;
	vtmp[1] = -0xffff;
	vtmp[2] = -1.0f;
	return vtmp;
  }
  else
	vtmp /= vtmp[2];
  
  vtmp[0] = std::max((float)(INT_MIN + 3.0f),
			 std::min((float)(INT_MAX - 3.0f),
				  vtmp[0]));
  vtmp[1] = std::max((float)(INT_MIN + 3.0f),
			 std::min((float)(INT_MAX - 3.0f),
				  vtmp[1]));
  
  return vtmp;
};
