#include "StdAfx.h"
#include "Camera.h"

#include <ios>
#include <fstream>

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
	m_projection = Mat::eye(3,4,CV_64F);
	m_intrinsics = Mat::eye(3,3,CV_64F);
	m_extrinsics = Mat::eye(3,4,CV_64F);

	//instinsics 
	for (int i=0; i<3; i++) 
		for (int j=0; j<3; j++) 
			ifstr >> m_intrinsics.at<double>(i,j);
	//extrinsics

	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) 
			ifstr >> m_extrinsics.at<double>(i,j);
		ifstr >> m_extrinsics.at<double>(i,3);
	}

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


	//updateCamera();
}

//void Camera::updateCamera(void) {
//  updateProjection();
//  
//  //----------------------------------------------------------------------
//  m_oaxis = m_projection[0][2];
//  m_oaxis[3] = 0.0;
//  const float ftmp = norm(m_oaxis);
//  m_oaxis[3] = m_projection[0][2][3];
//  m_oaxis /= ftmp;
//
//  m_center = getOpticalCenter();
//
//  m_zaxis = Vec3f(m_oaxis[0], m_oaxis[1], m_oaxis[2]);
//  m_xaxis = Vec3f(m_projection[0][0][0],
//		  m_projection[0][0][1],
//		  m_projection[0][0][2]);
//  m_yaxis = cross(m_zaxis, m_xaxis);
//  unitize(m_yaxis);
//  m_xaxis = cross(m_yaxis, m_zaxis);
//  
//  Vec4f xaxis = m_projection[0][0];  xaxis[3] = 0.0f;    
//  Vec4f yaxis = m_projection[0][1];  yaxis[3] = 0.0f;
//  float ftmp2 = (norm(xaxis) + norm(yaxis)) / 2.0f;
//  if (ftmp2 == 0.0f)
//	ftmp2 = 1.0f;
//  m_ipscale = ftmp2;
//}