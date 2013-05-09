// MVS.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include "MVS.h"

namespace MVS {
	Reconstructor::Reconstructor(void){
	}
	void Reconstructor::init(const Settings &s) {
		settings = s;
		df = DetectFeatures();
	}

	void Reconstructor::run() {
		namedWindow( "Display window", CV_WINDOW_AUTOSIZE );// Create a window for display.
		//imshow( "Display window", ps.getPhoto(4).image );  

		df.run(ps,16,16,4);

		Mat tmp = ps.photoList[15].image.clone();
		for (int i = 0; i<df.m_points[15].size(); i++) {
			if ((int)df.m_points[15][i].m_type == 0) {
				circle( tmp, cv::Point( (int)df.m_points[15][i].m_icoord[0],(int)df.m_points[15][i].m_icoord[1]), 2,  Scalar(0,255,0), -1, 8, 0 );
			}
			else {
				circle( tmp, cv::Point( (int)df.m_points[15][i].m_icoord[0],(int)df.m_points[15][i].m_icoord[1]), 2,  Scalar(0,0,255), -1, 8, 0 );
			}
		}
		imshow( "Display window", tmp );
		cout<<"df.m_points[15].size() = " << df.m_points[15].size();

		Mat g1, g2, m_result;
		GaussianBlur(ps.photoList[15].image,g1,Size(0,0),1.0f,0,0);
		GaussianBlur(ps.photoList[15].image,g2,Size(0,0),4.0f,0,0);
		m_result = g1 - g2;
	}


};