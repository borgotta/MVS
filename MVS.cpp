// MVS.cpp : Defines the entry point for the console application.
#include "stdafx.h"

#include "opencv2/imgproc/imgproc.hpp"

#include "Settings.h"
#include "DetectFeatures.h"
#include "Point.h"

#include <iostream>
#include <sstream>

using namespace std;
using namespace MVS;

int _tmain(int argc, char* argv[])
{
	//mainCalib("in_VID5.xml");
	string path;
	Settings settings;
	DetectFeatures df;

	if (argc >= 2) {
		path = argv[1];
	} else {
		path = "xml\\default.xml";
	}
	const string constpath = path;
	FileStorage fs(constpath, FileStorage::READ); // Read the settings
	if (!fs.isOpened())
	{
		cout << "Could not open the configuration file: \"" << path << "\"" << endl;
		return -1;
	}
	fs["Settings"] >> settings;
	fs.release(); 

	//loading images(photos)

	cout<<path<<"\n"<<settings.imageList.size()<<"\n";
	namedWindow( "Display window", CV_WINDOW_AUTOSIZE );// Create a window for display.
	//imshow( "Display window", settings.data.getPhoto(4).image );  

	df.run(settings.data,16,16,4);

	Mat tmp = settings.data.photoList[15].image.clone();
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
	GaussianBlur(settings.data.photoList[15].image,g1,Size(0,0),1.0f,0,0);
	GaussianBlur(settings.data.photoList[15].image,g2,Size(0,0),4.0f,0,0);
	m_result = g1 - g2;
	//imshow("Display window",m_result);

	waitKey(0);
	return 0;
}

