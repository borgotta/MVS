#include "StdAfx.h"
#include "Photo.h"

#include <opencv2/highgui/highgui.hpp>

using namespace MVS;

Photo::Photo(void)
{
}


Photo::~Photo(void)
{
}

void Photo::init(const std::string iname, const std::string cname, const int maxLevel){
	try {
		Camera::init(cname, maxLevel);
		image = imread(iname, CV_LOAD_IMAGE_COLOR);
		if(! image.data )                              // Check for invalid input
		{
			cout <<  "Could not open or find the image" << std::endl ;
			return;
		}
	}
	catch (Exception e) {
		cout << "Caught exception number:  " << e.code << endl;
	}
}