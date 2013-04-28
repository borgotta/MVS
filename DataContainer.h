#pragma once

#include "stdafx.h"
#include "Photo.h"

#include <iostream>
#include <sstream>

#include <opencv2/core/core.hpp>

using namespace cv;
using namespace std;

class DataContainer
{
public:
	DataContainer(void) {	}

	void init(vector<Mat> iL, vector<Mat> tL);

	int size();
	Mat getImage(const int i);
	Mat getPhoto(const int i);

	vector<Photo> photoList;
};

