#pragma once

#include "stdafx.h"

#include <iostream>
#include <sstream>

#include <opencv2/core/core.hpp>

using namespace cv;
using namespace std;

class DataContainer
{
public:
	DataContainer(vector<string> iL, vector<string> tL) : imageList(iL), txtList(tL) {
	}

	virtual int size();
	virtual Mat getImage(const int i);
	virtual Mat getProjMat(const int i);

	vector<string> imageList;
	vector<string> txtList;
};

