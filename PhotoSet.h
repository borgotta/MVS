#pragma once

#include "stdafx.h"
#include "Photo.h"

#include <iostream>
#include <sstream>

#include <opencv2/core/core.hpp>

using namespace cv;
using namespace std;
namespace MVS {
	class DataContainer
	{
	public:
		DataContainer(void) {	}

		void init(const vector<string> &iL, const vector<string> &tL);

		int size();
		Mat getImage(const int i);
		Photo getPhoto(const int i);

		vector<Photo> photoList;
	};
};