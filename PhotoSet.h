#pragma once

#include "stdafx.h"
#include "Photo.h"

#include <iostream>
#include <sstream>

#include <opencv2/core/core.hpp>

using namespace cv;
using namespace std;
namespace MVS {
	class PhotoSet
	{
	public:
		PhotoSet(void) {	}

		void init(const vector<string> &iL, const vector<string> &tL);
		int getWidth(const int index);
		int getHeight(const int index);
		int size();
		Mat getImage(const int i);
		Photo getPhoto(const int i);
		int checkAngles(const Vec4f& coord, const std::vector<int>& indexes,
			const float minAngle, const float maxAngle,
			const int tau) const;
		inline Vec3f project(const int index, const Vec4f& coord) const;

		// pairwise distance based on optical center and viewing direction
		void setDistances(void);
		std::vector<std::vector<float> > m_distances;
		vector<Photo> photoList;
	};

	Vec3f PhotoSet::project(const int index, const Vec4f& coord) const{
		return photoList[index].project(coord);
};

};