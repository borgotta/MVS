#include "stdafx.h"
#include <opencv2\imgproc\imgproc.hpp>
#include "Dog.h"

namespace MVS {
	void Dog::run(const Photo &photo, const float first_scale, const float second_scale, multiset<Point> &result) {
		
				Mat g1, g2, m_result;
				GaussianBlur(photo.image,g1,Size(0,0),first_scale,0,0);
				GaussianBlur(photo.image,g2,Size(0,0),second_scale,0,0);
				m_result = g1 - g2;
	}
};