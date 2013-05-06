//difference-of-gaussians 
#pragma once

#include "stdafx.h"
#include "Photo.h"
#include "Detector.h"

#include <set>

#include <opencv2/core/core.hpp>

using namespace cv;
using namespace std;

namespace MVS {
	class Dog : public Detector {
	public: 
		Dog(void);
		virtual ~Dog(void);

		static void run(const Photo &photo, const float first_scale, const float second_scale, multiset<Point> &result);

	protected:

		//some other assisting methods and functions
		
	};
};