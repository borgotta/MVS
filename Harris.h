//Harris
#pragma once

#include "stdafx.h"
#include "Photo.h"
#include "Detector.h"

#include <set>

#include <opencv2/core/core.hpp>

using namespace cv;
using namespace std;

namespace MVS {
	class Harris : public Detector {
	public: 
		Harris(void);
		virtual ~Harris(void);

		static void run(const Photo &photo, const float sigma, multiset<Point> &result);

	protected:
		
		//some other assisting methods and functions
		
	};
};