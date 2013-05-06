#pragma once

#include <set>
#include <opencv2\core\core.hpp>
#include "Point.h"

using namespace cv;
using namespace std;

namespace MVS {
	class Detector
	{
	public:
		Detector(void);
		virtual ~Detector(void);

		static void fillGrid(const Mat response, int h, int w, int type, multiset<Point> &result);
	protected:
		static float Detector::setThreshold(std::multiset<Point>& grid);
	};
};
