#include "StdAfx.h"
#include "Detector.h"

#include <iostream>

using namespace std;


namespace MVS {
	Detector::Detector(void)
	{	}
	Detector::~Detector(void)
	{	}
	float Detector::setThreshold(std::multiset<Point>& grid) {
		float ave = 0.0;
		float ave2 = 0.0;
		multiset<Point>::iterator begin = grid.begin();
		multiset<Point>::iterator end = grid.end();
		int count = 0;
		while (begin != end) {
			count++;
			ave += begin->m_response;
			ave2 += begin->m_response * begin->m_response;
			begin++;
		}
		if (count == 0)
			count = 1;
		ave /= count;      ave2 /= count;
		ave2 = sqrt(max(0.0f, ave2 - ave * ave));
		const float threshold = ave + ave2;

		//cout << ave << ' ' << ave2 << endl;

		return threshold;
	}

	void Detector::fillGrid(const Mat response, const int h, const int w, int type, multiset<Point> &result){
		const int maxPointsGrid = 4;
		int cols = (response.cols + h - 1) / w;
		int rows = (response.rows + w - 1) / h;

		//result grids
		vector<vector<multiset<Point> > > resultgrids;
		resultgrids.resize(rows);
		for (int y = 0; y < rows; ++y)
			resultgrids[y].resize(cols);

		//const int margin = 2;
		for (int y = 0; y < response.rows; ++y) {
			for (int x = 0; x < response.cols; ++x) {
				if (response.at<float>(y,x) == 0.0f)
					continue;

				const int x0 = min(x / w, cols - 1);
				const int y0 = min(y / h, rows - 1);

				if ((int)resultgrids[y0][x0].size() < maxPointsGrid ||
					resultgrids[y0][x0].begin()->m_response < response.at<float>(y,x)) {
						Point p;
						p.m_icoord = Vec3f(x, y, 1.0f);
						p.m_response = response.at<float>(y,x);
						p.m_type = type;

						resultgrids[y0][x0].insert(p);
						if (maxPointsGrid < (int)resultgrids[y0][x0].size())
							resultgrids[y0][x0].erase(resultgrids[y0][x0].begin ());
				}
			}
		}  

		for (int y = 0; y < rows; ++y)
			for (int x = 0; x < cols; ++x) {
				const float threshold = setThreshold(resultgrids[y][x]) + 0.0001f;//*1.00001;      
				multiset<Point>::iterator begin = resultgrids[y][x].begin();
				multiset<Point>::iterator end = resultgrids[y][x].end();
				while (begin != end) {
					if (threshold <= begin->m_response)
					result.insert(*begin);
					begin++;
				}
			}

	}
};
