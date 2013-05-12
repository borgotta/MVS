#include "StdAfx.h"
#include "Reconstructor.h"
#include "InitialMatcher.h"

namespace MVS {

	InitialMatcher::InitialMatcher(Reconstructor& r) : rec(r) {}
	void InitialMatcher::init(const std::vector<std::vector<Point> >& points) {
		m_ppoints.clear();
		m_ppoints.resize(rec.n_im);

		for (int index = 0; index < rec.n_im; ++index) {
			const int gheight = rec.po.m_gheights[index];
			const int gwidth = rec.po.m_gwidths[index];
			m_ppoints[index].resize(gwidth * gheight);
		}

		readPoints(points);
	}

	void InitialMatcher::readPoints(const std::vector<std::vector<Point> >& points) {
		for (int index = 0; index < rec.n_im; ++index) {
			for (int i = 0; i < (int)points[index].size(); ++i) {
				Point point(points[index][i]);
				point.m_itmp = index;
				const int ix = ((int)floor(point.m_icoord[0] + 0.5f)) / rec.csize;
				const int iy = ((int)floor(point.m_icoord[1] + 0.5f)) / rec.csize;
				const int index2 = iy * rec.po.m_gwidths[index] + ix;
				m_ppoints[index][index2].push_back(point);
			}
		}
	}
};