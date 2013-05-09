#include "StdAfx.h"
#include "InitialMatcher.h"

namespace MVS {

	InitialMatcher::InitialMatcher(Reconstructor& r) : rec(r) {}
	void InitialMatcher::init(const std::vector<std::vector<Point> >& points) {
		m_ppoints.clear();
		//m_ppoints.resize(m_fm.m_num);

		//for (int index = 0; index < m_fm.m_num; ++index) {
		//	const int gheight = m_fm.m_pos.m_gheights[index];
		//	const int gwidth = m_fm.m_pos.m_gwidths[index];
		//	m_ppoints[index].resize(gwidth * gheight);
		//}

		readPoints(points);
	}

	void InitialMatcher::readPoints(const std::vector<std::vector<Point> >& points) {
		for (int index = 0; index < rec.ps.size(); ++index) {
			for (int i = 0; i < (int)points[index].size(); ++i) {
				/*Ppoint ppoint(new Cpoint(points[index][i]));
				ppoint->m_itmp = index;
				const int ix = ((int)floor(ppoint->m_icoord[0] + 0.5f)) / m_fm.m_csize;
				const int iy = ((int)floor(ppoint->m_icoord[1] + 0.5f)) / m_fm.m_csize;
				const int index2 = iy * m_fm.m_pos.m_gwidths[index] + ix;
				m_ppoints[index][index2].push_back(ppoint);*/
			}
		}
	}
};