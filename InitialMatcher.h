#pragma once

#include "Reconstructor.h"

namespace MVS {
	class Reconstructor; //что бы это ни значило
	class InitialMatcher
	{
	public:
		InitialMatcher(Reconstructor &r);
		~InitialMatcher(void){};

		void init(const std::vector<std::vector<Point> >& points);
		void run();
	protected:
		void readPoints(const std::vector<std::vector<Point> >& points);
	public:
		Reconstructor& rec;
		// points in a grid. For each index, grid
		std::vector<std::vector<std::vector<Point> > > m_ppoints;

		//----------------------------------------------------------------------
		// thread related
		//----------------------------------------------------------------------
		void initialMatchThread(void);
		static void* initialMatchThreadTmp(void* arg);
	};
};
