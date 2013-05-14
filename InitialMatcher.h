#pragma once
#include "Point.h"
#include "Patch.h"

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
		int canAdd(const int index, const int x, const int y);  

		void initialMatch(const int index, const int id);
		void collectCells(const int index0, const int index1,
			const Point& p0, std::vector<Vec2i>& cells);

		void collectCandidates(const int index, const std::vector<int>& indexes,
			const Point& point, std::vector<Point>& vcp);

		int initialMatchSub(const int index0, const int index1,
			const int id, Patch& patch);

		void unproject(const int index0, const int index1,
			const Point& p0, const Point& p1,
			Vec4f& coord) const;
		void clear(void);
	public:
		// points in a grid. For each index, grid
		std::vector<std::vector<std::vector<Point> > > m_Points;

		//----------------------------------------------------------------------
		// thread related
		//----------------------------------------------------------------------
		void initialMatchThread(void);
		static void* initialMatchThreadTmp(void* arg);

		// Number of trials
		std::vector<int> m_scounts;
		// Number of failures in the prep
		std::vector<int> m_fcounts0;
		// Number of failures in the post processing
		std::vector<int> m_fcounts1;
		// Number passes
		std::vector<int> m_pcounts;

	protected:
		Reconstructor& rec;
	};
};
