#pragma once

#include <pthread.h>
#include <list>

#include "PhotoSet.h"
#include "Point.h"


using namespace cv;
using namespace std;

namespace MVS {
	class DetectFeatures {

	public:
		DetectFeatures(void);
		virtual ~DetectFeatures();

		void run(const PhotoSet& pss,
			const int num, const int csize,
			const int CPU = 1);

		std::vector<std::vector<Point> > m_points;

	protected:
		const PhotoSet* m_ppss;
		int m_csize;
		int m_level;

		//----------------------------------------------------------------------
		// thread related
		//----------------------------------------------------------------------  
		pthread_rwlock_t m_rwlock;
		int m_CPU;

		list<int> m_jobs;

		void runThread(void);
		static void* runThreadTmp(void*arg);
	};
};