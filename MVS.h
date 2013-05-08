#include "opencv2/imgproc/imgproc.hpp"

#include "Settings.h"
#include "DetectFeatures.h"
#include "Point.h"
#include "InitialMatcher.h"

#include <iostream>
#include <sstream>

using namespace std;

namespace MVS {
	class MVS {
	public:
		void init(const Settings &s);
		void run();

		//----------------------------------------------------------------------
		// For threads related
		//----------------------------------------------------------------------
		// General lock
		pthread_rwlock_t m_lock;
		// For each image
		std::vector<pthread_rwlock_t> m_imageLocks;
		std::vector<pthread_rwlock_t> m_countLocks;
		// count
		int m_count;
		// jobs
		std::list<int> m_jobs;
		// job unit
		int m_junit;

		//----------------------------------------------------------------------
		// Core members
		//----------------------------------------------------------------------
		// Patch organizer
		//CpatchOrganizerS m_pos;
		DetectFeatures df;
		Settings settings;
	};
};