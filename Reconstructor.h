#pragma once

#include "Settings.h"
#include "DetectFeatures.h"
#include "Point.h"
#include "InitialMatcher.h"
#include "PatchOrganizer.h"


#include <iostream>
#include <sstream>

using namespace std;

namespace MVS {
	class Reconstructor {
	public:
		Reconstructor(void);
		~Reconstructor(void);
		void init(Settings &s);

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
		PatchOrganizer po;
		DetectFeatures df;
		PhotoSet ps;
		Settings settings;
		InitialMatcher im;
		//number of images
		int n_im;
		//cell size
		int csize;
	};
};
