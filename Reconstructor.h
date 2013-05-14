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

		//number of images
		int n_im;
		//cell size
		int csize;
		//min number of images with the point visible
		int m_minImageNumThreshold;
		// nccThreshold - Normalized cross correlation
		float m_nccThreshold;

		// sequence Threshold
		int m_sequenceThreshold;

		// Threshold on filterQuad
		float m_quadThreshold;

		// Maximum number of images used in the optimization
		int m_tau;

		// If patches are dense or not, that is, if we use check(patch) after patch optimization
		int m_depth;

		//----------------------------------------------------------------------
		// Thresholds
		//----------------------------------------------------------------------
		// For first feature matching. Images within this angle are used in
		// matching.
		float m_angleThreshold0;
		// tigher angle
		float m_angleThreshold1;

		// Number of success generation from each seed point
		int m_countThreshold0;
		// Number of counts, expansion can be tried
		int m_countThreshold1;

		// Number of trials for each cell in seed
		int m_countThreshold2;

		// Parameter for isNeighbor in findemptyblocks
		float m_neighborThreshold;
		// Parameter for isNeighbor in filterOutside
		float m_neighborThreshold1;
		// Parameter for filterNeighbor
		float m_neighborThreshold2;

		// ncc threshold before optim
		float m_nccThresholdBefore;
		// Maximum angle of images must be at least as large as this
		float m_maxAngleThreshold;

		// visibility consistency threshold
		float m_visibleThreshold;
		float m_visibleThresholdLoose;

		// Epipolar distance in seed generation
		float m_epThreshold;

		//----------------------------------------------------------------------
		// Core members
		//----------------------------------------------------------------------
		// Patch organizer
		PatchOrganizer po;
		DetectFeatures df;
		PhotoSet ps;
		Settings settings;
		InitialMatcher im;

	};
};
