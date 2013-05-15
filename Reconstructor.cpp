// MVS.cpp : Defines the entry point for the console application.
#include "stdafx.h"

#include "Reconstructor.h"

namespace MVS {

	Reconstructor::Reconstructor(void) :im(*this), po(*this), optim(*this) {}
	Reconstructor::~Reconstructor(void) {
		pthread_rwlock_destroy(&m_lock);

		for (int image = 0; image < (int)m_imageLocks.size(); ++image)
			pthread_rwlock_destroy(&m_imageLocks[image]);
		for (int image = 0; image < (int)m_countLocks.size(); ++image)
			pthread_rwlock_destroy(&m_countLocks[image]);
	}

	void Reconstructor::init(Settings &s) {
		settings = s;
		settings.loadPhotos(ps);
		n_im = ps.size();
		ps.setDistances();

		csize = 5; //TODO добавить чтение из настроек
		m_nccThreshold = 0.7f; //TODO добавить чтение из настроек
		m_minImageNumThreshold = 3; //TODO добавить чтение из настроек
		m_sequenceThreshold = 3; //TODO добавить чтение из настроек (по умолчанию -1)
		m_wsize = 7;

		//init visdata
		
		m_visdata.resize(n_im);
		m_visdata2.resize(n_im);
		for (int y = 0; y < n_im; ++y) {
			m_visdata[y].resize(n_im);
			for (int x = 0; x < n_im; ++x)
				if (x == y)
					m_visdata[y][x] = 0;
				else {
					m_visdata[y][x] = 1;
					m_visdata2[y].push_back(x);
				}
		}

		m_junit = 100;
		// This initialization does not matter
		m_visibleThreshold = 0.0f;
		m_visibleThresholdLoose = 0.0f;

		//m_tau = max(option.m_minImageNum * 2, min(m_num, 5));
		m_tau = min(m_minImageNumThreshold * 2, n_im);

		m_depth = 0;



		//----------------------------------------------------------------------
		pthread_rwlock_init(&m_lock, NULL);
		m_imageLocks.resize(n_im);
		m_countLocks.resize(n_im);
		for (int image = 0; image < n_im; ++image) {
			pthread_rwlock_init(&m_imageLocks[image], NULL);
			pthread_rwlock_init(&m_countLocks[image], NULL);
		}

		//initialization 
		df = DetectFeatures();
		//TODO
		df.run(ps,16,16,4);
		po.init();
		im.init(df.m_points);

		//----------------------------------------------------------------------
		// Init thresholds
		m_angleThreshold0 = 60.0f * CV_PI / 180.0f;
		m_angleThreshold1 = 60.0f * CV_PI / 180.0f;

		m_countThreshold0 = 2;
		m_countThreshold1 = 4;
		m_countThreshold2 = 2;

		m_neighborThreshold = 0.5f;
		m_neighborThreshold1 = 1.0f;

		m_neighborThreshold2 = 1.0f;

		m_maxAngleThreshold = 10;

		m_nccThresholdBefore = m_nccThreshold - 0.3f;

		m_quadThreshold = 2.5f;

		m_epThreshold = 2.0f;

	}

	void Reconstructor::run() {
		namedWindow( "Display window", CV_WINDOW_AUTOSIZE );// Create a window for display.



		Mat tmp = ps.photoList[15].image.clone();
		for (int i = 0; i<df.m_points[15].size(); i++) {
			if ((int)df.m_points[15][i].m_type == 0) {
				circle( tmp, cv::Point( (int)df.m_points[15][i].m_icoord[0],(int)df.m_points[15][i].m_icoord[1]), 2,  Scalar(0,255,0), -1, 8, 0 );
			}
			else {
				circle( tmp, cv::Point( (int)df.m_points[15][i].m_icoord[0],(int)df.m_points[15][i].m_icoord[1]), 2,  Scalar(0,0,255), -1, 8, 0 );
			}
		}
		imshow( "Display window", tmp );
		cout<<"df.m_points[15].size() = " << df.m_points[15].size();

	}


};