// MVS.cpp : Defines the entry point for the console application.
#include "stdafx.h"

#include "Reconstructor.h"

namespace MVS {

	Reconstructor::Reconstructor(void) :im(*this), po(*this) {}
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

		csize = 5; //TODO добавить чтение из настроек
		n_im = ps.size();

		//----------------------------------------------------------------------
		pthread_rwlock_init(&m_lock, NULL);
		m_imageLocks.resize(n_im);
		m_countLocks.resize(n_im);
		for (int image = 0; image < n_im; ++image) {
			pthread_rwlock_init(&m_imageLocks[image], NULL);
			pthread_rwlock_init(&m_countLocks[image], NULL);
		}

		df = DetectFeatures();
		//TODO
		df.run(ps,16,16,4);
		po.init();
		im.init(df.m_points);


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