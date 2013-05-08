#include "stdafx.h"

#include "DetectFeatures.h"
#include "Harris.h"
#include "Dog.h"

#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <set>

//temp func
void temp1(const MVS::Photo &photo, const float k, multiset<MVS::Point> &result) {

		vector<Point2f> response; 
		Mat src_grey;
		cvtColor(photo.image,src_grey,CV_BGR2GRAY );
		goodFeaturesToTrack(src_grey,response,1000,0.00001,3,noArray(),3,true,0.06);
		vector<Point2f>::iterator begin = response.begin();
		vector<Point2f>::iterator end = response.end();
		while(begin != end) {
			MVS::Point p;
						p.m_icoord = Vec3f((*begin).x, (*begin).y, 1.0f);
						//p.m_response = *begin;
						p.m_type = 0;
			result.insert(p);
		begin++;
		}
}
namespace MVS {
	DetectFeatures::DetectFeatures(void) {
		pthread_rwlock_init(&m_rwlock, NULL);
	}

	DetectFeatures::~DetectFeatures() {
		pthread_rwlock_destroy(&m_rwlock);
	}

	void DetectFeatures::run(const PhotoSet& pss, const int num,
		const int csize, const int CPU) {
			m_ppss = &pss;
			m_csize = csize;
			m_CPU = CPU;

			m_points.clear();
			m_points.resize(num);

			//----------------------------------------------------------------------
			for (int index = 0; index < num; ++index)
				m_jobs.push_back(index);

			vector<pthread_t> threads(m_CPU);
			for (int i = 0; i < m_CPU; ++i)
				pthread_create(&threads[i], NULL, runThreadTmp, (void*)this);
			for (int i = 0; i < m_CPU; ++i)
				pthread_join(threads[i], NULL);
			//----------------------------------------------------------------------
			cerr << "done" << endl;
	}

	void* DetectFeatures::runThreadTmp(void* arg) {
		DetectFeatures* detectFeatures = (DetectFeatures*)arg;  
		detectFeatures->runThread();
		return NULL;
	}

	void DetectFeatures::runThread(void) {
		while (1) {
			int index = -1;
			pthread_rwlock_wrlock(&m_rwlock);
			if (!m_jobs.empty()) {
				index = m_jobs.front();
				m_jobs.pop_front();
			}
			pthread_rwlock_unlock(&m_rwlock);
			if (index == -1)
				break;

			const int image = index; //m_ppss->m_images[index];
			cerr << image << ' ' << flush;

			//?????????????  May need file lock, because targetting images
			//should not overlap among multiple processors.    
			/*	char buffer[1024];
			sprintf(buffer, "%smodels/%08d.affin%d", m_ppss->m_prefix.c_str(), image, m_level);
			ifstream ifstr;
			ifstr.open(buffer);
			if (ifstr.is_open()) {
			ifstr.close();
			continue;
			}
			ifstr.close();*/

			//----------------------------------------------------------------------
			// parameters
			// for harris
			const float sigma = 4.0f;
			// for DoG
			const float firstScale = 1.0f;    const float lastScale = 3.0f;

			//----------------------------------------------------------------------
			// Harris
			{
				//Charris harris;
				multiset<Point> result;
				Photo temp = m_ppss->photoList[index];

				/*Mat response = Mat::zeros( temp.image.size(), CV_32FC1 );
				Mat src_grey;
				cvtColor(temp.image,src_grey,CV_BGR2GRAY );
				cornerHarris(src_grey, response, 3, 3, sigma, BORDER_DEFAULT);
				/// Normalizing
				Mat response_norm;
				normalize( response, response_norm, 0, 255, NORM_MINMAX, CV_32FC1, Mat() );
				Mat showM;
				convertScaleAbs( response_norm, showM );
				imshow("Display window",showM);*/
				//pthread_rwlock_wrlock(&m_rwlock);
				temp1(temp,0.06,result);
				//Harris::run(temp,0.06f,result);

				/*harris.run(m_ppss->m_photos[index].getImage(m_level),
				m_ppss->m_photos[index].Cimage::getMask(m_level),
				m_ppss->m_photos[index].Cimage::getEdge(m_level),
				m_ppss->m_photos[index].getWidth(m_level),
				m_ppss->m_photos[index].getHeight(m_level), m_csize, sigma, result);*/

				multiset<Point>::reverse_iterator rbegin = result.rbegin();
				while (rbegin != result.rend()) {
					m_points[index].push_back(*rbegin);
					rbegin++;
				}
			}

			//----------------------------------------------------------------------
			// DoG
			{
				
				/*multiset<Point> result;
				Photo temp = m_ppss->photoList[index];
				Dog::run(temp,1.0f, 4.0f,result);
				

				multiset<Point>::reverse_iterator rbegin = result.rbegin();      
				while (rbegin != result.rend()) {
					m_points[index].push_back(*rbegin);
					rbegin++;
				}*/
			}


		}
	}
};

