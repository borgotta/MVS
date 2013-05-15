#include "stdafx.h"
#include "PhotoSet.h"
#include "MatLib.h"


using namespace MVS;

Photo PhotoSet::getPhoto(int i) {
	return photoList[i];
}
Mat PhotoSet::getImage(int i) {
	return photoList[i].image;
}
void PhotoSet::init(const vector<string> &iL, const vector<string> &tL) {
	int size = 0;
	if (iL.size() > tL.size()) {
		size = tL.size();
	} else size = iL.size();
	//main filling
	photoList = vector<Photo>(size);
	for (int i = 0; i < size; i++) {
		Photo temp;
		temp.init(iL[i], tL[i], 0);
		photoList[i] = temp;
	}
}
int PhotoSet::size(void) {
	return photoList.size();
}
int PhotoSet::getWidth(const int index){
	return photoList[index].image.cols;
}
int PhotoSet::getHeight(const int index) {
	return photoList[index].image.rows;
}
void PhotoSet::setDistances(void) {
	m_distances.resize(photoList.size());
	float avedis = 0.0f;
	int denom = 0;
	for (int i = 0; i < photoList.size(); ++i) {
		m_distances[i].resize(photoList.size());
		for (int j = 0; j < photoList.size(); ++j) {
			if (i == j)
				m_distances[i][j] = 0.0f;
			else {
				const float ftmp = norm(photoList[i].m_center - photoList[j].m_center);
				m_distances[i][j] = ftmp;
				avedis += ftmp;
				denom++;
			}
		}
	}
	if (denom == 0)
		return;

	avedis /= denom;
	if (avedis == 0.0f) {
		cerr << "All the optical centers are identical..?" << endl;
		exit (1);
	}

	// plus angle difference
	for (int i = 0; i < photoList.size(); ++i) {
		Vec4f ray0 = photoList[i].m_oaxis;
		ray0[3] = 0.0f;
		for (int j = 0; j < photoList.size(); ++j) {
			Vec4f ray1 = photoList[j].m_oaxis;
			ray1[3] = 0.0f;

			m_distances[i][j] /= avedis;
			const float margin = cos(10.0f * CV_PI / 180.0f);
			const float dis = max(0.0f, 1.0f - mult(ray0 , ray1) - margin);
			m_distances[i][j] += dis;
		}
	}
}
int PhotoSet::checkAngles(const Vec4f& coord, const std::vector<int>& indexes,
	const float minAngle, const float maxAngle,
	const int tau) const{
		  int count = 0;
  
  vector<Vec4f> rays;  rays.resize((int)indexes.size());
  for (int i = 0; i < (int)indexes.size(); ++i) {
	const int index = indexes[i];
	rays[i] = m_photos[index].m_center - coord;
	unitize(rays[i]);
  }
  
  for (int i = 0; i < (int)indexes.size(); ++i) {
	for (int j = i+1; j < (int)indexes.size(); ++j) {
	  const float dot = max(-1.0f, min(1.0f, rays[i] * rays[j]));
	  const float angle = acos(dot);
	  if (minAngle < angle && angle < maxAngle)
		++count;
	}
  }

  //if (count < num * (num - 1) / 2)
  if (count < 1)
	return 1;
  else
	return 0;
}