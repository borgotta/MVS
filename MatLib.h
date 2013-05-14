#pragma once

#include <cmath>
#include <opencv2\core\core.hpp>
#include "Camera.h"

using namespace std;
using namespace cv;

static void gray2rgb(const float gray, float& r, float& g, float& b);
static void setF(const MVS::Camera& lhs, const MVS::Camera& rhs,
	  Mat& F);

inline float mult(Vec4f x, Vec4f y) {
	return x[0]*y[0] + x[1]*y[1] + x[2]*y[2] + x[3]*y[3];
}

inline float mult(Vec3f x, Vec3f y) {
	return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

inline float getLength(Vec4f x) {
	return sqrt(mult(x,x));
}

inline float getLength(Vec3f x) {
	return sqrt(mult(x,x));
}

inline void unitize(Vec4f& x) {
	const float inverseLength = 1.0f / getLength(x); 
	x[0] *= inverseLength;
	x[1] *= inverseLength;
	x[2] *= inverseLength;
	x[3] *= inverseLength;
}

inline void unitize(Vec3f& x) {
	const float inverseLength = 1.0f / getLength(x); 
	x[0] *= inverseLength;
	x[1] *= inverseLength;
	x[2] *= inverseLength;
}

inline Vec3f cross(const Vec3f& u, const Vec3f& v) {
		return Vec3f( u[1]*v[2] - v[1]*u[2],
		-u[0]*v[2] + v[0]*u[2],
		 u[0]*v[1] - v[0]*u[1] );
}

void gray2rgb(const float gray, float& r, float& g, float& b) {
  if (gray < 0.5) {
	r = 0.0f;
	g = 2.0f * gray;
	b = 1.0f - g;
  }
  else {
	r = (gray - 0.5f) * 2.0f;
	g = 1.0f - r;
	b = 0.0f;
  }
}
void setF(const MVS::Camera& lhs, const MVS::Camera& rhs,
	Mat& F) {
		
		const vector<float> v00 = lhs.m_projection.row(0);
		const vector<float> v01 = lhs.m_projection.row(1);
		const vector<float> v02 = lhs.m_projection.row(2);

		const vector<float> v10 = rhs.m_projection.row(0);
		const vector<float> v11 = rhs.m_projection.row(1);
		const vector<float> v12 = rhs.m_projection.row(2);

		Mat m00 = Mat(0,4,CV_32F);
		m00.push_back(v01);
		m00.push_back(v02);
		m00.push_back(v11);
		m00.push_back(v12);
		Mat m01 = Mat(0,4,CV_32F);
		m01.push_back(v01);
		m01.push_back(v02);
		m01.push_back(v12);
		m01.push_back(v10);
		Mat m02 = Mat(0,4,CV_32F);
		m02.push_back(v01);
		m02.push_back(v02);
		m02.push_back(v10);
		m02.push_back(v11);

		Mat m10 = Mat(0,4,CV_32F);
		m10.push_back(v02);
		m10.push_back(v00);
		m10.push_back(v11);
		m10.push_back(v12);
		Mat m11 = Mat(0,4,CV_32F);
		m11.push_back(v02);
		m11.push_back(v00);
		m11.push_back(v12);
		m11.push_back(v10);
		Mat m12 = Mat(0,4,CV_32F);
		m12.push_back(v02);
		m12.push_back(v00);
		m12.push_back(v10);
		m12.push_back(v11);

		Mat m20 = Mat(4,4,CV_32F);
		m20.push_back(v00);
		m20.push_back(v01);
		m20.push_back(v11);
		m20.push_back(v12);
		Mat m21 = Mat(4,4,CV_32F);
		m21.push_back(v00);
		m21.push_back(v01);
		m21.push_back(v12);
		m21.push_back(v10);
		Mat m22 = Mat(4,4,CV_32F);
		m22.push_back(v00);
		m22.push_back(v01);
		m22.push_back(v10);
		m22.push_back(v11);

		F.at<float>(0,0) = (float)determinant(m00);
		F.at<float>(0,1) = (float)determinant(m01);
		F.at<float>(0,2) = (float)determinant(m02);

		F.at<float>(1,0) = (float)determinant(m10);
		F.at<float>(1,1) = (float)determinant(m11);
		F.at<float>(1,2) = (float)determinant(m12);

		F.at<float>(2,0) = (float)determinant(m20);
		F.at<float>(2,1) = (float)determinant(m21);
		F.at<float>(2,2) = (float)determinant(m22);
}