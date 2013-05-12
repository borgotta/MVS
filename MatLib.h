#pragma once

#include <cmath>
#include <opencv2\core\core.hpp>

using namespace std;
using namespace cv;

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