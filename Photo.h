#pragma once

#include "Camera.h"

class Photo :
	public Camera
{
public:
	Photo(void);
	virtual ~Photo(void);

	int getWidth();
	int getHight();

	//Furukawa's methods
	inline Vec3f getColor(const float fx, const float fy, const int level) const;  
	inline Vec3f getColor(const Vec4f& coord, const int level) const;
	inline int getMask(const Vec4f& coord, const int level) const;
	inline int getEdge(const Vec4f& coord, const int level) const;

	Mat image;
};

