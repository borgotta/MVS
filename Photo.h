#pragma once

#include "Camera.h"

class Photo :
	public Camera
{
public:
	Photo(void);
	virtual ~Photo(void);

	void init(const std::string iname, const std::string cname, const int maxLevel);

	int getWidth();
	int getHight();

	//Furukawa's methods
	inline Vec3f getColor(const float fx, const float fy, const int level) const;  
	inline Vec3f getColor(const Vec4f& coord, const int level) const;
	inline int getMask(const Vec4f& coord, const int level) const;
	inline int getEdge(const Vec4f& coord, const int level) const;

	//ordered number
	int n;

	Mat image;
};

