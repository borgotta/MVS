#pragma once

#include "Camera.h"
#include <opencv2\calib3d\calib3d.hpp>

namespace MVS {
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
		inline Vec3f getColor(const float fx, const float fy) const;  
		inline Vec3f getColor(const Vec4f& coord) const;
		//inline int getMask(const Vec4f& coord) const;
		//inline int getEdge(const Vec4f& coord) const;

		//ordered number
		int n;

		Mat image;
	};

	Vec3f Photo::getColor(const Vec4f& coord) const {
		vector<Point3f> dst;
		vector<Point2f> p2d;
		convertPointsHomogeneous(coord, dst);
		projectPoints(dst,rvec,tvec,m_intrinsics,Mat(),p2d);
		return getColor(p2d[0].x, p2d[0].y);
	}
	Vec3f Photo::getColor(float fx, float fy) const {
		return image.at<Vec3f>((int)(floor(fy+0.5f)),(int)(floor(fx+0.5f)));
	}
};
