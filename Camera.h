#pragma once

#include "stdafx.h"

#include <iostream>
#include <sstream>

#include <opencv2/core/core.hpp>

using namespace cv;
using namespace std;

namespace MVS {
	class Camera
	{
	public:
		Camera(void);
		virtual ~Camera(void);

		// Update projection matrices from intrinsics and extrinsics
		void updateProjection(void);
		// Update all the camera related parameters. 
		void updateCamera(void);

		virtual void init(const std::string cname, const int maxLevel);
		void write(const std::string file);

		inline Vec3f project(const Vec4f& coord) const;
		//inline Vec3f mult(const Vec4f& coord, const int level) const;

		//----------------------------------------------------------------------
		// txt file name
		std::string m_cname;  
		// Optical center
		Vec4f m_center;
		// Optical axis
		Vec4f m_oaxis;

		float m_ipscale;
		// 3x4 projection matrix
		//std::vector<std::vector<Vec4f> > m_projection;
		Mat m_projection;

		// intrinsic and extrinsic camera parameters. Compact form.
		Mat m_intrinsics;;
		Mat m_extrinsics;

		Vec3f tvec;
		Vec3f rvec;

		Vec3f m_xaxis;
		Vec3f m_yaxis;
		Vec3f m_zaxis;



	protected:
		int m_maxLevel;

		float m_axesScale;

		Vec4f getOpticalCenter(void) const;
	};
};