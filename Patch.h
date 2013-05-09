#pragma once

#include <vector>
#include <iostream>
#include <opencv2/core/core.hpp>

using namespace cv;
namespace MVS {
	//Furukawa's class
	class Patch {
	public:
		Patch(void) {
			m_ncc = -1.0;
			m_timages = 0;
			m_fix = 0;
			// dflag is initialized only once. if failed in one direction, we
			// never try that.
			m_dflag = 0;
		}

		//----------------------------------------------------------------------
		// saved information
		// 3D coordinates of the center of the patch
		Vec4f m_coord;
		// patch outward normal vector
		Vec4f m_normal;

		// associated image ids. first image id is the reference one. images
		// can be non-targetting image.
		std::vector<int> m_images;
		std::vector<Vec2i> m_grids;

		// visible images. m_vimages must be targetting images.
		std::vector<int> m_vimages;
		std::vector<Vec2i> m_vgrids;

		//----------------------------------------------------------------------
		inline float score(const float threshold) const{
			return std::max(0.0f, m_ncc - threshold) * (int)m_images.size();
		}
		inline float score2(const float threshold) const{
			return std::max(0.0f, m_ncc - threshold) * m_timages;
		}

		// average ncc
		float m_ncc;
		// number of targetting images in m_images
		int m_timages;

		// flat for expansion
		// 0: not yet tested
		// 1: done
		int m_flag;

		// for directional flag
		unsigned char m_dflag;

		// fixed patch or not
		char m_fix;

		// id number in m_ppatches
		int m_id;

		// scaling factor corresponding to one pixel difference
		float m_dscale;
		float m_ascale;

		float m_tmp;
	};
};