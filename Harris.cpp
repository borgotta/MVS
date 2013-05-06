#include "stdafx.h"
#include "Harris.h"
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/highgui/highgui.hpp>

namespace MVS {
	void Harris::run(const Photo &photo, const float sigma, multiset<Point> &result) {

		Mat response = Mat::zeros( photo.image.size(), CV_32FC1 );
		Mat src_grey;
		cvtColor(photo.image,src_grey,CV_BGR2GRAY );
		cornerHarris(src_grey, response, 3, 3, sigma, BORDER_DEFAULT);
		/// Normalizing
		Mat response_norm;
		normalize( response, response_norm, 0, 255, NORM_MINMAX, CV_32FC1, Mat() );

		//----------------------------------------------------------------------
		//suppress non local max
		Mat vvftmp = response_norm;
		for (int y = 1; y < response_norm.rows - 1; ++y) {
			for (int x = 1; x < response_norm.cols - 1; ++x) {
				if (response_norm.at<float>(y,x) < response_norm.at<float>(y,x+1) ||
					response_norm.at<float>(y,x) < response_norm.at<float>(y,x-1) ||
					response_norm.at<float>(y,x) < response_norm.at<float>(y+1,x) ||
					response_norm.at<float>(y,x) < response_norm.at<float>(y-1,x))
					vvftmp.at<float>(y,x) = 0.0;
			}
		}

		response_norm = vvftmp;
		//vvftmp.swap(m_response);
		Mat dst_norm_scaled;
		//convertScaleAbs( response_norm, dst_norm_scaled );
		//imshow("Display window",response_norm);
		fillGrid(response_norm, 32, 32, result);
	}
};