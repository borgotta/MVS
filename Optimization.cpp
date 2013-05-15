#include "StdAfx.h"
#include <algorithm>
#include <numeric>
#include <cstdio>
#include "Optimization.h"
#include "Reconstructor.h"
#include <gsl_deriv.h>
#include "lmmin.h"

namespace MVS {
	Optimization::Optimization(Reconstructor& r) : rec(r) {}
	Optimization::~Optimization(void)	{}

	void Optimization::init(void) {
		m_vect0T.resize(rec.settings.CPU);
		m_centersT.resize(rec.settings.CPU);
		m_raysT.resize(rec.settings.CPU);
		m_indexesT.resize(rec.settings.CPU);
		m_dscalesT.resize(rec.settings.CPU);
		m_ascalesT.resize(rec.settings.CPU);
		m_paramsT.resize(rec.settings.CPU);

		m_texsT.resize(rec.settings.CPU);
		m_weightsT.resize(rec.settings.CPU);

		for (int c = 0; c < rec.settings.CPU; ++c) {
			m_texsT[c].resize(rec.n_im);
			m_weightsT[c].resize(rec.n_im);
			for (int j = 0; j < rec.m_tau; ++j)
				m_texsT[c][j].resize(3 * rec.m_wsize * rec.m_wsize);
		}

		setAxesScales();
	}
	void Optimization::setAxesScales(void) { 
		m_xaxes.resize(rec.n_im);
		m_yaxes.resize(rec.n_im);
		m_zaxes.resize(rec.n_im);
		for (int index = 0; index < rec.n_im; ++index) {
			m_zaxes[index] = Vec3f(rec.ps.photoList[index].m_oaxis[0],
				rec.ps.photoList[index].m_oaxis[1],
				rec.ps.photoList[index].m_oaxis[2]);
			m_xaxes[index] = Vec3f(rec.ps.photoList[index].m_projection.at<float>(0,0),
				rec.ps.photoList[index].m_projection.at<float>(0,1),
				rec.ps.photoList[index].m_projection.at<float>(0,2));
			m_yaxes[index] = cross(m_zaxes[index], m_xaxes[index]);
			unitize(m_yaxes[index]);
			m_xaxes[index] = cross(m_yaxes[index], m_zaxes[index]);
		}

		m_ipscales.resize(rec.n_im);
		for (int index = 0; index < rec.n_im; ++index) {
			const Vec4f xaxe(m_xaxes[index][0], m_xaxes[index][1], m_xaxes[index][2], 0.0);
			const Vec4f yaxe(m_yaxes[index][0], m_yaxes[index][1], m_yaxes[index][2], 0.0);

			Vec4f fx_temp;
			Vec4f fy_temp;
			rec.ps.photoList[index].m_projection.row(0).copyTo(fx_temp);
			rec.ps.photoList[index].m_projection.row(1).copyTo(fy_temp);
			const float fx = mult(xaxe, fx_temp) ;
			const float fy = mult(yaxe, fy_temp) ;
			m_ipscales[index] = fx + fy;
		}  
	}

	int Optimization::preProcess(Patch& patch, const int id, const int seed) {
		addImages(patch);

		// Here define reference images, and sort images.
		// Something similar to constraintImages is done inside.
		constraintImages(patch, rec.m_nccThresholdBefore, id);

		// Fix the reference image and sort the other  m_tau - 1 images.
		sortImages(patch);

		// Pierre Moulon (it avoid crash in some case)
		if( (int)patch.m_images.size() > 0)
		{
			// setSscales should be here to avoid noisy output
			rec.po.setScales(patch);
		}

		// Check minimum number of images
		if ((int)patch.m_images.size() < rec.m_minImageNumThreshold)
			return 1;

		const int flag =
			rec.ps.checkAngles(patch.m_coord, patch.m_images,
			rec.m_maxAngleThreshold,
			rec.m_angleThreshold1,
			rec.m_minImageNumThreshold);

		if (flag) {
			patch.m_images.clear();
			return 1;
		}

		return 0;
	}

	void Optimization::collectImages(const int index, std::vector<int>& indexes) const {
		// Find images with constraints m_angleThreshold, m_visdata,
		// m_sequenceThreshold, m_targets. Results are sorted by
		// CphotoSet::m_distances.
		indexes.clear();
		Vec4f ray0 = rec.ps.photoList[index].m_oaxis;
		ray0[3] = 0.0f;

		vector<Vec2f> candidates;
		// Search for only related images
		for (int i = 0; i < (int)rec.m_visdata2[index].size(); ++i) {
			const int indextmp = rec.m_visdata2[index][i];

			//if (rec.m_tnum <= indextmp)
			//continue;
			if (rec.m_sequenceThreshold != -1 &&
				rec.m_sequenceThreshold < abs(index - indextmp))
				continue;

			Vec4f ray1 = rec.ps.photoList[indextmp].m_oaxis;
			ray1[3] = 0.0f;
			float tmp = mult(ray0,ray1);
			if (tmp < cos(rec.m_angleThreshold0))
				continue;

			candidates.push_back(Vec2f(rec.ps.m_distances[index][indextmp], indextmp));
		}

		sort(candidates.begin(), candidates.end());//TODO не знаю работает или нет //Svec2cmp<float>());
		for (int i = 0; i < min(rec.m_tau, (int)candidates.size()); ++i)
			indexes.push_back((int)candidates[i][1]);
	}

	void Optimization::sortImages(Patch& patch) const{
  const int newm = 1;
  if (newm == 1) {
	const float threshold = 1.0f - cos(10.0 * CV_PI / 180.0);
	vector<int> indexes, indexes2;
	vector<float> units, units2;
	vector<Vec4f> rays, rays2;
	
	computeUnits(patch, indexes, units, rays);
	
	patch.m_images.clear();
	if (indexes.size() < 2)
	  return;
	
	units[0] = 0.0f;
	
	while (!indexes.empty()) {
	  vector<float>::iterator ite = min_element(units.begin(), units.end());
	  const int index = ite - units.begin();
	  
	  patch.m_images.push_back(indexes[index]);
	  
	  // Remove other images within 5 degrees
	  indexes2.clear();    units2.clear();
	  rays2.clear();
	  for (int j = 0; j < (int)rays.size(); ++j) {
		if (j == index)
		  continue;
		
		indexes2.push_back(indexes[j]);
		rays2.push_back(rays[j]);
		const float ftmp = min(threshold, max(threshold / 2.0f, 1.0f - (mult(rays[index], rays[j]))));
		
		units2.push_back(units[j] * (threshold / ftmp));
	  }
	  indexes2.swap(indexes);
	  units2.swap(units);
	  rays2.swap(rays);
	}
  }
  else {
	//----------------------------------------------------------------------
	//Sort and grab the best m_tau images. All the other images don't
	//matter.  First image is the reference and fixed
	const float threshold = cos(5.0 * CV_PI / 180.0);
	vector<int> indexes, indexes2;
	vector<float> units, units2;
	vector<Vec4f> rays, rays2;
	
	computeUnits(patch, indexes, units, rays);
	
	patch.m_images.clear();
	if (indexes.size() < 2)
	  return;
	
	units[0] = 0.0f;
	
	while (!indexes.empty()) {
	  //for (int i = 0; i < size; ++i) {
	  vector<float>::iterator ite = min_element(units.begin(), units.end());
	  const int index = ite - units.begin();
	  
	  patch.m_images.push_back(indexes[index]);    
	  
	  // Remove other images within 5 degrees
	  indexes2.clear();    units2.clear();
	  rays2.clear();
	  for (int j = 0; j < (int)rays.size(); ++j) {
		if ((mult(rays[index],rays[j])) < threshold) {
		  indexes2.push_back(indexes[j]);
		  units2.push_back(units[j]);
		  rays2.push_back(rays[j]);
		}
	  }
	  indexes2.swap(indexes);
	  units2.swap(units);
	  rays2.swap(rays);
	}
  }
}
	
void Optimization::computeUnits(const Patch& patch,
						  std::vector<float>& units) const{
  const int size = (int)patch.m_images.size();
  units.resize(size);

  vector<int>::const_iterator bimage = patch.m_images.begin();
  vector<int>::const_iterator eimage = patch.m_images.end();

  vector<float>::iterator bfine = units.begin();

  while (bimage != eimage) {
	*bfine = INT_MAX/2;
	
	*bfine = getUnit(*bimage, patch.m_coord);
	Vec4f ray = rec.ps.photoList[*bimage].m_center - patch.m_coord;
	unitize(ray);
	const float denom = mult(ray,patch.m_normal);
	if (0.0 < denom)
	  *bfine /= denom;
	else
	  *bfine = INT_MAX/2;
	
	++bimage;
	++bfine;
  }
}

void Optimization::computeUnits(const Patch& patch,
						  std::vector<int>& indexes,
						  std::vector<float>& units,
						  std::vector<Vec4f>& rays) const{
  vector<int>::const_iterator bimage = patch.m_images.begin();
  vector<int>::const_iterator eimage = patch.m_images.end();
  
  while (bimage != eimage) {
	Vec4f ray = rec.ps.photoList[*bimage].m_center - patch.m_coord;
	unitize(ray);
	const float dot = mult(ray, patch.m_normal);
	if (dot <= 0.0f) {
	  ++bimage;
	  continue;
	}

	const float scale = getUnit(*bimage, patch.m_coord);
	const float fine = scale / dot;

	indexes.push_back(*bimage);
	units.push_back(fine);
	rays.push_back(ray);    
	++bimage;
  }
}

void Optimization::refinePatch(Patch& patch, const int id,
						 const int time) {
  //refinePatchBFGS(patch, id, 1000, 1);

  // Try the faster version, if that fails try the slower one
  if(!refinePatchBFGS2(patch, id, 1000, 1)) {
	refinePatchBFGS(patch, id, 1000, 1);
  }

  // WORKED REALLY WELL
  if (patch.m_images.empty())
	return;
}

//----------------------------------------------------------------------
// BFGS functions
//----------------------------------------------------------------------
void Optimization::my_f_lm(const double *par, int m_dat, const void *data, double *fvec, int *info) {
  double xs[3] = {par[0], par[1], par[2]};
  const int id = *((int*)data);

  //const float sigma = CV_PI / 16.0f;//4.0f * CV_PI;
  //const float sigma2 = 2 * sigma * sigma;
  const float angle1 = xs[1] * m_one->m_ascalesT[id];
  const float angle2 = xs[2] * m_one->m_ascalesT[id];

  double ret = 0.0;

  if (angle1 <= - CV_PI / 2.0f || CV_PI / 2.0f <= angle1 ||
	  angle2 <= - CV_PI / 2.0f || CV_PI / 2.0f <= angle2) {
	ret = 2.0f;

	fvec[0] = ret;
	fvec[1] = ret;
	fvec[2] = ret;

	return;
 }
  
  //?????
  const double bias = 0.0f;//2.0 - exp(- angle1 * angle1 / sigma2) - exp(- angle2 * angle2 / sigma2);
  
  Vec4f coord, normal;
  m_one->decode(coord, normal, xs, id);
  
  const int index = m_one->m_indexesT[id][0];
  Vec4f pxaxis, pyaxis;
  m_one->getPAxes(index, coord, normal, pxaxis, pyaxis);
  
  const int size = min(m_one->rec.m_tau, (int)m_one->m_indexesT[id].size());
  const int mininum = min(m_one->rec.m_minImageNumThreshold, size);

  for (int i = 0; i < size; ++i) {
	int flag;
	flag = m_one->grabTex(coord, pxaxis, pyaxis, normal, m_one->m_indexesT[id][i],
						  m_one->rec.m_wsize, m_one->m_texsT[id][i]);

	if (flag == 0)
	  m_one->normalize(m_one->m_texsT[id][i]);
  }

  const int pairwise = 0;
  if (pairwise) {
	double ans = 0.0f;
	int denom = 0;
	for (int i = 0; i < size; ++i) {
	  for (int j = i+1; j < size; ++j) {
		if (m_one->m_texsT[id][i].empty() || m_one->m_texsT[id][j].empty())
		  continue;
		
		ans += robustincc(1.0 - m_one->dot(m_one->m_texsT[id][i], m_one->m_texsT[id][j]));
		denom++;
	  }
	}
	if (denom <
		//m_one->rec.m_minImageNumThreshold *
		//(m_one->rec.m_minImageNumThreshold - 1) / 2)
		mininum * (mininum - 1) / 2)
	  ret = 2.0f;
	else
	  ret = ans / denom + bias;
  }
  else {
	if (m_one->m_texsT[id][0].empty())
	  ret = 2.0f;
	  
	double ans = 0.0f;
	int denom = 0;
	for (int i = 1; i < size; ++i) {
	  if (m_one->m_texsT[id][i].empty())
		continue;
	  ans +=
		robustincc(1.0 - m_one->dot(m_one->m_texsT[id][0], m_one->m_texsT[id][i]));
	  denom++;
	}
	//if (denom < m_one->rec.m_minImageNumThreshold - 1)
	if (denom < mininum - 1)
	  ret = 2.0f;
	else
	  ret = ans / denom + bias;
  }

  fvec[0] = ret;
  fvec[1] = ret;
  fvec[2] = ret;
}
double Optimization::my_f(const gsl_vector *v, void *params) {
  double xs[3] = {gsl_vector_get(v, 0),
				  gsl_vector_get(v, 1),
				  gsl_vector_get(v, 2)};
  const int id = *((int*)params);

  //const float sigma = CV_PI / 16.0f;//4.0f * CV_PI;
  //const float sigma2 = 2 * sigma * sigma;
  const float angle1 = xs[1] * m_one->m_ascalesT[id];
  const float angle2 = xs[2] * m_one->m_ascalesT[id];

  if (angle1 <= - CV_PI / 2.0f || CV_PI / 2.0f <= angle1 ||
	  angle2 <= - CV_PI / 2.0f || CV_PI / 2.0f <= angle2)
	return 2.0f;      
  
  //?????
  const double bias = 0.0f;//2.0 - exp(- angle1 * angle1 / sigma2) - exp(- angle2 * angle2 / sigma2);
  
  Vec4f coord, normal;
  m_one->decode(coord, normal, xs, id);
  
  const int index = m_one->m_indexesT[id][0];
  Vec4f pxaxis, pyaxis;
  m_one->getPAxes(index, coord, normal, pxaxis, pyaxis);
  
  const int size = min(m_one->rec.m_tau, (int)m_one->m_indexesT[id].size());
  const int mininum = min(m_one->rec.m_minImageNumThreshold, size);

  for (int i = 0; i < size; ++i) {
	int flag;
	flag = m_one->grabTex(coord, pxaxis, pyaxis, normal, m_one->m_indexesT[id][i],
						  m_one->rec.m_wsize, m_one->m_texsT[id][i]);

	if (flag == 0)
	  m_one->normalize(m_one->m_texsT[id][i]);
  }

  const int pairwise = 0;
  if (pairwise) {
	double ans = 0.0f;
	int denom = 0;
	for (int i = 0; i < size; ++i) {
	  for (int j = i+1; j < size; ++j) {
		if (m_one->m_texsT[id][i].empty() || m_one->m_texsT[id][j].empty())
		  continue;
		
		ans += robustincc(1.0 - m_one->dot(m_one->m_texsT[id][i], m_one->m_texsT[id][j]));
		denom++;
	  }
	}
	if (denom <
		//m_one->rec.m_minImageNumThreshold *
		//(m_one->rec.m_minImageNumThreshold - 1) / 2)
		mininum * (mininum - 1) / 2)
	  return 2.0f;
	else
	  return ans / denom + bias;
  }
  else {
	if (m_one->m_texsT[id][0].empty())
	  return 2.0f;
	  
	double ans = 0.0f;
	int denom = 0;
	for (int i = 1; i < size; ++i) {
	  if (m_one->m_texsT[id][i].empty())
		continue;
	  ans +=
		robustincc(1.0 - m_one->dot(m_one->m_texsT[id][0], m_one->m_texsT[id][i]));
	  denom++;
	}
	//if (denom < m_one->rec.m_minImageNumThreshold - 1)
	if (denom < mininum - 1)
	  return 2.0f;
	else
	  return ans / denom + bias;
  }
}

/* The gradient of f, df = (df/dx, df/dy). */
void Optimization::my_df(const gsl_vector *v, void *params,
				   gsl_vector *df) {
  int which = 1;
  if (which == 0) {
	  const double x0 = gsl_vector_get(v, 0);
  const double x1 = gsl_vector_get(v, 1);
  const double x2 = gsl_vector_get(v, 2);

  double dfdx0, dfdx1, dfdx2;

  //????
  const double step = 0.4f;//0.4f;
  gsl_vector* x = gsl_vector_alloc (3);
  //-------------------------------------------------------
  gsl_vector_set(x, 0, x0 + step);
  gsl_vector_set(x, 1, x1);
  gsl_vector_set(x, 2, x2);
  dfdx0 = my_f(x, params);
  gsl_vector_set(x, 0, x0 - step);
  dfdx0 -= my_f(x, params);
  dfdx0 /= 2.0 * step;
  //-------------------------------------------------------
  gsl_vector_set(x, 0, x0);
  gsl_vector_set(x, 1, x1 + step);
  dfdx1 = my_f(x, params);
  gsl_vector_set(x, 1, x1 - step);
  dfdx1 -= my_f(x, params);
  dfdx1 /= 2.0 * step;
  //-------------------------------------------------------
  gsl_vector_set(x, 1, x1);
  gsl_vector_set(x, 2, x2 + step);
  dfdx2 = my_f(x, params);
  gsl_vector_set(x, 2, x2 - step);
  dfdx2 -= my_f(x, params);
  dfdx2 /= 2.0 * step;
  //-------------------------------------------------------
  gsl_vector_set(df, 0, dfdx0);
  gsl_vector_set(df, 1, dfdx1);
  gsl_vector_set(df, 2, dfdx2);

  gsl_vector_free(x);
  }
  else {
	//???check h is 0.5 or 1.0?
	const float h = 0.5f;//0.5f;//1.0f;
	const int id = *((int*)params);

  const double x0 = gsl_vector_get(v, 0);
  const double x1 = gsl_vector_get(v, 1);
  const double x2 = gsl_vector_get(v, 2);

  m_one->m_paramsT[id][0] = x0;
  m_one->m_paramsT[id][1] = x1;
  m_one->m_paramsT[id][2] = x2;

  //???check debug abserr here how big are they
  gsl_function func;
  func.params = params;
  
  double result, abserr;
  func.function = &my_f0;
  gsl_deriv_central(&func, x0, h, &result, &abserr);
  gsl_vector_set(df, 0, result);
  
  func.function = &my_f1;
  gsl_deriv_central(&func, x1, h, &result, &abserr);
  gsl_vector_set(df, 1, result);
  
  func.function = &my_f2;  
  gsl_deriv_central(&func, x2, h, &result, &abserr);
  
  gsl_vector_set(df, 2, result);
  }
}
	 
/* Compute both f and df together. */
void Optimization::my_fdf(const gsl_vector *x, void *params, 
					double *f, gsl_vector *df) {  
  *f = my_f(x, params); 
  my_df(x, params, df);
}

double Optimization::my_f0(double x, void* params) {
  const int id = *((int*)params);

  gsl_vector* v = gsl_vector_alloc (3);
  gsl_vector_set(v, 0, x);
  gsl_vector_set(v, 1, m_one->m_paramsT[id][1]);
  gsl_vector_set(v, 2, m_one->m_paramsT[id][2]);

  const double score = my_f(v, params);
  gsl_vector_free(v);
  
  return score;
}

double Optimization::my_f1(double x, void* params) {
  const int id = *((int*)params);

  gsl_vector* v = gsl_vector_alloc (3);
  gsl_vector_set(v, 0, m_one->m_paramsT[id][0]);
  gsl_vector_set(v, 1, x);
  gsl_vector_set(v, 2, m_one->m_paramsT[id][2]);

  const double score = my_f(v, params);
  gsl_vector_free(v);
  
  return score;
}

double Optimization::my_f2(double x, void* params) {
  const int id = *((int*)params);

  gsl_vector* v = gsl_vector_alloc (3);
  gsl_vector_set(v, 0, m_one->m_paramsT[id][0]);
  gsl_vector_set(v, 1, m_one->m_paramsT[id][1]);
  gsl_vector_set(v, 2, x);
  
  const double score = my_f(v, params);
  gsl_vector_free(v);
  
  return score;
}

//----------------------------------------------------------------------
// BFGS functions SSD
//----------------------------------------------------------------------
double Optimization::my_f_ssd(const gsl_vector *v, void *params) {
  double xs[3] = {gsl_vector_get(v, 0),
				  gsl_vector_get(v, 1),
				  gsl_vector_get(v, 2)};
  const int id = *((int*)params);

  //const float sigma = 4.0f * CV_PI;
  //const float sigma2 = 2 * sigma * sigma;
  const float angle1 = xs[1] * m_one->m_ascalesT[id];
  const float angle2 = xs[2] * m_one->m_ascalesT[id];

  if (angle1 <= - CV_PI / 2.0f || CV_PI / 2.0f <= angle1 ||
	  angle2 <= - CV_PI / 2.0f || CV_PI / 2.0f <= angle2)
	return 2.0f;      
  
  const double bias = 0.0f;//2.0 - exp(- angle1 * angle1 / sigma2) - exp(- angle2 * angle2 / sigma2);
  
  Vec4f coord, normal;
  m_one->decode(coord, normal, xs, id);
  
  const int index = m_one->m_indexesT[id][0];
  Vec4f pxaxis, pyaxis;
  m_one->getPAxes(index, coord, normal, pxaxis, pyaxis);
  
  const int size = min(m_one->rec.m_tau, (int)m_one->m_indexesT[id].size());
  
  for (int i = 0; i < size; ++i) {
	int flag;
	flag = m_one->grabTex(coord, pxaxis, pyaxis, normal, m_one->m_indexesT[id][i],
						  m_one->rec.m_wsize, m_one->m_texsT[id][i]);
  }
  m_one->normalize(m_one->m_texsT[id], size);
  
  const int pairwise = 0;
  if (pairwise) {
	double ans = 0.0f;
	int denom = 0;
	for (int i = 0; i < size; ++i) {
	  for (int j = i+1; j < size; ++j) {
		if (m_one->m_texsT[id][i].empty() || m_one->m_texsT[id][j].empty())
		  continue;
		
		ans += m_one->ssd(m_one->m_texsT[id][i], m_one->m_texsT[id][j]);
		denom++;
	  }
	}
	if (denom <
		m_one->rec.m_minImageNumThreshold *
		(m_one->rec.m_minImageNumThreshold - 1) / 2)
	  return 2.0f;
	else
	  return ans / denom + bias;
  }
  else {
	if (m_one->m_texsT[id][0].empty())
	  return 2.0f;
	  
	double ans = 0.0f;
	int denom = 0;
	for (int i = 1; i < size; ++i) {
	  if (m_one->m_texsT[id][i].empty())
		continue;
	  ans += m_one->ssd(m_one->m_texsT[id][0], m_one->m_texsT[id][i]);
	  denom++;
	}
	if (denom < m_one->rec.m_minImageNumThreshold - 1)
	  return 2.0f;
	else
	  return ans / denom + bias;
  }
}

/* The gradient of f, df = (df/dx, df/dy). */
void Optimization::my_df_ssd(const gsl_vector *v, void *params,
					   gsl_vector *df) {
  const float h = 0.5f;
  const int id = *((int*)params);

  const double x0 = gsl_vector_get(v, 0);
  const double x1 = gsl_vector_get(v, 1);
  const double x2 = gsl_vector_get(v, 2);

  m_one->m_paramsT[id][0] = x0;
  m_one->m_paramsT[id][1] = x1;
  m_one->m_paramsT[id][2] = x2;

  gsl_function func;
  func.params = params;
  
  //???check debug abserr here how big are they
  double result, abserr;
  func.function = &my_f_ssd0;
  gsl_deriv_central(&func, x0, h, &result, &abserr);
  gsl_vector_set(df, 0, result);

  func.function = &my_f_ssd1;
  gsl_deriv_central(&func, x1, h, &result, &abserr);
  gsl_vector_set(df, 1, result);

  func.function = &my_f_ssd2;  
  gsl_deriv_central(&func, x2, h, &result, &abserr);
  gsl_vector_set(df, 2, result);
}
	 
/* Compute both f and df together. */
void Optimization::my_fdf_ssd(const gsl_vector *x, void *params, 
						double *f, gsl_vector *df) {  
  *f = my_f_ssd(x, params); 
  my_df_ssd(x, params, df);
}

double Optimization::my_f_ssd0(double x, void* params) {
  const int id = *((int*)params);

  gsl_vector* v = gsl_vector_alloc (3);
  gsl_vector_set(v, 0, x);
  gsl_vector_set(v, 1, m_one->m_paramsT[id][1]);
  gsl_vector_set(v, 2, m_one->m_paramsT[id][2]);

  const double score = my_f_ssd(v, params);
  gsl_vector_free(v);
  
  return score;
}

double Optimization::my_f_ssd1(double x, void* params) {
  const int id = *((int*)params);

  gsl_vector* v = gsl_vector_alloc (3);
  gsl_vector_set(v, 0, m_one->m_paramsT[id][0]);
  gsl_vector_set(v, 1, x);
  gsl_vector_set(v, 2, m_one->m_paramsT[id][2]);

  const double score = my_f_ssd(v, params);
  gsl_vector_free(v);
  
  return score;
}

double Optimization::my_f_ssd2(double x, void* params) {
  const int id = *((int*)params);

  gsl_vector* v = gsl_vector_alloc (3);
  gsl_vector_set(v, 0, m_one->m_paramsT[id][0]);
  gsl_vector_set(v, 1, m_one->m_paramsT[id][1]);
  gsl_vector_set(v, 2, x);
  
  const double score = my_f_ssd(v, params);
  gsl_vector_free(v);
  
  return score;
}
//----------------------------------------------------------------------
// Depth
double Optimization::my_f_depth(const gsl_vector *v, void *params) {
  const int id = *((int*)params);
  double xs[3] = {m_one->m_paramsT[id][0],
				  m_one->m_paramsT[id][1],
				  m_one->m_paramsT[id][2]};
  xs[0] = gsl_vector_get(v, 0);
  
  Vec4f coord, normal;
  m_one->decode(coord, normal, xs, id);
  
  const int index = m_one->m_indexesT[id][0];
  Vec4f pxaxis, pyaxis;
  m_one->getPAxes(index, coord, normal, pxaxis, pyaxis);
  
  const int size = min(m_one->rec.m_tau, (int)m_one->m_indexesT[id].size());

  for (int i = 0; i < size; ++i) {
	int flag;
	flag = m_one->grabTex(coord, pxaxis, pyaxis, normal, m_one->m_indexesT[id][i],
						  m_one->rec.m_wsize, m_one->m_texsT[id][i]);

	if (flag == 0)
	  m_one->normalize(m_one->m_texsT[id][i]);
  }

  const int pairwise = 0;
  if (pairwise) {
	double ans = 0.0f;
	int denom = 0;
	for (int i = 0; i < size; ++i) {
	  for (int j = i+1; j < size; ++j) {
		if (m_one->m_texsT[id][i].empty() || m_one->m_texsT[id][j].empty())
		  continue;
		
		ans += robustincc(1.0 - m_one->dot(m_one->m_texsT[id][i], m_one->m_texsT[id][j]));
		denom++;
	  }
	}
	if (denom <
		m_one->rec.m_minImageNumThreshold *
		(m_one->rec.m_minImageNumThreshold - 1) / 2)
	  return 2.0f;
	else
	  return ans / denom;
  }
  else {
	if (m_one->m_texsT[id][0].empty())
	  return 2.0f;
	  
	double ans = 0.0f;
	int denom = 0;
	for (int i = 1; i < size; ++i) {
	  if (m_one->m_texsT[id][i].empty())
		continue;
	  ans +=
		robustincc(1.0 - m_one->dot(m_one->m_texsT[id][0], m_one->m_texsT[id][i]));
	  denom++;
	}
	if (denom < m_one->rec.m_minImageNumThreshold - 1)
	  return 2.0f;
	else
	  return ans / denom;
  }
}

/* The gradient of f, df = (df/dx, df/dy). */
void Optimization::my_df_depth(const gsl_vector *v, void *params,
						 gsl_vector *df) {
  int which = 0;
  if (which == 0) {
	const double x0 = gsl_vector_get(v, 0);

	double dfdx0;
  
	const double step = 0.05f;
	gsl_vector* x = gsl_vector_alloc (1);
	//-------------------------------------------------------
	gsl_vector_set(x, 0, x0 + step);
	dfdx0 = my_f_depth(x, params);
	gsl_vector_set(x, 0, x0 - step);
	dfdx0 -= my_f_depth(x, params);
	dfdx0 /= 2.0 * step;

	gsl_vector_set(df, 0, dfdx0);
	
	gsl_vector_free(x);
  }
  else {
	//???check h is 0.5 or 1.0?
	const float h = 0.5f;//0.5f;//1.0f;
	
	const double x0 = gsl_vector_get(v, 0);

	//???check debug abserr here how big are they
	gsl_function func;
	func.params = params;
	
	double result, abserr;
	func.function = &my_f0_depth;
	gsl_deriv_central(&func, x0, h, &result, &abserr);
	gsl_vector_set(df, 0, result);
  }
}
	 
/* Compute both f and df together. */
void Optimization::my_fdf_depth(const gsl_vector *x, void *params, 
						  double *f, gsl_vector *df) {  
  *f = my_f_depth(x, params); 
  my_df_depth(x, params, df);
}

double Optimization::my_f0_depth(double x, void* params) {
  gsl_vector* v = gsl_vector_alloc (1);
  gsl_vector_set(v, 0, x);

  const double score = my_f_depth(v, params);
  gsl_vector_free(v);
  
  return score;
}

void Optimization::refinePatchBFGS(Patch& patch, const int id,
							 const int time, const int ncc) {
  int status;
  const gsl_multimin_fminimizer_type *T;
  gsl_multimin_fminimizer *s;
  
  gsl_multimin_function my_func;
  int idtmp = id;
  my_func.n = 3;
  
  if (ncc)
	my_func.f = &my_f;
  else
	my_func.f = &my_f_ssd;
  my_func.params = (void *)&idtmp;
  
  m_centersT[id] = patch.m_coord;
  m_raysT[id] = patch.m_coord -
	rec.ps.photoList[patch.m_images[0]].m_center;
  unitize(m_raysT[id]);
  m_indexesT[id] = patch.m_images;
  
  m_dscalesT[id] = patch.m_dscale;
  m_ascalesT[id] = CV_PI / 48.0f;//patch.m_ascale;
  
  computeUnits(patch, m_weightsT[id]);
  for (int i = 1; i < (int)m_weightsT[id].size(); ++i)
	m_weightsT[id][i] = min(1.0f, m_weightsT[id][0] / m_weightsT[id][i]);  
  m_weightsT[id][0] = 1.0f;
  
  double p[3];
  encode(patch.m_coord, patch.m_normal, p, id);
  
  gsl_vector* x = gsl_vector_alloc (3);
  gsl_vector_set(x, 0, p[0]);
  gsl_vector_set(x, 1, p[1]);
  gsl_vector_set(x, 2, p[2]);
  
  gsl_vector* ss = gsl_vector_alloc (3);
  gsl_vector_set_all(ss, 1.5);

  T = gsl_multimin_fminimizer_nmsimplex;
  s = gsl_multimin_fminimizer_alloc (T, 3);
  gsl_multimin_fminimizer_set (s, &my_func, x, ss);
  
  int iter = 0;
  do {
	++iter;
	status = gsl_multimin_fminimizer_iterate (s);
	
	if (status)
	  break;
	
	double size = gsl_multimin_fminimizer_size(s);
	//status = gsl_multimin_test_size (size, 1e-2);
	status = gsl_multimin_test_size (size, 1e-3);
  } while (status == GSL_CONTINUE && iter < time);
  p[0] = gsl_vector_get(s->x, 0);
  p[1] = gsl_vector_get(s->x, 1);
  p[2] = gsl_vector_get(s->x, 2);
  
  if (status == GSL_SUCCESS) {
	decode(patch.m_coord, patch.m_normal, p, id);
	
	patch.m_ncc = 1.0 -
	  unrobustincc(computeINCC(patch.m_coord,
							   patch.m_normal, patch.m_images, id, 1));
  }
  else
	patch.m_images.clear();   
  
  ++m_status[status + 2];
  
  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (x);
  gsl_vector_free (ss);
}

bool Optimization::refinePatchBFGS2(Patch& patch, const int id,
							 const int time, const int ncc) {
  int idtmp = id;
  
  m_centersT[id] = patch.m_coord;
  m_raysT[id] = patch.m_coord - rec.ps.photoList[patch.m_images[0]].m_center;
  unitize(m_raysT[id]);
  m_indexesT[id] = patch.m_images;
  
  m_dscalesT[id] = patch.m_dscale;
  m_ascalesT[id] = CV_PI / 48.0f;//patch.m_ascale;
  
  computeUnits(patch, m_weightsT[id]);
  for (int i = 1; i < (int)m_weightsT[id].size(); ++i)
	m_weightsT[id][i] = min(1.0f, m_weightsT[id][0] / m_weightsT[id][i]);  
  m_weightsT[id][0] = 1.0f;
  
  double p[3];
  encode(patch.m_coord, patch.m_normal, p, id);
  
  double x[3] = {p[0], p[1], p[2]};
  
  lm_control_struct control = lm_control_float;
  control.epsilon = 1e-5; // default step size is too small for floats
  control.printflags = 0;

  lm_status_struct status;

  // this function requires data >= params, so the later 3 is a fudge
  lmmin(3, x, 3, (void *)&idtmp, my_f_lm, &control, &status, lm_printout_std);

  p[0] = x[0];
  p[1] = x[1];
  p[2] = x[2];

  // status.info 0 to 3 are "good", the rest are bad
  if (status.info >= 0 && status.info <= 3) {
	decode(patch.m_coord, patch.m_normal, p, id);
	
	patch.m_ncc = 1.0 -
	  unrobustincc(computeINCC(patch.m_coord,
							   patch.m_normal, patch.m_images, id, 1));
  }
  else {
	return false;
  }

  return true;
}

void Optimization::encode(const Vec4f& coord,
			double* const vect, const int id) const {
  vect[0] = (mult((coord - m_centersT[id]), m_raysT[id])) / m_dscalesT[id];
}

void Optimization::encode(const Vec4f& coord, const Vec4f& normal,
			double* const vect, const int id) const {
  encode(coord, vect, id);
  
  const int image = m_indexesT[id][0];
  const float fx = m_xaxes[image] * proj(normal); // projects from 4D to 3D, divide by last value
  const float fy = m_yaxes[image] * proj(normal);
  const float fz = m_zaxes[image] * proj(normal);

  vect[2] = asin(max(-1.0f, min(1.0f, fy)));
  const float cosb = cos(vect[2]);

  if (cosb == 0.0)
	vect[1] = 0.0;
  else {
	const float sina = fx / cosb;
	const float cosa = - fz / cosb;
	vect[1] = acos(max(-1.0f, min(1.0f, cosa)));
	if (sina < 0.0)
	  vect[1] = - vect[1];
  }

  vect[1] = vect[1] / m_ascalesT[id];
  vect[2] = vect[2] / m_ascalesT[id];  
}

void Optimization::decode(Vec4f& coord, Vec4f& normal,
			const double* const vect, const int id) const {
  decode(coord, vect, id);
  const int image = m_indexesT[id][0];

  const float angle1 = vect[1] * m_ascalesT[id];
  const float angle2 = vect[2] * m_ascalesT[id];

  const float fx = sin(angle1) * cos(angle2);
  const float fy = sin(angle2);
  const float fz = - cos(angle1) * cos(angle2);

  Vec3f ftmp = m_xaxes[image] * fx + m_yaxes[image] * fy + m_zaxes[image] * fz;
  normal = Vec4f(ftmp[0], ftmp[1], ftmp[2], 0.0f);
}

void Optimization::decode(Vec4f& coord, const double* const vect, const int id) const {
  coord = m_centersT[id] + m_dscalesT[id] * vect[0] * m_raysT[id];
}

void Optimization::setINCCs(const Patch& patch,
					  std::vector<float> & inccs,
					  const std::vector<int>& indexes,
					  const int id, const int robust) {
  const int index = indexes[0];
  Vec4f pxaxis, pyaxis;
  getPAxes(index, patch.m_coord, patch.m_normal, pxaxis, pyaxis);
  
  vector<vector<float> >& texs = m_texsT[id];
  
  const int size = (int)indexes.size();
  for (int i = 0; i < size; ++i) {
	const int flag = grabTex(patch.m_coord, pxaxis, pyaxis, patch.m_normal,
							 indexes[i], rec.m_wsize, texs[i]);
	if (flag == 0) {
	  normalize(texs[i]);
	}
  }
  
  inccs.resize(size);
  if (texs[0].empty()) {
	fill(inccs.begin(), inccs.end(), 2.0f);
	return;
  }
  
  for (int i = 0; i < size; ++i) {
	if (i == 0)
	  inccs[i] = 0.0f;
	else if (!texs[i].empty()) {
	  if (robust == 0)
		inccs[i] = 1.0f - dot(texs[0], texs[i]);
	  else
		inccs[i] = robustincc(1.0f - dot(texs[0], texs[i]));
	}
	else
	  inccs[i] = 2.0f;
  }
}

void Optimization::setINCCs(const Patch& patch,
					  std::vector<std::vector<float> >& inccs,
					  const std::vector<int>& indexes,
					  const int id, const int robust) {
  const int index = indexes[0];
  Vec4f pxaxis, pyaxis;
  getPAxes(index, patch.m_coord, patch.m_normal, pxaxis, pyaxis);
  
  vector<vector<float> >& texs = m_texsT[id];

  const int size = (int)indexes.size();
  for (int i = 0; i < size; ++i) {    
	const int flag = grabTex(patch.m_coord, pxaxis, pyaxis, patch.m_normal,
							 indexes[i], rec.m_wsize, texs[i]);
	
	if (flag == 0)
	  normalize(texs[i]);
  }
  
  inccs.resize(size);
  for (int i = 0; i < size; ++i)
	inccs[i].resize(size);
	
  for (int i = 0; i < size; ++i) {
	inccs[i][i] = 0.0f;
	for (int j = i+1; j < size; ++j) {
	  if (!texs[i].empty() && !texs[j].empty()) {
		if (robust == 0)
		  inccs[j][i] = inccs[i][j] = 1.0f - dot(texs[i], texs[j]);
		else
		  inccs[j][i] = inccs[i][j] = robustincc(1.0f - dot(texs[i], texs[j]));
	  }
	  else
		inccs[j][i] = inccs[i][j] = 2.0f;
	}
  }
}

int Optimization::grabSafe(const int index, const int size, const Vec3f& center,
					 const Vec3f& dx, const Vec3f& dy, const int level) const {
  const int margin = size / 2;

  const Vec3f tl = center - dx * margin - dy * margin;
  const Vec3f tr = center + dx * margin - dy * margin;

  const Vec3f bl = center - dx * margin + dy * margin;
  const Vec3f br = center + dx * margin + dy * margin;

  const float minx = min(tl[0], min(tr[0], min(bl[0], br[0])));
  const float maxx = max(tl[0], max(tr[0], max(bl[0], br[0])));
  const float miny = min(tl[1], min(tr[1], min(bl[1], br[1])));
  const float maxy = max(tl[1], max(tr[1], max(bl[1], br[1])));

  // 1 should be enough
  const int margin2 = 3;
  // ??? may need to change if we change interpolation method
  if (minx < margin2 ||
	  rec.ps.getWidth(index) - 1 - margin2 <= maxx ||
	  miny < margin2 ||
	  rec.ps.getHeight(index) - 1 - margin2 <= maxy)
	return 0;
  return 1;
}

// My own optimisaton
float MyPow2(int x)
{
	const float answers[] = {0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

	return answers[x + 4];
}

static float Log2 = log(2.0f);

int Optimization::grabTex(const Vec4f& coord, const Vec4f& pxaxis, const Vec4f& pyaxis,
			const Vec4f& pzaxis, const int index, const int size,
					std::vector<float>& tex) const {
  tex.clear();

  Vec4f ray = rec.ps.photoList[index].m_center - coord;
  unitize(ray);
  const float weight = max(0.0f, mult(ray, pzaxis));

  //???????
  //if (weight < cos(rec.m_angleThreshold0))
  if (weight < cos(rec.m_angleThreshold1))
	return 1;

  const int margin = size / 2;

  Vec3f center = rec.ps.project(index, coord);  
  Vec3f dx = rec.ps.project(index, coord + pxaxis) - center;
  Vec3f dy = rec.ps.project(index, coord + pyaxis) - center;
  
  const float ratio = (norm(dx) + norm(dy)) / 2.0f;
  //int leveldif = (int)floor(log(ratio) / log(2.0f) + 0.5f);
  int leveldif = (int)floor(log(ratio) / Log2 + 0.5f);

  // Upper limit is 2
  //leveldif = max(-rec.m_level, min(2, leveldif));
  leveldif = max(-0, min(2, leveldif));

  //const float scale = pow(2.0f, (float)leveldif);

  const float scale = MyPow2(leveldif);
  //const int newlevel = rec.m_level + leveldif;
  const int newlevel = 0 + leveldif;

  center /= scale;  dx /= scale;  dy /= scale;
  
  if (grabSafe(index, size, center, dx, dy, newlevel) == 0)
	return 1;
  
  Vec3f left = center - dx * margin - dy * margin;

  tex.resize(3 * size * size);
  float* texp = &tex[0] - 1;
  for (int y = 0; y < size; ++y) {
	Vec3f vftmp = left;
	left += dy;
	for (int x = 0; x < size; ++x) {
		Vec3f color = rec.ps.photoList[index].getColor( vftmp[0], vftmp[1]);
	  *(++texp) = color[0];
	  *(++texp) = color[1];
	  *(++texp) = color[2];
	  vftmp += dx;
	}
  }
  
  return 0;
}

double Optimization::computeINCC(const Vec4f& coord, const Vec4f& normal,
			   const std::vector<int>& indexes, const int id,
						   const int robust) {
  if ((int)indexes.size() < 2)
	return 2.0;

  const int index = indexes[0];
  Vec4f pxaxis, pyaxis;
  getPAxes(index, coord, normal, pxaxis, pyaxis);
  
  return computeINCC(coord, normal, indexes, pxaxis, pyaxis, id, robust);
}

double Optimization::computeINCC(const Vec4f& coord, const Vec4f& normal,
			   const std::vector<int>& indexes, const Vec4f& pxaxis,
			   const Vec4f& pyaxis, const int id,
						   const int robust) {
  if ((int)indexes.size() < 2)
	return 2.0;

  const int size = min(rec.m_tau, (int)indexes.size());
  vector<vector<float> >& texs = m_texsT[id];

  for (int i = 0; i < size; ++i) {
	int flag;
	flag = grabTex(coord, pxaxis, pyaxis, normal,
				   indexes[i], rec.m_wsize, texs[i]);
	
	if (flag == 0)
	  normalize(texs[i]);
  }
  
  if (texs[0].empty())
	return 2.0;

  double score = 0.0;

  // pure pairwise of reference based
#ifdef PMVS_PAIRNCC
  float totalweight = 0.0;
  for (int i = 0; i < size; ++i) {
	for (int j = i+1; j < size; ++j) {
	  if (!texs[i].empty() && !texs[j].empty()) {
		const float ftmp = m_weightsT[id][i] * m_weightsT[id][j];
		totalweight += ftmp;
		if (robust)
		  score += robustincc(1.0 - dot(texs[i], texs[j])) * ftmp;
		else
		  score += (1.0 - dot(texs[i], texs[j])) * ftmp;
	  }	
	}
  }
  
  if (totalweight == 0.0)
	score = 2.0;
  else
	score /= totalweight;
#else
  float totalweight = 0.0;
  for (int i = 1; i < size; ++i) {
	if (!texs[i].empty()) {
	  totalweight += m_weightsT[id][i];
	  if (robust)
		score += robustincc(1.0 - dot(texs[0], texs[i])) * m_weightsT[id][i];
	  else
		score += (1.0 - dot(texs[0], texs[i])) * m_weightsT[id][i];
	}
  }
  if (totalweight == 0.0)
	score = 2.0;
  else
	score /= totalweight;
#endif  

  return score;
}

void Optimization::lfunc(double* p, double* hx, int m, int n, void* adata) {
  int iflag;
  m_one->func(n, m, p, hx, &iflag, adata);
}

void Optimization::func(int m, int n, double* x, double* fvec, int* iflag, void* arg) {
  const int id = *((int*)arg);
  double xs[3] = {x[0], x[1], x[2]};
  
  for (int i = 0; i < m; ++i)
	fvec[i] = 2.0;

  Vec4f coord, normal;
  decode(coord, normal, xs, id);

  const int index = m_indexesT[id][0];
  Vec4f pxaxis, pyaxis;
  getPAxes(index, coord, normal, pxaxis, pyaxis);
  
  const int size = min(rec.m_tau, (int)m_indexesT[id].size());
  
  for (int i = 0; i < size; ++i) {
	int flag;
	flag = grabTex(coord, pxaxis, pyaxis, normal, m_indexesT[id][i],
		   rec.m_wsize, m_texsT[id][i]);

	if (flag == 0)
	  normalize(m_texsT[id][i]);
  }
  
  int count = -1;
  for (int i = 0; i < size; ++i) {
	for (int j = i+1; j < size; ++j) {
	  count++;
	  if (m_texsT[id][i].empty() || m_texsT[id][j].empty())
	continue;
	  
	  fvec[count] = robustincc(1.0 - dot(m_texsT[id][i], m_texsT[id][j]));
	}
  }
}

// Normalize only scale for each image
void Optimization::normalize(std::vector<std::vector<float> >& texs,
					   const int size) {
  // compute average rgb
  Vec3f ave;
  int denom = 0;
  
  vector<Vec3f> rgbs;
  rgbs.resize(size);
  for (int i = 0; i < size; ++i) {
	if (texs[i].empty())
	  continue;
	
	int count = 0;
	while (count < (int)texs[i].size()) {
	  rgbs[i][0] += texs[i][count++];
	  rgbs[i][1] += texs[i][count++];
	  rgbs[i][2] += texs[i][count++];
	}
	rgbs[i] /= (int)texs[i].size() / 3;

	ave += rgbs[i];
	++denom;
  }

  // overall average
  if (denom == 0)
	return;

  ave /= denom;

  // Scale all the colors
  for (int i = 0; i < size; ++i) {
	if (texs[i].empty())
	  continue;
	int count = 0;
	// compute scale
	Vec3f scale;
	for (int j = 0; j < 3; ++j)
	  if (rgbs[i][j] != 0.0f)
		scale[j] = ave[j] / rgbs[i][j];
	
	while (count < (int)texs[i].size()) {
	  texs[i][count++] *= scale[0];
	  texs[i][count++] *= scale[1];
	  texs[i][count++] *= scale[2];
	}
  }
}

void Optimization::normalize(std::vector<float>& tex) {
  const int size = (int)tex.size();
  const int size3 = size / 3;
  Vec3f ave;

  float* texp = &tex[0] - 1;
  for (int i = 0; i < size3; ++i) {
	ave[0] += *(++texp);
	ave[1] += *(++texp);
	ave[2] += *(++texp);
  }
  
  ave /= size3;

  float ave2 = 0.0;
  texp = &tex[0] - 1;
  for (int i = 0; i < size3; ++i) {
	const float f0 = ave[0] - *(++texp);
	const float f1 = ave[1] - *(++texp);
	const float f2 = ave[2] - *(++texp);
	
	ave2 += f0 * f0 + f1 * f1 + f2 * f2;
  }
  
  ave2 = sqrt(ave2 / size);
  
  if (ave2 == 0.0f)
	ave2 = 1.0f;
  
  texp = &tex[0] - 1;
  for (int i = 0; i < size3; ++i) {
	*(++texp) -= ave[0];    *texp /= ave2;
	*(++texp) -= ave[1];    *texp /= ave2;
	*(++texp) -= ave[2];    *texp /= ave2;
  }
}

float Optimization::dot(const std::vector<float>& tex0,
		  const std::vector<float>& tex1) const{
#ifndef PMVS_WNCC
  // Pierre Moulon (use classic access to array, windows STL do not like begin()-1)
  const int size = (int)tex0.size();
  float ans = 0.0f;
  for (int i = 0; i < size; ++i) {
	ans += tex0[i] * tex1[i];
  }  
  return ans / size;
#else
  const int size = (int)tex0.size();
  vector<float>::const_iterator i0 = tex0.begin();
  vector<float>::const_iterator i1 = tex1.begin();
  float ans = 0.0f;
  for (int i = 0; i < size; ++i, ++i0, ++i1) {
	ans += (*i0) * (*i1) * m_template[i];
  }
  return ans;
#endif
}

float Optimization::ssd(const std::vector<float>& tex0,
		  const std::vector<float>& tex1) const{
  const float scale = 0.01;

#ifndef PMVS_WNCC
  // Pierre Moulon (use classic access to array, windows STL do not like begin()-1)
  const int size = (int)tex0.size();
  float ans = 0.0f;
  for(int i=0; i < size; ++i) {
	ans += fabs( tex0[i] - tex1[i] );
  }
  
  return scale * ans / size;
#else
  const int size = (int)tex0.size();
  vector<float>::const_iterator i0 = tex0.begin();
  vector<float>::const_iterator i1 = tex1.begin();
  float ans = 0.0f;
  for (int i = 0; i < size; ++i, ++i0, ++i1) {
	const float ftmp = fabs((*i0) - (*i1));
	//ans += (*i0) * (*i1) * m_template[i];
	ans += ftmp * m_template[i];
  }
  return scale * ans;
#endif
}

float Optimization::getUnit(const int index, const Vec4f& coord) const {
  const float fz = norm(coord - rec.ps.photoList[index].m_center);
  const float ftmp = m_ipscales[index];
  if (ftmp == 0.0)
	return 1.0;

  //return 2.0 * fz * (0x0001 << rec.m_level) / ftmp;
  return 2.0 * fz * (0x0001 << 0) / ftmp;
}

// get x and y axis to collect textures given reference image and normal
void Optimization::getPAxes(const int index, const Vec4f& coord, const Vec4f& normal,
			  Vec4f& pxaxis, Vec4f& pyaxis) const{  
  // yasu changed here for fpmvs
  const float pscale = getUnit(index, coord);

  Vec3f normal3(normal[0], normal[1], normal[2]);
  Vec3f yaxis3 = cross(normal3, m_xaxes[index]);
  unitize(yaxis3);
  Vec3f xaxis3 = cross(yaxis3, normal3);
  pxaxis[0] = xaxis3[0];  pxaxis[1] = xaxis3[1];  pxaxis[2] = xaxis3[2];  pxaxis[3] = 0.0;
  pyaxis[0] = yaxis3[0];  pyaxis[1] = yaxis3[1];  pyaxis[2] = yaxis3[2];  pyaxis[3] = 0.0;

  pxaxis *= pscale;
  pyaxis *= pscale;
  const float xdis = norm(rec.ps.project(index, coord + pxaxis) -
						  rec.ps.project(index, coord));
  const float ydis = norm(rec.ps.project(index, coord + pyaxis) -
						  rec.ps.project(index, coord));
  pxaxis /= xdis;
  pyaxis /= ydis;
}

void Optimization::setWeightsT(const Patch& patch, const int id) {
  computeUnits(patch, m_weightsT[id]);
  for (int i = 1; i < (int)m_weightsT[id].size(); ++i)
	m_weightsT[id][i] = min(1.0f, m_weightsT[id][0] / m_weightsT[id][i]);  
  m_weightsT[id][0] = 1.0f;
}

};