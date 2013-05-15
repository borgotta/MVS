#include "StdAfx.h"
#include "Reconstructor.h"
#include "InitialMatcher.h"
#include <numeric>

namespace MVS {

	InitialMatcher::InitialMatcher(Reconstructor& r) : rec(r) {}

	void InitialMatcher::init(const std::vector<std::vector<Point> >& points) {
		m_Points.clear();
		m_Points.resize(rec.n_im);

		for (int index = 0; index < rec.n_im; ++index) {
			const int gheight = rec.po.m_gheights[index];
			const int gwidth = rec.po.m_gwidths[index];
			m_Points[index].resize(gwidth * gheight);
		}

		readPoints(points);
		cerr<<"InitialMatcher initialized"<<endl;
	}

	void InitialMatcher::readPoints(const std::vector<std::vector<Point> >& points) {
		for (int index = 0; index < rec.n_im; ++index) {
			for (int i = 0; i < (int)points[index].size(); ++i) {
				Point point(points[index][i]);
				point.m_itmp = index;
				const int ix = ((int)floor(point.m_icoord[0] + 0.5f)) / rec.csize;
				const int iy = ((int)floor(point.m_icoord[1] + 0.5f)) / rec.csize;
				const int index2 = iy * rec.po.m_gwidths[index] + ix;
				m_Points[index][index2].push_back(point);
			}
		}
	}



	void InitialMatcher::run(void) {
  rec.m_count = 0;
  rec.m_jobs.clear();
  m_scounts.resize(rec.settings.CPU);
  m_fcounts0.resize(rec.settings.CPU);
  m_fcounts1.resize(rec.settings.CPU);
  m_pcounts.resize(rec.settings.CPU);
  fill(m_scounts.begin(), m_scounts.end(), 0);
  fill(m_fcounts0.begin(), m_fcounts0.end(), 0);
  fill(m_fcounts1.begin(), m_fcounts1.end(), 0);
  fill(m_pcounts.begin(), m_pcounts.end(), 0);
  
  vector<int> vitmp;
  for (int i = 0; i < rec.n_im; ++i)
	vitmp.push_back(i);

  random_shuffle(vitmp.begin(), vitmp.end());
  rec.m_jobs.insert(rec.m_jobs.end(), vitmp.begin(), vitmp.end());

  cerr << "adding seeds " << endl;
  
  rec.po.clearCounts(); //todo

  // If there already exists a patch, don't use
  for (int index = 0; index < (int)rec.n_im; ++index) {
	for (int j = 0; j < (int)rec.po.m_pgrids[index].size(); ++j) {
	  if (!rec.po.m_pgrids[index][j].empty())
		rec.po.m_counts[index][j] = rec.m_countThreshold2;
	}
  }

  time_t tv;
  time(&tv);
  time_t curtime = tv;
  vector<pthread_t> threads(rec.settings.CPU);
  for (int i = 0; i < rec.settings.CPU; ++i)
	pthread_create(&threads[i], NULL, initialMatchThreadTmp, (void*)this);
  for (int i = 0; i < rec.settings.CPU; ++i)
	pthread_join(threads[i], NULL);
  //----------------------------------------------------------------------
  cerr << "done" << endl;
  time(&tv);
  cerr << "---- Initial: " << (tv - curtime)/CLOCKS_PER_SEC << " secs ----" << endl;

  const int trial = accumulate(m_scounts.begin(), m_scounts.end(), 0);
  const int fail0 = accumulate(m_fcounts0.begin(), m_fcounts0.end(), 0);
  const int fail1 = accumulate(m_fcounts1.begin(), m_fcounts1.end(), 0);
  const int pass = accumulate(m_pcounts.begin(), m_pcounts.end(), 0);
  cerr << "Total pass fail0 fail1 refinepatch: "
	   << trial << ' ' << pass << ' '
	   << fail0 << ' ' << fail1 << ' ' << pass + fail1 << endl;
  cerr << "Total pass fail0 fail1 refinepatch: "
	   << 100 * trial / (float)trial << ' '
	   << 100 * pass / (float)trial << ' '
	   << 100 * fail0 / (float)trial << ' '
	   << 100 * fail1 / (float)trial << ' '
	   << 100 * (pass + fail1) / (float)trial << endl;
}

void InitialMatcher::initialMatchThread(void) {
  pthread_rwlock_wrlock(&rec.m_lock);
  const int id = rec.m_count++;
  pthread_rwlock_unlock(&rec.m_lock);

  while (1) {
	int index = -1;
	pthread_rwlock_wrlock(&rec.m_lock);
	if (!rec.m_jobs.empty()) {
	  index = rec.m_jobs.front();
	  rec.m_jobs.pop_front();
	}
	pthread_rwlock_unlock(&rec.m_lock);
	if (index == -1)
	  break;

	initialMatch(index, id);
 }
}

void* InitialMatcher::initialMatchThreadTmp(void* arg) {
  ((InitialMatcher*)arg)->initialMatchThread();
  return NULL;
}

void InitialMatcher::clear(void) {
  vector<vector<vector<Point> > >().swap(m_Points);
}

void InitialMatcher::initialMatch(const int index, const int id) {
  vector<int> indexes;
  rec.optim.collectImages(index, indexes);

  if (rec.m_tau < (int)indexes.size())
	indexes.resize(rec.m_tau);
  
  if (indexes.empty())
	return;  

  int totalcount = 0;
  //======================================================================
  // for each feature point, starting from the optical center, keep on
  // matching until we find candidateThreshold patches
  const int gheight = rec.po.m_gheights[index];
  const int gwidth = rec.po.m_gwidths[index];

  int index2 = -1;
  for (int y = 0; y < gheight; ++y) {
	for (int x = 0; x < gwidth; ++x) {
	  ++index2;
	  if (!canAdd(index, x, y))
	continue;

	  for (int p = 0; p < (int)m_Points[index][index2].size(); ++p) {
	// collect features that satisfies epipolar geometry
	// constraints and sort them according to the differences of
	// distances between two cameras.
	vector<Point> vcp;
	collectCandidates(index, indexes,
						  m_Points[index][index2][p], vcp);//*m_Points[index][index2][p]
		
	int count = 0;
	Patch bestpatch;
	//======================================================================
	for (int i = 0; i < (int)vcp.size(); ++i) {
	  Patch patch;
	  patch.m_coord = vcp[i].m_coord;
	  patch.m_normal =
		  rec.ps.photoList[index].m_center - patch.m_coord;

	  unitize(patch.m_normal);
	  patch.m_normal[3] = 0.0;
	  patch.m_flag = 0;

		  ++rec.po.m_counts[index][index2];
		  const int ix = ((int)floor(vcp[i].m_icoord[0] + 0.5f)) / rec.csize;
		  const int iy = ((int)floor(vcp[i].m_icoord[1] + 0.5f)) / rec.csize;
		  const int index3 = iy * rec.po.m_gwidths[vcp[i].m_itmp] + ix;
		  if (vcp[i].m_itmp < rec.n_im)
			++rec.po.m_counts[vcp[i].m_itmp][index3];
		  
	  const int flag = initialMatchSub(index, vcp[i].m_itmp, id, patch);
	  if (flag == 0) {
		++count;
		if (bestpatch.score(rec.m_nccThreshold) <
				patch.score(rec.m_nccThreshold))
		  bestpatch = patch;
		if (rec.m_countThreshold0 <= count)
		  break;
	  }
		}
	if (count != 0) {
	  Patch ppatch(bestpatch);
	  rec.po.addPatch(ppatch);
	  ++totalcount;
		  break;
	}
	  }
	}
  }
  cerr << '(' << index << ',' << totalcount << ')' << flush;
}

void InitialMatcher::collectCells(const int index0, const int index1,
						 const Point& p0, std::vector<Vec2i>& cells) {
  Vec3f point(p0.m_icoord[0], p0.m_icoord[1], p0.m_icoord[2]);
#ifdef DEBUG
  if (p0.m_icoord[2] != 1.0f) {
	cerr << "Impossible in collectCells" << endl;    exit (1);
  }
#endif
  
  Mat F(3,3,CV_32F); //Mat3
  setF(rec.ps.photoList[index0], rec.ps.photoList[index1],
			  F);
  const int gwidth = rec.po.m_gwidths[index1];
  const int gheight = rec.po.m_gheights[index1];
  
  Mat F_t;
  transpose(F,F_t);
  Vec3f line;
  Mat m_line =  F_t * Mat(point).t(); //may be bugged
  m_line.row(0).copyTo(line); //too

  if (line[0] == 0.0 && line[1] == 0.0) {
	cerr << "Point right on top of the epipole?"
		 << index0 << ' ' << index1 << endl;
	return;
  }
  // vertical
  if (fabs(line[0]) > fabs(line[1])) {
	for (int y = 0; y < gheight; ++y) {
	  const float fy = (y + 0.5) * rec.csize - 0.5f;
	  float fx = (- line[1] * fy - line[2]) / line[0];
	  fx = max((float)(INT_MIN + 3.0f), std::min((float)(INT_MAX - 3.0f), fx));
	  
	  const int ix = ((int)floor(fx + 0.5f)) / rec.csize;
	  if (0 <= ix && ix < gwidth)
		cells.push_back(Vec2i(ix, y));
	  if (0 <= ix - 1 && ix - 1 < gwidth)
		cells.push_back(Vec2i(ix - 1, y));
	  if (0 <= ix + 1 && ix + 1 < gwidth)
		cells.push_back(Vec2i(ix + 1, y));
	}
  }
  else {
	for (int x = 0; x < gwidth; ++x) {
	  const float fx = (x + 0.5) * rec.csize - 0.5f;
	  float fy = (- line[0] * fx - line[2]) / line[1];
	  fy = max((float)(INT_MIN + 3.0f), std::min((float)(INT_MAX - 3.0f), fy));
	  
	  const int iy = ((int)floor(fy + 0.5f)) / rec.csize;
	  if (0 <= iy && iy < gheight)
		cells.push_back(Vec2i(x, iy));
	  if (0 <= iy - 1 && iy - 1 < gheight)
		cells.push_back(Vec2i(x, iy - 1));
	  if (0 <= iy + 1 && iy + 1 < gheight)
		cells.push_back(Vec2i(x, iy + 1));
	}
  }
}

// make sorted array of feature points in images, that satisfy the
// epipolar geometry coming from point in image
void InitialMatcher::collectCandidates(const int index, const std::vector<int>& indexes,
							  const Point& point, std::vector<Point>& vcp) {
  const Vec3f p0(point.m_icoord[0], point.m_icoord[1], 1.0);
  for (int i = 0; i < (int)indexes.size(); ++i) {        
	const int indexid = indexes[i];
	
	vector<Vec2i> cells;
	collectCells(index, indexid, point, cells);
	Mat F(3,3,CV_32F); //Mat3
	setF(rec.ps.photoList[index], rec.ps.photoList[indexid],
				F);
	
	for (int i = 0; i < (int)cells.size(); ++i) {
	  const int x = cells[i][0];      const int y = cells[i][1];
	  if (!canAdd(indexid, x, y))
	continue;
	  const int index2 = y * rec.po.m_gwidths[indexid] + x;

	  vector<Point>::iterator begin = m_Points[indexid][index2].begin();
	  vector<Point>::iterator end = m_Points[indexid][index2].end();
	  while (begin != end) {
		Point& rhs = *begin; //**
		// ? use type to reject candidates?
		if (point.m_type != rhs.m_type) {
		  ++begin;
		  continue;
		}
		  
		const Vec3f p1(rhs.m_icoord[0], rhs.m_icoord[1], 1.0);
		if (rec.m_epThreshold <= computeEPD(F, p0, p1)) {
		  ++begin;          
		  continue;
		}
		vcp.push_back(*begin);
		++begin;
	  }
	}
  }
  
  // set distances to m_response
  vector<Point> vcptmp;
  for (int i = 0; i < (int)vcp.size(); ++i) {
	unproject(index, vcp[i].m_itmp, point, vcp[i], vcp[i].m_coord);
	
	if ((mult(rec.ps.photoList[index].m_projection.row(2) , // m_... было другим
		vcp[i].m_coord)) <= 0.0)
	  continue;

	//if (m_fm.m_pss.getMask(vcp[i].m_coord, m_fm.m_level) == 0 ||
	//	m_fm.insideBimages(vcp[i].m_coord) == 0)
	//  continue;

	//??? from the closest
	vcp[i].m_response =
		fabs(norm(vcp[i].m_coord - rec.ps.photoList[index].m_center) -
		   norm(vcp[i].m_coord - rec.ps.photoList[vcp[i].m_itmp].m_center));
	
	vcptmp.push_back(vcp[i]);
  }
  vcptmp.swap(vcp);
  sort(vcp.begin(), vcp.end());
}

int InitialMatcher::canAdd(const int index, const int x, const int y) {
  /*if (!m_fm.m_pss.getMask(index, m_fm.m_csize * x, m_fm.m_csize * y, m_fm.m_level))
	return 0;*/

  const int index2 = y * rec.po.m_gwidths[index] + x;

  if (rec.n_im <= index)
	return 1;
  
  // Check if m_pgrids already contains something
  if (!rec.po.m_pgrids[index][index2].empty())
	return 0;

  //??? critical
  if (rec.m_countThreshold2 <= rec.po.m_counts[index][index2])
	return 0;
  
  return 1;
}

void InitialMatcher::unproject(const int index0, const int index1,
					  const Point& p0, const Point& p1,
					  Vec4f& coord) const{
  Mat A(4,4,CV_32F); //Mat4
  A.at<float>(0,0) =
	rec.ps.photoList[index0].m_projection.at<float>(0,0) -
	p0.m_icoord[0] * rec.ps.photoList[index0].m_projection.at<float>(2,0);
  A.at<float>(0,1) =
	rec.ps.photoList[index0].m_projection.at<float>(0,1) -
	p0.m_icoord[0] * rec.ps.photoList[index0].m_projection.at<float>(2,1);
  A.at<float>(0,2) =
	rec.ps.photoList[index0].m_projection.at<float>(0,2) -
	p0.m_icoord[0] * rec.ps.photoList[index0].m_projection.at<float>(2,2);
  A.at<float>(1,0) =
	rec.ps.photoList[index0].m_projection.at<float>(1,0) -
	p0.m_icoord[1] * rec.ps.photoList[index0].m_projection.at<float>(2,0);
  A.at<float>(1,1) =
	rec.ps.photoList[index0].m_projection.at<float>(1,1) -
	p0.m_icoord[1] * rec.ps.photoList[index0].m_projection.at<float>(2,1);
  A.at<float>(1,2) =
	rec.ps.photoList[index0].m_projection.at<float>(1,2) -
	p0.m_icoord[1] * rec.ps.photoList[index0].m_projection.at<float>(2,2);
  A.at<float>(2,0) =
	rec.ps.photoList[index1].m_projection.at<float>(0,0) -
	p1.m_icoord[0] * rec.ps.photoList[index1].m_projection.at<float>(2,0);
  A.at<float>(2,1) =
	rec.ps.photoList[index1].m_projection.at<float>(0,1) -
	p1.m_icoord[0] * rec.ps.photoList[index1].m_projection.at<float>(2,1);
  A.at<float>(2,2) =
	rec.ps.photoList[index1].m_projection.at<float>(0,2) -
	p1.m_icoord[0] * rec.ps.photoList[index1].m_projection.at<float>(2,2);
  A.at<float>(3,0) =
	rec.ps.photoList[index1].m_projection.at<float>(1,0) -
	p1.m_icoord[1] * rec.ps.photoList[index1].m_projection.at<float>(2,0);
  A.at<float>(3,1) =
	rec.ps.photoList[index1].m_projection.at<float>(1,1) -
	p1.m_icoord[1] * rec.ps.photoList[index1].m_projection.at<float>(2,1);
  A.at<float>(3,2) =
	rec.ps.photoList[index1].m_projection.at<float>(1,2) -
	p1.m_icoord[1] * rec.ps.photoList[index1].m_projection.at<float>(2,2);

  Vec4f b;
  b[0] =
	p0.m_icoord[0] * rec.ps.photoList[index0].m_projection.at<float>(2,3) -
	rec.ps.photoList[index0].m_projection.at<float>(0,3);
  b[1] =
	p0.m_icoord[1] * rec.ps.photoList[index0].m_projection.at<float>(2,3) -
	rec.ps.photoList[index0].m_projection.at<float>(1,3);
  b[2] =
	p1.m_icoord[0] * rec.ps.photoList[index1].m_projection.at<float>(2,3) -
	rec.ps.photoList[index1].m_projection.at<float>(0,3);
  b[3] =
	p1.m_icoord[1] * rec.ps.photoList[index1].m_projection.at<float>(2,3) -
	rec.ps.photoList[index1].m_projection.at<float>(1,3);

  Mat AT;
  transpose(A, AT);
  Mat ATA = AT * A;
  Vec4f ATb;
  Mat m_ATb = AT * (Mat(b).t());//TODO may be bugged
  m_ATb.row(0).copyTo(ATb);

  Mat ATA3(3,3,CV_32F); //3x3

  for (int y = 0; y < 3; ++y)
	for (int x = 0; x < 3; ++x)
	  ATA3.at<float>(y,x) = ATA.at<float>(y,x);
  Vec3f ATb3;
  for (int y = 0; y < 3; ++y)
	ATb3[y] = ATb[y];
  
  Mat iATA3; //3x3
  invert(ATA3, iATA3);
  Vec3f ans;
  Mat m_ans = iATA3 * (Mat(ATb3).t());
  m_ans.row(0).copyTo(ans); //TODO may be bugged
  for (int y = 0; y < 3; ++y)
	coord[y] = ans[y];
  coord[3] = 1.0f;
}		       

// starting with (index, indexs), set visible images by looking at correlation.
int InitialMatcher::initialMatchSub(const int index0, const int index1,
						   const int id, Patch& patch) {
  //----------------------------------------------------------------------
  patch.m_images.clear();
  patch.m_images.push_back(index0);
  patch.m_images.push_back(index1);

  ++m_scounts[id];

  //----------------------------------------------------------------------
  // We know that patch.m_coord is inside bimages and inside mask
  if (rec.optim.preProcess(patch, id, 1)) {
	++m_fcounts0[id];
	return 1;
  }
  
  //----------------------------------------------------------------------  
  rec.optim.refinePatch(patch, id, 100);

  //----------------------------------------------------------------------
  if (rec.optim.postProcess(patch, id, 1)) {
	++m_fcounts1[id];
	return 1;
  }
  
  ++m_pcounts[id];
  //----------------------------------------------------------------------
  return 0;
}


};