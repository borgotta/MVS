#include "StdAfx.h"
#include "PatchOrganizer.h"
#include "Reconstructor.h"

#include <iostream>
#include <fstream>


namespace MVS {

	Patch PatchOrganizer::m_MAXDEPTH;
	Patch PatchOrganizer::m_BACKGROUND;

	PatchOrganizer::PatchOrganizer(Reconstructor& r) : rec(r)
	{	}

	PatchOrganizer::~PatchOrganizer(void)
	{	}

	void PatchOrganizer::init(void) {
		m_pgrids.clear();   m_pgrids.resize(rec.n_im);
		m_vpgrids.clear();  m_vpgrids.resize(rec.n_im);
		m_dpgrids.clear();  m_dpgrids.resize(rec.n_im);
		m_counts.clear();   m_counts.resize(rec.n_im);

		m_gwidths.clear();  m_gwidths.resize(rec.n_im);
		m_gheights.clear(); m_gheights.resize(rec.n_im);
		for (int index = 0; index < rec.n_im; ++index) {
			const int gwidth = (rec.ps.getWidth(index)
				+ rec.csize - 1) / rec.csize;
			const int gheight = (rec.ps.getHeight(index)
				+ rec.csize - 1) / rec.csize;
			m_gwidths[index] = gwidth;
			m_gheights[index] = gheight;

			if (index < rec.n_im) {
				m_pgrids[index].resize(gwidth * gheight);
				m_vpgrids[index].resize(gwidth * gheight);
				m_dpgrids[index].resize(gwidth * gheight);
				m_counts[index].resize(gwidth * gheight);
				fill(m_dpgrids[index].begin(), m_dpgrids[index].end(), m_MAXDEPTH);
			}
		}
		cerr<<"PatchOrganizer initialized"<<endl;
	}

	void PatchOrganizer::clearCounts(void) {
		for (int index = 0; index < rec.n_im; ++index) {
			vector<unsigned char>::iterator begin = m_counts[index].begin();
			vector<unsigned char>::iterator end = m_counts[index].end();
			while (begin != end) {
				*begin = (unsigned char)0;
				++begin;
			}
		}
	}

	void PatchOrganizer::addPatch(Patch& ppatch) {
		// First handle m_vimages
		vector<int>::iterator bimage = ppatch.m_images.begin();
		vector<int>::iterator eimage = ppatch.m_images.end();
		vector<Vec2i>::iterator bgrid = ppatch.m_grids.begin();
		while (bimage != eimage) {
			const int index = *bimage;
			if (rec.n_im <= index) {
				++bimage;      ++bgrid;
				continue;
			}

			const int index2 = (*bgrid)[1] * m_gwidths[index] + (*bgrid)[0];
			pthread_rwlock_wrlock(&rec.m_imageLocks[index]);
			m_pgrids[index][index2].push_back(ppatch);
			pthread_rwlock_unlock(&rec.m_imageLocks[index]);
			++bimage;
			++bgrid;
		}

		// If depth, set vimages
		if (rec.m_depth == 0)
			return;

		bimage = ppatch.m_vimages.begin();
		eimage = ppatch.m_vimages.end();
		bgrid = ppatch.m_vgrids.begin();

		while (bimage != eimage) {
			const int index = *bimage;
			const int index2 = (*bgrid)[1] * m_gwidths[index] + (*bgrid)[0];
			pthread_rwlock_wrlock(&rec.m_imageLocks[index]);
			m_vpgrids[index][index2].push_back(ppatch);
			pthread_rwlock_unlock(&rec.m_imageLocks[index]);
			++bimage;
			++bgrid;
		}

		updateDepthMaps(ppatch);
	}
	void PatchOrganizer::removePatch(const Patch& ppatch) {
		for (int i = 0; i < (int)ppatch.m_images.size(); ++i) {
			const int image = ppatch.m_images[i];
			if (rec.n_im <= image)
				continue;

			const int& ix = ppatch.m_grids[i][0];
			const int& iy = ppatch.m_grids[i][1];
			const int index = iy * m_gwidths[image] + ix;
			m_pgrids[image][index].erase(remove(m_pgrids[image][index].begin(),
				m_pgrids[image][index].end(),
				ppatch),
				m_pgrids[image][index].end());
		}

		for (int i = 0; i < (int)ppatch.m_vimages.size(); ++i) {
			const int image = ppatch.m_vimages[i];
#ifdef DEBUG
			if (rec.n_im <= image) {
				cerr << "Impossible in removePatch. m_vimages must be targetting images" << endl;
				exit (1);
			}
#endif

			const int& ix = ppatch.m_vgrids[i][0];
			const int& iy = ppatch.m_vgrids[i][1];
			const int index = iy * m_gwidths[image] + ix;
			m_vpgrids[image][index].erase(remove(m_vpgrids[image][index].begin(),
				m_vpgrids[image][index].end(),
				ppatch),
				m_vpgrids[image][index].end());
		}
	}

	void PatchOrganizer::setScales(Patch& patch) const {
  const float unit = rec.optim.getUnit(patch.m_images[0], patch.m_coord);
  const float unit2 = 2.0f * unit;
  Vec4f ray = patch.m_coord - rec.ps.photoList[patch.m_images[0]].m_center;
  unitize(ray);

  const int inum = min(rec.m_tau, (int)patch.m_images.size());
  
  // First compute, how many pixel difference per unit along vertical
  //for (int i = 1; i < (int)patch.m_images.size(); ++i) {
  for (int i = 1; i < inum; ++i) {
    Vec3f diff = rec.ps.project(patch.m_images[i], patch.m_coord) -
      rec.ps.project(patch.m_images[i], patch.m_coord - ray * unit2);
	
    patch.m_dscale += norm(diff);
  }

  // set m_dscale to the vertical distance where average pixel move is half pixel
  //patch.m_dscale /= (int)patch.m_images.size() - 1;
  patch.m_dscale /= inum - 1;
  patch.m_dscale = unit2 / patch.m_dscale;
  
  patch.m_ascale = atan(patch.m_dscale / (unit * rec.m_wsize / 2.0f));
}

	void PatchOrganizer::writePLY(const std::vector<Patch>& patches,
		const std::string filename) {
			ofstream ofstr;
			ofstr.open(filename.c_str());
			ofstr << "ply" << endl
				<< "format ascii 1.0" << endl
				<< "element vertex " << (int)patches.size() << endl
				<< "property float x" << endl
				<< "property float y" << endl
				<< "property float z" << endl
				<< "property float nx" << endl
				<< "property float ny" << endl
				<< "property float nz" << endl
				<< "property uchar diffuse_red" << endl
				<< "property uchar diffuse_green" << endl
				<< "property uchar diffuse_blue" << endl
				<< "end_header" << endl;

			vector<Patch>::const_iterator bpatch = patches.begin();
			vector<Patch>::const_iterator bend = patches.end();

			while (bpatch != bend) {
				// Get color
				Vec3i color;

				const int mode = 0;
				// 0: color from images
				// 1: fix
				// 2: angle
				if (mode == 0) {
					int denom = 0;
					Vec3f colorf;
					for (int i = 0; i < (int)(bpatch)->m_images.size(); ++i) {
						const int image = (bpatch)->m_images[i];
						colorf += rec.ps.getPhoto(image).getColor(bpatch->m_coord);
						denom++;
					}
					colorf /= denom;
					color[0] = min(255,(int)floor(colorf[0] + 0.5f));
					color[1] = min(255,(int)floor(colorf[1] + 0.5f));
					color[2] = min(255,(int)floor(colorf[2] + 0.5f));
				}
				else if (mode == 1) {
					if ((bpatch)->m_tmp == 1.0f) {
						color[0] = 255;
						color[1] = 0;
						color[2] = 0;
					}
					else {
						color[0] = 255;
						color[1] = 255;
						color[2] = 255;
					}
				}
				else if (mode == 2) {
					float angle = 0.0f;
					vector<int>::const_iterator bimage = (bpatch)->m_images.begin();
					vector<int>::const_iterator eimage = (bpatch)->m_images.end();

					while (bimage != eimage) {
						const int index = *bimage;
						Vec4f ray = rec.ps.photoList[index].m_center - (bpatch)->m_coord;
						ray[3] = 0.0f;

						unitize(ray);

						angle += acos(mult(ray, (bpatch)->m_normal));
						++bimage;
					}

					angle = angle / (CV_PI / 2.0f);
					float r, g, b;
					gray2rgb(angle, r, g, b);
					color[0] = (int)(r * 255.0f);
					color[1] = (int)(g * 255.0f);
					color[2] = (int)(b * 255.0f);
				}

				ofstr << (bpatch)->m_coord[0] << ' '
					<< (bpatch)->m_coord[1] << ' '
					<< (bpatch)->m_coord[2] << ' '
					<< (bpatch)->m_normal[0] << ' '
					<< (bpatch)->m_normal[1] << ' '
					<< (bpatch)->m_normal[2] << ' '
					<< color[0] << ' ' << color[1] << ' ' << color[2] << endl;
				++bpatch;
			}
			ofstr.close();  
	}

	void PatchOrganizer::writePLY(const std::vector<Patch>& patches,
		const std::string filename,
		const std::vector<Vec3i>& colors) {
			ofstream ofstr;
			ofstr.open(filename.c_str());
			ofstr << "ply" << endl
				<< "format ascii 1.0" << endl
				<< "element vertex " << (int)patches.size() << endl
				<< "property float x" << endl
				<< "property float y" << endl
				<< "property float z" << endl
				<< "property float nx" << endl
				<< "property float ny" << endl
				<< "property float nz" << endl
				<< "property uchar diffuse_red" << endl
				<< "property uchar diffuse_green" << endl
				<< "property uchar diffuse_blue" << endl
				<< "end_header" << endl;

			vector<Patch>::const_iterator bpatch = patches.begin();
			vector<Patch>::const_iterator bend = patches.end();
			vector<Vec3i>::const_iterator colorb = colors.begin();

			while (bpatch != bend) {
				ofstr << (bpatch)->m_coord[0] << ' '
					<< (bpatch)->m_coord[1] << ' '
					<< (bpatch)->m_coord[2] << ' '
					<< (bpatch)->m_normal[0] << ' '
					<< (bpatch)->m_normal[1] << ' '
					<< (bpatch)->m_normal[2] << ' '
					<< *colorb << endl;
				++bpatch;
				++colorb;
			}
			ofstr.close();  
	}
}