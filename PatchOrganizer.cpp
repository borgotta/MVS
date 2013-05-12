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
					vector<int>::iterator bimage = (bpatch)->m_images.begin();
					vector<int>::iterator eimage = (bpatch)->m_images.end();

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
					Image::Cimage::gray2rgb(angle, r, g, b);
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
				ofstr << (*bpatch)->m_coord[0] << ' '
					<< (*bpatch)->m_coord[1] << ' '
					<< (*bpatch)->m_coord[2] << ' '
					<< (*bpatch)->m_normal[0] << ' '
					<< (*bpatch)->m_normal[1] << ' '
					<< (*bpatch)->m_normal[2] << ' '
					<< *colorb << endl;
				++bpatch;
				++colorb;
			}
			ofstr.close();  
	}
}