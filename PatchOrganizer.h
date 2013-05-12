#pragma once
#include "MatLib.h"
#include "Patch.h"
#include <queue>



namespace MVS {
	

	class Reconstructor;
	class PatchOrganizer
	{
	public:
		PatchOrganizer(Reconstructor& r);
		~PatchOrganizer(void);
		void init(void);
		void writePLY(const std::vector<Patch>& patches,
			const std::string filename);
		void writePLY(const std::vector<Patch>& patches,
			const std::string filename,
			const std::vector<Vec3i>& colors);
		//----------------------------------------------------------------------
		// Widths of grids
		std::vector<int> m_gwidths;
		std::vector<int> m_gheights;

		//----------------------------------------------------------------------
		// image, grid
		std::vector<std::vector<std::vector<Patch> > > m_pgrids;  
		// image, grid
		std::vector<std::vector<std::vector<Patch> > > m_vpgrids;
		// Closest patch
		std::vector<std::vector<Patch> > m_dpgrids;

		// all the patches in the current level of m_pgrids 
		std::vector<Patch> m_ppatches;

		// Check how many times patch optimization was performed for expansion
		std::vector<std::vector<unsigned char> > m_counts;

		static Patch m_MAXDEPTH;
		static Patch m_BACKGROUND;
	protected:
		Reconstructor& rec;
	};

};