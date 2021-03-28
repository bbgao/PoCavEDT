/*
 * MolecularSurface.h
 *
 *  Created on: Nov 11, 2013
 *      Author: daberdaku
 */
/**
 * This header defines the MolecularSurface object, representing
 * the surface of a molecule, with the related methods for its
 * calculation.
 */
#ifndef MOLECULARSURFACE_H_
#define MOLECULARSURFACE_H_

#include "../exceptions/ParsingPDBException.h"
#include "../Geometry/Cavity.h"
#include "../Geometry/Pocket.h"
#include "../Geometry/point3D.h"
#include "../Geometry/voxelGrid.h"
#include "../PDB/atom.h"
#include "../utils/numerical_limits.h"
#include <list>
#include <queue>
#include <vector>

using namespace std;

class MolecularSurface {
	friend class Molecule;
public:
	voxelGrid *le_surface, *le_cpk_model, *le_cavities, *pe_surface, *pe_cavities, *pockets, *ses;
	uint16_t * distanceMap;
	float resolution;
	vector<Pocket> pocket_list;
	vector<Cavity> le_cavity_list, pe_cavity_list;


	MolecularSurface(vector<atom> const & atoms, float ligandRadius,
			float probeRadius, float resolution, uint16_t length,
			uint16_t width, uint16_t height, point3D translation);

	virtual ~MolecularSurface();

	static void boundingBox(vector<atom> const & atomsInModel,
			float probeRadius, float resolution, uint16_t & length,
			uint16_t & width, uint16_t & height, point3D & translation,
			float & max_atm_radius, float & min_atm_radius);
	void computePESandPockets(float minPocketDepth);
	void computeLESandCavities();

	void computeSES();

	void outputPCDModel(voxelGrid const * grid, std::string const & filename);
	void outputOpenDXModel(voxelGrid const * grid, string const & filename);
	void outputVTKStructuredPointsModel(voxelGrid const * grid, string const & filename);
	void outputVTKPolyDataModel(voxelGrid const * grid, string const & filename);

private:
	void buildCPKModel(voxelGrid ** cpk_model, bool ligExclSurf);
	void buildSurface(voxelGrid const *cpk_model, voxelGrid ** surface);
	void extractCavities(voxelGrid const * cpk_model, voxelGrid ** cavities, vector<Cavity> & cavity_list);
	void extractPockets(voxelGrid const * pe_cpk_model, float minPocketDepth, vector<Pocket> & pocket_list);


	point3D ptran;
	vector<atom> const * atoms;
	float ligandRadius, probeRadius;
	uint16_t pheight, pwidth, plength;



};
#endif /* MOLECULARSURFACE_H_ */
