/*
 * MolecularSurface.cpp
 *
 *  Created on: Dec 23, 2013
 *      Author: daberdaku
 */
/**
 * Implementation of the MolecularSurface class. This class
 * defines the molecular surface object, with all the methods
 * needed for its calculation.
 */
#include "../Geometry/HierarchicalQueue.h"
#include "../Geometry/partialEDMap.h"
#include "../Geometry/SeedFill3D/rapid3DOuterSurfaceExtract.h"
#include "../Geometry/SeedFill3D/rapid3DPocketExtract.h"
#include "../Geometry/SeedFill3D/rapid3DSeedFill.h"
#include "../Geometry/SeedFill3D/rapid3DSurfaceExtract.h"
#include "../Geometry/Sphere/DrawSphere.h"

#include "../utils/disclaimer.h"

#include "../utils/makeDirectory.h"
#include "MolecularSurface.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

/**
 * Constructor of the class. Initializes the voxel grid and other data structures.
 *
 * \param atoms				List containing the atoms of the molecule.
 * \param probeRadius		Desired probe radius of the solvent rolling sphere.
 * \param resolution		resolution^3 = number of voxels per Å^3
 * 							false: use original atomic radius values.
 * \param length			Length of the voxel grid.
 * \param width				Width of the voxel grid.
 * \param height			Height of the voxel grid.
 * \param translation		Translation vector to the grid's coordinate system
 * 							(already scaled).
 *
 */
MolecularSurface::MolecularSurface(std::vector<atom> const & atoms,
		float ligandRadius, float probeRadius, float resolution,
		uint16_t length, uint16_t width, uint16_t height, point3D translation) :
		atoms(&atoms), ligandRadius(ligandRadius), probeRadius(probeRadius),
		resolution(resolution), plength(length), pwidth(width), pheight(height),
		ptran(translation), le_surface(NULL), pe_surface(NULL), distanceMap(NULL),
		le_cavities(NULL), pe_cavities(NULL), le_cpk_model(NULL), pockets(NULL), ses(NULL) {
	if (atoms.empty())
		throw invalid_argument("MolecularSurface::MolecularSurface() - The molecule has no atoms!");
} /* MolecularSurface() */

/**
 * Destructor of the class.
 */
MolecularSurface::~MolecularSurface() {
	if (le_surface != NULL)
		delete le_surface;
	if (pe_surface != NULL)
		delete pe_surface;
	if (distanceMap != NULL)
		delete[] distanceMap;
	if (le_cpk_model != NULL)
		delete le_cpk_model;
	if (le_cavities != NULL)
		delete le_cavities;
	if (pe_cavities != NULL)
		delete pe_cavities;
	if (pockets != NULL)
		delete pockets;
	if (ses != NULL)
		delete ses;
} /* ~MolecularSurface() */

/** Calculates a bounding box containing the molecule. This method
 * calculates the maximal and minimal coordinates reached by the
 * molecule by checking all the coordinates of the atoms composing it.
 *
 * \param atomsInModel		List containing all the model's atoms.
 * \param probeRadius		Desired probe radius of the solvent rolling sphere.
 * \param resolution		resolution^3 = number of voxels per Å^3
 * \param length			(return value) the length of the voxel grid.
 * \param width				(return value) the width of the voxel grid.
 * \param height			(return value) the height of the voxel grid.
 * \param translation		(return value) the translation vector to the
 * 							grid's coordinate system (already scaled).
 * \param max_atm_radii		(return value) the maximum atomic radii in
 * 							atomsInModel
 */
void MolecularSurface::boundingBox(std::vector<atom> const & atomsInModel,
		float probeRadius, float resolution, uint16_t & length,
		uint16_t & width, uint16_t & height, point3D & translation,
		float & max_atm_radius, float & min_atm_radius) {
	if (atomsInModel.empty())
		throw ParsingPDBException("No atoms in PDB model.",
				"MolecularSurface::boundingBox", "Probably due to corrupt PDB file.");
	point3D minp(max_float, max_float, max_float);
	point3D maxp(min_float, min_float, min_float);
	for (auto const & atm : atomsInModel) {
		if (atm.x < minp.x)
			minp.x = atm.x;
		if (atm.y < minp.y)
			minp.y = atm.y;
		if (atm.z < minp.z)
			minp.z = atm.z;
		if (atm.x > maxp.x)
			maxp.x = atm.x;
		if (atm.y > maxp.y)
			maxp.y = atm.y;
		if (atm.z > maxp.z)
			maxp.z = atm.z;
		if (max_atm_radius < atm.radius)
			max_atm_radius = atm.radius;
		if (min_atm_radius > atm.radius && atm.radius > 0)
			min_atm_radius = atm.radius;
	}
	float d = max_atm_radius + probeRadius + 0.5;

	maxp += point3D(d, d, d);
	minp -= point3D(d, d, d);

	/* transformation values */
	translation = -minp;

	/* bounding box dimensions */
	double boxLength = ceil(resolution * (maxp.x - minp.x));
	double boxWidth  = ceil(resolution * (maxp.y - minp.y));
	double boxHeight = ceil(resolution * (maxp.z - minp.z));

	if (boxLength <= UINT16_MAX && boxWidth <= UINT16_MAX && boxHeight <= UINT16_MAX) {
		length = static_cast<uint16_t>(boxLength);
		width  = static_cast<uint16_t>(boxWidth);
		height = static_cast<uint16_t>(boxHeight);
	} else {
		std::stringstream ss;
		ss << "MolecularSurface::boundingBox() - ";
		ss << "The bounding box's dimensions exceed the maximum value of " << UINT16_MAX << ". ";
		ss << "Try setting a lower \"resolution\" value.";
		throw std::invalid_argument(ss.str());
	}
} /* boundingBox() */

/** Fills the voxels in the grid occupied by the molecule (protein).
 * This method implements a space-filling algorithm which is the
 * preliminary step for our grid-based macro-molecular surface
 * generation.
 *
 * In chemistry, a space-filling model, also known as a calotte model,
 * is a type of three-dimensional molecular model where the atoms are
 * represented by spheres whose radii are proportional to the radii of
 * the atoms and whose center-to-center distances are proportional to
 * the distances between the atomic nuclei, all in the same scale.
 * Atoms of different chemical elements are usually represented by
 * spheres of different colors.
 *
 * Calotte models are distinguished from other 3D representations,
 * such as the ball-and-stick and skeletal models, by the use of "full
 * size" balls for the atoms. They are useful for visualizing the
 * effective shape and relative dimensions of the molecule, in particular
 * the region of space occupied by it. On the other hand, calotte models
 * do not show explicitly the chemical bonds between the atoms, nor the
 * structure of the molecule beyond the first layer of atoms.
 *
 * Space-filling models are also called CPK models after the chemists
 * Robert Corey, Linus Pauling and Walter Koltun, who pioneered their use.
 */
void MolecularSurface::buildCPKModel(voxelGrid ** cpk_model, bool ligExclSurf) {
	if (*cpk_model != NULL)
		delete *cpk_model;
	*cpk_model = new voxelGrid(plength, pwidth, pheight);
	int cx, cy, cz; /**< Discretized coordinates of the atom's center. */
	int radius; /**< radius */

	/* For every atom in our list, calculate the voxels it occupies. */
	for (auto const & atm : *atoms) {
		if (atm.radius <= 0.0)
			continue;
		/* Translate and discretize the coordinates */
		cx = static_cast<int>(round((atm.x + ptran.x) * resolution));
		cy = static_cast<int>(round((atm.y + ptran.y) * resolution));
		cz = static_cast<int>(round((atm.z + ptran.z) * resolution));
		if (ligExclSurf)
			radius = static_cast<int>(round((atm.radius + ligandRadius) * resolution));
		else
			radius = static_cast<int>(round((atm.radius + probeRadius) * resolution));

		DrawBall(**cpk_model, cx, cy, cz, radius);
	}
} /* buildCPKModel() */


/**
 * Build the boundary of the solid created with the space-filling algorithm.
 * A voxel belongs to the boundary if it has at least one neighbor which is not
 * occupied by any atom.
 */
void MolecularSurface::buildSurface(voxelGrid const *cpk_model, voxelGrid ** surface) {
	if (*surface != NULL)
		delete *surface;

	*surface = new voxelGrid(plength, pwidth, pheight); /**< already initialized to all zeros */
	int cx, cy, cz;

	voxelGrid temp(*cpk_model);
	/*
	 * the model could have disjoint parts,
	 * all atom centers must be used as seeds
	 */
	for (auto const & atm : *atoms) {
		if (atm.radius <= 0.0)
			continue;
		/* Translate and discretize the coordinates */
		cx = static_cast<int>(round((atm.x + ptran.x) * resolution));
		cy = static_cast<int>(round((atm.y + ptran.y) * resolution));
		cz = static_cast<int>(round((atm.z + ptran.z) * resolution));
		if (temp.getVoxel(cx, cy, cz)) {
			rapid3DSurfaceExtract::surfaceExtract3D(*cpk_model, temp, **surface, voxel(cx, cy, cz));
		}
	}
} /* buildSurface() */


void MolecularSurface::extractCavities(voxelGrid const * cpk_model, voxelGrid ** cavities, vector<Cavity> & cavity_list) {
	voxelGrid * temp = new voxelGrid(*cpk_model);
	if (*cavities != NULL)
		delete *cavities;
	rapid3DSeedFill::seedFill3D(*temp, voxel(0, 0, 0));
	*cavities = &(temp->invert());


	voxelGrid tmp = voxelGrid(**cavities);
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (tmp.getVoxel(i, j, k)) {
					voxelGrid currentCavity = voxelGrid(plength, pwidth, pheight);
					rapid3DPocketExtract::poCavExtract3D(**cavities, tmp, currentCavity, voxel(i, j, k));
					Cavity cCavity = Cavity(currentCavity, resolution, ptran);
					cavity_list.push_back(cCavity);
				}
			}
		}
	}

	sort(cavity_list.begin(), cavity_list.end(),
			[](Cavity const & a, Cavity const & b){return (a.getVolume() > b.getVolume());});
}


/** Calculates a partial Euclidean Distance Map of the voxelized
 * representation of the molecular surface. Starting from the Solvent-
 * Accessible Surface of the molecule, it is possible to obtain the
 * molecular surface (or Solvent-Excluded Surface) applying the EDT
 * (Euclidean Distance Transform). The SES is then extracted from the
 * partial Euclidean Distance Map, as it's voxels have a particular
 * distance value, i.e. greater or equal to the solvent/probe sphere
 * radius. The calculation starts from the boundary (surface) voxels,
 * and procedes towards the inside of the molecule, one shell at a time.
 * The initial shell is obviously the surface of the molecule. The
 * second shell is composed by voxels internal to the molecule, that
 * are neighbors to the first shell's voxels, and don't belong to the
 * first shell. The third shell will be composed by neighbors of the
 * second shell's voxels, internal to the molecule, and that don't
 * belong to the first or second shell, and so on. The distance map
 * values are propagated from the molecular surface towards the inside
 * of the molecule one shell at a time.
 */
void MolecularSurface::computeLESandCavities() {
	cout << "Calculating Ligand Accessible CPK Model." << endl;
	buildCPKModel(&le_cpk_model, true);
	cout << "Building Ligand Accessible surface." << endl;
	voxelGrid * sa_surface = NULL;
	buildSurface(le_cpk_model, &sa_surface);

	cout << "Calculating the Euclidean Distance Transform." << endl;
	voxel * nearestSurfaceVoxel = new voxel[plength * pwidth * pheight];
	if (distanceMap != NULL)
		delete[] distanceMap;
	distanceMap = new uint16_t[plength * pwidth * pheight];
	uint16_t min_sdist = (uint16_t)(pow(ligandRadius * resolution, 2.0) + 0.5);
	uint16_t numberOfQueues = static_cast<uint16_t>(ceil(pow((ligandRadius * resolution) + 0.5, 2.0)));

	HierarchicalQueue* HQ1 = new HierarchicalQueue(numberOfQueues);
	HierarchicalQueue* HQ2 = new HierarchicalQueue(numberOfQueues);

	/* for all voxels */
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				/* create the first shell by pushing the surface voxels in the list */
				if (sa_surface->getVoxel(i, j, k)) {
					HQ1->push(voxel(i, j, k), 0);
					/* a surface voxel has zero distance from itself */
					nearestSurfaceVoxel[(i * pwidth + j) * pheight + k] = voxel(i, j, k);
					distanceMap[(i * pwidth + j) * pheight + k] = 0;
				} else {
					distanceMap[(i * pwidth + j) * pheight + k] = UINT16_MAX;
				}
			}
		}
	}

	delete sa_surface;
	sa_surface = NULL;

	while (!HQ1->empty()) {
		/*current voxel*/
		voxel cv = HQ1->front();
		HQ1->pop();
		voxel nearestSurfVox = nearestSurfaceVoxel[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		uint16_t squaredDistance = distanceMap[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		voxel nb; // neighbour
		bool isEnd = true;
		for (int i = 0; i < 26; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (le_cpk_model->getVoxel(nb)) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				if (newSDistance < distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz]
						&& newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = nearestSurfVox;
					HQ1->push(nb, newSDistance);
					isEnd = false;
				}
			}
		}
		if (isEnd && squaredDistance >= 24) {
			HQ2->push(cv, squaredDistance);
		}
	}

	delete HQ1;
	HQ1 = NULL;

	while (!HQ2->empty()) {
		/*current voxel*/
		voxel cv = HQ2->front();
		HQ2->pop();
		voxel nearestSurfVox = nearestSurfaceVoxel[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		// neighbor
		voxel nb;
		for (int i = 0; i < 124; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (le_cpk_model->getVoxel(nb)) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				if (newSDistance < distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] && newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = nearestSurfVox;
					HQ2->push(nb, newSDistance);
				}
			}
		}
	}
	cout << "Calculating Ligand Excluded CPK Model." << endl;
	/* update the volumetric model */
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (distanceMap[(i * pwidth + j) * pheight + k] < min_sdist) {
					le_cpk_model->clearVoxel(i, j, k);
				}
			}
		}
	}

	delete[] nearestSurfaceVoxel;
	nearestSurfaceVoxel = NULL;

	delete HQ2;
	HQ2 = NULL;
	cout << "Calculating Ligand Excluded Surface." << endl;
	buildSurface(le_cpk_model, &le_surface);
	cout << "Extracting Ligand Excluded Cavities." << endl;
	extractCavities(le_cpk_model, &le_cavities, le_cavity_list);
} /*  */

void MolecularSurface::extractPockets(voxelGrid const * pe_cpk_model, float minPocketDepth, vector<Pocket> & pocket_list) {
	if (pockets != NULL)
		delete pockets;
	pockets = new voxelGrid(plength, pwidth, pheight);
	uint16_t min_sdist = (uint16_t)(pow((probeRadius + minPocketDepth) * resolution, 2.0) + 0.5);

	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (pe_cpk_model->getVoxel(i, j, k) && !le_cpk_model->getVoxel(i, j, k) && !le_cavities->getVoxel(i, j, k))
					if (distanceMap[(i * pwidth + j) * pheight + k] >= min_sdist /*&& distanceMap[(i * pwidth + j) * pheight + k] < UINT16_MAX*/)
						pockets->setVoxel(i, j, k);
			}
		}
	}

	voxelGrid temp = voxelGrid(*pockets);
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (temp.getVoxel(i, j, k)) {
					voxelGrid currentPocket = voxelGrid(plength, pwidth, pheight);
					rapid3DPocketExtract::poCavExtract3D(*pockets, temp, currentPocket, voxel(i, j, k));
					Pocket cPocket = Pocket(currentPocket, distanceMap, resolution, ptran);
					pocket_list.push_back(cPocket);
				}
			}
		}
	}

	sort(pocket_list.begin(), pocket_list.end(),
			[](Pocket const & a, Pocket const & b){return (a.getScore() > b.getScore());});
}



void MolecularSurface::computePESandPockets(float minPocketDepth) {
	cout << "Calculating Probe-sphere Accessible CPK Model." << endl;
	voxelGrid *cpk_model = NULL;
	buildCPKModel(&cpk_model, false);
	cout << "Building Probe-sphere Accessible surface." << endl;
	voxelGrid * la_surface = NULL;
	buildSurface(cpk_model, &la_surface);

	cout << "Calculating the Euclidean Distance Transform." << endl;
	voxel * nearestSurfaceVoxel = new voxel[plength * pwidth * pheight];
	if (distanceMap != NULL)
		delete[] distanceMap;
	distanceMap = new uint16_t[plength * pwidth * pheight];
	uint16_t min_sdist = (uint16_t)(pow(probeRadius * resolution, 2.0) + 0.5);
	uint16_t numberOfQueues = plength * plength;

	HierarchicalQueue* HQ1 = new HierarchicalQueue(numberOfQueues);
	HierarchicalQueue* HQ2 = new HierarchicalQueue(numberOfQueues);

	/* for all voxels */
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				/* create the first shell by pushing the surface voxels in the list */
				if (la_surface->getVoxel(i, j, k)) {
					HQ1->push(voxel(i, j, k), 0);
					/* a surface voxel has zero distance from itself */
					nearestSurfaceVoxel[(i * pwidth + j) * pheight + k] = voxel(i, j, k);
					distanceMap[(i * pwidth + j) * pheight + k] = 0;
				} else {
					distanceMap[(i * pwidth + j) * pheight + k] = UINT16_MAX;
				}
			}
		}
	}
	delete la_surface;
	la_surface = NULL;

	while (!HQ1->empty()) {
		/*current voxel*/
		voxel cv = HQ1->front();
		HQ1->pop();
		voxel nearestSurfVox = nearestSurfaceVoxel[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		uint16_t squaredDistance = distanceMap[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		voxel nb; // neighbour
		bool isEnd = true;
		for (int i = 0; i < 26; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (cpk_model->getVoxel(nb) && !(le_cpk_model->getVoxel(nb))) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				if (newSDistance < distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz]
						&& newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = nearestSurfVox;
					HQ1->push(nb, newSDistance);
					isEnd = false;
				}
			}
		}
		if (isEnd && squaredDistance >= 24) {
			HQ2->push(cv, squaredDistance);
		}
	}

	delete HQ1;
	HQ1 = NULL;

	while (!HQ2->empty()) {
		/*current voxel*/
		voxel cv = HQ2->front();
		HQ2->pop();
		voxel nearestSurfVox = nearestSurfaceVoxel[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		// neighbor
		voxel nb;
		for (int i = 0; i < 124; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (cpk_model->getVoxel(nb) && !le_cpk_model->getVoxel(nb)) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				if (newSDistance < distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz]
						&& newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = nearestSurfVox;
					HQ2->push(nb, newSDistance);
				}
			}
		}
	}
	cout << "Calculating Probe-sphere Excluded CPK Model." << endl;
	/* update the volumetric model */
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (distanceMap[(i * pwidth + j) * pheight + k] < min_sdist) {
					cpk_model->clearVoxel(i, j, k);
				}
			}
		}
	}

	delete[] nearestSurfaceVoxel;
	nearestSurfaceVoxel = NULL;

	delete HQ2;
	HQ2 = NULL;

	cout << "Calculating Probe-sphere Excluded Surface." << endl;
	buildSurface(cpk_model, &pe_surface);
	cout << "Extracting Probe-sphere Excluded Cavities." << endl;
	extractCavities(cpk_model, &pe_cavities, pe_cavity_list);
	cout << "Extracting potential pocket sites." << endl;
	extractPockets(cpk_model, minPocketDepth, pocket_list);
	delete cpk_model;
	cpk_model = NULL;
} /*  */

void MolecularSurface::computeSES() {
	cout << "Calculating Solvent Accessible CPK Model." << endl;
	voxelGrid *cpk_model = NULL;
	buildCPKModel(&cpk_model, false);
	cout << "Building Solvent Accessible surface." << endl;
	voxelGrid * sa_surface = NULL;
	buildSurface(cpk_model, &sa_surface);

	cout << "Calculating the Euclidean Distance Transform." << endl;
	voxel * nearestSurfaceVoxel = new voxel[plength * pwidth * pheight];
	if (distanceMap != NULL)
		delete[] distanceMap;
	distanceMap = new uint16_t[plength * pwidth * pheight];
	uint16_t min_sdist = (uint16_t)(pow(probeRadius * resolution, 2.0) + 0.5);
	uint16_t numberOfQueues = plength * plength;

	HierarchicalQueue* HQ1 = new HierarchicalQueue(numberOfQueues);
	HierarchicalQueue* HQ2 = new HierarchicalQueue(numberOfQueues);

	/* for all voxels */
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				/* create the first shell by pushing the surface voxels in the list */
				if (sa_surface->getVoxel(i, j, k)) {
					HQ1->push(voxel(i, j, k), 0);
					/* a surface voxel has zero distance from itself */
					nearestSurfaceVoxel[(i * pwidth + j) * pheight + k] = voxel(i, j, k);
					distanceMap[(i * pwidth + j) * pheight + k] = 0;
				} else {
					distanceMap[(i * pwidth + j) * pheight + k] = UINT16_MAX;
				}
			}
		}
	}
	delete sa_surface;
	sa_surface = NULL;

	while (!HQ1->empty()) {
		/*current voxel*/
		voxel cv = HQ1->front();
		HQ1->pop();
		voxel nearestSurfVox = nearestSurfaceVoxel[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		uint16_t squaredDistance = distanceMap[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		voxel nb; // neighbour
		bool isEnd = true;
		for (int i = 0; i < 26; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (cpk_model->getVoxel(nb)) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				if (newSDistance < distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz]
						&& newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = nearestSurfVox;
					HQ1->push(nb, newSDistance);
					isEnd = false;
				}
			}
		}
		if (isEnd && squaredDistance >= 24) {
			HQ2->push(cv, squaredDistance);
		}
	}

	delete HQ1;
	HQ1 = NULL;

	while (!HQ2->empty()) {
		/*current voxel*/
		voxel cv = HQ2->front();
		HQ2->pop();
		voxel nearestSurfVox = nearestSurfaceVoxel[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		// neighbor
		voxel nb;
		for (int i = 0; i < 124; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (cpk_model->getVoxel(nb)) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				if (newSDistance < distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz]
						&& newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = nearestSurfVox;
					HQ2->push(nb, newSDistance);
				}
			}
		}
	}
	cout << "Calculating Solvent Excluded CPK Model." << endl;
	/* update the volumetric model */
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (distanceMap[(i * pwidth + j) * pheight + k] < min_sdist) {
					cpk_model->clearVoxel(i, j, k);
				}
			}
		}
	}

	delete[] nearestSurfaceVoxel;
	nearestSurfaceVoxel = NULL;

	delete HQ2;
	HQ2 = NULL;

	cout << "Calculating Solvent Excluded Surface." << endl;
	buildSurface(cpk_model, &ses);
	delete cpk_model;
	cpk_model = NULL;
} /*  */

/** Prints the 3D voxelized representation to file using the PCD
 * (Point Cloud Data) file format.
 *
 * Each PCD file contains a header that identifies and declares
 * certain properties of the point cloud data stored in the file.
 * The header of a PCD must be encoded in ASCII.
 * Storing point cloud data in both a simple ascii form with each
 * point on a line, space or tab separated, without any other
 * characters on it, as well as in a binary dump format, allows
 * us to have the best of both worlds: simplicity and speed,
 * depending on the underlying application. The ascii format
 * allows users to open up point cloud files and plot them using
 * standard software tools like gnuplot or manipulate them using
 * tools like sed, awk, etc.
 *
 * For a detailed description of the PCD (Point Cloud Data) file
 * format specification see:
 * http://pointclouds.org/documentation/tutorials/pcd_file_format.php
 *
 * \param filename	Name of the output file. The '.pcd' extension
 * 					is added automatically.
 *
 * \throws ofstream::failure
 */
void MolecularSurface::outputPCDModel(voxelGrid const * grid, std::string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open("./output/" + filename + ".pcd");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	file_stream  << std::setprecision(5);
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (grid->getVoxel(i, j, k))
					surfaceVoxels.push(voxel(i, j, k));
			}
		}
	}
	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "# created with " << PROGRAM_NAME << "\n"
			<< "# version " << PROGRAM_VERSION << "\n"
			<< "VERSION .7\n" << "FIELDS x y z\n" << "SIZE 4 4 4\n"
			<< "TYPE F F F\n" << "COUNT 1 1 1\n" << "WIDTH "
			<< surfaceVoxels.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << surfaceVoxels.size()
			<< "\n" << "DATA ascii";

	while (!surfaceVoxels.empty()) {
		file_stream << "\n" << surfaceVoxels.front().ix / resolution - ptran.x
				<< " " << surfaceVoxels.front().iy / resolution - ptran.y
				<< " " << surfaceVoxels.front().iz / resolution - ptran.z;
		surfaceVoxels.pop();
	}
	file_stream.close();
} /* outputSurfacePCDModel() */



void MolecularSurface::outputOpenDXModel(voxelGrid const * grid, std::string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open("./output/" + filename + ".dx");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	size_t gridSize = plength * pwidth * pheight;

	file_stream  << std::setprecision(5);
	/* File header */
	file_stream << "# created with " << PROGRAM_NAME << "\n";
	file_stream << "# version " << PROGRAM_VERSION << "\n";
	file_stream << "object 1 class gridpositions counts " << plength << " "	<< pwidth << " " << pheight << "\n";
	file_stream << "origin " << -ptran.x << " " << -ptran.y << " " << -ptran.z << "\n";
	file_stream << "delta " << 1/resolution << " " << 0.0f << " " << 0.0f << "\n";
	file_stream << "delta " << 0.0f << " " << 1/resolution << " " << 0.0f << "\n";
	file_stream << "delta " << 0.0f << " " << 0.0f << " " << 1/resolution << "\n";
	file_stream << "object 2 class gridpositions counts " << plength << " "	<< pwidth << " " << pheight << "\n";
	file_stream << "object 3 class array type double rank 0 items "<< gridSize << " data follows\n";
	/* data */
	size_t c = 1;
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (grid->getVoxel(i, j, k))
					file_stream << 1.0f;
				else
					file_stream << 0.0f;
				if (c % 3 == 0 || c == gridSize) {
					file_stream << "\n";
				}
				else {
					file_stream << " ";
				}
				++c;
			}
		}
	}
	file_stream << "attribute \"dep\" string \"positions\"\n";
	file_stream	<< "object \"regular positions regular connections\" class field\n";
	file_stream << "component \"positions\" value 1\n";
	file_stream << "component \"connections\" value 2\n";
	file_stream << "component \"data\" value 3\n";
	file_stream.close();
} /* outputSurfaceOpenDXModel() */
void MolecularSurface::outputVTKStructuredPointsModel(voxelGrid const * grid, std::string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open("./output/" + filename + ".vtk");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}

	file_stream  << std::setprecision(5);
	/* File header */
	file_stream << "# vtk DataFile Version 3.0\n";
	file_stream << filename << ".vtk created with " << PROGRAM_NAME << " - version " << PROGRAM_VERSION << ".\n";
	file_stream << "ASCII\n";
	file_stream << "DATASET STRUCTURED_POINTS\n";
	file_stream << "DIMENSIONS " << plength << " " << pwidth << " " << pheight << "\n";
	file_stream << "ORIGIN " << -ptran.x << " " << -ptran.y << " " << -ptran.z << "\n";
	file_stream << "SPACING " << 1/resolution<< " " << 1/resolution<< " " << 1/resolution<< "\n";
	file_stream << "POINT_DATA " << plength * pwidth * pheight << "\n";
	file_stream << "SCALARS " << filename << " int 1\n";
	file_stream << "LOOKUP_TABLE default";

	/* data */
	for (int k = 0; k < pheight; ++k) {
		for (int j = 0; j < pwidth; ++j) {
			for (int i = 0; i < plength; ++i) {
				file_stream << "\n";
				if (grid->getVoxel(i, j, k))
					file_stream << "1";
				else
					file_stream << "0";
			}
		}
	}

	file_stream.close();
}
void MolecularSurface::outputVTKPolyDataModel(voxelGrid const * grid, std::string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open("./output/" + filename + ".vtk");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	file_stream  << std::setprecision(5);
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (grid->getVoxel(i, j, k))
					surfaceVoxels.push(voxel(i, j, k));
			}
		}
	}
	/* File header */
	file_stream << "# vtk DataFile Version 3.0\n";
	file_stream << filename << ".vtk created with " << PROGRAM_NAME << " - version " << PROGRAM_VERSION <<".\n";
	file_stream << "ASCII\n";
	file_stream << "DATASET POLYDATA\n";
	file_stream << "POINTS " << surfaceVoxels.size() << " float";

	/* data */
	while (!surfaceVoxels.empty()) {
		file_stream << "\n" << surfaceVoxels.front().ix / resolution - ptran.x
				<< " " << surfaceVoxels.front().iy / resolution - ptran.y
				<< " " << surfaceVoxels.front().iz / resolution - ptran.z;
		surfaceVoxels.pop();
	}

	file_stream.close();
} /* outputSurfaceVTKModel() */

