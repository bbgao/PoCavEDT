/*
 * Molecule.cpp
 *
 *  Created on: 07/set/2014
 *      Author: sebastian
 */

#include "../utils/disclaimer.h"
#include "../utils/doubleCompare.h"
#include "../utils/makeDirectory.h"
#include "Molecule.h"
#include <string>
#include <unordered_map>
float Molecule::max_atm_radius = min_float;
float Molecule::min_atm_radius = max_float;
using namespace std;

Molecule::Molecule(float ligandRadius, float probeRadius, float resolution,
		float minPocketDepth, int topKPockets, int topKSECavities, int topKLECavities,
		bool hydrogen, bool hetatm,
		string const & inname, string const & outname,
		string const & inname_radii) :
	centroid(point3D(0, 0, 0)), inname_molecule(inname), resolution(resolution) {

	cout << "Loading PDB file: "<< inname <<"\n";
	pdbModel = new PDBModel(inname, inname_radii, hydrogen, hetatm);
	surface = NULL;

	atoms = pdbModel->atomsInModel;

	computeCentroid();

    auto const t_start = std::chrono::high_resolution_clock::now();
	cout << "Initializing parameters for " << inname << ".\n";
	cout << "Number of atoms in model: " << atoms.size() << "\n";
	MolecularSurface::boundingBox(atoms, probeRadius, resolution, length,
			width, height, translation, max_atm_radius, min_atm_radius);


	cout << "Calculating ligand excluded surface.\n";
	if (surface != NULL)
		delete surface;
	surface = new MolecularSurface(atoms, ligandRadius, probeRadius,
			resolution, length, width, height, translation);
	surface->computeLESandCavities();
	surface->computePESandPockets(minPocketDepth);
	cout << "Computation time:\t"
			<< elapsedTime(t_start, std::chrono::high_resolution_clock::now())
			<< "\n";
	num_pockets = ((surface->pocket_list.size() > topKPockets && topKPockets >= 0) ? topKPockets : (surface->pocket_list.size()));
	num_le_cavities = ((surface->le_cavity_list.size() > topKSECavities && topKSECavities >= 0) ? topKSECavities : (surface->le_cavity_list.size()));
	num_pe_cavities = ((surface->pe_cavity_list.size() > topKLECavities && topKLECavities >= 0) ? topKLECavities : (surface->pe_cavity_list.size()));


#ifndef NO_OUTPUT_TEST
	cout << "Output surface PCD model.\n";
	surface->outputPCDModel(surface->le_surface, outname + "_le");
	surface->outputPCDModel(surface->pe_surface, outname + "_pe");
	if (num_le_cavities != 0)
		surface->outputPCDModel(surface->le_cavities, outname + "_le_cavities");
	if (num_pe_cavities != 0)
		surface->outputPCDModel(surface->pe_cavities, outname + "_pe_cavities");
	if (num_pockets != 0)
		surface->outputPCDModel(surface->pockets, outname + "_pockets");

	for (int cc = 0; cc < num_pockets; ++cc) {
		surface->pocket_list[cc].outputPocketPCDModel(outname + "_pocket_" + to_string(cc), resolution, translation);
	}
	for (int cc = 0; cc < num_le_cavities; ++cc) {
		surface->le_cavity_list[cc].outputCavityPCDModel(outname + "_le_cavity_" + to_string(cc), resolution, translation);
	}
	for (int cc = 0; cc < num_pe_cavities; ++cc) {
		surface->pe_cavity_list[cc].outputCavityPCDModel(outname + "_pe_cavity_" + to_string(cc), resolution, translation);
	}
#endif
}


Molecule::Molecule(float solventRadius, float resolution, bool hydrogen,
		bool hetatm,
		string const & inname, string const & outname,
		string const & inname_radii) :
	centroid(point3D(0, 0, 0)), inname_molecule(inname), resolution(resolution),
	num_le_cavities(0), num_pe_cavities(0), num_pockets(0) {

	cout << "Loading PDB file: "<< inname <<"\n";
	pdbModel = new PDBModel(inname, inname_radii, hydrogen, hetatm);
	surface = NULL;

	atoms = pdbModel->atomsInModel;
	atoms_tree = kdtree_atom(&atoms);

    auto const t_start = std::chrono::high_resolution_clock::now();
	cout << "Initializing parameters for " << inname << ".\n";
	cout << "Number of atoms in model: " << atoms.size() << "\n";
	MolecularSurface::boundingBox(atoms, solventRadius, resolution, length,
			width, height, translation, max_atm_radius, min_atm_radius);


	cout << "Calculating ligand's solvent excluded surface.\n";
	if (surface != NULL)
		delete surface;
	surface = new MolecularSurface(atoms, solventRadius, solventRadius,
			resolution, length, width, height, translation);
	surface->computeSES();

	cout << "Computation time:\t"
			<< elapsedTime(t_start, std::chrono::high_resolution_clock::now())
			<< "\n";

#ifndef NO_OUTPUT_TEST
	cout << "Output surface PCD model.\n";
	surface->outputPCDModel(surface->ses, outname + "_ses");
#endif
}

void Molecule::getMoleculeDimensions(float & length, float & width, float & height) {
	if (pdbModel->atomsInModel.empty())
		throw ParsingPDBException("No atoms in PDB model.",
				"Molecule::getMoleculeDimensions", "Probably due to corrupt PDB file.");

	vector<atom> atoms;
	poseNormalization(pdbModel->atomsInModel, atoms);

	point3D minp(max_float, max_float, max_float);
	point3D maxp(min_float, min_float, min_float);

	float max_atm_radius = min_float;
	float min_atm_radius = max_float;

	for (auto const & atm : atoms) {
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
		if (min_atm_radius > atm.radius)
			min_atm_radius = atm.radius;
	}
	float d = max_atm_radius + 0.5;

	maxp += point3D(d, d, d);
	minp -= point3D(d, d, d);


	length = maxp.x - minp.x;
	width  = maxp.y - minp.y;
	height = maxp.z - minp.z;

}

float Molecule::distanceToClosestAtom(point3D const & p) {
	float s_distance;
	size_t n = atoms_tree.nnSearch(p, s_distance);
	return sqrt(s_distance);
}

Molecule::~Molecule() {
	if (surface != NULL)
		delete surface;
	if (pdbModel != NULL)
		delete pdbModel;
}

