/*
 * Molecule.h
 *
 *  Created on: 07/set/2014
 *      Author: sebastian
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "../Geometry/eig3/eig3.h"
#include "../kdtree/kdtree.h"
#include "../MolecularSurface/MolecularSurface.h"
#include "../PDB/PDBModel.h"
#include "../utils/elapsedTime.h"
#include <iostream>
#include <string.h>

using namespace std;

class Molecule {
public:
	MolecularSurface * surface;
	PDBModel * pdbModel;
	uint16_t length, width, height; // bounding box dimensions
	point3D translation; // translation vector
	vector<atom> atoms;
	float resolution;
	string inname_molecule; // input filename

	static float max_atm_radius;
	static float min_atm_radius;

	kdtree_atom atoms_tree;

	point3D centroid;
	double R[3][3];
	double RT[3][3];

	int num_pockets, num_le_cavities, num_pe_cavities;

	Molecule(float ligandRadius, float probeRadius, float resolution,
			float minPocketDepth, int topKPockets, int topKSECavities, int topKLECavities,
			bool hydrogen, bool hetatm,
			string const & inname, string const & outname,
			string const & inname_radii);

	Molecule(float solventRadius, float resolution, bool hydrogen,
			bool hetatm,
			string const & inname, string const & outname,
			string const & inname_radii);

	virtual ~Molecule();

	void getMoleculeDimensions(float & length, float & width, float & height);

	float distanceToClosestAtom(point3D const & p);


private :

	inline void computeCentroid() {
		for (auto const & atm : atoms)
			centroid += point3D(atm.x, atm.y, atm.z);
		centroid /= atoms.size();
	}

	void poseNormalization(vector<atom> const & input_atoms, vector<atom> & output_atoms) {
		/**
		 * Calculate the centroid of the atom coordinates
		 */
		size_t N = input_atoms.size();

		for (auto const & atm : input_atoms) {
			centroid += point3D(atm.x, atm.y, atm.z) / N;
		}
		/**
		 * Subtract the mean from each of the data dimensions
		 */
		vector<double> cX, cY, cZ;

		for (auto & atm : input_atoms) {
			cX.push_back(atm.x - centroid.x);
			cY.push_back(atm.y - centroid.y);
			cZ.push_back(atm.z - centroid.z);
		}

		/**
		 * Calculate the covariance matrix
		 */
		double covXX = 0, covXY = 0, covXZ = 0, covYY = 0, covYZ = 0, covZZ = 0;
		for (int ii = 0; ii < N; ++ii) {
			covXX += cX[ii] * cX[ii] / N;
			covXY += cX[ii] * cY[ii] / N;
			covXZ += cX[ii] * cZ[ii] / N;
			covYY += cY[ii] * cY[ii] / N;
			covYZ += cY[ii] * cZ[ii] / N;
			covZZ += cZ[ii] * cZ[ii] / N;
		}

		/**
		 * The covariance matrix is symmetric
		 */
		double C[3][3] = {
				{covXX, covXY, covXZ},
				{covXY, covYY, covYZ},
				{covXZ, covYZ, covZZ}
		};
		/**
		 * Eigenvalues
		 */
		double lambda[3];
		/**
		 * Eigenvectors
		 */
		eigen_decomposition(C, RT, lambda);

//		cout << "Pose normalization\n";
//		cout << "Eigenvalues:\t" << lambda[0] << "\t" << lambda[1] << "\t" << lambda[2] <<endl;

		/**
		 * Rotation matrix
		 */
	    transpose(RT, R);
//	    cout << "Rotation matrix (Eigenvectors in rows): \n";
//		cout << "\t" << R[0][0] << "\t" << R[0][1] << "\t" << R[0][2] <<endl;
//		cout << "\t" << R[1][0] << "\t" << R[1][1] << "\t" << R[1][2] <<endl;
//		cout << "\t" << R[2][0] << "\t" << R[2][1] << "\t" << R[2][2] <<endl;

		/**
		 * Apply the transformation
		 */
		output_atoms = input_atoms;

		for (auto & atm : output_atoms) {
			atm.transform(R, -centroid);
		}
	}
};

#endif /* MOLECULE_H_ */
