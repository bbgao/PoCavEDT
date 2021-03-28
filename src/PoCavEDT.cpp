/*
 * PoCavEDT.cpp
 *
 *  Created on: 9/jan/2018
 *      Author: Sebastian Daberdaku
 */
#include "Molecule/Molecule.h"

#include "CommandLineParser/CommandLineParser.h"
#include "utils/disclaimer.h"
#include "utils/makeDirectory.h"
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <utility>
#include "exceptions/ParsingPDBException.h"
using namespace std;

int main (int argc, char* argv[]) {
//---------------------------variables and parameters------------------------------------
	float ligandRadius; // ligand-probe sphere radius
	float probeRadius; // probe-sphere radius

	float resolution; // resolution^3 = #voxels/Å^3
	string inname_protein; // input filename
	string outname_protein; // output filename
	string inname_radii; // input atom radii

	string inname_ligand;
	string outname_ligand;
	float solventRadius;


	float minPocketDepth;
	int topKPockets;
	int topKLECavities;
	int topKPECavities;

	bool help = false; // if true, print the help message
	bool version = false; // if true, print the program version
	bool license = false; // if true, print the license information
	bool evaluate_prediction;
	bool hydrogen;
	bool hetatm;
    auto const t_start = std::chrono::high_resolution_clock::now();

	if (!CommandLineParser::parseCommandLine(argc, argv,
			ligandRadius, probeRadius, solventRadius, resolution,
			inname_protein, outname_protein, inname_ligand, outname_ligand, inname_radii,
			minPocketDepth, topKPockets, topKLECavities, topKPECavities,
			hydrogen, hetatm, evaluate_prediction,
			help, version, license)) {
		return EXIT_FAILURE;
	}
	if (help || version || license) {
		return EXIT_SUCCESS;
	}
//-----------------------------print config-------------------------------------------

	PROGRAM_INFO
	/* summary of the parsed parameters */
	cout << "The specification is: \n";
	cout << "input protein filename: " << inname_protein << "\n";
	cout << "output protein filename: " << outname_protein << "\n";
	if (inname_ligand != "") {
		cout << "input ligand filename: " << inname_ligand << "\n";
		cout << "output ligand filename: " << outname_ligand << "\n";
		cout << "solvent-sphere radius:\t" << solventRadius << "Å, \n";
	}
	cout << "ligand-probe radius:\t" << ligandRadius << "Å, \n";
	cout << "probe-sphere radius:\t" << probeRadius << "Å, \n";
	cout << "minimum pocket depth:\t" << minPocketDepth << "Å, \n";
	cout << "number of top candidate pockets:\t" << topKPockets << "\n";
	cout << "number of top candidate ligand-excluded cavities:\t" << topKLECavities << "\n";
	cout << "number of top candidate probe-sphere-excluded cavities:\t" << topKPECavities << "\n";

	cout << "resolution:\t" << pow((double) resolution, 3.0) << " voxels per Å³, \n";
	if (evaluate_prediction)
		cout << "evaluate prediction: yes\n";
	else
		cout << "evaluate prediction: no\n";
	if (hetatm)
		cout << "include HETATM records: yes\n";
	else
		cout << "include HETATM records: no\n";
	if (hydrogen)
		cout << "include hydrogen atoms: yes\n";
	else
		cout << "include hydrogen atoms: no\n";
	cout << "atomic radii: " << inname_radii << "\n";
	cout << "**************************************************\n";

//-----------------------------computation--------------------------------------------
	try {
		Molecule * ligand = NULL;
		if (inname_ligand != "") {

			ligand = new Molecule(solventRadius, resolution, hydrogen, hetatm,
					inname_ligand, outname_ligand, inname_radii);
			float ligand_length, ligand_width, ligand_height;
			ligand->getMoleculeDimensions(ligand_length, ligand_width, ligand_height);
			ligandRadius = (Molecule::min_atm_radius >= 1.4) ? Molecule::min_atm_radius : 1.4;
			probeRadius = (ligand_height/2.0 >= 4.0) ? ligand_height/2.0 : 4.0;
			minPocketDepth = ligand_height/3.0;

			cout << "ligand length: " <<ligand_length<< "Å, \n";
			cout << "ligand width: " <<ligand_width<< "Å, \n";
			cout << "ligand height: " <<ligand_height<< "Å, \n";

			cout << "new ligand-probe radius:\t" << ligandRadius << "Å, \n";
			cout << "new probe-sphere radius:\t" << probeRadius << "Å, \n";
			cout << "new minimum pocket depth:\t" << minPocketDepth << "Å, \n";

		}

		Molecule protein(ligandRadius, probeRadius, resolution, minPocketDepth, topKPockets, topKLECavities, topKPECavities,
				hydrogen, hetatm, inname_protein, outname_protein, inname_radii);


		if (evaluate_prediction) {
			for (int cc = 0; cc < protein.num_pockets; ++cc) {
				point3D pcenter = protein.surface->pocket_list[cc].getCentroid();
				point3D wcenter = protein.surface->pocket_list[cc].getWeightedCentroid();
				float d = ligand->distanceToClosestAtom(pcenter);
				float w = ligand->distanceToClosestAtom(wcenter);
				float v = protein.surface->pocket_list[cc].getVolume();
				cout << "Pocket[" << cc << "] - centroid distance from ligand: " << d << "Å.\n";
				cout << "Pocket[" << cc << "] - weighted centroid distance from ligand: " << w << "Å.\n";
				cout << "Pocket[" << cc << "] - volume in voxels: " << v << ".\n";
				cout << "Pocket[" << cc << "] - centroid distance from protein centroid: " << pcenter.distance(protein.centroid) << "Å.\n";
				cout << "Pocket[" << cc << "] - weighted centroid distance from protein centroid: " << wcenter.distance(protein.centroid) << "Å.\n";


			}
			for (int cc = 0; cc < protein.num_le_cavities; ++cc) {
				point3D ccenter = protein.surface->le_cavity_list[cc].getCentroid();
				float d = ligand->distanceToClosestAtom(ccenter);
				float v = protein.surface->le_cavity_list[cc].getVolume();
				cout << "Cavity[" << cc << "] - centroid distance from ligand: " << d << "Å.\n";
				cout << "Cavity[" << cc << "] - volume in voxels: " << v << ".\n";
				cout << "Cavity[" << cc << "] - centroid distance from protein centroid: " << ccenter.distance(protein.centroid) << "Å.\n";
			}
		}

		if (ligand != NULL)
			delete ligand;
//-----------------------------conclusion --------------------------------------------
		cout << "**************************************************\n";
		cout << "Total calculation time:\t" << elapsedTime(t_start, std::chrono::high_resolution_clock::now()) << "\n";
		cout << "**************************************************\n";

	} catch (ParsingPDBException const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (fstream::failure const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (out_of_range const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (invalid_argument const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (logic_error const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (runtime_error const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
