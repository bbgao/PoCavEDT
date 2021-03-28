/*
 * CommandLineParser.cpp
 *
 *  Created on: Sep 6, 2014
 *      Author: sebastian
 */

#include "../utils/basename.h"
#include "../utils/disclaimer.h"
#include "CommandLineParser.h"
#include "CustomValidators.h"

/**
 * This method parses the input command line options.
 * @param argc		input argc reference, i.e. the first argument of the main method
 * @param argv		input argv reference, i.e. the second argument of the main method
 *
 * @return 			true if all arguments are parsed correctly, false otherwise.
 */
bool CommandLineParser::parseCommandLine(int const & argc, char const * const argv[],
		float & ligandRadius, float & probeRadius, float & solventRadius,
		float & resolution, string & inname_protein,
		string & outname_protein, string & inname_ligand,
		string & outname_ligand, string & inname_radii,
		float & minPocketDepth, int & topKPockets, int & topKLECavities,
		int & topKPECavities, bool & hydrogen, bool & hetatm,
		bool & evaluate_prediction, bool & help,
		bool & version, bool & license) {

	//------------------------------- command line options --------------------------------------
	boost::program_options::options_description description("PoCavEDT - Protein pocket/cavity detection with Euclidean Distance Transform. Usage");
	description.add_options()("help,h", "Display this help message.")
	/*
	 * This token is used to specify the input filename. It is mandatory.
	 * The pqr command line token has no option name. The command line tokens which have
	 * no option name are called "positional options" in the boost::program_options library.
	 */
	("protein", boost::program_options::value<input_PDB_filename>()->required(), "Input protein PDB file.")
	("ligand", boost::program_options::value<input_PDB_filename>(), "Input ligand PDB file (optional). If a candidate ligand is provided as input, the 'ligandRadius', 'probeRadius' and 'minPocketDepth' parameters will be automatically selected by the program (user input will be ignored).")
	/*
	 * This token is used to specify the output filename.
	 */
	("protein_output", boost::program_options::value<filename>(), "Protein output filename. If not specified, the input protein filename (with no extension) will be used.")
	("ligand_output", boost::program_options::value<filename>(), "Ligand output filename. If not specified, the input ligand filename (with no extension) will be used.")

	/*
	 * This token is used to specify the input atom radii filename.
	 */
	("atom_radii", boost::program_options::value<filename>(), "File containing the radius information of each atom. If not specified, the default CHARMM22 radius values will be used.")
	/*
	 * The -p (--probe_radius) flag is used to specify the probe radius. Default is 1.4. Optional.
	 */
	("ligandRadius,l", boost::program_options::value<probe_radius>()->default_value(1.4), "Ligand-sphere radius (in Å), positive floating point number (default is 1.4).")
	/*
	 * The -p (--probe_radius) flag is used to specify the probe radius. Default is 1.4. Optional.
	 */
	("probeRadius,R", boost::program_options::value<probe_radius>()->default_value(4.0), "Probe-sphere radius (in Å), positive floating point number (default is 4.0). Must be greater than the ligand-sphere radius.")
	("solventRadius,s", boost::program_options::value<probe_radius>()->default_value(1.4), "Solvent-sphere radius (in Å), positive floating point number (default is 1.4). Only used to compute the ligand's SES.")
	/*
	 * The -r (--resolution) flag is used to specify the resolution of the voxel grid. Default is 4. Optional.
	 */
	("resolution,r", boost::program_options::value<resolution_param>()->default_value(4.0), "Resolution factor, positive floating point (default is 4.0). This value's cube determines the number of voxels per Å³.")

	("minPocketDepth,d", boost::program_options::value<probe_radius>(), "Minimum depth from the Probe-excluded surface of a given pocket (in Å), positive floating point number (by default it is half the probe-sphere radius).")
	("topKPockets,p", boost::program_options::value<positive_int>()->default_value(3), "Number of top candidate pockets to consider, integer (default is 3). If set to -1, the program outputs all candidates.")
	("topKLECavities,c", boost::program_options::value<positive_int>()->default_value(3), "Number of top candidate ligand-excluded cavities to consider, integer (default is 3). If set to -1, the program outputs all candidates.")
	("topKPECavities,k", boost::program_options::value<positive_int>()->default_value(0), "Number of top candidate probe-excluded cavities to consider, integer (default is 0). If set to -1, the program outputs all candidates.")

	/*
	 * The --hetatm flag is used to include HETATM records in the surface computation.
	 */
	("hetatm", "Include HETATM records in the surface computation.")
	/*
	 * The --hetatm flag is used to include HETATM records in the surface computation.
	 */
	("hydrogen", "Include hydrogen atoms in the surface computation.\nNOTE: X-ray crystallography cannot resolve hydrogen atoms \
in most protein crystals, so in most PDB files, hydrogen atoms are absent. Sometimes hydrogens are added by modeling. \
Hydrogens are always present in PDB files resulting from NMR analysis, and are usually present in theoretical models.")

	("evaluate_prediction", "Can be used when a candidate ligand is provided. If selected, the program will return the distance from the centers of the predicted top k pockets to the clostest ligand atom.")
	/*
	 * The (--license) flag is used to view the program's license
	 * information.
	 */
	("license", "View license information.")
	/*
	 * The -v (--version) flag is used to view the program's version
	 * information.
	 */
	("version,v", "Display the version number");
	/*
	 * The input filename must be declared as positional.
	 */
	boost::program_options::positional_options_description p;
	p.add("protein", 1);
	p.add("ligand", 2);

	boost::program_options::variables_map vm;
	try {
		//--------------------------------parsing command line options------------------------------
		/*
		 * And it is finally specified when parsing the command line.
		 */
		store(boost::program_options::command_line_parser(argc, argv).options(description).positional(p).run(), vm);

		if (vm.count("help")) {
			cout << description;
			help = true;
			return true;
		}
		if (vm.count("version")) {
			cout << "Program: " << argv[0] << ", version: " << PROGRAM_VERSION << "\n";
			version = true;
			return true;
		}
		if (vm.count("license")) {
			DISCLAIMER
			license = true;
			return true;
		}
		if (vm.count("evaluate_prediction") && vm.count("ligand"))
			evaluate_prediction = true;
		else
			evaluate_prediction = false;
		if (vm.count("hetatm"))
			hetatm = true;
		else
			hetatm = false;
		if (vm.count("hydrogen"))
			hydrogen = true;
		else
			hydrogen = false;
		/*
		 * notify throws exceptions so we call it after the above checks
		 */
		notify(vm);

		/* initializing variables */
		inname_protein = vm["protein"].as<input_PDB_filename>().filename;
		if (vm.count("protein_output")) {
			outname_protein = vm["protein_output"].as<filename>().fname;
		} else {
			int lastindex = inname_protein.find_last_of(".");
			outname_protein = basename(inname_protein.substr(0, lastindex));
		}

		if (vm.count("ligand")) {
			inname_ligand = vm["ligand"].as<input_PDB_filename>().filename;
			if (vm.count("ligand_output")) {
				outname_ligand = vm["ligand_output"].as<filename>().fname;
			} else {
				int lastindex = inname_ligand.find_last_of(".");
				outname_ligand = basename(inname_ligand.substr(0, lastindex));
			}
		} else {
			inname_ligand = "";
			outname_ligand  ="";
		}

		if (vm.count("atom_radii")) {
			inname_radii = vm["atom_radii"].as<filename>().fname;
		} else {
			inname_radii = "CHARMM22";
		}
		ligandRadius = vm["ligandRadius"].as<probe_radius>().p;
		probeRadius = vm["probeRadius"].as<probe_radius>().p;
		solventRadius = vm["solventRadius"].as<probe_radius>().p;

		topKPockets = vm["topKPockets"].as<positive_int>().k;
		topKLECavities = vm["topKLECavities"].as<positive_int>().k;
		topKPECavities = vm["topKPECavities"].as<positive_int>().k;

		if (vm.count("minPocketDepth")) {
			minPocketDepth = vm["minPocketDepth"].as<probe_radius>().p;
		} else {
			minPocketDepth = 0.5 * probeRadius;
		}

		if (probeRadius <= ligandRadius)
			throw boost::program_options::error("The ligand-probe radius must be greater than the solvent-probe radius!");

		resolution = vm["resolution"].as<resolution_param>().r;

	} catch (boost::program_options::error const & e) {
		cerr << "error: " << e.what() << "\n";
		cerr << description << "\n";
		return false;
	}
	return true;
};


