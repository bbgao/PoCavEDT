/*
 * CommandLineParser.h
 *
 *  Created on: Sep 6, 2014
 *      Author: sebastian
 */

#ifndef COMMANDLINEPARSER_H_
#define COMMANDLINEPARSER_H_

#include <string>

using namespace std;

/**
 * This class provides a simple interface for accessing the command line arguments.
 */
class CommandLineParser {
public:
	static bool parseCommandLine(int const & argc, char const * const argv[],
			float & ligandRadius, float & probeRadius, float & solventRadius,
			float & resolution, string & inname_protein,
			string & outname_protein, string & inname_ligand,
			string & outname_ligand, string & inname_radii,
			float & minPocketDepth, int & topKPockets, int & topKLECavities,
			int & topKPECavities, bool & hydrogen, bool & hetatm,
			bool & evaluate_prediction, bool & help,
			bool & version, bool & license);
};

#endif /* COMMANDLINEPARSER_H_ */
