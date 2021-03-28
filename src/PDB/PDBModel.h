/**
 * In this header we define the PDBModel class and relative methods.
 */

#ifndef PDBMODEL_H_
#define PDBMODEL_H_

#include "../Geometry/point3D.h"
#include "../PDB/atom.h"
#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <vector>

typedef enum lineType {
	ATOM, HETATM, OTHER_LN, TER, REMARK, END
} lineType;

namespace std {
template<>
struct hash<pair<string,string> > {
	size_t operator()(pair<string,string> const & p) const {
		size_t result = hash<string>()(p.first);
		boost::hash_combine(result, p.second);
		return result;
	}
};

template<>
struct equal_to<pair<string,string> > {
	bool operator()(pair<string,string> const & p, pair<string,string> const & q) const {
		return (p.first == q.first) && (p.second == q.second);
	}
};
}

class PDBModel {
public:
	size_t line_num;
	std::vector<atom> atomsInModel;
	std::string filename;
	std::string header;
	bool hydrogen;
	bool hetatm;

	static std::unordered_map<pair<string,string>, float> atom_radii;
	static std::unordered_map<string, float> single_atom_radii;

	PDBModel(string const & filename, string const & inname_radii, bool hydrogen, bool hetatm);
	PDBModel(PDBModel const & model);
	/**
	 * Copy assignment operator
	 */
	PDBModel & operator=(PDBModel const & model);

private:
	lineType getLineType(string const & line);
	void loadFile(string const & filename);
	void loadAtomRadii(string const & inname_radii);
	void parseLine(string const & line);
};
#endif /* PDBMODEL_H_ */
