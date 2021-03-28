/*
 * CustomValidators.h
 *
 *  Created on: Jan 11, 2014
 *      Author: sebastian
 */
/**
 * Here we overload the validate function in order to accept only valid
 * command line input parameters.
 *
 * For a detailed description on overloading the default validator see:
 * http://www.boost.org/doc/libs/1_55_0/doc/html/program_options/howto.html#idp163429032
 */
#ifndef CUSTOMVALIDATORS_H_
#define CUSTOMVALIDATORS_H_

#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <string>

using namespace std;

struct input_PDB_filename {
	input_PDB_filename(string const & val) :
			filename(val) {	}
	string filename;
	friend ostream& operator <<(ostream& s, const input_PDB_filename& idpbf) {
		s << idpbf.filename;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * PDB input filenames.
 *
 * The function takes four parameters. The first is the storage for
 * the value, and in this case is either empty or contains an instance
 * of the magic_number class. The second is the list of strings found
 * in the next occurrence of the option. The remaining two parameters
 * are needed to workaround the lack of partial template specialization
 * and partial function template ordering on some compilers.
 *
 * The function first checks that we don't try to assign to the same
 * option twice. Then it checks that only a single string was passed
 * in. Next the string is verified. If that test is passed, the parsed
 * value is stored into the v variable.
 */

void validate(boost::any& v, vector<string> const& values,
		input_PDB_filename* /* target_type */, int) {
	static boost::regex r("[^[.NUL.]]+(?:\\.pdb|\\.PDB)");
	// Make sure no previous assignment to 'v' was made.
	boost::program_options::validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = boost::program_options::validators::get_single_string(values);
	// Do regex match
	if (boost::regex_match(s, r))
		v = boost::any(input_PDB_filename(s));
	 else
		throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
};
struct filename {
	filename(string const & outname) :
		fname(outname) {	}
	string fname;
	friend ostream& operator <<(ostream& s, const filename& out) {
		s << out.fname;
		return s;
	}
};
void validate(boost::any& v, vector<string> const& values,
		filename* /* target_type */, int) {
	static boost::regex r("[^[.NUL.]]+");
	// Make sure no previous assignment to 'v' was made.
	boost::program_options::validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = boost::program_options::validators::get_single_string(values);
	// Do regex match
	if (regex_match(s, r))
		v = boost::any(filename(s));
	 else
		throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
};

struct probe_radius {
public:
	probe_radius(float p) :
			p(p) { }
	float p;
	friend ostream& operator <<(ostream& s, const probe_radius& pr) {
		s << pr.p;
		return s;
	}
};
/**
 * Here we overload the validate function in order to accept only valid
 * probe radius parameters.
 */

void validate(boost::any& v, vector<string> const& values,
		probe_radius* /* target_type */, int) {
	static boost::regex r("[+-]?[0-9]*\\.?[0-9]*");

	// Make sure no previous assignment to 'v' was made.
	boost::program_options::validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = boost::program_options::validators::get_single_string(values);

	if (regex_match(s, r)) {
		float temp = boost::lexical_cast<float>(s);
		if (temp > 0)
			v = boost::any(probe_radius(temp));
		else
			throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
	}
	else
		throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
};

struct resolution_param {
public:
	resolution_param(float r) :
			r(r) { }
	float r;
	friend ostream& operator <<(ostream& s, const resolution_param& rp) {
		s << rp.r;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * resolution values.
 */

void validate(boost::any& v, vector<string> const& values,
		resolution_param* /* target_type */, int) {
	static boost::regex r("[+-]?[0-9]*\\.?[0-9]*");

	// Make sure no previous assignment to 'v' was made.
	boost::program_options::validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = boost::program_options::validators::get_single_string(values);

	if (regex_match(s, r)) {
		float temp = boost::lexical_cast<float>(s);
		if (temp > 0)
			v = boost::any(resolution_param(temp));
		else
			throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
	}
	else
		throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
};

struct out_type {
public:
	out_type(int st) :
			st(st) { }
	int st;
	friend ostream& operator <<(ostream& s, const out_type& type) {
		s << type.st;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * resolution values.
 */

void validate(boost::any& v, vector<string> const& values,
		out_type* /* target_type */, int) {
	static boost::regex r("[1234]");

	// Make sure no previous assignment to 'v' was made.
	boost::program_options::validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = boost::program_options::validators::get_single_string(values);

	if (regex_match(s, r)) {
		v = boost::any(out_type(boost::lexical_cast<int>(s)));
	}
	else
		throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
};


struct positive_int {
public:
	positive_int(int k) :
			k(k) { }
	int k;
	friend ostream& operator <<(ostream& s, const positive_int& type) {
		s << type.k;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * resolution values.
 */

void validate(boost::any& v, vector<string> const& values,
		positive_int* /* target_type */, int) {
	static boost::regex r("[0-9]+|-1");

	// Make sure no previous assignment to 'v' was made.
	boost::program_options::validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = boost::program_options::validators::get_single_string(values);

	if (regex_match(s, r)) {
		v = boost::any(positive_int(boost::lexical_cast<int>(s)));
	}
	else
		throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
};

#endif /* CUSTOMVALIDATORS_H_ */
