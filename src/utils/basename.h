
#ifndef UTILS_BASENAME_H_
#define UTILS_BASENAME_H_

#include <algorithm>
#include <string>

using namespace std;

struct MatchPathSeparator {
	bool operator()(char ch) const {
		return ch == '/';
	}
};

string basename(string const& pathname) {
	return string(find_if(pathname.rbegin(), pathname.rend(), MatchPathSeparator()).base(), pathname.end());
}

#endif /* UTILS_TRIM_H_ */
