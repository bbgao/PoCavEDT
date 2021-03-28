/*
 * range.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: sebastian
 */


#include "range.h"


range::range() : xl(0), xr(0), y(0) {}

range::range(range const & r) : xl(r.xl), xr(r.xr), y(r.y) {}

range::range(int l, int r, int y) : xl(l), xr(r), y(y) {}

range & range::operator=(range const & r) {
	if (this != &r) { // protect against invalid self-assignment
		this->xl = r.xl;
		this->xr = r.xr;
		this->y = r.y;
	}
	// by convention, always return *this
	return *this;
}
