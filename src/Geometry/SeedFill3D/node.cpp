/*
 * node.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: sebastian
 */


#include "node.h"

node::node() :
		xl(0), xr(0), y(0), direction_y(1) { }

node::node(node const & n) :
		xl(n.xl), xr(n.xr), y(n.y), direction_y(n.direction_y) { }

node::node(int xl, int xr, int y, int direction_y) :
		xl(xl), xr(xr), y(y), direction_y(direction_y) { }

/**
 * Copy assignment operator
 */
node & node::operator=(node const & n) {
	if (this != &n) { // protect against invalid self-assignment
		this->xl = n.xl;
		this->xr = n.xr;
		this->y = n.y;
		this->direction_y = n.direction_y;
	}
	// by convention, always return *this
	return *this;
}
