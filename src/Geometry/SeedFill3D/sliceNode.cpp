/*
 * sliceNode.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: sebastian
 */

#include "sliceNode.h"

sliceNode::sliceNode() : z(0), direction_z(0) { };
sliceNode::sliceNode(sliceNode const & sn) :
		z(sn.z), seedList(sn.seedList), direction_z(sn.direction_z) { };
sliceNode::sliceNode(IdList const & seedList, int z, int direction_z) :
		z(z), seedList(seedList), direction_z(direction_z) { };
/**
 * Copy assignment operator
 */
sliceNode & sliceNode::operator=(sliceNode const & sn) {
	if (this != &sn) { // protect against invalid self-assignment
		this->z = sn.z;
		this->seedList = sn.seedList;
		this->direction_z = sn.direction_z;
	}
	// by convention, always return *this
	return *this;
}

