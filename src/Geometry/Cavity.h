/*
 * Cavity.h
 *
 *  Created on: Jun 24, 2018
 *      Author: sebastian
 */

#include "voxel.h"
#include "voxelGrid.h"
#include <vector>
#include "point3D.h"

#ifndef CAVITY_H_
#define CAVITY_H_


class Cavity {

public:
	Cavity(voxelGrid const & cavity_grid, float resolution, point3D const & ptran);
	Cavity(Cavity const & p);
	double getVolume() const;
	Cavity & operator=(Cavity const & p);
	void outputCavityPCDModel(std::string const & filename, float resolution, point3D const & ptran) const;
	point3D getCentroid() const;

private:
	double volume;
	point3D centroid;
	std::vector<voxel> cavityVoxels;
};


#endif /* CAVITY_H_ */
