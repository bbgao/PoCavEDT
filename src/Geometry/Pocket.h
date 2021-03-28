/*
 * Pocket.h
 *
 *  Created on: Jun 23, 2018
 *      Author: sebastian
 */

#include "voxel.h"
#include "voxelGrid.h"
#include <vector>
#include "point3D.h"
#ifndef POCKET_H_
#define POCKET_H_


class Pocket {

public:
	Pocket(voxelGrid const & pocket_grid, uint16_t const * distanceMap, float resolution, point3D const & ptran);
	Pocket(Pocket const & p);
	Pocket & operator=(Pocket const & p);
	double getScore() const;
	double getVolume() const;
	point3D getCentroid() const;
	point3D getWeightedCentroid() const;

	void outputPocketPCDModel(std::string const & filename, float resolution, point3D const & ptran) const;

private:
	double score, volume;
	point3D centroid, weighted_centroid;
	std::vector<voxel> pocketVoxels;
};


#endif /* POCKET_H_ */
