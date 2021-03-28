/*
 * Cavity.cpp
 *
 *  Created on: Jun 24, 2018
 *      Author: sebastian
 */



#include "Cavity.h"
#include <iomanip>
#include <iostream>
#include "../utils/disclaimer.h"
#include "../utils/makeDirectory.h"

Cavity::Cavity(voxelGrid const & cavity_grid, float resolution, point3D const & ptran) {
	volume = 0;
	centroid = point3D(0, 0, 0);
	size_t n_voxels = 0;
	for (int ii = 0; ii < cavity_grid.length; ++ii) {
		for (int jj = 0; jj < cavity_grid.width; ++jj) {
			for (int kk = 0; kk < cavity_grid.height; ++kk) {
				if (cavity_grid.getVoxel(ii, jj, kk)) {
					cavityVoxels.push_back(voxel(ii, jj, kk));
					point3D cp(ii / resolution - ptran.x, jj / resolution - ptran.y, kk / resolution - ptran.z);
					centroid += cp;
					volume += 1;
				}
			}
		}
	}
	centroid /= volume;
}

Cavity::Cavity(Cavity const & c) : volume(c.volume), centroid(c.centroid), cavityVoxels(c.cavityVoxels) {}

double Cavity::getVolume() const {
	return volume;
}

point3D Cavity::getCentroid() const { return centroid;}

Cavity & Cavity::operator=(Cavity const & c) {
	this->volume = c.volume;
	this->centroid = c.centroid;
	this->cavityVoxels = c.cavityVoxels;
	return *this;
}
void Cavity::outputCavityPCDModel(std::string const & filename, float resolution, point3D const & ptran) const{
	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open("./output/" + filename + ".pcd");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	file_stream  << std::setprecision(5);

	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "# created with " << PROGRAM_NAME << "\n"
			<< "# version " << PROGRAM_VERSION << "\n"
			<< "VERSION .7\n" << "FIELDS x y z\n" << "SIZE 4 4 4\n"
			<< "TYPE F F F\n" << "COUNT 1 1 1\n" << "WIDTH "
			<< cavityVoxels.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << cavityVoxels.size()
			<< "\n" << "DATA ascii";

	for(auto const & v : cavityVoxels) {
		file_stream << "\n" << v.ix / resolution - ptran.x
				<< " " << v.iy / resolution - ptran.y
				<< " " << v.iz / resolution - ptran.z;
	}
	file_stream.close();
}
