

#include "Pocket.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include "../utils/disclaimer.h"
#include "../utils/makeDirectory.h"

using namespace std;

Pocket::Pocket(voxelGrid const & pocket_grid, uint16_t const * distanceMap, float resolution, point3D const & ptran) {
	score = 0;
	volume = 0;
	centroid = point3D(0, 0, 0);
	weighted_centroid = point3D(0, 0, 0);

	for (int ii = 0; ii < pocket_grid.length; ++ii) {
		for (int jj = 0; jj < pocket_grid.width; ++jj) {
			for (int kk = 0; kk < pocket_grid.height; ++kk) {
				if (pocket_grid.getVoxel(ii, jj, kk)) {
					pocketVoxels.push_back(voxel(ii, jj, kk));
//					float d = sqrt(distanceMap[(ii * pocket_grid.width + jj) * pocket_grid.height + kk]);
					float d = distanceMap[(ii * pocket_grid.width + jj) * pocket_grid.height + kk];
					score += d;
					point3D cp(ii / resolution - ptran.x, jj / resolution - ptran.y, kk / resolution - ptran.z);
					++volume;
					centroid += cp;
					weighted_centroid += d * cp;
				}
			}
		}
	}

	centroid /= volume;
	weighted_centroid /= score;
}
Pocket::Pocket(Pocket const & p) :
		score(p.score), volume(p.volume), pocketVoxels(p.pocketVoxels),
		centroid(p.centroid), weighted_centroid(p.weighted_centroid) { }

Pocket & Pocket::operator=(Pocket const & p) {
	this->score = p.score;
	this->volume = p.volume;
	this->pocketVoxels = p.pocketVoxels;
	this->centroid = p.centroid;
	this->weighted_centroid = p.weighted_centroid;
	return *this;
}
double Pocket::getScore() const{
	return score;
}
double Pocket::getVolume() const{
	return volume;
}
void Pocket::outputPocketPCDModel(string const & filename, float resolution, point3D const & ptran) const {
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
			<< pocketVoxels.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << pocketVoxels.size()
			<< "\n" << "DATA ascii";

	for(auto const & v : pocketVoxels) {
		file_stream << "\n" << v.ix / resolution - ptran.x
				<< " " << v.iy / resolution - ptran.y
				<< " " << v.iz / resolution - ptran.z;
	}
	file_stream.close();
}

point3D Pocket::getCentroid() const { return centroid;}
point3D Pocket::getWeightedCentroid() const { return weighted_centroid; }

