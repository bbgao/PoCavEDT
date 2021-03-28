/*
 * EDT.h
 *
 *  Created on: 20 gen 2018
 *      Author: dsl
 */

#ifndef EDT_EDT_H_
#define EDT_EDT_H_


#include "../voxel.h"
#include "../voxelGrid.h"
#include "../voxel_offset.h"
#include "../HierarchicalQueue.h"

inline bool inCube(voxel const & cc, uint16_t d, voxel const & v) {
	return (abs(cc.ix - v.ix) <= d && abs(cc.iy - v.iy) <= d && abs(cc.iz - v.iz) <= d);
}


inline bool inCube(int cx, int cy, int cz, uint16_t d, voxel const & v) {
	return (abs(cx - v.ix) <= d && abs(cy - v.iy) <= d && abs(cz - v.iz) <= d);
}

inline void computeEDT(voxelGrid const * cpk_model, voxelGrid const * sa_surface, uint16_t numberOfQueues,
		uint16_t min_sdist, voxel const & cc, uint16_t d, uint16_t pwidth, uint16_t pheight,
		uint16_t * distanceMap, voxel * nearestSurfaceVoxel) {

	HierarchicalQueue* HQ1 = new HierarchicalQueue(numberOfQueues);
	HierarchicalQueue* HQ2 = new HierarchicalQueue(numberOfQueues);

	/* for all voxels */
	for (int i = cc.ix - d; i <= cc.ix + d; ++i) {
		for (int j = cc.iy - d; j <= cc.iy + d; ++j) {
			for (int k = cc.iz - d; k <= cc.iz + d; ++k) {
				/* create the first shell by pushing the surface voxels in the list */
				if (sa_surface->getVoxel(i, j, k)) {
					HQ1->push(voxel(i, j, k), 0);
					/* a surface voxel has zero distance from itself */
					nearestSurfaceVoxel[(i * pwidth + j) * pheight + k] = voxel(i, j, k);
					distanceMap[(i * pwidth + j) * pheight + k] = 0;
				}
			}
		}
	}

	while (!HQ1->empty()) {
		/*current voxel*/
		voxel cv = HQ1->front();
		HQ1->pop();
		voxel nearestSurfVox = nearestSurfaceVoxel[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		uint16_t squaredDistance = distanceMap[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		voxel nb; // neighbour
		bool isEnd = true;
		for (int i = 0; i < 26; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (cpk_model->getVoxel(nb) && inCube(cc, d, nb)) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				if (newSDistance < distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz]
						&& newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = nearestSurfVox;
					HQ1->push(nb, newSDistance);
					isEnd = false;
				}
			}
		}
		if (isEnd && squaredDistance >= 24) {
			HQ2->push(cv, squaredDistance);
		}
	}

	delete HQ1;
	HQ1 = NULL;

	while (!HQ2->empty()) {
		/*current voxel*/
		voxel cv = HQ2->front();
		HQ2->pop();
		voxel nearestSurfVox = nearestSurfaceVoxel[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		// neighbor
		voxel nb;
		for (int i = 0; i < 124; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (cpk_model->getVoxel(nb) && inCube(cc, d, nb)) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				if (newSDistance < distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] && newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = nearestSurfVox;
					HQ2->push(nb, newSDistance);
				}
			}
		}
	}
	delete HQ2;
	HQ2 = NULL;
}

inline void computeEDT(voxelGrid const * cpk_model,
		voxelGrid const * sa_surface, uint16_t numberOfQueues,
		uint16_t min_sdist, atom const & atm, point3D const & ptran,
		float resolution, float probeRadius, uint16_t pwidth, uint16_t pheight,
		uint16_t * distanceMap, voxel * nearestSurfaceVoxel) {

	HierarchicalQueue* HQ1 = new HierarchicalQueue(numberOfQueues);
	HierarchicalQueue* HQ2 = new HierarchicalQueue(numberOfQueues);

	int cx = static_cast<int>(round((atm.x + ptran.x) * resolution));
	int cy = static_cast<int>(round((atm.y + ptran.y) * resolution));
	int cz = static_cast<int>(round((atm.z + ptran.z) * resolution));
	uint16_t d = static_cast<int>(round((atm.radius + probeRadius) * resolution));

	/* for all voxels */
	for (int i = cx - d; i <= cx + d; ++i) {
		for (int j = cy - d; j <= cy + d; ++j) {
			for (int k = cz - d; k <= cz + d; ++k) {
				/* create the first shell by pushing the surface voxels in the list */
				if (sa_surface->getVoxel(i, j, k)) {
					HQ1->push(voxel(i, j, k), 0);
					/* a surface voxel has zero distance from itself */
					nearestSurfaceVoxel[(i * pwidth + j) * pheight + k] = voxel(i, j, k);
					distanceMap[(i * pwidth + j) * pheight + k] = 0;
				}
			}
		}
	}

	while (!HQ1->empty()) {
		/*current voxel*/
		voxel cv = HQ1->front();
		HQ1->pop();
		voxel nearestSurfVox = nearestSurfaceVoxel[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		uint16_t squaredDistance = distanceMap[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		voxel nb; // neighbour
		bool isEnd = true;
		for (int i = 0; i < 26; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (cpk_model->getVoxel(nb) && inCube(cx, cy, cz, d, nb)) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				if (newSDistance < distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz]
						&& newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = nearestSurfVox;
					HQ1->push(nb, newSDistance);
					isEnd = false;
				}
			}
		}
		if (isEnd && squaredDistance >= 24) {
			HQ2->push(cv, squaredDistance);
		}
	}

	delete HQ1;
	HQ1 = NULL;

	while (!HQ2->empty()) {
		/*current voxel*/
		voxel cv = HQ2->front();
		HQ2->pop();
		voxel nearestSurfVox = nearestSurfaceVoxel[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		// neighbor
		voxel nb;
		for (int i = 0; i < 124; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (cpk_model->getVoxel(nb) && inCube(cx, cy, cz, d, nb)) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				if (newSDistance < distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] && newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = nearestSurfVox;
					HQ2->push(nb, newSDistance);
				}
			}
		}
	}
	delete HQ2;
	HQ2 = NULL;
}

#endif /* EDT_EDT_H_ */
