/*
 * EuclideanDistanceTransform.h
 *
 *  Created on: Jan 19, 2018
 *      Author: sebastian
 */

#ifndef EDT_EUCLIDEANDISTANCETRANSFORM_H_
#define EDT_EUCLIDEANDISTANCETRANSFORM_H_


#include "../voxel.h"
#include "../voxelGrid.h"
#include "../HierarchicalQueue.h"


#include "math.h"


void computeEDT(voxelGrid const * sa_surface, voxelGrid const * cpk_model,
		voxel const & cc, uint16_t clength, uint16_t min_sdist,
		uint16_t numberOfQueues, uint16_t * distanceMap) {
	voxel * nearestSurfaceVoxel = new voxel[clength * clength * clength];

	HierarchicalQueue* HQ1 = new HierarchicalQueue(numberOfQueues);
	HierarchicalQueue* HQ2 = new HierarchicalQueue(numberOfQueues);

	/* for all voxels */

	int d = clength / 2;


	for (int i = -d; i <= d; ++i) {
		for (int j = -d; j <= d; ++j) {
			for (int k = -d; k <= d; ++k) {
				/* create the first shell by pushing the surface voxels in the list */
				if (sa_surface->getVoxel(i, j, k)) {
					HQ1->push(voxel(i, j, k), 0);
					/* a surface voxel has zero distance from itself */
					nearestSurfaceVoxel[(i * clength + j) * clength + k] = voxel(i, j, k);
					distanceMap[(i * clength + j) * clength + k] = 0;
				} else {
					distanceMap[(i * clength + j) * clength + k] = UINT16_MAX;
				}
			}
		}
	}

	delete sa_surface;
	sa_surface = NULL;

	while (!HQ1->empty()) {
		/*current voxel*/
		voxel cv = HQ1->front();
		HQ1->pop();
		voxel nearestSurfVox = nearestSurfaceVoxel[(cv.ix * clength + cv.iy) * clength + cv.iz];
		uint16_t squaredDistance = distanceMap[(cv.ix * clength + cv.iy) * clength + cv.iz];
		voxel nb; // neighbour
		bool isEnd = true;
		for (int i = 0; i < 26; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (cpk_model->getVoxel(nb)) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				if (newSDistance < distanceMap[(nb.ix * clength + nb.iy) * clength + nb.iz]
						&& newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * clength + nb.iy) * clength + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(nb.ix * clength + nb.iy) * clength + nb.iz] = nearestSurfVox;
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
		voxel nearestSurfVox = nearestSurfaceVoxel[(cv.ix * clength + cv.iy) * clength + cv.iz];
		// neighbor
		voxel nb;
		for (int i = 0; i < 124; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (se_cpk_model->getVoxel(nb)) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				if (newSDistance < distanceMap[(nb.ix * clength + nb.iy) * clength + nb.iz] && newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * clength + nb.iy) * clength + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(nb.ix * clength + nb.iy) * clength + nb.iz] = nearestSurfVox;
					HQ2->push(nb, newSDistance);
				}
			}
		}
	}
	cout << "Calculating Solvent Excluded CPK Model." << endl;
	/* update the volumetric model */
	for (int i = 0; i < clength; ++i) {
		for (int j = 0; j < clength; ++j) {
			for (int k = 0; k < clength; ++k) {
				if (distanceMap[(i * clength + j) * clength + k] < min_sdist) {
					se_cpk_model->clearVoxel(i, j, k);
				}
			}
		}
	}

	delete[] nearestSurfaceVoxel;
	nearestSurfaceVoxel = NULL;

	delete HQ2;
}




#endif /* EDT_EUCLIDEANDISTANCETRANSFORM_H_ */
