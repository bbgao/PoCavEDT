/*
 * rapid3DSurfaceExtract.h
 *
 *  Created on: Feb 5, 2015
 *      Author: sebastian
 */
/**
 * Here we implement a method that extracts the molecular surface from
 * a volumetric model, using the 3D seed filling-algorithm proposed in:
 * Wei-Wei Yu, Fei He, Ping Xi, A rapid 3D seed-filling algorithm based on scan slice,
 * Computers & Graphics, Volume 34, Issue 4, August 2010, Pages 449-459, ISSN 0097-8493,
 * http://dx.doi.org/10.1016/j.cag.2010.05.005.
 * http://www.sciencedirect.com/science/article/pii/S0097849310000725
 */
#ifndef SEEDFILL3D_RAPID3DSURFACEEXTRACT_H_
#define SEEDFILL3D_RAPID3DSURFACEEXTRACT_H_

#include "../voxelGrid.h"
#include "node.h"
#include "sliceNode.h"
#include <stack>

namespace rapid3DSurfaceExtract {

/**
 * Clears a whole range in temp, from xl to xr, on line y of slice z.
 * Also checks if any of the voxels in the range belongs to the surface.
 * @param image		target image
 * @param temp		current 3D image
 * @param output	output voxelGrid
 * @param xl		range lower bound
 * @param xr		range upper bound
 * @param y			range line
 * @param z			image slice
 */
static inline void clearRange(voxelGrid const & cpkModel, voxelGrid & temp, voxelGrid & output,
		int const & xl, int const & xr, int const & y, int const & z) {
	for (int i = xl; i <= xr; i++) {
		voxel nb, v(i, y, z);
		for (int ii = 0; ii < 4; ++ii) {
			nb = v + voxel_offset::neighbours[ii];
			if (!cpkModel.getVoxel(nb))
				output.setVoxel(v);
		}
		temp.clearVoxel(v);
	}
	if (xl > 0)
		output.setVoxel(xl, y, z);
	if (xr < output.length - 1)
		output.setVoxel(xr, y, z);
};

/**
 * Returns the filled and valid range discovered scanning the target image from seed
 * @param temp		current 3D image
 * @param seed		the seed
 * OUTPUT:
 * @param xl 		the left range delimiter
 * @param xr 		the right range delimiter
 */
static inline void obtainFilledValidRange(voxelGrid const & temp, voxel const & seed,
		int & xl, int & xr) {
	xr = seed.ix;
	while ((xr < temp.length) && temp.getVoxel(xr, seed.iy, seed.iz)) {
		++xr;
	}
	--xr;
	xl = seed.ix;
	while ((xl > -1) && temp.getVoxel(xl, seed.iy, seed.iz)) {
		--xl;
	}
	++xl;
};

/**
 * Checks if the current range [xpl, xpr] contains any valid seeds. If a seed is found
 * the corresponding range is extracted, filled and inserted in dList.
 * @param imageData	the current voxelGrid
 * @param temp		temporary voxelGrid
 * @param output	output voxelGrid
 * @param xpl		lower bound of the range to be checked
 * @param xpr		upper bound of the range to be checked
 * @param y			line of the current range
 * @param z			slice of the current range
 * OUTPUT:
 * @param dList		list of the filled ranges in the current slice
 * @param xl		lower bound of the newly found range
 * @param xr		upper bound of the newly found range
 * @return			true if a valid seed is found, false otherwise
 */
static inline bool checkAndClearRange(voxelGrid const & cpkModel, voxelGrid & temp, voxelGrid & output, int xpl, int xpr, int y, int z,
		IdList & dList, int & xl, int & xr) {
	for (int i = xpl; i <= xpr; ++i) {
		if (temp.getVoxel(i, y, z)) { // if a valid seed is found
			// extract the unfilled range,and write the filled seeds ID into dList
			// search the valid and unfilled range ([xl,xr]) from [xpl, xpr]
			obtainFilledValidRange(temp, voxel(i, y, z), xl, xr);
			//extract and fill the range [xl,xr], mark the seeds and write the ID of them into dList
			clearRange(cpkModel, temp, output, xl, xr, y, z);
			dList.push_back(range(xl, xr, y));
			return true;
		}
	}
	return false;
};

/**
 * This method flood-fills a single slice of the input voxelGrid, starting from voxel seed.
 * The slice is identified by the z coordinate of seed. The ID of each seed is extracted from
 * dList.
 * @param imageData	target voxelGrid
 * @param temp		temporary voxelGrid
 * @param output	output voxelGrid
 * @param seed		starting voxel for the flood-filling procedure
 * OUTPUT:
 * @param seedList	list of the filled ranges in the current slice
 */
static inline void surfaceExtract2D(voxelGrid const & cpkModel, voxelGrid & temp, voxelGrid & output, voxel const & seed, IdList & seedList) {
	// empty the seed list
	seedList.clear();
	// initialize the empty 2D stack
	stack<node> stack2D;
	// obtain the unfilled and valid range ([xl,xr]) by using seed
	int xl, xr; // range [xl, xr]
	obtainFilledValidRange(temp, seed, xl, xr);
	// extract and fill the range [xl,xr],
	clearRange(cpkModel, temp, output,  xl, xr, seed.iy, seed.iz);
	// mark the seeds and write the ID of them into dList
	seedList.push_back(range(xl, xr, seed.iy));
	// create two nodes with range ([xl, xr]), push them into 2D stack, along two opposite directions
	if (seed.iy > 0)
		stack2D.push(node(xl, xr, seed.iy, -1));
	if (seed.iy < temp.width - 1)
		stack2D.push(node(xl, xr, seed.iy, 1));
	// obtain the length of this range
	while (!stack2D.empty()) { // if the 2D stack is not empty
		//pop a node and then set the search range ([xpl,xpr]) onto the next scan line
		node c(stack2D.top()); // current node
		stack2D.pop();
		int xpl = c.xl, xpr = c.xr, y = c.y + c.direction_y;

		while (checkAndClearRange(cpkModel, temp, output, xpl, xpr, y, seed.iz, seedList, xl, xr)) { //if the range [xpl, xpr] is valid
			if (xr < xpr - 1) { //it maybe exist other valid and unfilled ranges on the same scan line
				for (int i = xr + 2; i <= xpr; ++i) { // search the leftmost seed of the next valid range, which locates at the same scan line
					if (temp.getVoxel(i, y, seed.iz)) { /* a valid seed exists */
						stack2D.push(node(xr + 2, xpr, y - c.direction_y, c.direction_y));
						break;
					}
				}
			}
			//execute the necessary rollback operation
			if (xl < xpl - 1) {	//the rollback operation may occur on the left side of xpl
			// rollback to the previous scan line
				int new_y = y - c.direction_y;
				// search a valid and unfilled seed on the left side of xpl
				for (int i = xl; i <= xpl - 2; ++i) {
					if (temp.getVoxel(i, new_y, seed.iz)) { /* a valid seed exists */
						//the search range is valid, push it into 2D stack
						// push a new node into 2D stack, the range is [xl, xpl-2], the direction is opposite
						stack2D.push(node(xl, xpl - 2, y, -c.direction_y));
						break;
					}
				}
			}
			if (xr > xpr + 1) {	//the rollback operation may occur on the right side of xpr
			// rollback to the previous scan line
				int new_y = y - c.direction_y;
				// search a valid and unfilled seed on the right side of xpr
				for (int i = xpr + 2; i <= xr; ++i) {
					// if (a valid seed exists){//the search range is valid, push it into 2D stack
					if (temp.getVoxel(i, new_y, seed.iz)) { /* a valid seed exists */
						// the search range is valid, push it into 2D stack
						// push a new node into 2D stack, the range is [xpr+ 2, xr], the direction is opposite
						stack2D.push(node(xpr + 2, xr, y, -c.direction_y));
						break;
					}
				}
			}
			// continue the loop from this range [xl, xr], along the same direction as before
			// obtain the new search range by assigning the range [xl, xr] to [xpl, xpr]
			xpl = xl;
			xpr = xr;
			y += c.direction_y;
			if (y < 0 || y >= temp.width)
				break;
		}
	}
};

/**
 * This method projects the seed list of the current slice onto slice z and returns true if
 * a valid seed is found in the projected area, false otherwise.
 * @param temp		the current voxelGrid
 * @param seed_list	the list of seeds in the current slice
 * @param z			the neighbor slice
 * @param leap_var	the leaping variable is used to leap over plenty of invalid
 * 					voxels in order to accelerate search process for the first
 * 					valid seed
 * OUTPUT:
 * @param sx 		x coordinate of the first valid seed found
 * @param sy		y coordinate of the first valid seed found
 * @return			true if a new seed is found, false otherwise
 */
static inline bool searchNeighborSlice(voxelGrid const & temp, IdList & seed_list,
		int z, int leap_var, int & sx, int & sy) {
	IdList::iterator id, end;
	for (id = seed_list.begin(), end = seed_list.end(); id != end; ++id) {
		for (int ix = id->xl; ix <= id->xr; ix += leap_var) {
			if (temp.getVoxel(ix, id->y, z)) {
				sx = ix;
				sy = id->y;
				/*
				 * Erase previous seeds, excluding the current one!
				 */
				seed_list.erase(seed_list.begin(), id);
				/*
				 * Resize the current seed, as all voxels from xl to the current ix are occupied
				 */
				id->xl = ix;
				return true;
			}
		}
	}
	seed_list.clear();
	return false;
};

/**
 * Here we implement a rapid 3D seed-filling algorithm, to extract the portion
 * of surface of the object in imageData which is connected to the given seed.
 * @param imageData input image
 * @param temp		temporary voxel grid, initially a copy of imageData
 * @param output	output voxelGrid
 * @param seed		starting voxel
 * @param leap_var	the leaping variable is used to leap over plenty of invalid
 * 					voxels in order to accelerate search process for the first
 * 					valid seed
 */
static inline void surfaceExtract3D(voxelGrid const & imageData, voxelGrid & temp, voxelGrid & output, voxel const & seed, int leap_var = 1) { // the pseudocode of proposed algorithm — 2009
	if (!imageData.getVoxel(seed))
		return;
	//create empty 3D stack
	stack<sliceNode> stack3D;

	IdList currentSeeds;
	//extract and fill the initial unfilled region by using an initial seed (seed)
	surfaceExtract2D(imageData, temp, output, seed, currentSeeds);
	//create two nodes with this initial region(dList), push them into 3D stack, along two opposite directions
	if (seed.iz < imageData.height - 1)
		stack3D.push(sliceNode(currentSeeds, seed.iz, +1));
	if (seed.iz > 0)
		stack3D.push(sliceNode(currentSeeds, seed.iz, -1));
	while (!stack3D.empty()) {
		//pop a node from 3D stack
		sliceNode n(stack3D.top());
		stack3D.pop();
		//set the search region (PassList) by assigning the region (dList) to it
		currentSeeds = n.seedList; /*
							 * std::list::operator= c++11 Assigns new contents to the container,
							 * replacing its current contents, and modifying its size accordingly.
							 */
		int new_z = n.z + n.direction_z;
		//obtain a search region (Plist_mapped) by mapping the region (PassList) onto the next slice
		//and then check the validation of this mapped region (Plist_mapped)
		int sx, sy; // new seed coordinates
		while(searchNeighborSlice(temp, currentSeeds, new_z, leap_var, sx, sy)){ //if the region is valid
			//set a value to the leaping variable
			IdList newSeeds;
			//search a valid seed in the mapped region (Plist_mapped)
			//extract and fill this new unfilled region (newSeeds here represents the first new region)
			surfaceExtract2D(imageData, temp, output, voxel(sx, sy, new_z), newSeeds);
			//continue to search other new unfilled regions in the region (Plist_mapped)
			while (searchNeighborSlice(temp, currentSeeds, new_z, leap_var, sx, sy)){
				//extract and fill these new unfilled regions
				IdList otherSeeds;
				surfaceExtract2D(imageData, temp, output, voxel(sx, sy, new_z), otherSeeds);
				//create two nodes with each region, push them into 3D stack, along two opposite directions
				if (new_z < imageData.height - 1)
					stack3D.push(sliceNode(otherSeeds, new_z, +1));
				if (new_z > 0)
					stack3D.push(sliceNode(otherSeeds, new_z, -1));
			}
			//execute the rollback operation
			//map the region (CurrentList) onto the previous slice
			currentSeeds = newSeeds;
			while (searchNeighborSlice(temp, newSeeds, new_z - n.direction_z, leap_var, sx, sy)) {
				IdList otherSeeds;
				//extract and fill these new unfilled regions
				surfaceExtract2D(imageData, temp, output, voxel(sx, sy, new_z - n.direction_z), otherSeeds);
				//create two nodes with each region, push them into 3D stack, along two opposite directions
				if (new_z - n.direction_z < imageData.height - 1)
					stack3D.push(sliceNode(otherSeeds, new_z - n.direction_z, +1));
				if (new_z - n.direction_z > 0)
					stack3D.push(sliceNode(otherSeeds, new_z - n.direction_z, -1));
			}
			//continue the loop from the region (CurrentList), along the same direction as before
			//obtain the new search region by assigning the region (CurrentList) to (PassList)
			new_z += n.direction_z;
			if (new_z < 0 || new_z == imageData.height)
				break;
		}
	}
};
}
#endif /* SEEDFILL3D_RAPID3DSURFACEEXTRACT_H_ */
