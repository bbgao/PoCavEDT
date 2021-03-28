/*
 * range.h
 *
 *  Created on: Sep 24, 2015
 *      Author: sebastian
 */

#ifndef RANGE_H_
#define RANGE_H_


typedef struct range {
	range();
	range(range const & r);
	range(int l, int r, int y);
	range & operator=(range const & r);

	int xl, xr, y;
} range;


#endif /* RANGE_H_ */
