/*
 * doubleCompare.h
 *
 *  Created on: Oct 7, 2014
 *      Author: sebastian
 */

#ifndef UTILS_DOUBLECOMPARE_H_
#define UTILS_DOUBLECOMPARE_H_

#include <math.h>

/**
 * Simple method for comparing two double-precision floating point numbers.
 * @param l first number
 * @param r second number
 * @return true if l == r, false otherwise.
 */
static inline bool doubleCompare(double const & l, double const & r) {
	return (fabs(l - r) < 1.0e-10);
};
/**
 * Simple method for comparing single-precision floating point numbers.
 * @param l first number
 * @param r second number
 * @return true if l == r, false otherwise.
 */
static inline bool floatCompare(float const & l, float const & r) {
	return (fabs(l - r) < 1.0e-5);
};


#endif /* UTILS_DOUBLECOMPARE_H_ */
