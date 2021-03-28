/*
 * elapsedTime.h
 *
 *  Created on: 25/feb/2015
 *      Author: sebastian
 */

#ifndef UTILS_ELAPSEDTIME_H_
#define UTILS_ELAPSEDTIME_H_


#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>
/**
 * Simple method that calculates the time difference between the start and end time points
 * and returns the value in hh:mm:ss.fff format.
 * @param start		starting time point
 * @param end 		ending time point
 * @return string containing the elapsed time between the two time points in hh:mm:ss.fff format
 */
using namespace std;
template<typename T>
static string elapsedTime(std::chrono::time_point<T> const & start, std::chrono::time_point<T> const & end) {
	int t_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	int hh = t_ms / (1000 * 60 * 60);
	int mm = t_ms / (1000 * 60) - 60 * hh;
	double ss = t_ms / 1000.0 - 60 * (mm + 60 * hh);
	stringstream s;
	s << fixed << setprecision(3);
	if (hh < 10)
		s << "0";
	s << hh << ":";
	if (mm < 10)
		s << "0";
	s<< mm << ":";
	if (ss < 10)
		s << "0";
	s<< ss;
	return s.str();
};



#endif /* UTILS_ELAPSEDTIME_H_ */
