/*
 * HierarchicalQueue.cpp
 *
 *  Created on: Sep 6, 2014
 *      Author: sebastian
 */

/**
 * Implementation of the HierarchicalQueue data structure.
 *
 */
#include "HierarchicalQueue.h"
#include <stdexcept>
/**
 * Constructor.
 * @param numberOfQueues	number of queues composing the HierarchicalQueue
 */
HierarchicalQueue::HierarchicalQueue(uint16_t numberOfQueues) :
		numberOfQueues(numberOfQueues), min_id(numberOfQueues - 1), max_id(0) {
	queues = new std::queue<voxel>[numberOfQueues];
}
/**
 * Destructor.
 */
HierarchicalQueue::~HierarchicalQueue() {
	delete[] queues;
}
/**
 * Pops the highest priority element from the HierarchicalQueue.
 * Does nothing if the HierarchicalQueue is empty.
 */
void HierarchicalQueue::pop() {
	queues[min_id].pop();
	if (queues[min_id].empty()) {
		bool found = false;
		int ii = min_id + 1;
		while (!found && ii <= max_id) {
			if (!queues[ii].empty()) {
				min_id = ii;
				found = true;
			}
			else
				++ii;
		}
		if (!found) {
			min_id = numberOfQueues - 1;
			max_id = 0;
		}
	}
}
/**
 * Returns the number of elements currently in the HierarchicalQueue.
 */
size_t HierarchicalQueue::size() {
	size_t num = 0;
	for (int i = min_id; i <= max_id; ++i) {
		num += queues[i].size();
	}
	return num;
}
/**
 * Returns true if the HierarchicalQueue is empty (i.e. has no elements),
 * false otherwise.
 */
bool HierarchicalQueue::empty() {
	return (min_id > max_id);
}
/**
 * Returns a read/write reference to the highest priority element
 * of the HierarchicalQueue.
 */
voxel& HierarchicalQueue::front() {
	return queues[min_id].front();
}
/**
 * Returns a const reference to the highest priority element
 * of the HierarchicalQueue.
 */
const voxel& HierarchicalQueue::front() const{
	return queues[min_id].front();
}
/**
 * Inserts a new element in the HierarchicalQueue. The element is inserted at
 * the end of the queue with idx index, after its current last element.
 * @param vox	the voxel to be inserted
 * @param idx	the priority of the voxel to be inserted
 */
void HierarchicalQueue::push(voxel const & vox, uint16_t idx) {
#if defined RANGECHECK_TEST
	if (idx >= numberOfQueues)
		throw std::out_of_range("HierarchicalQueue::push() - Index is out of bounds!");
#endif
	queues[idx].push(vox);
	if (idx < min_id)
		min_id = idx;
	if (idx > max_id)
		max_id = idx;
}
/**
 * Returns a read/write reference to the queue with priority idx of the current
 * HierarchicalQueue.
 */
std::queue<voxel>& HierarchicalQueue::operator[] (size_t idx) {
#if defined RANGECHECK_TEST
	if (idx >= numberOfQueues)
		throw std::out_of_range("HierarchicalQueue::operator[] - Index is out of bounds!");
#endif
	return queues[idx];
}
/**
 * Returns a const reference to the queue with priority idx of the current
 * HierarchicalQueue.
 */
const std::queue<voxel>& HierarchicalQueue::operator[] (size_t idx) const {
#if defined RANGECHECK_TEST
	if (idx >= numberOfQueues)
		throw std::out_of_range("HierarchicalQueue::operator[] - Index is out of bounds!");
#endif
	return queues[idx];
}
