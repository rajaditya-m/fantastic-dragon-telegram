#ifndef BBox_h
#define BBox_h
#include <algorithm>
#include <cstdint>
#include "Vector3.h"

struct BBox {
	Vector3 min, max, extent;
	BBox() { }
	BBox(const Vector3& min, const Vector3& max)
		: min(min), max(max) {
		extent = max - min;
	}
	BBox(const Vector3& p)
		: min(p), max(p) {
		extent = max - min;
	}

	void expandToInclude(const Vector3& p) {
		min = ::min(min, p);
		max = ::max(max, p);
		extent = max - min;
	}
	void expandToInclude(const BBox& b) {
		min = ::min(min, b.min);
		max = ::max(max, b.max);
		extent = max - min;
	}
	uint32_t maxDimension() const {
		uint32_t result = 0;
		if (extent.y > extent.x) result = 1;
		if (extent.z > extent.y) result = 2;
		return result;
	}

	bool overlap(const BBox& b) const {
		bool overlap = ((max[0] >= b.min[0]) && (min[0] <= b.max[0]));
		overlap = overlap && ((max[1] >= b.min[1]) && (min[1] <= b.max[1]));
		overlap = overlap && ((max[2] >= b.min[2]) && (min[2] <= b.max[2]));
		return overlap;
	}
};

#endif


