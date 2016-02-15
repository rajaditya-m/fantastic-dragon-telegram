#ifndef BVHPrimitive_h_
#define BVHPrimitive_h_
#include "BBox.h"

// base class for building BVH
class BVHBasePrimitive {
public:
	BVHBasePrimitive(){}
	virtual ~BVHBasePrimitive() {}

	//! Return a bounding box for this object
	virtual BBox getBBox() const = 0;

	//! Return the centroid for this object. (Used in BVH Sorting)
	virtual Vector3 getCentroid() const = 0;
};


class BVHPointPrimitive : public BVHBasePrimitive
{
public:
	BVHPointPrimitive(){
	}
	BVHPointPrimitive(float *points, int idx,  float t = 0.05) : points_(points), pIdex_(idx), pointSize_(t){
	}
	BVHPointPrimitive(const BVHPointPrimitive &p){
		points_ = p.getMesh();
		std::vector<int> index = p.getIndex();
		pIdex_ = index[0];
		pointSize_ = p.getPointSize();
	}
	virtual ~BVHPointPrimitive(){}

	virtual BBox getBBox() const{
		Vector3 position(points_[3 * pIdex_ + 0], points_[3 * pIdex_ + 1], points_[3 * pIdex_ + 2]);
		return BBox(position - Vector3(pointSize_, pointSize_, pointSize_), position + Vector3(pointSize_, pointSize_, pointSize_));
	}

	virtual Vector3 getCentroid() const{
		Vector3 position(points_[3 * pIdex_ + 0], points_[3 * pIdex_ + 1], points_[3 * pIdex_ + 2]);
		return position;
	}

	float* getMesh() const{
		return points_;
	}

	std::vector<int> getIndex() const{
		std::vector<int> idx;
		idx.push_back(pIdex_);
		return idx;
	}
	
	float getPointSize() const{
		return pointSize_;
	}

private:
	float *points_;
	int pIdex_;
	float pointSize_;
};

class BVHTrianglePrimitive : public BVHBasePrimitive
{
public:
	BVHTrianglePrimitive(){
	}
	BVHTrianglePrimitive(float *points, int idx0, int idx1, int idx2, float t = 0.05) :points_(points), idx0_(idx0), idx1_(idx1), idx2_(idx2), pointSize_(t){
	}
	BVHTrianglePrimitive(const BVHTrianglePrimitive &p){
		points_ = p.getMesh();
		std::vector<int> index = p.getIndex();
		idx0_ = index[0];
		idx1_ = index[1];
		idx2_ = index[2];
	}
	virtual ~BVHTrianglePrimitive(){}
	
	virtual BBox getBBox() const{
		Vector3 t0(points_[3 * idx0_ + 0], points_[3 * idx0_ + 1], points_[3 * idx0_ + 2]);
		Vector3 t1(points_[3 * idx1_ + 0], points_[3 * idx1_ + 1], points_[3 * idx1_ + 2]);
		Vector3 t2(points_[3 * idx2_ + 0], points_[3 * idx2_ + 1], points_[3 * idx2_ + 2]);

		BBox box0(t0 - Vector3(pointSize_, pointSize_, pointSize_), t0 + Vector3(pointSize_, pointSize_, pointSize_));
		BBox box1(t1 - Vector3(pointSize_, pointSize_, pointSize_), t1 + Vector3(pointSize_, pointSize_, pointSize_));
		BBox box2(t2 - Vector3(pointSize_, pointSize_, pointSize_), t2 + Vector3(pointSize_, pointSize_, pointSize_));

		BBox bbox(box0);
		bbox.expandToInclude(box1);
		bbox.expandToInclude(box2);
		return bbox;
	}

	virtual Vector3 getCentroid() const{
		Vector3 t0(points_[3 * idx0_ + 0], points_[3 * idx0_ + 1], points_[3 * idx0_ + 2]);
		Vector3 t1(points_[3 * idx1_ + 0], points_[3 * idx1_ + 1], points_[3 * idx1_ + 2]);
		Vector3 t2(points_[3 * idx2_ + 0], points_[3 * idx2_ + 1], points_[3 * idx2_ + 2]);

		Vector3 centroid = (t0+t1+t2)/3.0;
		return centroid;
	}

	float* getMesh() const{
		return points_;
	}

	std::vector<int> getIndex() const{
		std::vector<int> idx;
		idx.push_back(idx0_);
		idx.push_back(idx1_);
		idx.push_back(idx2_);
		return idx;
	}

private:
	float *points_;
	float pointSize_;
	int idx0_, idx1_, idx2_;
};



#endif
