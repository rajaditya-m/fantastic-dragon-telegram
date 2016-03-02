#ifndef BVHPrimitive_h_
#define BVHPrimitive_h_
#include "BBox.h"
#include "Types.h"

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
	BVHPointPrimitive(){}
	BVHPointPrimitive(ScalarType *points, int idx,  ScalarType t = 0.05) : points_(points), pIdex_(idx), pointSize_(t){
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

	ScalarType* getMesh() const{
		return points_;
	}

	std::vector<int> getIndex() const{
		std::vector<int> idx;
		idx.push_back(pIdex_);
		return idx;
	}
	
	ScalarType getPointSize() const{
		return pointSize_;
	}

private:
	ScalarType *points_;
	int pIdex_;
	ScalarType pointSize_;
};

class BVHEdgePrimitive : public BVHBasePrimitive
{
public:
	BVHEdgePrimitive(){}
	BVHEdgePrimitive(ScalarType *points, int eIdx, int idx0, int idx1, ScalarType t = 0.05) :points_(points), eIdx_(eIdx), idx0_(idx0), idx1_(idx1), pointSize_(t){}
	BVHEdgePrimitive(const BVHEdgePrimitive &p){
		points_ = p.getMesh();
		eIdx_ = p.getEIndex();
		std::vector<int> index = p.getIndices();
		idx0_ = index[0];
		idx1_ = index[1];
	}
	virtual ~BVHEdgePrimitive(){}

	virtual BBox getBBox() const{
		Vector3 e0(points_[3 * idx0_ + 0], points_[3 * idx0_ + 1], points_[3 * idx0_ + 2]);
		Vector3 e1(points_[3 * idx1_ + 0], points_[3 * idx1_ + 1], points_[3 * idx1_ + 2]);

		BBox box0(e0 - Vector3(pointSize_, pointSize_, pointSize_), e0 + Vector3(pointSize_, pointSize_, pointSize_));
		BBox box1(e1 - Vector3(pointSize_, pointSize_, pointSize_), e1 + Vector3(pointSize_, pointSize_, pointSize_));

		BBox bbox(box0);
		bbox.expandToInclude(box1);
		return bbox;
	}

	virtual Vector3 getCentroid() const{
		Vector3 e0(points_[3 * idx0_ + 0], points_[3 * idx0_ + 1], points_[3 * idx0_ + 2]);
		Vector3 e1(points_[3 * idx1_ + 0], points_[3 * idx1_ + 1], points_[3 * idx1_ + 2]);

		Vector3 centroid = (e0 + e1) / 3.0;
		return centroid;
	}

	ScalarType* getMesh() const{
		return points_;
	}

	int getEIndex() const{
		return eIdx_;
	}

	std::vector<int> getIndices() const{
		std::vector<int> idx;
		idx.push_back(idx0_);
		idx.push_back(idx1_);
		return idx;
	}

private:
	ScalarType *points_;
	ScalarType pointSize_;
	int eIdx_;
	int idx0_, idx1_;
};

class BVHTrianglePrimitive : public BVHBasePrimitive
{
public:
	BVHTrianglePrimitive(){}
	BVHTrianglePrimitive(ScalarType *points, int tIdx, int idx0, int idx1, int idx2, ScalarType t = 0.05) :points_(points), tIdx_(tIdx),idx0_(idx0), idx1_(idx1), idx2_(idx2), pointSize_(t){
	}
	BVHTrianglePrimitive(const BVHTrianglePrimitive &p){
		points_ = p.getMesh();
		tIdx_ = p.getTIndex();
		std::vector<int> index = p.getIndices();
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

	ScalarType* getMesh() const{
		return points_;
	}

	int getTIndex() const {
		return tIdx_;
	}

	std::vector<int> getIndices() const{
		std::vector<int> idx;
		idx.push_back(idx0_);
		idx.push_back(idx1_);
		idx.push_back(idx2_);
		return idx;
	}

private:
	ScalarType *points_;
	ScalarType pointSize_;
	int tIdx_;
	int idx0_, idx1_, idx2_;
};



#endif
