#pragma once

#include "cloth_data.h"
#include "body_data.h"
#include "DISTANCE.H"
#include "BVH.h"


class CollisionEngine
{
public:
	 CollisionEngine();
	~CollisionEngine();

	//Collision Resolution between body and cloth 
	void resolveClothBodyCollision(Cloth_Data* cloth,Body_Data* body);
	void resolveClothClothCollision(Cloth_Data* cloth);


	void run(float *vPos, int numPoints, float *tPos, int *triIndex, int numTriangles);
	void run(float *vPos_old, float *vPos, float *vVel, int numPoints, float *tPos, int *triIndex, int numTriangles, float timeStep);
	void resolveOnePoint(float *vPos, float *tPos, int pidx, std::vector<BVHTrianglePrimitive> &candidates);
	void resolveOnePoint(float *vPos_old, float *vPos, float *vVel, float *tPos, int pidx, std::vector<BVHTrianglePrimitive> &candidates, float timeStep);


private:
	bool hasBVH;
	BVH<BVHPointPrimitive> *pointsBVH_;
	BVH<BVHTrianglePrimitive> *triBVH_;
	std::vector<std::pair<BVHPointPrimitive, BVHTrianglePrimitive> > candicates_;

	float cullingDist_;
	float dampingCoeff_;
	float frictionCoeff_;
};

