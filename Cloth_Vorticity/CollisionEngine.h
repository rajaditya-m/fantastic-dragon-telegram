#pragma once

#include "cloth_data.h"
#include "body_data.h"
#include "DISTANCE.H"
#include "INTERSECTION.h"
#include "BVH.h"


class CollisionEngine
{
public:
	 CollisionEngine();
	~CollisionEngine();

	//Collision Resolution between body and cloth 
	void resolveClothBodyCollision(Cloth_Data* cloth,Body_Data* body);
	void resolveClothClothCollision(Cloth_Data* cloth);


	void run(float *vPos_old, float *vPos, float *vVel, int numPoints, float *tPos, int *triIndex, int numTriangles, float timeStep);
	void resolveOnePoint(float *vPos_old, float *vPos, float *vVel, float *tPos, int pidx, std::vector<BVHTrianglePrimitive> &candidates, float timeStep);

	void run(ScalarType *ePos_old_1, ScalarType *ePos_1, ScalarType *eVel_1, int *eIndices_1, int numEdges_1,
		ScalarType *ePos_2, int *eIndices_2, int *tIndices_2, int *ET_2, int numEdges_2, ScalarType timeStep);
	void resoloveOneEdge(ScalarType *ePos_old_1, ScalarType *ePos_1, ScalarType *eVel_1, ScalarType *ePos_2, int *eIndices_1, int *tIndices_2, int e, 
		int *ET_2, std::vector<BVHEdgePrimitive> &candidates, ScalarType timeStep);

private:
	bool hasTriBVH_;
	bool hasEdgeBVH_;
	BVH<BVHTrianglePrimitive> *triBVH_;
	BVH<BVHEdgePrimitive>	*edgeBVH_;


	float cullingDist_;
	float dampingCoeff_;
	float frictionCoeff_;
};

