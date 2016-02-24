#include "CollisionEngine.h"

CollisionEngine::CollisionEngine()
{
	cullingDist_ = 0.05;
	dampingCoeff_ = 0;
	frictionCoeff_ = 1;
	hasBVH = false;

	pointsBVH_ = nullptr;
	triBVH_ = nullptr;
}


CollisionEngine::~CollisionEngine()
{
	delete pointsBVH_;
	delete triBVH_;
}

void CollisionEngine::resolveClothBodyCollision(Cloth_Data* cloth,Body_Data* body)
{
	int lastFrameId = cloth->getMesh()->get_number_frames();

	int numPoints = cloth->getMesh()->get_number_vertices();
	float time_step = cloth->get_property_obj()->get_timestep();
	float* clothPos = cloth->getPositions();
	float* clothOldPos = cloth->getMesh()->getPositions();
	float* clothVel = cloth->getVelocities();

	int numTriangles = body->getMesh()->get_number_triangles();
	float* bodyPos = body->getMesh()->getPositions();
	std::vector<int> triIndex = body->getMesh()->getTriangles();

	//run(clothPos,numPoints, bodyPos, &triIndex[0] ,numTriangles);
	run(clothOldPos,clothPos,clothVel,numPoints,bodyPos,&triIndex[0],numTriangles,time_step);
	cloth->setPositions(clothPos);
	cloth->setVelocities(clothVel);
}


void CollisionEngine::run(float *vPos, int numPoints, float *tPos, int *triIndex, int numTriangles){

	if (!hasBVH){	
		std::vector<BVHPointPrimitive *> *primitives = new std::vector<BVHPointPrimitive *>(numPoints, nullptr);
		for (size_t i = 0; i < numPoints; i++){
			(*primitives)[i] = (new BVHPointPrimitive(vPos, i, 0.05));
		}
		pointsBVH_ = new BVH<BVHPointPrimitive>(primitives);		
		
		std::vector<BVHTrianglePrimitive *> *primitives_2 = new std::vector<BVHTrianglePrimitive *>(numTriangles, nullptr);
		for (size_t i = 0; i < numTriangles; i++){
			(*primitives_2)[i] = (new BVHTrianglePrimitive(tPos, triIndex[3 * i + 0], triIndex[3 * i + 1], triIndex[3 * i + 2], 0.05));
		}
		triBVH_ = new BVH<BVHTrianglePrimitive>(primitives_2);

		hasBVH = true;
	}
	else{
		pointsBVH_->update();
		triBVH_->update();
	}

	for (size_t i = 0; i < numPoints; i++){
		std::vector<BVHTrianglePrimitive> candidates;
		BBox b = BBox(Vector3(vPos[3 * i + 0], vPos[3 * i + 1], vPos[3 * i + 2]));
		bool intersect = triBVH_->testIntersection(b, candidates);
		if (intersect){
			resolveOnePoint(vPos, tPos, i, candidates);
			//std::cout<<"collistion detected...................................."<<std::endl;
			//getchar();
		}
	}

}

void CollisionEngine::run(float *vPos_old, float *vPos, float *vVel, int numPoints, float *tPos, int *triIndex, int numTriangles, float timeStep){

	if (!hasBVH){
		std::vector<BVHPointPrimitive *> *primitives = new std::vector<BVHPointPrimitive *>(numPoints, nullptr);
		for (size_t i = 0; i < numPoints; i++){
			(*primitives)[i] = (new BVHPointPrimitive(vPos, i, 0.05));
		}
		pointsBVH_ = new BVH<BVHPointPrimitive>(primitives);		
		
		std::vector<BVHTrianglePrimitive *> *primitives_2 = new std::vector<BVHTrianglePrimitive *>(numTriangles, nullptr);
		for (size_t i = 0; i < numTriangles; i++){
			(*primitives_2)[i] = (new BVHTrianglePrimitive(tPos, triIndex[3 * i + 0], triIndex[3 * i + 1], triIndex[3 * i + 2], 0.05));
		}
		triBVH_ = new BVH<BVHTrianglePrimitive>(primitives_2);

		hasBVH = true;
		hasBVH = true;
	}
	else{
		pointsBVH_->update();
		triBVH_->update();
	}

	for (size_t i = 0; i < numPoints; i++){
		std::vector<BVHTrianglePrimitive> candidates;
		BBox b = BBox(Vector3(vPos[3 * i + 0], vPos[3 * i + 1], vPos[3 * i + 2]));
		bool intersect = triBVH_->testIntersection(b, candidates);
		if (intersect){
			resolveOnePoint(vPos_old, vPos, vVel, tPos, i, candidates, timeStep);
		}
		else{
			for (size_t j = 0; j < 3; j++){
				vVel[3 * i + j] = (vPos[3 * i + j] - vPos_old[3 * i + j]) / timeStep;
			}
		}
	}
}

void CollisionEngine::resolveOnePoint(float *vPos, float *tPos, int pidx, std::vector<BVHTrianglePrimitive> &candidates){

	float min_d = MY_INFINITE;
	float min_xni[3], min_N[3];
	float *min_xi, *min_xa, *min_xb, *min_xc;

	for (size_t i = 0; i < candidates.size(); i++){
		BVHTrianglePrimitive &t = candidates[i];
		std::vector<int> tidx = t.getIndex();

		float *xi = vPos + 3 * pidx;
		float *xa = tPos + 3 * tidx[0];
		float *xb = tPos + 3 * tidx[1];
		float *xc = tPos + 3 * tidx[2];

		float ba, bb, bc;
		float d;
		float xni[3], N[3];

		d = Squared_VT_Distance<float>(xi, xa, xb, xc, ba, bb, bc, xni);
		d = sqrt(d);

		if (min_d > d){
			min_d = d;
			for (size_t i = 0; i < 3; i++){
				min_xni[i] = xni[i];
			}
			min_xi = xi;
			min_xa = xa;
			min_xb = xb;
			min_xc = xc;
		}
	}

	//test in or out
	Normalize(min_xni);
	float xba[3], xca[3], xia[3];
	for (int n = 0; n < 3; n++) {
		xba[n] = min_xb[n] - min_xa[n];
		xca[n] = min_xc[n] - min_xa[n];
	}
	CROSS(xba, xca, min_N);
	Normalize(min_N);

	bool inside = (DOT(min_N, min_xni) < 0);

	if ((!inside) && (min_d > cullingDist_)){
		return;
	}

	//position correction
	float penetration_depth = (inside ? cullingDist_ + min_d : cullingDist_ - min_d);
	for (size_t i = 0; i < 3; i++){
		min_xi[i] += penetration_depth*min_N[i];
	}

}

void CollisionEngine::resolveOnePoint(float *vPos_old, float *vPos, float *vVel, float *tPos, int pidx, std::vector<BVHTrianglePrimitive> &candidates, float timeStep){

	float min_d = MY_INFINITE;
	float min_xni[3], min_N[3];
	float *min_xi, *min_xa, *min_xb, *min_xc;

	for (size_t i = 0; i < candidates.size(); i++){
		BVHTrianglePrimitive &t = candidates[i];
		std::vector<int> tidx = t.getIndex();

		float *xi = vPos + 3 * pidx;
		float *xa = tPos + 3 * tidx[0];
		float *xb = tPos + 3 * tidx[1];
		float *xc = tPos + 3 * tidx[2];

		float ba, bb, bc;
		float d;
		float xni[3], N[3];

		d = Squared_VT_Distance<float>(xi, xa, xb, xc, ba, bb, bc, xni);
		d = sqrt(d);

		if (min_d > d){
			min_d = d;
			for (size_t i = 0; i < 3; i++){
				min_xni[i] = xni[i];
			}
			min_xi = xi;
			min_xa = xa;
			min_xb = xb;
			min_xc = xc;
		}
	}

	//test in or out
	Normalize(min_xni);
	float xba[3], xca[3], xia[3];
	for (int n = 0; n < 3; n++) {
		xba[n] = min_xb[n] - min_xa[n];
		xca[n] = min_xc[n] - min_xa[n];
	}
	CROSS(xba, xca, min_N);
	Normalize(min_N);

	bool inside = (DOT(min_N, min_xni) < 0);

	if ((!inside) && (min_d > cullingDist_)){
		return;
	}

	//position correction
	float penetration_depth = (inside ? cullingDist_ + min_d : cullingDist_ - min_d);
	for (size_t i = 0; i < 3; i++){
		min_xi[i] += penetration_depth*min_N[i];
	}

	//velocity update
	for (size_t i = 0; i < 3; i++){
		vVel[3 * pidx + i] = (vPos[3 * pidx + i] - vPos_old[3 * pidx + i]) / timeStep;
	}

	//velocity correction
	float *vi = vVel + 3 * pidx;
	float vel_N = DOT(vi, min_N);
	if (vel_N < 0){
		for (size_t i = 0; i < 3; i++){
			//extract tangent velocity
			vVel[3 * pidx + i] -= (-vel_N)*(-min_N[i]);
			vVel[3 * pidx + i] *= frictionCoeff_;
			vVel[3 * pidx + i] += (-vel_N)*(-min_N[i])*dampingCoeff_;
		}
	}
}


