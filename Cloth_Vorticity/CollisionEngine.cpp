#include "CollisionEngine.h"

CollisionEngine::CollisionEngine()
{
	cullingDist_ = 0.05;
	dampingCoeff_ = 0;
	frictionCoeff_ = 1;
	hasTriBVH_ = false;
	hasEdgeBVH_ = false;

	triBVH_ = nullptr;
	edgeBVH_ = nullptr;
}


CollisionEngine::~CollisionEngine()
{
	delete triBVH_;
	delete edgeBVH_;
}

void CollisionEngine::resolveClothBodyCollision(Cloth_Data* cloth,Body_Data* body)
{
	int lastFrameId = cloth->getMesh()->get_number_frames();

	int numPoints = cloth->getMesh()->get_number_vertices();
	float time_step = cloth->get_property_obj()->get_timestep();
	float* clothPos = cloth->getPositions();
	float* clothOldPos = cloth->getMesh()->getPositions();
	float* clothVel = cloth->getVelocities();
	std::vector<int> clothEdges = cloth->getMesh()->getEdges();

	int numTriangles = body->getMesh()->get_number_triangles();
	float* bodyPos = body->getMesh()->getPositions();
	std::vector<int> triIndex = body->getMesh()->getTriangles();
	std::vector<int> bodyEdges = body->getMesh()->getEdges();
	std::vector<int> bodyET = body->getMesh()->getET();

	run(clothOldPos,clothPos,clothVel,numPoints,bodyPos,&triIndex[0],numTriangles,time_step);
	run(clothOldPos,clothPos,clothVel,clothEdges.data(),clothEdges.size()/2,
		bodyPos, bodyEdges.data(),triIndex.data(), bodyET.data(), bodyEdges.size()/2,time_step);


	cloth->setPositions(clothPos);
	cloth->setVelocities(clothVel);
}


void CollisionEngine::run(float *vPos_old, float *vPos, float *vVel, int numPoints, float *tPos, int *triIndex, int numTriangles, float timeStep){

	if (!hasTriBVH_){	
		std::vector<BVHTrianglePrimitive *> *primitives_2 = new std::vector<BVHTrianglePrimitive *>(numTriangles, nullptr);
		for (size_t i = 0; i < numTriangles; i++){
			(*primitives_2)[i] = new BVHTrianglePrimitive(tPos,(int)i,triIndex[3 * i + 0], triIndex[3 * i + 1], triIndex[3 * i + 2]);

		}
		triBVH_ = new BVH<BVHTrianglePrimitive>(primitives_2);
		hasTriBVH_ = true;
	}
	else{
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

void CollisionEngine::resolveOnePoint(float *vPos_old, float *vPos, float *vVel, float *tPos, int pidx, std::vector<BVHTrianglePrimitive> &candidates, float timeStep){

	float min_d = MY_INFINITE;
	float min_xni[3], min_N[3];
	float *min_xi, *min_xa, *min_xb, *min_xc;

	for (size_t i = 0; i < candidates.size(); i++){
		BVHTrianglePrimitive &t = candidates[i];
		std::vector<int> tidx = t.getIndices();

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

	//vPos_old corrction so that it is consistent
	for (size_t i = 0; i < 3; i++){
		vPos_old[3 * pidx + i] = vPos[3 * pidx + i] - vVel[3 * pidx + i] * timeStep;
	}
}


void CollisionEngine::run(ScalarType *ePos_old_1, ScalarType *ePos_1, ScalarType *eVel_1, int *eIndices_1, int numEdges_1,
						  ScalarType *ePos_2, int *eIndices_2, int *tIndices_2, int *ET_2, int numEdges_2, ScalarType timeStep){

		if (!hasEdgeBVH_){
			std::vector<BVHEdgePrimitive *> *primitives = new std::vector<BVHEdgePrimitive *>();
			for (size_t i = 0; i < numEdges_2; i++){
				primitives->push_back(new BVHEdgePrimitive(ePos_2, i, eIndices_2[2 * i + 0], eIndices_2[2 * i + 1]));
			}
			edgeBVH_ = new BVH<BVHEdgePrimitive>(primitives);
			hasEdgeBVH_ = true;
		}
		else{
			edgeBVH_->update();
		}

		for (size_t i = 0; i < numEdges_1; i++){
			int eIndex0 = eIndices_1[2 * i + 0];
			int eIndex1 = eIndices_1[2 * i + 1];
			Vector3 ePos0 = Vector3(ePos_1[3 * eIndex0 + 0], ePos_1[3 * eIndex0 + 1], ePos_1[3 * eIndex0 + 2]);
			Vector3 ePos1 = Vector3(ePos_1[3 * eIndex1 + 0], ePos_1[3 * eIndex1 + 1], ePos_1[3 * eIndex1 + 2]);
			BBox b = BBox(ePos0);
			b.expandToInclude(ePos1);

			std::vector<BVHEdgePrimitive> candidates;
			bool intersect = edgeBVH_->testIntersection(b, candidates);
			if (intersect){
				resoloveOneEdge(ePos_old_1, ePos_1, eVel_1, ePos_2, eIndices_1, tIndices_2, i, ET_2, candidates, timeStep);
			}
			else{
				for (size_t j = 0; j < 3; j++){
					eVel_1[3 * eIndex0 + j] = (ePos_1[3 * eIndex0 + j] - ePos_old_1[3 * eIndex0 + j]) / timeStep;
					eVel_1[3 * eIndex1 + j] = (ePos_1[3 * eIndex1 + j] - ePos_old_1[3 * eIndex1 + j]) / timeStep;
				}
			}
		}
}


void CollisionEngine::resoloveOneEdge(ScalarType *ePos_old_1, ScalarType *ePos_1, ScalarType *eVel_1, ScalarType *ePos_2, int *eIndices_1, int *tIndices_2, int e, 
									  int *ET_2, std::vector<BVHEdgePrimitive> &candidates, ScalarType timeStep){
		
		int eIndex0 = eIndices_1[2 * e + 0];
		int eIndex1 = eIndices_1[2 * e + 1];

		ScalarType min_d = MY_INFINITE;
		ScalarType min_xni[3], min_N[3];
		ScalarType *min_xi, *min_xj, *min_xa, *min_xb;
		int min_e;

		for (size_t i = 0; i < candidates.size(); i++){
			BVHEdgePrimitive &e = candidates[i];
			std::vector<int> eIndices = e.getIndices();
			int eIdx = e.getEIndex();

			ScalarType *xi = ePos_1 + 3 * eIndex0;
			ScalarType *xj = ePos_1 + 3 * eIndex1;
			ScalarType *xa = ePos_2 + 3 * eIndices[0];
			ScalarType *xb = ePos_2 + 3 * eIndices[1];

			ScalarType br, bs;
			ScalarType d;
			ScalarType xni[3];

			d = Squared_EE_Distance<ScalarType>(xi, xj, xa, xb, br, bs, xni);
			d = sqrt(d);
			if (min_d > d){
				min_d = d;
				for (size_t i = 0; i < 3; i++)
					min_xni[i] = xni[i];
				min_xi = xi;
				min_xj = xj;
				min_xa = xa;
				min_xb = xb;	
				min_e = eIdx;
			}
		}

		int t1_min_e = ET_2[2 * min_e + 0];
		int t2_min_e = ET_2[2 * min_e + 1];

		//test in or out
		bool inside = false;

		ScalarType xba[3], xca[3], xN[3];
		if (t1_min_e != -1){
			int a = tIndices_2[3 * t1_min_e + 0];
			int b = tIndices_2[3 * t1_min_e + 1];
			int c = tIndices_2[3 * t1_min_e + 2];
			for (int i = 0; i < 3; ++i){
				xba[i] = ePos_2[3 * b + i] - ePos_2[3 * a + i];
				xca[i] = ePos_2[3 * c + i] - ePos_2[3 * a + i];
			}

			ScalarType x[3];
			inside = Triangle_Edge_Intersection(&ePos_2[3 * a], &ePos_2[3 * b], &ePos_2[3 * c], min_xi, min_xj, x);

		}else if (t2_min_e != -1) {
			int a = tIndices_2[3 * t2_min_e + 0];
			int b = tIndices_2[3 * t2_min_e + 1];
			int c = tIndices_2[3 * t2_min_e + 2];
			for (int i = 0; i < 3; ++i){
				xba[i] = ePos_2[3 * b + i] - ePos_2[3 * a + i];
				xca[i] = ePos_2[3 * c + i] - ePos_2[3 * a + i];
			}

			ScalarType x[3];
			inside = Triangle_Edge_Intersection(&ePos_2[3 * a], &ePos_2[3 * b], &ePos_2[3 * c], min_xi, min_xj, x);
		}

		CROSS(xba, xca, xN);
		Normalize(xN);
		Normalize(min_xni);
		//inside = DOT(min_N, min_xni) < 0;

		if ((!inside) && (min_d > cullingDist_)){
			//std::cout << "outside: NO EE Collision, distance = " << min_d << std::endl;
			//getchar();
			return;
		}

		//position correction
		ScalarType penetration_depth = (inside ? cullingDist_ + min_d : cullingDist_ - min_d);
		for (size_t i = 0; i < 3; i++){
			min_xi[i] += penetration_depth*min_xni[i];
			min_xj[i] += penetration_depth*min_xni[i];
		}

		//velocity update
		for (size_t i = 0; i < 3; i++){
			eVel_1[3 * eIndex0 + i] = (ePos_1[3 * eIndex0 + i] - ePos_old_1[3 * eIndex0 + i]) / timeStep;
			eVel_1[3 * eIndex1 + i] = (ePos_1[3 * eIndex1 + i] - ePos_old_1[3 * eIndex1 + i]) / timeStep;

		}

		//velocity correction
		ScalarType *vi = eVel_1 + 3 * eIndex0;
		ScalarType *vj = eVel_1 + 3 * eIndex1;
		ScalarType vel_N_i = DOT(vi, min_xni);
		ScalarType vel_N_j = DOT(vj, min_xni);
		if (vel_N_i < 0){
			for (size_t i = 0; i < 3; i++){
				//extract tangent velocity
				eVel_1[3 * eIndex0 + i] -= (-vel_N_i)*(-min_xni[i]);
				eVel_1[3 * eIndex0 + i] *= frictionCoeff_;
				eVel_1[3 * eIndex0 + i] += (-vel_N_i)*(-min_xni[i])*dampingCoeff_;

			}
		}
		if (vel_N_j < 0){
			for (size_t i = 0; i < 3; i++){
				//extract tangent velocity
				eVel_1[3 * eIndex1 + i] -= (-vel_N_j)*(-min_xni[i]);
				eVel_1[3 * eIndex1 + i] *= frictionCoeff_;
				eVel_1[3 * eIndex1 + i] += (-vel_N_j)*(-min_xni[i])*dampingCoeff_;

			}
		}

		//vPos_old corrction so that it is consistent
		for (size_t i = 0; i < 3; i++){
			ePos_old_1[3 * eIndex0 + i] = ePos_1[3 * eIndex0 + i] - eVel_1[3 * eIndex0 + i] * timeStep;
			ePos_old_1[3 * eIndex1 + i] = ePos_1[3 * eIndex1 + i] - eVel_1[3 * eIndex1 + i] * timeStep;
		}

		//inside ? std::cout << "warning: inside," : std::cout << "warning: outside,";
		//std::cout << "EE Collision: distance =  " << min_d << std::endl;
		//getchar();
}