#include "ImplicitMassSpringSolver.h"

#define USECONJUGATEGRADIENTSOLVER 1
#define SINC(X) (sin(X)/X)

ImplicitMassSpringSolver::ImplicitMassSpringSolver(void)
	{
	}


ImplicitMassSpringSolver::~ImplicitMassSpringSolver(void)
	{
		 delete[] massMatrixDiagonal_;
		 delete LHSMatrix_;
		 delete constantEnergyHessian_;
		 delete mKeQ_;
		 delete mKdQ_;
		 delete mhSquaredKeQ_;
		 delete mhKeQ_;
		 delete mhKdQ_;
		 delete[] constrainedVerts_;
	}

void ImplicitMassSpringSolver::initialize(Cloth_Data *cloth) {
	//Generate all the springs 
	generateAllSprings(cloth);

	//Initialize the sparse matrix strructure 
	initializeSparseMatrixFromOutline(cloth);

	//Initialize the bending hessian matrix 
	generateConstantEnergyHessian(cloth);
	
	int numVertices = cloth->getMesh()->get_number_vertices();
	massMatrixDiagonal_ = new double[3*numVertices];
	OMP_FOR
	for(int p=0; p<numVertices;p++) {
		float vertex_mass = cloth->get_vertex_mass(p);
		massMatrixDiagonal_[3*p+0] = vertex_mass;
		massMatrixDiagonal_[3*p+1] = vertex_mass;
		massMatrixDiagonal_[3*p+2] = vertex_mass;
	}

	RHSVector_ = Eigen::VectorXd::Zero(numVertices*3);

	//Constrained verts
	numConstrainedVerts_ = 2;
	constrainedVerts_ = new int[6];
	constrainedVerts_[3] = 7800;
	constrainedVerts_[4] = 7801;
	constrainedVerts_[5] = 7802;
	constrainedVerts_[0] = 150;
	constrainedVerts_[1] = 151;
	constrainedVerts_[2] = 152;
	/* constrainedVerts_[3] = 9;
	constrainedVerts_[4] = 10;
	constrainedVerts_[5] = 11;
	constrainedVerts_[0] = 6;
	constrainedVerts_[1] = 7;
	constrainedVerts_[2] = 8;*/

}

void ImplicitMassSpringSolver::generateAllSprings(Cloth_Data* cloth) {
	std::vector<Edge> edgeVector = cloth->getMesh()->getEdgeVector();
	int numEdges = cloth->getMesh()->getEdgeVectorSize();
	for(int i=0;i<numEdges;i++) {
		Edge curEdge = edgeVector[i];
		int start = curEdge.get_start_vertex();
		int end = curEdge.get_end_vertex();
		stretchSprings_.push_back(std::pair<int,int>(start,end));
		Eigen::Vector3d sPos = cloth->getMesh()->get_point_data(start);
		Eigen::Vector3d ePos = cloth->getMesh()->get_point_data(end);
		double len = (ePos-sPos).norm();
		restLenStretchSprings_.push_back(len);

		//Check if this qualifies as a bending spriong candidate 
		int adjTri1 = curEdge.get_tri_1();
		int adjTri2 = curEdge.get_tri_2();
		if(adjTri1 != -1 && adjTri2 != -1) {
			//Get the oppsing edges 
			Triangles t1 = cloth->getMesh()->get_triangle(adjTri1);
			Triangles t2 = cloth->getMesh()->get_triangle(adjTri2);

			std::vector<int> v1;
			v1.push_back(t1.a);
			v1.push_back(t1.b);
			v1.push_back(t1.c);

			std::vector<int> v2;
			v2.push_back(t2.a);
			v2.push_back(t2.b);
			v2.push_back(t2.c);

			std::sort(v1.begin(),v1.end());
			std::sort(v2.begin(),v2.end());
			
			int oppV1;
			int oppV2;

			if(v1[0]==start && v1[1]==end) {
				oppV1 = v1[2];
			} else if(v1[1]==start && v1[2]==end) {
				oppV1 = v1[0];
			} else if(v1[0]==start && v1[2]==end) {
				oppV1 = v1[1];
			} else {
				std::cout << "[ERROR] This should never be reached\n";
			}

			if(v2[0]==start && v2[1]==end) {
				oppV2 = v2[2];
			} else if(v2[1]==start && v2[2]==end) {
				oppV2 = v2[0];
			} else if(v2[0]==start && v2[2]==end) {
				oppV2 = v2[1];
			} else {
				std::cout << "[ERROR] This should never be reached\n";
			}

			bendingSprings_.push_back(std::pair<int,int>(oppV1,oppV2));
			Eigen::Vector3d sPos = cloth->getMesh()->get_point_data(oppV1);
			Eigen::Vector3d ePos = cloth->getMesh()->get_point_data(oppV2);
			double len = (ePos-sPos).norm();
			restLenbendingSprings_.push_back(len);
		}
	}
	
	//Add all of them to the single spring vector 
	allSprings_.reserve(bendingSprings_.size() + stretchSprings_.size());
	allSprings_.insert(allSprings_.end(),stretchSprings_.begin(),stretchSprings_.end());
	//allSprings_.insert(allSprings_.end(),bendingSprings_.begin(),bendingSprings_.end());

	restLenAllSprings_.reserve(restLenbendingSprings_.size() + restLenStretchSprings_.size());
	restLenAllSprings_.insert(restLenAllSprings_.end(),restLenStretchSprings_.begin(),restLenStretchSprings_.end());
	//restLenAllSprings_.insert(restLenAllSprings_.end(),restLenbendingSprings_.begin(),restLenbendingSprings_.end());

}

void ImplicitMassSpringSolver::initializeSparseMatrixFromOutline(Cloth_Data *cloth) {
	int numVertices = cloth->getMesh()->get_number_vertices();
	int numTriangles = cloth->getMesh()->get_number_triangles();
	SparseMatrixOutline *spMatOutline = new SparseMatrixOutline(numVertices*3);

	double *zero3x3 = new double[9];
	memset(zero3x3,0,sizeof(double)*9);

	//Initialize the diagonal 3x3 blocks
	for(int p=0; p<numVertices;p++) {
		spMatOutline->AddBlock3x3Entry(p,p,zero3x3);
	}
 
	//Initialize the springs 
	int sizeStretchSprings = stretchSprings_.size();
	for(int i=0;i<sizeStretchSprings;i++) {
		int a1 = stretchSprings_[i].first;
		int a2 = stretchSprings_[i].second;
		spMatOutline->AddBlock3x3Entry(a1,a2,zero3x3);
		spMatOutline->AddBlock3x3Entry(a2,a1,zero3x3);
	}

	int sizeBendSprings = bendingSprings_.size();
	for(int i=0;i<sizeBendSprings;i++) {
		int a1 = bendingSprings_[i].first;
		int a2 = bendingSprings_[i].second;
		spMatOutline->AddBlock3x3Entry(a1,a2,zero3x3);
		spMatOutline->AddBlock3x3Entry(a2,a1,zero3x3);
	}


	//Initialize the actual sparse matrix 
	LHSMatrix_= new SparseMatrix(spMatOutline);

	//Delete the outline 
	delete (spMatOutline);
	delete[] zero3x3;

	//Create the diagonal indices
	LHSMatrix_->BuildDiagonalIndices();
}

void ImplicitMassSpringSolver::generateConstantEnergyHessian(Cloth_Data *cloth) {
	int numVertices = cloth->getMesh()->get_number_vertices();
	int numTriangles = cloth->getMesh()->get_number_triangles();
	std::vector<Edge> edgeVector = cloth->getMesh()->getEdgeVector();
	int numEdges = cloth->getMesh()->getEdgeVectorSize();

	SparseMatrixOutline *spMatOutline = new SparseMatrixOutline(numVertices*3);

	double *fakeIdentity = new double[9];
	memset(fakeIdentity,0,sizeof(double)*9);
	double ke = cloth->get_property_obj()->getQuadBendStiffness();
	double kd = cloth->get_property_obj()->getQuadBendDamping();
	double h = cloth->get_property_obj()->get_timestep();

	for(int i=0;i<numEdges;i++) {

		Edge curEdge = edgeVector[i];

		int i0 = curEdge.start;
		int i1 = curEdge.end;

		//Check if this qualifies as a bending edge candidate
		int adjTri1 = curEdge.get_tri_1();
		int adjTri2 = curEdge.get_tri_2();
		if(adjTri1 != -1 && adjTri2 != -1) {

			//Get the oppsing edges 
			Triangles t1 = cloth->getMesh()->get_triangle(adjTri1);
			Triangles t2 = cloth->getMesh()->get_triangle(adjTri2);

			std::vector<int> v1;
			v1.push_back(t1.a);
			v1.push_back(t1.b);
			v1.push_back(t1.c);

			std::vector<int> v2;
			v2.push_back(t2.a);
			v2.push_back(t2.b);
			v2.push_back(t2.c);

			std::sort(v1.begin(),v1.end());
			std::sort(v2.begin(),v2.end());
			
			int i2;
			int i3;

			if(v1[0]==i0 && v1[1]==i1) {
				i2 = v1[2];
			} else if(v1[1]==i0 && v1[2]==i1) {
				i2 = v1[0];
			} else if(v1[0]==i0 && v1[2]==i1) {
				i2 = v1[1];
			} else {
				std::cout << "[ERROR] This should never be reached\n";
			}

			if(v2[0]==i0 && v2[1]==i1) {
				i3 = v2[2];
			} else if(v2[1]==i0 && v2[2]==i1) {
				i3 = v2[0];
			} else if(v2[0]==i0 && v2[2]==i1) {
				i3 = v2[1];
			} else {
				std::cout << "[ERROR] This should never be reached\n";
			}

			//All the vertices are ready 
			Eigen::Vector3d x0 = cloth->getMesh()->get_point_data(i0);
			Eigen::Vector3d x1 = cloth->getMesh()->get_point_data(i1);
			Eigen::Vector3d x2 = cloth->getMesh()->get_point_data(i2);
			Eigen::Vector3d x3 = cloth->getMesh()->get_point_data(i3);

			Eigen::Vector3d e0 = x1-x0;
			Eigen::Vector3d e1 = x2-x0;
			Eigen::Vector3d e2 = x3-x0;
			Eigen::Vector3d e3 = x2-x1;
			Eigen::Vector3d e4 = x3-x1;

			double c01 = cotTheta( e0, e1);
			double c02 = cotTheta( e0, e2);
			double c03 = cotTheta(-e0, e3);
			double c04 = cotTheta(-e0, e4);

			double K0[] = {c03+c04, c01+c02, -c01-c03, -c02-c04};

			double A0 = 0.5 * (e0.cross(e1)).norm();
			double A1 = 0.5 * (e0.cross(e2)).norm();

			double coef = -3. / (2.*(A0+A1));
		
			double Q[4][4];

			for (int i=0; i<4; ++i) {
				for (int j=0; j<i; ++j) {
					Q[i][j] = Q[j][i] = coef * K0[i] * K0[j];
				}
				Q[i][i] = coef * K0[i] * K0[i];
			}

			//Now we will add them as index matrices to the result
			int idx[] = {i0,i1,i2,i3};
			for(int i=0;i<4;i++) {
				for(int j=0;j<4;j++) {
					fakeIdentity[0] = fakeIdentity[4] = fakeIdentity[8] = Q[i][j];
					spMatOutline->AddBlock3x3Entry(idx[i],idx[j],fakeIdentity);
				}
			}
		}
	}

	//Create the actual SparseMatrix 
	constantEnergyHessian_ = new SparseMatrix(spMatOutline);

	//create its miniions 
	
	SparseMatrix temp = ke*(*constantEnergyHessian_);
	mKeQ_ = new SparseMatrix(temp);
	temp = kd*(*constantEnergyHessian_);
	mKdQ_ = new SparseMatrix(temp);
	temp = ke*(-h)*h*(*constantEnergyHessian_);
	mhSquaredKeQ_ = new SparseMatrix(temp);
	temp = ke*h*(*constantEnergyHessian_);
	mhKeQ_ = new SparseMatrix(temp);
	temp = kd*(-h)*(*constantEnergyHessian_);
	mhKdQ_ = new SparseMatrix(temp);

	//Delete the residuals
	delete (spMatOutline);
	delete[] fakeIdentity;
}

void ImplicitMassSpringSolver::advance_time_step(Cloth_Data* cloth) {
	//Set the last frame information
	lastFrameId_ = (cloth->getMesh()->get_number_frames()) - 1;
	
	std::cout << lastFrameId_ << "\n";

	//Prepare the LHS matrix for usage
	LHSMatrix_->ResetToZero();
	LHSMatrix_->AddDiagonalMatrix(massMatrixDiagonal_);
	
	//Prepare the RHS Vector for usage
	RHSVector_.setZero();

	//Add the physics components
	addAllSprings(cloth);
	addQuadraticBending(cloth);
	addGravityComponents(cloth);

	//Solve and report
	finalizeAndSolve(cloth);
}

void ImplicitMassSpringSolver::addGravityComponents(Cloth_Data* cloth) {
	int num_vertices = cloth->getMesh()->get_number_vertices();
	OMP_FOR
	for(int p=0;p<num_vertices;p++)
	{
		float gravForce = cloth->get_vertex_mass(p) * 9.8f;
		RHSVector_[3*p+1] -= gravForce; 
	}
}

void ImplicitMassSpringSolver::addShearSprings(Cloth_Data* cloth) {
	int numShearSprings = stretchSprings_.size();
	double k = 100.0;
	double kd = cloth->get_property_obj()->get_damping_param();
	std::cout << kd << "\n";
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	float h = cloth->get_property_obj()->get_timestep();

	for(int idx =0; idx< numShearSprings;idx++) {
		int i = stretchSprings_[idx].first;
		int j = stretchSprings_[idx].second;
		double r = restLenStretchSprings_[idx];

		Eigen::Vector3d xi = cloth->getMesh()->get_point_data(i,lastFrameId_);
		Eigen::Vector3d xj = cloth->getMesh()->get_point_data(j,lastFrameId_);

		Eigen::Vector3d vi = cloth->get_velocity(i);
		Eigen::Vector3d vj = cloth->get_velocity(j);

		Eigen::Vector3d xij = xi - xj;
		Eigen::Vector3d vij = vi - vj;
		Eigen::Vector3d xij_Hat = xij.normalized(); 

		double xij_Norm = xij.norm();

		//Spring Forces
		Eigen::Vector3d f_i = -k*(xij_Norm - r) * xij_Hat;
		Eigen::Vector3d f_j = -f_i;

		//Damping forces
		Eigen::Vector3d fd_i = (-kd)*(xij_Hat)*(vij.dot(xij_Hat));
		Eigen::Vector3d fd_j = -fd_i;

		//Forces on the RHS
		/*RHSVector_[3*i+0] += f_i[0];
		RHSVector_[3*i+1] += f_i[1];
		RHSVector_[3*i+2] += f_i[2];
		RHSVector_[3*j+0] += f_j[0];
		RHSVector_[3*j+1] += f_j[1];
		RHSVector_[3*j+2] += f_j[2];*/

		RHSVector_[3*i+0] += (f_i[0]+fd_i[0]);
		RHSVector_[3*i+1] += (f_i[1]+fd_i[1]);
		RHSVector_[3*i+2] += (f_i[2]+fd_i[2]);
		RHSVector_[3*j+0] += (f_j[0]+fd_j[0]);
		RHSVector_[3*j+1] += (f_j[1]+fd_j[1]);
		RHSVector_[3*j+2] += (f_j[2]+fd_j[2]);

		//Jacobians on the LHS
		Eigen::Matrix3d outerPdk = xij_Hat*xij_Hat.transpose();
		Eigen::Matrix3d J_Fi_xi = outerPdk + (1.0-(r/xij_Norm))*(I - outerPdk);
		J_Fi_xi *= (-k);
		Eigen::Matrix3d J_Fj_xi = -J_Fi_xi;

		//Damping Jacobians
		Eigen::Matrix3d J_D = (-kd*I);
		Eigen::Matrix3d negJ_D = -J_D;;

		LHSMatrix_->Add3x3Entry(3*i,3*i,J_Fi_xi.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*j,3*j,J_Fi_xi.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*j,3*i,J_Fj_xi.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*i,3*j,J_Fj_xi.data(),(-h*h));

		LHSMatrix_->Add3x3Entry(3*i,3*i,J_D.data(),(-h));
		LHSMatrix_->Add3x3Entry(3*j,3*j,J_D.data(),(-h));
		LHSMatrix_->Add3x3Entry(3*j,3*i,negJ_D.data(),(-h));
		LHSMatrix_->Add3x3Entry(3*i,3*j,negJ_D.data(),(-h));

		Eigen::Vector3d df_dx_i = h*(J_Fi_xi * vi + J_Fj_xi * vj);
		Eigen::Vector3d df_dx_j = h*(J_Fj_xi * vi + J_Fi_xi * vj);

		RHSVector_[3*i+0] += df_dx_i[0];
		RHSVector_[3*i+1] += df_dx_i[1];
		RHSVector_[3*i+2] += df_dx_i[2];

		RHSVector_[3*j+0] += df_dx_j[0];
		RHSVector_[3*j+1] += df_dx_j[1];
		RHSVector_[3*j+2] += df_dx_j[2];

	}
}

void ImplicitMassSpringSolver::addBendingSprings(Cloth_Data* cloth) {
	int numbendingSprings = allSprings_.size();
	double k = cloth->get_property_obj()->getKStiffness();
	double kd = cloth->get_property_obj()->get_damping_param();
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	float h = cloth->get_property_obj()->get_timestep();

	for(int idx =0; idx< numbendingSprings;idx++) {
		int i = bendingSprings_[idx].first;
		int j = bendingSprings_[idx].second;
		double r = restLenbendingSprings_[idx];

		Eigen::Vector3d xi = cloth->getMesh()->get_point_data(i,lastFrameId_);
		Eigen::Vector3d xj = cloth->getMesh()->get_point_data(j,lastFrameId_);

		Eigen::Vector3d vi = cloth->get_velocity(i);
		Eigen::Vector3d vj = cloth->get_velocity(j);

		Eigen::Vector3d xij = xi - xj;
		Eigen::Vector3d vij = vi - vj;
		double xij_Norm = xij.norm();
		Eigen::Vector3d xij_Hat = xij.normalized(); 

		Eigen::Vector3d f_i = -k*(xij_Norm - r) * xij_Hat;
		Eigen::Vector3d f_j = -f_i;

		//Damping forces
		Eigen::Vector3d fd_i = (-kd)*(xij_Hat)*(vij.dot(xij_Hat));
		Eigen::Vector3d fd_j = -fd_i;


		//Forces on the RHS
		RHSVector_[3*i+0] += (f_i[0]+fd_i[0]);
		RHSVector_[3*i+1] += (f_i[1]+fd_i[1]);
		RHSVector_[3*i+2] += (f_i[2]+fd_i[2]);
		RHSVector_[3*j+0] += (f_j[0]+fd_j[0]);
		RHSVector_[3*j+1] += (f_j[1]+fd_j[1]);
		RHSVector_[3*j+2] += (f_j[2]+fd_j[2]);
		/*RHSVector_[3*i+0] += f_i[0];
		RHSVector_[3*i+1] += f_i[1];
		RHSVector_[3*i+2] += f_i[2];
		RHSVector_[3*j+0] += f_j[0];
		RHSVector_[3*j+1] += f_j[1];
		RHSVector_[3*j+2] += f_j[2];*/

		//Jacobians on the LHS
		Eigen::Matrix3d outerPdk = xij_Hat*xij_Hat.transpose();
		Eigen::Matrix3d J_Fi_xi = outerPdk + (1.0-(r/xij_Norm))*(I - outerPdk);
		J_Fi_xi *= (-k);
		Eigen::Matrix3d J_Fj_xi = -J_Fi_xi;

		//Damping Jacobians
		Eigen::Matrix3d J_D = (-kd*I);
		Eigen::Matrix3d negJ_D = -J_D;;

		LHSMatrix_->Add3x3Entry(3*i,3*i,J_Fi_xi.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*j,3*j,J_Fi_xi.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*j,3*i,J_Fj_xi.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*i,3*j,J_Fj_xi.data(),(-h*h));

		LHSMatrix_->Add3x3Entry(3*i,3*i,J_D.data(),(-h));
		LHSMatrix_->Add3x3Entry(3*j,3*j,J_D.data(),(-h));
		LHSMatrix_->Add3x3Entry(3*j,3*i,negJ_D.data(),(-h));
		LHSMatrix_->Add3x3Entry(3*i,3*j,negJ_D.data(),(-h));

		Eigen::Vector3d df_dx_i = h*(J_Fi_xi * vi + J_Fj_xi * vj);
		Eigen::Vector3d df_dx_j = h*(J_Fj_xi * vi + J_Fi_xi * vj);

		RHSVector_[3*i+0] += df_dx_i[0];
		RHSVector_[3*i+1] += df_dx_i[1];
		RHSVector_[3*i+2] += df_dx_i[2];

		RHSVector_[3*j+0] += df_dx_j[0];
		RHSVector_[3*j+1] += df_dx_j[1];
		RHSVector_[3*j+2] += df_dx_j[2];

	}
}

void ImplicitMassSpringSolver::addAllSprings(Cloth_Data *cloth) {
	int numSprings = allSprings_.size();
	double k = cloth->get_property_obj()->getKStiffness();
	double kd = cloth->get_property_obj()->get_damping_param();
	double kb = 0.1;
	double cb = k;
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	float h = cloth->get_property_obj()->get_timestep();

	//OMP_FOR
	for(int idx =0; idx< numSprings;idx++) {
		int i = allSprings_[idx].first;
		int j = allSprings_[idx].second;
		double r = restLenAllSprings_[idx];

		Eigen::Vector3d xi = cloth->getMesh()->get_point_data(i,lastFrameId_);
		Eigen::Vector3d xj = cloth->getMesh()->get_point_data(j,lastFrameId_);

		Eigen::Vector3d vi = cloth->get_velocity(i);
		Eigen::Vector3d vj = cloth->get_velocity(j);

		Eigen::Vector3d xij = xi - xj;
		Eigen::Vector3d vij = vi - vj;
		double xij_Norm = xij.norm();
		Eigen::Vector3d xij_Hat = xij.normalized(); 

		//Now we model the type 1  
		if(xij_Norm > r) {
			//Force 
			Eigen::Vector3d f_i = (-k)*(xij_Norm - r) * xij_Hat;
			Eigen::Vector3d f_j = -f_i;

			//Jacobians on the LHS
			Eigen::Matrix3d outerPdk = xij_Hat*xij_Hat.transpose();
			Eigen::Matrix3d J_Fi_xi = outerPdk + (1.0-(r/xij_Norm))*(I - outerPdk);
			J_Fi_xi *= -k;
			Eigen::Matrix3d J_Fj_xi = -J_Fi_xi;

			//Multiplicands for the right hand side
			Eigen::Vector3d df_dx_i = h*(J_Fi_xi * vi + J_Fj_xi * vj);
			Eigen::Vector3d df_dx_j = h*(J_Fj_xi * vi + J_Fi_xi * vj);

			//Add them appropriately 
			//Forces on the RHS
			RHSVector_[3*i+0] += (f_i[0]+df_dx_i[0]);
			RHSVector_[3*i+1] += (f_i[1]+df_dx_i[1]);
			RHSVector_[3*i+2] += (f_i[2]+df_dx_i[2]);
			RHSVector_[3*j+0] += (f_j[0]+df_dx_j[0]);
			RHSVector_[3*j+1] += (f_j[1]+df_dx_j[1]);
			RHSVector_[3*j+2] += (f_j[2]+df_dx_j[2]);

			//Then the jacobians on the left
			LHSMatrix_->Add3x3Entry(3*i,3*i,J_Fi_xi.data(),(-h*h));
			LHSMatrix_->Add3x3Entry(3*j,3*j,J_Fi_xi.data(),(-h*h));
			LHSMatrix_->Add3x3Entry(3*j,3*i,J_Fj_xi.data(),(-h*h));
			LHSMatrix_->Add3x3Entry(3*i,3*j,J_Fj_xi.data(),(-h*h));
		}
		else { //Type 2 interactions
			double ratio = xij_Norm/r;
			double inverseSincRatio = inverseSinc(ratio);
			double derInverseSincRatio = derivativeInverseSinc(ratio,inverseSincRatio);

			double kappa = (2.0*inverseSincRatio)/r;
			double ratio2 = kappa*r*0.5;
			double fb = (kb * kappa * kappa) /(cos(ratio2)-SINC(ratio2));
			double quant = cb*(xij_Norm-r);
			double fbStar;
			double fbStarJacobian;
			if(fb < quant) {
				fbStar = quant;
				fbStarJacobian = cb;
			}
			else {
				fbStar = fb;
				//This is the part which is slightly complex 
				double a = (8.0*kb*inverseSincRatio*derInverseSincRatio);
				double b = (r*r*r)*(cos(inverseSincRatio)-ratio);
				double c = (4*kb*inverseSincRatio*inverseSincRatio)*((sin(inverseSincRatio)*derInverseSincRatio)+(derInverseSincRatio*derivativeSinc(inverseSincRatio)));
				double d = (r*r*r)*(cos(inverseSincRatio)-ratio)*(cos(inverseSincRatio)-ratio);
				fbStarJacobian = (a/b)+(c/d);
			}
			if(fbStarJacobian<0) {
				std::cout << "[INFO] fb_Star Jacobian is negative.\n";
			}

			Eigen::Vector3d f_i = -fbStar * xij_Hat;
			Eigen::Vector3d f_j = -f_i;

			Eigen::Matrix3d outerPdk = xij_Hat*xij_Hat.transpose();
			Eigen::Matrix3d J_Fi_xi = -fbStarJacobian*outerPdk;
			Eigen::Matrix3d J_Fj_xi = -J_Fi_xi;

			//Multiplicands for the right hand side
			Eigen::Vector3d df_dx_i = h*(J_Fi_xi * vi + J_Fj_xi * vj);
			Eigen::Vector3d df_dx_j = h*(J_Fj_xi * vi + J_Fi_xi * vj);

			//Add them appropriately 
			//Forces on the RHS
			RHSVector_[3*i+0] += (f_i[0]+df_dx_i[0]);
			RHSVector_[3*i+1] += (f_i[1]+df_dx_i[1]);
			RHSVector_[3*i+2] += (f_i[2]+df_dx_i[2]);
			RHSVector_[3*j+0] += (f_j[0]+df_dx_j[0]);
			RHSVector_[3*j+1] += (f_j[1]+df_dx_j[1]);
			RHSVector_[3*j+2] += (f_j[2]+df_dx_j[2]);

			//Then the jacobians on the left
			LHSMatrix_->Add3x3Entry(3*i,3*i,J_Fi_xi.data(),(-h*h));
			LHSMatrix_->Add3x3Entry(3*j,3*j,J_Fi_xi.data(),(-h*h));
			LHSMatrix_->Add3x3Entry(3*j,3*i,J_Fj_xi.data(),(-h*h));
			LHSMatrix_->Add3x3Entry(3*i,3*j,J_Fj_xi.data(),(-h*h));

		}

		//Damping forces
		Eigen::Vector3d fd_i = (-kd * vij);
		Eigen::Vector3d fd_j = -fd_i;

		//Damping Jacobians
		Eigen::Matrix3d J_D = (kd*I);
		Eigen::Matrix3d negJ_D = -J_D;

		//Forces on the RHS
		RHSVector_[3*i+0] += fd_i[0];
		RHSVector_[3*i+1] += fd_i[1];
		RHSVector_[3*i+2] += fd_i[2];
		RHSVector_[3*j+0] += fd_j[0];
		RHSVector_[3*j+1] += fd_j[1];
		RHSVector_[3*j+2] += fd_j[2];

		LHSMatrix_->Add3x3Entry(3*i,3*i,J_D.data(),(-h));
		LHSMatrix_->Add3x3Entry(3*j,3*j,J_D.data(),(-h));
		LHSMatrix_->Add3x3Entry(3*j,3*i,negJ_D.data(),(-h));
		LHSMatrix_->Add3x3Entry(3*i,3*j,negJ_D.data(),(-h));
	}

}

void ImplicitMassSpringSolver::addQuadraticBending(Cloth_Data *cloth) {
	// Assemble the force and velocity 
	int numVertices = cloth->getMesh()->get_number_vertices();
	double *x = new double[numVertices*3];
	double *v = new double[numVertices*3];
	std::vector<Eigen::Vector3d> posVec = cloth->getMesh()->getPositionVector(lastFrameId_);
	memcpy(x,&posVec[0],sizeof(double)*3*numVertices);
	std::vector<Eigen::Vector3d> velVec = cloth->get_velocity_vector();
	memcpy(v,&velVec[0],sizeof(double)*3*numVertices);
	
	//Find the elastic force and add it to LHS 
	double *fe = new double[numVertices*3];
	mKeQ_->MultiplyVector(x,fe);
	Eigen::Map<Eigen::VectorXd> fe_Map(fe,numVertices*3);
	RHSVector_ += fe_Map;
	
	//Find the elastic jacobian and add it to LHS
	(*LHSMatrix_) += (*mhSquaredKeQ_);

	//Find the product of velocity*h*df/dx 
	double *pdk = new double[numVertices*3];
	mhKeQ_->MultiplyVector(v,pdk);
	Eigen::Map<Eigen::VectorXd> pdk_Map(pdk,numVertices*3);
	RHSVector_ += pdk_Map;

	//Find the damping force and add it to the system
	double *fd = new double[numVertices*3];
	mKdQ_->MultiplyVector(v,fd);
	Eigen::Map<Eigen::VectorXd> fd_Map(fd,numVertices*3);
	RHSVector_ += fd_Map;

	//Add the damping jacobian to the  LHS mAtrix 
	(*LHSMatrix_) += (*mhKdQ_);

	//Delete this stuff 
	delete[] x;
	delete[] v;
	delete[] fe;
	delete[] pdk;
	delete[] fd;
}

void ImplicitMassSpringSolver::finalizeAndSolve(Cloth_Data* cloth)  {
	float time_step = cloth->get_property_obj()->get_timestep();
	int num_vertices = cloth->getMesh()->get_number_vertices();
	
	//Add the UI Force if any 
	int lastVClicked = cloth->getLastClickedVertex();
	Eigen::Vector3d uiForce = cloth->getCurrentUIForce();
	if(lastVClicked!=-1) {
		RHSVector_[3*lastVClicked+0] += uiForce[0];
		RHSVector_[3*lastVClicked+1] += uiForce[1];
		RHSVector_[3*lastVClicked+2] += uiForce[2];
	}

	RHSVector_ *= time_step;

	
	//Setup the constrained version of the whole thingy 
	int n3 = num_vertices*3;
	int numConstrainedDOFs = n3 - (3*numConstrainedVerts_);
	//Make a copy of the LHSMatrix_ and remove rows and columns from it
	SparseMatrix *tempSpMatCopy = new SparseMatrix(*LHSMatrix_);
	tempSpMatCopy->RemoveRowsColumns(numConstrainedVerts_*3,constrainedVerts_);
	double suum = tempSpMatCopy->SumEntries();
	std::cout << "Checksum:" << suum << "\n";
	//Make a copy of RHS Vector and remove rows and columns from it 
	double* rhsConstrained = new double[numConstrainedDOFs];
	RemoveRows(num_vertices*3,rhsConstrained,RHSVector_.data(),numConstrainedVerts_*3,constrainedVerts_);
	//Make a correct size row vector to hold the reult
	double *resultConstrained = new double[numConstrainedDOFs];
	memset(resultConstrained,0,sizeof(double)*numConstrainedDOFs);

	if(USECONJUGATEGRADIENTSOLVER) {
		 PardisoSolver solver(tempSpMatCopy,7,0,0,0);
		 solver.ComputeCholeskyDecomposition(tempSpMatCopy);
		 int retVal = solver.SolveLinearSystem(resultConstrained,rhsConstrained);
	}
	else {
		//Solve this using the Gauss Seidel Iterations
		int numGaussSeidelIterations = 800;
		for(int iter=0; iter<numGaussSeidelIterations; iter++)
			tempSpMatCopy->DoOneGaussSeidelIteration(resultConstrained, rhsConstrained);
		tempSpMatCopy->CheckLinearSystemSolution(resultConstrained,rhsConstrained); 
	}

	//Expand the rows now to fill it 
	double *delV = new double[num_vertices*3];
	InsertRows(num_vertices*3,resultConstrained,delV,numConstrainedVerts_*3,constrainedVerts_);
	//Free the stuff we allocated
	delete(tempSpMatCopy);
	delete[] rhsConstrained;
	delete[] resultConstrained;

	//For now just leapfrog this one 
	OMP_FOR
	for(int p=0;p<num_vertices;p++)
	{
		Eigen::Vector3d old_pos = cloth->getMesh()->get_point_data(p,lastFrameId_);
		Eigen::Vector3d old_vel = cloth->get_velocity(p);

		Eigen::Vector3d elemDelV(delV[3*p+0],delV[3*p+1],delV[3*p+2]);
		//DBG__Vector_Dump(elemDelV,p,"DElv");
		Eigen::Vector3d new_vel = old_vel +  elemDelV;
		Eigen::Vector3d new_pos = old_pos + (time_step * new_vel);
		cloth->set_next_step_pos(new_pos,p);
		cloth->set_next_step_velocity(new_vel,p);
	}


	delete[] delV;
}

void ImplicitMassSpringSolver::resetParameters() {
	RHSVector_.setZero();
	LHSMatrix_->ResetToZero();
}

