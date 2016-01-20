#include "ImplicitHyperElasticFEMSolver.h"

#ifndef ELT
#define ELT(numRows,i,j) (((long)j)*((long)numRows)+((long)i))
#endif

ImplicitHyperElasticFEMSolver::ImplicitHyperElasticFEMSolver(void)
{
	dDSdU.resize(54,0);

	dDSdU[compute6x9TensorIndex(0, 0, 0, 0)] = -1.0;
	dDSdU[compute6x9TensorIndex(1, 0, 0, 1)] = -1.0;
	dDSdU[compute6x9TensorIndex(2, 0, 0, 2)] = -1.0;
	dDSdU[compute6x9TensorIndex(0, 1, 1, 0)] = -1.0;
	dDSdU[compute6x9TensorIndex(1, 1, 1, 1)] = -1.0;
	dDSdU[compute6x9TensorIndex(2, 1, 1, 2)] = -1.0;

	dDSdU[compute6x9TensorIndex(0, 0, 2, 0)] = 1.0;
	dDSdU[compute6x9TensorIndex(0, 1, 2, 0)] = 1.0;
	dDSdU[compute6x9TensorIndex(1, 0, 2, 1)] = 1.0;
	dDSdU[compute6x9TensorIndex(1, 1, 2, 1)] = 1.0;
	dDSdU[compute6x9TensorIndex(2, 0, 2, 2)] = 1.0;
	dDSdU[compute6x9TensorIndex(2, 1, 2, 2)] = 1.0;

	error_ = 0;
	errorQuotient_ = 0;
	numIters_ = 0;
}

void ImplicitHyperElasticFEMSolver::initializeSparseMatrixFromOutline(Cloth_Data* cloth)
{
	int numVertices = cloth->getMesh()->get_number_vertices();
	int numElementVertices = 3;
	std::vector<int> vertices(numElementVertices);

	// build the non-zero locations of the tangent stiffness matrix
	SparseMatrixOutline * emptyMatrix = new SparseMatrixOutline(3 * numVertices);
	int numElements = cloth->getMesh()->get_number_triangles();

	for (int el = 0; el < numElements; el++) {

		Triangles tri = cloth->getMesh()->get_triangle(el);
		vertices[0] = tri.a;
		vertices[1] = tri.b;
		vertices[2] = tri.c;

		for (int i = 0; i < numElementVertices; i++)
			for (int j = 0; j < numElementVertices; j++) {
				for (int k = 0; k < 3; k++)
					for (int l = 0; l < 3; l++) {
						// only add the entry if both vertices are free (non-fixed)
						// the corresponding elt is in row 3*i+k, column 3*j+l
						emptyMatrix->AddEntry( 3 * vertices[i] + k, 3 * vertices[j] + l, 0.0 );
					}
			}
	}

#ifdef ADDBENDINGFORCES
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
					emptyMatrix->AddBlock3x3Entry(idx[i],idx[j],fakeIdentity);
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

#endif
	tangentStiffnessMatrix_ = new SparseMatrix(emptyMatrix);
	delete(emptyMatrix);

	rayleighDampingMatrix_ = new SparseMatrix(*tangentStiffnessMatrix_);
}

void ImplicitHyperElasticFEMSolver::initializeMassMatrixFromOutline(Cloth_Data* cloth)
{
	int numVertices = cloth->getMesh()->get_number_vertices();
	SparseMatrixOutline * emptyMatrix = new SparseMatrixOutline(3 * numVertices);
	for(int p=0;p<numVertices;p++)
	{
		float vertexMass = cloth->get_vertex_mass(p);
		for(int d=0;d<3;d++)
		{
			emptyMatrix->AddEntry(3*p+d,3*p+d,vertexMass);
		}
	}
	massMatrix_ = new SparseMatrix(emptyMatrix);
	delete emptyMatrix;
}

//Some notes that are of importance here 
// The structure of the matrix u which denotes displacement (3x3) is as follows: 
//   | uax uay uaz |        | u00 u01 u02 |
//   | ubx uby ubz | =====  | u10 u11 u12 | 
//   | ucx ucy ucz |        | u20 u21 u22 |
// The strucure of the matrix F which denotes the deformation gradient (3x2) is as follows :
//   | F00 F01 |
//   | F10 F11 |
//   | F20 F21 |
// The strucure of the matrix dF\dU which is a 3x2x3x3 second order tensor is as follows (here it is reprented as a ROW-MAJOR 6x9 matrix)
//   | dF_00\du_00 dF_00\du_01 dF_00\du_02 dF_00\du_10 dF_00\du_11 dF_00\du_12 dF_00\du_20 dF_00\du_21 dF_00\du_22 |
//   | dF_01\du_00 dF_01\du_01 dF_01\du_02 dF_01\du_10 dF_01\du_11 dF_01\du_12 dF_01\du_20 dF_01\du_21 dF_01\du_22 |
//   | dF_10\du_00 dF_10\du_01 dF_10\du_02 dF_10\du_10 dF_10\du_11 dF_10\du_12 dF_10\du_20 dF_10\du_21 dF_10\du_22 |
//   | dF_11\du_00 dF_11\du_01 dF_11\du_02 dF_11\du_10 dF_11\du_11 dF_11\du_12 dF_11\du_20 dF_11\du_21 dF_11\du_22 |
//   | dF_20\du_00 dF_20\du_01 dF_20\du_02 dF_20\du_10 dF_20\du_11 dF_20\du_12 dF_20\du_20 dF_20\du_21 dF_20\du_22 |
//   | dF_21\du_00 dF_21\du_01 dF_21\du_02 dF_21\du_10 dF_21\du_11 dF_21\du_12 dF_21\du_20 dF_21\du_21 dF_21\du_22 |
//  To calculate the index in this matrix use the function 
//  compute6x9tensor(i,j,m,n) where (i,j) denotes the F index (ROW INDEX) and (m,n) denotes the u index (COLUMN INDEX) 
void ImplicitHyperElasticFEMSolver::computedFdU(Cloth_Data* cloth) {
	int numElements = cloth->getMesh()->get_number_triangles();
	for (int el = 0; el < numElements; el++) {
		double * dFdU = &dFdUs[54 * el];
		Eigen::Matrix2d dmInv = cloth->getDmInv(el);
		for (int index = 0; index < 54; index++) {
			int n = index % 3;
			int m = (int)(index / 3) % 3;
			int j = (int)(index / 9) % 2;
			int i = (int)(index / 18) % 3;
			double result = 0.0;
			for (int k = 0; k < 2; k++)
				result += dDSdU[compute6x9TensorIndex(i, k, m, n)] * dmInv(k,j);
			dFdU[compute6x9TensorIndex(i, j, m, n)] = result;
		}
	}
}

ImplicitHyperElasticFEMSolver::~ImplicitHyperElasticFEMSolver(void)
{
	delete tangentStiffnessMatrix_;
	delete massMatrix_;
	delete rayleighDampingMatrix_;
}

void ImplicitHyperElasticFEMSolver::initialize(Cloth_Data* cloth)
{
	int numVertices = cloth->getMesh()->get_number_vertices();
	int numTriangles = cloth->getMesh()->get_number_triangles();

	massMatrixDiagonal_ = new double[3*numVertices];
	//OMP_FOR
	for(int p=0; p<numVertices;p++) {
		float vertex_mass = cloth->get_vertex_mass(p);
		massMatrixDiagonal_[3*p+0] = vertex_mass;
		massMatrixDiagonal_[3*p+1] = vertex_mass;
		massMatrixDiagonal_[3*p+2] = vertex_mass;
	}


	//Constrained verts
	numConstrainedVerts_ = 2;
	constrainedVerts_ = new int[6];
	constrainedVerts_[3] = 3*CORNERV2+0;
	constrainedVerts_[4] = 3*CORNERV2+1;
	constrainedVerts_[5] = 3*CORNERV2+2;
	constrainedVerts_[0] = 3*CORNERV1+0;
	constrainedVerts_[1] = 3*CORNERV1+1;
	constrainedVerts_[2] = 3*CORNERV1+2;

	//Add the complemented stuff 
	dFdUs.resize(54*numTriangles,0);
	computedFdU(cloth);

	//Compute the tangent Stiffness Matrix Structure 
	initializeSparseMatrixFromOutline(cloth);

	//Compute the mass matrix structure 
	initializeMassMatrixFromOutline(cloth);

	//Establish substructures
	rayleighDampingMatrix_->BuildSubMatrixIndices(*massMatrix_);
	tangentStiffnessMatrix_->BuildSubMatrixIndices(*massMatrix_);
#ifdef ADDBENDINGFORCES
	tangentStiffnessMatrix_->BuildSubMatrixIndices(*mKeQ_,1);
	tangentStiffnessMatrix_->BuildSubMatrixIndices(*mKdQ_,2);
#endif

	//Build the correct internal and external forces vector 
	internalForce_ =  Eigen::VectorXd::Zero(numVertices*3);
	externalForce_ =  Eigen::VectorXd::Zero(numVertices*3);
	residuals_ = std::vector<double>(numVertices*3,0.0);
	velocities_ = Eigen::VectorXd::Zero(numVertices*3);

	q_  =  Eigen::VectorXd::Zero(numVertices*3);
	q1_  =  Eigen::VectorXd::Zero(numVertices*3);
	qVel_  =  Eigen::VectorXd::Zero(numVertices*3);
	qVel1_  =  Eigen::VectorXd::Zero(numVertices*3);

#ifdef IMPLICITNEWMARKSOLVER
	qAccel_.resize(numVertices*3,0);
	qAccel1_.resize(numVertices*3,0);
#endif

	oldPosition_.resize(numVertices*3,0.0);

	int numConstrainedDOFs = (numVertices*3) - (3*numConstrainedVerts_);
	constrainedResiduals_ = std::vector<double>(numConstrainedDOFs);

#ifdef IMPLICITNEWMARKSOLVER
	initializeAlphas(cloth);
#endif

}

void ImplicitHyperElasticFEMSolver::loadQAndQVels(Cloth_Data* cloth)
{

	int numVertices = cloth->getMesh()->get_number_vertices();
	//OMP_FOR
	for(int i=0;i<numVertices;i++)
	{
		Eigen::Vector3d v = cloth->get_velocity(i);
		qVel_[3*i+0] = v[0];
		qVel_[3*i+1] = v[1];
		qVel_[3*i+2] = v[2];

		Eigen::Vector3d pv0 = cloth->getMesh()->get_point_data(i,0);
		Eigen::Vector3d pvn = cloth->getMesh()->get_point_data(i,lastFrameId_);

		oldPosition_[3*i+0] = pvn.x();
		oldPosition_[3*i+1] = pvn.y();
		oldPosition_[3*i+2] = pvn.z();

		q_[3*i+0] = pvn[0]-pv0[0];
		q_[3*i+1] = pvn[1]-pv0[1];
		q_[3*i+2] = pvn[2]-pv0[2];

	}
}

#ifdef IMPLICITNEWMARKSOLVER
void ImplicitHyperElasticFEMSolver::computeInitialAccelerations(Cloth_Data* cloth)
{

	int numVertices = cloth->getMesh()->get_number_vertices();
	int n3 = numVertices*3;
	int numConstrainedDOFs = n3 - (3*numConstrainedVerts_);

	tangentStiffnessMatrix_->ResetToZero();
	rayleighDampingMatrix_->ResetToZero();

	internalForce_.setZero();
	externalForce_.setZero();

	//Populate the matrices
	addShearComponents(cloth);
#ifdef ADDBENDINGFORCES
	addQuadraticBending(cloth);
#endif

	double dampingStiffnessCoeff = cloth->get_property_obj()->getRayleighDampingCoeff();
	double dampingMassCoeff = cloth->get_property_obj()->getRayleighMassCoeff();

	//tangentStiffnessMatrix_->AddSubMatrix(1.0,*mKeQ_,1);

	*rayleighDampingMatrix_ = dampingStiffnessCoeff * (*tangentStiffnessMatrix_);
  rayleighDampingMatrix_->AddSubMatrix(dampingMassCoeff, *massMatrix_);

	std::vector<double> buffer(n3,0);
	rayleighDampingMatrix_->MultiplyVector(qVel_.data(), &buffer[0]);

	for(int i=0; i<numVertices*3; i++)
    buffer[i] = -buffer[i] - internalForce_[i];
	double* bufferConstrained = new double[numConstrainedDOFs];
	RemoveRows(n3,bufferConstrained,&buffer[0],numConstrainedVerts_*3,constrainedVerts_);
	
	//Use tangent stiffness matrix also 
	 // use tangentStiffnessMatrix as the buffer place
  tangentStiffnessMatrix_->ResetToZero();
  tangentStiffnessMatrix_->AddSubMatrix(1.0, *massMatrix_);

	//Make a copy of the LHSMatrix_ and remove rows and columns from it
	SparseMatrix *tempSpMatCopy = new SparseMatrix(*tangentStiffnessMatrix_);
	tempSpMatCopy->RemoveRowsColumns(numConstrainedVerts_*3,constrainedVerts_);
	
	//Make a correct size row vector to hold the reult
	double *resultConstrained = new double[numConstrainedDOFs];
	memset(resultConstrained,0,sizeof(double)*numConstrainedDOFs);

#if defined USE_PARDISO_SOLVER
	PardisoSolver solver(tempSpMatCopy,7,0,0,0);
	solver.ComputeCholeskyDecomposition(tempSpMatCopy);
	int retVal = solver.SolveLinearSystem(resultConstrained,bufferConstrained);
#elif defined USE_CG_SOLVER
	CGSolver cgsolver(tempSpMatCopy);
	cgsolver.SolveLinearSystemWithJacobiPreconditioner(resultConstrained, bufferConstrained, 1e-6, 10000);
#else
	int numGaussSeidelIterations = 8000;
	for(int iter=0; iter<numGaussSeidelIterations; iter++)
		tempSpMatCopy->DoOneGaussSeidelIteration(resultConstrained, bufferConstrained);
#endif
	//Expand the rows now to fill it 
	InsertRows(n3,resultConstrained,&qAccel_[0],numConstrainedVerts_*3,constrainedVerts_);

	delete[] resultConstrained;

}
#endif

/*void ImplicitHyperElasticFEMSolver::advance_time_step(Cloth_Data* cloth)
{
	//Set the last frame information
	lastFrameId_ = (cloth->getMesh()->get_number_frames()) - 1;

	//Preapre the components as usual 
	tangentStiffnessMatrix_->ResetToZero();
	rayleighDampingMatrix_->ResetToZero();

	//Prepare the RHS Vector for usage
	force_.setZero();
	internalForce_.setZero();
	externalForce_.setZero();
	residuals_.setZero();
	velocities_.setZero();
	

	//Prepare the parameters
	error_ = 0;
	errorQuotient_ = 0;
	numIters_ = 0;

	//Add the physics components
	addShearComponents(cloth);
	addGravityComponents(cloth);

	//Solve and report
	//finalizeAndSolve(cloth);
	finalizeAndSolve_Explicit(cloth);
}*/

void ImplicitHyperElasticFEMSolver::updateQAndQ1s(Cloth_Data* cloth)
{
	int numVertices = cloth->getMesh()->get_number_vertices();
	
	for(int i=0;i<numVertices*3;i++)
	{
		q1_[i] = q_[i];
		qVel1_[i] = qVel_[i];

#ifdef IMPLICITNEWMARKSOLVER
		qAccel1_[i] = qAccel_[i];

		qAccel_[i] = alpha1_ * (q_[i] - q1_[i]) - alpha2_ * qVel1_[i] - alpha3_ * qAccel1_[i];
    qVel_[i] = alpha4_ * (q_[i] - q1_[i]) + alpha5_ * qVel1_[i] + alpha6_ * qAccel1_[i];
#endif

	}
}

void ImplicitHyperElasticFEMSolver::advance_time_step(Cloth_Data* cloth)
{
	//Set the last frame information
	lastFrameId_ = (cloth->getMesh()->get_number_frames()) - 1;

	//Prepare the parameters
	error_ = 0;
	errorQuotient_ = 0;
	numIters_ = 0;

	//Initialize q and q vels
	loadQAndQVels(cloth);

#ifdef IMPLICITNEWMARKSOLVER
	computeInitialAccelerations(cloth);
#endif
	updateQAndQ1s(cloth);

	externalForce_.setZero();
	int lastVClicked = cloth->getLastClickedVertex();
	Eigen::Vector3d uiForce = 1.0*(cloth->getCurrentUIForce());
	if(lastVClicked!=-1) {
		externalForce_[3*lastVClicked+0] += uiForce[0];
		externalForce_[3*lastVClicked+1] += uiForce[1];
		externalForce_[3*lastVClicked+2] += uiForce[2];
	}
	addGravityComponents(cloth);


	do
	{

		//Preapre the components as usual 
		tangentStiffnessMatrix_->ResetToZero();
		rayleighDampingMatrix_->ResetToZero();

		//Prepare the RHS Vector for usage
		internalForce_.setZero();
		std::fill(residuals_.begin(),residuals_.end(),0.0);

		//Add the physics components
		addShearComponents(cloth);

#ifdef ADDBENDINGFORCES
		addQuadraticBending(cloth);
#endif

		//Now yoiu prepare the numerical system for solution
		prepareNumericalEquations(cloth);



		//Check for convergence
		bool converged = checkConvergence(cloth);

		if(converged)
		{
			//Send this as final values to the solver 
			finalizeStep(cloth);
			std::cout << "Num Iters: " << numIters_ << "\n";
			break;
		}
		else
		{
			// Perform the solver step 
			solverStep(cloth);
		}

		numIters_++;

	}while(numIters_ < 10000);

	if(numIters_ == 10000)
	{
		std::cout << "[ERROR] Solver has failed to converge in 1000 iterations. Reduce timestep.";
	}

}

int ImplicitHyperElasticFEMSolver::compute6x9TensorIndex(int i, int j, int m, int n) {
	int rowIndex = 2 * i + j;
	int columnIndex = 3 * m + n;
	if(9*rowIndex+columnIndex>=54)
	{
		std::cout << i << " " << j << " " << m << " " << n << " " << (9*rowIndex+columnIndex) << "\n";
		return 0;
	}
	return (9 * rowIndex + columnIndex);
}

/*The Matrix Order is as follows:
      | P_00 P_01 |        | F_00 F_01 |
if P= | P_10 P_11 | and F= | F_10 F_11 |
      | P_20 P_21 |        | F_20 F_21 |

| dP_00/dF_00  dP_00/dF_01 dP_00/dF_10 dP_00/dF_11 dP_00/dF_20 dP_00/dF_21 |
| dP_01/dF_00  dP_01/dF_01 dP_01/dF_10 dP_01/dF_11 dP_01/dF_20 dP_01/dF_21 |
| dP_10/dF_00  dP_10/dF_01 dP_10/dF_10 dP_10/dF_11 dP_10/dF_20 dP_10/dF_21 |
| dP_11/dF_00  dP_11/dF_01 dP_11/dF_10 dP_11/dF_11 dP_11/dF_20 dP_11/dF_21 |
| dP_20/dF_00  dP_20/dF_01 dP_20/dF_10 dP_20/dF_11 dP_20/dF_20 dP_20/dF_21 |
| dP_21/dF_00  dP_21/dF_01 dP_21/dF_10 dP_21/dF_11 dP_21/dF_20 dP_21/dF_21 |
*/
//DPDF is assumed to be row-major
int ImplicitHyperElasticFEMSolver::compute6x6TensorIndex(int i, int j, int m, int n)
{
	int rowIndex = i * 2 + j;
	int colIndex = m * 2 + n;
#ifdef ROWMAJOR6x6TENSOR
	return (6*rowIndex+colIndex);
#else
	return (6*colIndex+rowIndex);
#endif
}

/*
void ImplicitHyperElasticFEMSolver::computeDPDFFrustrating_OLD(Eigen::MatrixXd &F,std::vector<double> &gradients, std::vector<double> &hessians,std::vector<double> &DPDF)
{
//First some convinienet shorthands because whynot 
double x1 = F(0,0); double x2 = F(0,1);
double y1 = F(1,0); double y2 = F(1,1);
double z1 = F(2,0); double z2 = F(2,1);

double x1_2 = x1*x1; double x2_2 = x2*x2;
double y1_2 = y1*y1; double y2_2 = y2*y2;
double z1_2 = z1*z1; double z2_2 = z2*z2;

DPDF.resize(36,0);

//A total of 36 numbers 
DPDF[compute6x6TensorIndex(0,0,0,0)] = (2.0*y2_2 + 2.0*z2_2)*gradients[2] + (8.0*x1_2 + 4.0*x2_2 + 4.0*(x1_2 + y1_2 + z1_2))*gradients[1] + 2.0*gradients[0] 
+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 -2.0*x2*z1*z2 + 2.0*x1*z2_2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
+ (2.0*x1)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

DPDF[compute6x6TensorIndex(0,1,0,0)] = (-2.0*y1*y2 -2.0*z1*z2)*gradients[2] + (4.0*x1*x2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
+ (2.0*x1)* ((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

DPDF[compute6x6TensorIndex(1,0,0,0)] = -2.0*x2*y2*gradients[2] + (8.0*x1*y1 + 4.0*x2*y2)*gradients[1] + 
+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
+ (2.0*x1)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

DPDF[compute6x6TensorIndex(1,1,0,0)] = (-2.0*x2*y1 + 4.0*x1*y2)*gradients[2] + (4.0*x2*y1)*gradients[1] 
+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
+ (2.0*x1)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

DPDF[compute6x6TensorIndex(2,0,0,0)] = -2.0*x2*z2*gradients[2] + (8.0*x1*z1 + 4.0*x2*z2)*gradients[1] 
+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
+ (2.0*x1)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

DPDF[compute6x6TensorIndex(2,1,0,0)] = (-2.0*x2*z1 + 4.0*x1*z2)*gradients[2] + (4.0*x2*z1)*gradients[1] 
+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 -2.0*x2*z1*z2 + 2.0*x1*z2_2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] +2.0*z2*hessians[2])
+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] +2.0*z2*hessians[1])
+ (2.0*x1)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] +2.0*z2*hessians[0]);

DPDF[compute6x6TensorIndex(0,0,1,0)] = -2.0*x2*y2*gradients[2] + (8.0*x1*y1 + 4.0*x2*y2)*gradients[1]
+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
+ (2*y1)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

DPDF[compute6x6TensorIndex(0,1,1,0)] = (4.0*x2*y1 -2.0*x1*y2)*gradients[2] + (4.0*x1*y2)*gradients[1]
+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
+ (2.0*y1)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

DPDF[compute6x6TensorIndex(1,0,1,0)] = (2.0*x2_2 + 2.0*z2_2)*gradients[2] + (8.0*y1_2 + 4.0*y2_2 + 4.0*(x1_2 + y1_2 + z1_2))*gradients[1] + 2.0*gradients[0]
+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
+ (2.0*y1)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

DPDF[compute6x6TensorIndex(1,1,1,0)] = (-2.0*x1*x2 - 2.0*z1*z2)*gradients[2] + (4.0*y1*y2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1]
+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
+ (2.0*y1)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

DPDF[compute6x6TensorIndex(2,0,1,0)] = (-2.0*y2*z2)*gradients[2] + (8.0*y1*z1 + 4.0*y2*z2)*gradients[1]
+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] +2.0*z1*hessians[2])
+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] +2.0*z1*hessians[1])
+ (2.0*y1)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] +2.0*z1*hessians[0]);

DPDF[compute6x6TensorIndex(2,1,1,0)] = (-2.0*y2*z1 + 4.0*y1*z2)*gradients[2] + (4.0*y2*z1)*gradients[1]
+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
+ (2.0*y1)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

DPDF[compute6x6TensorIndex(0,0,2,0)] = (-2.0*x2*z2)*gradients[2] + (8.0*x1*z1 + 4.0*x2*z2)*gradients[1]
+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
+ (2.0*z1)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

DPDF[compute6x6TensorIndex(0,1,2,0)] = (4.0*x2*z1 - 2.0*x1*z2)*gradients[2] + (4.0*x1*z2)*gradients[1] 
+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
+ (2.0*z1)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

DPDF[compute6x6TensorIndex(1,0,2,0)] = (-2.0*y2*z2)*gradients[2] + (8.0*y1*z1 + 4.0*y2*z2)*gradients[1]
+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
+ (2.0*z1)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

DPDF[compute6x6TensorIndex(1,1,2,0)] = (4.0*y2*z1 - 2.0*y1*z2)*gradients[2] + (4.0*y1*z2)*gradients[1] 
+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
+ (2.0*z1)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

DPDF[compute6x6TensorIndex(2,0,2,0)] = (2.0*x2_2 + 2.0*y2_2)*gradients[2] + (8.0*z1_2 + 4.0*(x1_2 + y1_2 + z1_2) + 4.0*z2_2)*gradients[1] + 2.0*gradients[0]
+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
+ (2.0*z1)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

DPDF[compute6x6TensorIndex(2,1,2,0)] = (-2.0*x1*x2 - 2.0*y1*y2)*gradients[2] + (4.0*z1*z2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
+ (2.0*z1)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

DPDF[compute6x6TensorIndex(0,0,0,1)] = (-2.0*y1*y2 - 2.0*z1*z2)*gradients[2] + (4.0*x1*x2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
+ (2.0*x2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

DPDF[compute6x6TensorIndex(0,1,0,1)] = (2.0*y1_2 + 2.0*z1_2)*gradients[2] + (4.0*x1_2 + 8.0*x2_2 + 4.0*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*gradients[0]
+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
+ (2.0*x2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

DPDF[compute6x6TensorIndex(1,0,0,1)] = (4.0*x2*y1 - 2.0*x1*y2)*gradients[2] + (4.0*x1*y2)*gradients[1]
+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
+ (2.0*x2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

DPDF[compute6x6TensorIndex(1,1,0,1)] = (-2.0*x1*y1)*gradients[2] + (4.0*x1*y1 + 8.0*x2*y2)*gradients[1]
+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
+ (2.0*x2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

DPDF[compute6x6TensorIndex(2,0,0,1)] = (4.0*x2*z1 - 2.0*x1*z2)*gradients[2] + (4.0*x1*z2)*gradients[1] 
+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] +  2.0*z1*hessians[2])
+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] +  2.0*z1*hessians[1])
+ (2.0*x2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] +  2.0*z1*hessians[0]);

DPDF[compute6x6TensorIndex(2,1,0,1)] = (-2.0*x1*z1)*gradients[2] + (4.0*x1*z1 + 8.0*x2*z2)*gradients[1]
+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
+ (2.0*x2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

DPDF[compute6x6TensorIndex(0,0,1,1)] = (-2.0*x2*y1 + 4.0*x1*y2)*gradients[2] + (4.0*x2*y1)*gradients[1] 
+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1* z2))*hessians[4] + 2.0*x1*hessians[2])
+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1* z2))*hessians[3] + 2.0*x1*hessians[1])
+ (2.0*y2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1* z2))*hessians[1] + 2.0*x1*hessians[0]);

DPDF[compute6x6TensorIndex(0,1,1,1)] = (-2.0*x1*y1)*gradients[2] + (4.0*x1*y1 + 8.0*x2*y2)*gradients[1] 
+  (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
+ (2.0*y2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

DPDF[compute6x6TensorIndex(1,0,1,1)] = (-2.0*x1*x2 - 2.0*z1*z2)*gradients[2] + (4.0*y1*y2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
+ (2.0*y2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

DPDF[compute6x6TensorIndex(1,1,1,1)] = (2.0*x1_2 + 2.0*z1_2)*gradients[2] + (4.0*y1_2 + 8.0*y2_2 + 4.0*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*gradients[0]
+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
+ (2.0*y2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

DPDF[compute6x6TensorIndex(2,0,1,1)] = (4.0*y2*z1 - 2.0*y1*z2)*gradients[2] + (4.0*y1*z2)*gradients[1] + 
+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
+ (2.0*y2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

DPDF[compute6x6TensorIndex(2,1,1,1)] = (-2.0*y1*z1)*gradients[2] + (4.0*y1*z1 + 8.0*y2*z2)*gradients[1]
+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
+ (2.0*y2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

DPDF[compute6x6TensorIndex(0,0,2,1)] = (-2.0*x2*z1 + 4.0*x1*z2)*gradients[2] + (4.0*x2*z1)*gradients[1]
+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
+ (2.0*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

DPDF[compute6x6TensorIndex(0,1,2,1)] = (-2.0*x1*z1)*gradients[2] + (4.0*x1*z1 + 8.0*x2*z2)*gradients[1]
+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
+ (2.0*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

DPDF[compute6x6TensorIndex(1,0,2,1)] = (-2.0*y2*z1 + 4.0*y1*z2)*gradients[2] + (4.0*y2*z1)*gradients[1]
+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
+ (2.0*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

DPDF[compute6x6TensorIndex(1,1,2,1)] = (-2.0*y1*z1)*gradients[2] + (4.0*y1*z1 + 8.0*y2*z2)*gradients[1]
+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
+ (2.0*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

DPDF[compute6x6TensorIndex(2,0,2,1)] = (-2.0*x1*x2 - 2.0*y1*y2)*gradients[2] + (4.0*z1*z2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
+ (2.0*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

DPDF[compute6x6TensorIndex(2,1,2,1)] = (2.0*x1_2 + 2.0*y1_2)*gradients[2] + (4.0*z1_2 + 8.0*z2_2 + 4.0*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*gradients[0]
+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
+ (2.0*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

//DBG__Matrix_Dump(&DPDF[0],6,6);
}
*/

void ImplicitHyperElasticFEMSolver::computeDPDF(Eigen::MatrixXd &F,std::vector<double> &gradients, std::vector<double> &hessians,std::vector<double> &DPDF)
{
	//First some convinienet shorthands because whynot 
	double x1 = F(0,0); double x2 = F(0,1);
	double y1 = F(1,0); double y2 = F(1,1);
	double z1 = F(2,0); double z2 = F(2,1);

	double x1_2 = x1*x1; double x2_2 = x2*x2;
	double y1_2 = y1*y1; double y2_2 = y2*y2;
	double z1_2 = z1*z1; double z2_2 = z2*z2;

	DPDF.resize(36,0);

	//A total of 36 numbers 
	DPDF[compute6x6TensorIndex(0,0,0,0)] = (2.0*y2_2 + 2.0*z2_2)*gradients[2] + (8.0*x1_2 + 4.0*x2_2 + 4.0*(x1_2 + y1_2 + z1_2))*gradients[1] + 2.0*gradients[0] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 -2.0*x2*z1*z2 + 2.0*x1*z2_2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
		+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
		+ (2.0*x1)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,0,1)] = (-2.0*y1*y2 -2.0*z1*z2)*gradients[2] + (4.0*x1*x2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
		+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
		+ (2.0*x1)* ((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,1,0)] = -2.0*x2*y2*gradients[2] + (8.0*x1*y1 + 4.0*x2*y2)*gradients[1] + 
		+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
		+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
		+ (2.0*x1)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,1,1)] = (-2.0*x2*y1 + 4.0*x1*y2)*gradients[2] + (4.0*x2*y1)*gradients[1] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
		+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
		+ (2.0*x1)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,2,0)] = -2.0*x2*z2*gradients[2] + (8.0*x1*z1 + 4.0*x2*z2)*gradients[1] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
		+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
		+ (2.0*x1)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,2,1)] = (-2.0*x2*z1 + 4.0*x1*z2)*gradients[2] + (4.0*x2*z1)*gradients[1] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 -2.0*x2*z1*z2 + 2.0*x1*z2_2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] +2.0*z2*hessians[2])
		+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] +2.0*z2*hessians[1])
		+ (2.0*x1)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] +2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,0,0)] = -2.0*x2*y2*gradients[2] + (8.0*x1*y1 + 4.0*x2*y2)*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
		+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
		+ (2*y1)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,0,1)] = (4.0*x2*y1 -2.0*x1*y2)*gradients[2] + (4.0*x1*y2)*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
		+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
		+ (2.0*y1)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,1,0)] = (2.0*x2_2 + 2.0*z2_2)*gradients[2] + (8.0*y1_2 + 4.0*y2_2 + 4.0*(x1_2 + y1_2 + z1_2))*gradients[1] + 2.0*gradients[0]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
		+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
		+ (2.0*y1)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,1,1)] = (-2.0*x1*x2 - 2.0*z1*z2)*gradients[2] + (4.0*y1*y2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
		+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
		+ (2.0*y1)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,2,0)] = (-2.0*y2*z2)*gradients[2] + (8.0*y1*z1 + 4.0*y2*z2)*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] +2.0*z1*hessians[2])
		+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] +2.0*z1*hessians[1])
		+ (2.0*y1)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] +2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,2,1)] = (-2.0*y2*z1 + 4.0*y1*z2)*gradients[2] + (4.0*y2*z1)*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
		+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
		+ (2.0*y1)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,0,0)] = (-2.0*x2*z2)*gradients[2] + (8.0*x1*z1 + 4.0*x2*z2)*gradients[1]
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
		+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
		+ (2.0*z1)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,0,1)] = (4.0*x2*z1 - 2.0*x1*z2)*gradients[2] + (4.0*x1*z2)*gradients[1] 
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
		+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
		+ (2.0*z1)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,1,0)] = (-2.0*y2*z2)*gradients[2] + (8.0*y1*z1 + 4.0*y2*z2)*gradients[1]
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
		+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
		+ (2.0*z1)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,1,1)] = (4.0*y2*z1 - 2.0*y1*z2)*gradients[2] + (4.0*y1*z2)*gradients[1] 
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
		+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
		+ (2.0*z1)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,2,0)] = (2.0*x2_2 + 2.0*y2_2)*gradients[2] + (8.0*z1_2 + 4.0*(x1_2 + y1_2 + z1_2) + 4.0*z2_2)*gradients[1] + 2.0*gradients[0]
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
		+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
		+ (2.0*z1)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,2,1)] = (-2.0*x1*x2 - 2.0*y1*y2)*gradients[2] + (4.0*z1*z2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
		+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
		+ (2.0*z1)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,0,0)] = (-2.0*y1*y2 - 2.0*z1*z2)*gradients[2] + (4.0*x1*x2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
		+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
		+ (2.0*x2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,0,1)] = (2.0*y1_2 + 2.0*z1_2)*gradients[2] + (4.0*x1_2 + 8.0*x2_2 + 4.0*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*gradients[0]
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
		+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
		+ (2.0*x2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,1,0)] = (4.0*x2*y1 - 2.0*x1*y2)*gradients[2] + (4.0*x1*y2)*gradients[1]
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
		+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
		+ (2.0*x2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,1,1)] = (-2.0*x1*y1)*gradients[2] + (4.0*x1*y1 + 8.0*x2*y2)*gradients[1]
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
		+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
		+ (2.0*x2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,2,0)] = (4.0*x2*z1 - 2.0*x1*z2)*gradients[2] + (4.0*x1*z2)*gradients[1] 
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] +  2.0*z1*hessians[2])
		+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] +  2.0*z1*hessians[1])
		+ (2.0*x2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] +  2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,2,1)] = (-2.0*x1*z1)*gradients[2] + (4.0*x1*z1 + 8.0*x2*z2)*gradients[1]
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
		+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
		+ (2.0*x2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,0,0)] = (-2.0*x2*y1 + 4.0*x1*y2)*gradients[2] + (4.0*x2*y1)*gradients[1] 
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1* z2))*hessians[4] + 2.0*x1*hessians[2])
		+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1* z2))*hessians[3] + 2.0*x1*hessians[1])
		+ (2.0*y2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1* z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,0,1)] = (-2.0*x1*y1)*gradients[2] + (4.0*x1*y1 + 8.0*x2*y2)*gradients[1] 
	+  (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
		+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
		+ (2.0*y2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,1,0)] = (-2.0*x1*x2 - 2.0*z1*z2)*gradients[2] + (4.0*y1*y2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
		+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
		+ (2.0*y2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,1,1)] = (2.0*x1_2 + 2.0*z1_2)*gradients[2] + (4.0*y1_2 + 8.0*y2_2 + 4.0*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*gradients[0]
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
		+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
		+ (2.0*y2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,2,0)] = (4.0*y2*z1 - 2.0*y1*z2)*gradients[2] + (4.0*y1*z2)*gradients[1] + 
		+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
		+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
		+ (2.0*y2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,2,1)] = (-2.0*y1*z1)*gradients[2] + (4.0*y1*z1 + 8.0*y2*z2)*gradients[1]
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
		+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
		+ (2.0*y2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,0,0)] = (-2.0*x2*z1 + 4.0*x1*z2)*gradients[2] + (4.0*x2*z1)*gradients[1]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
		+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
		+ (2.0*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,0,1)] = (-2.0*x1*z1)*gradients[2] + (4.0*x1*z1 + 8.0*x2*z2)*gradients[1]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
		+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
		+ (2.0*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,1,0)] = (-2.0*y2*z1 + 4.0*y1*z2)*gradients[2] + (4.0*y2*z1)*gradients[1]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
		+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
		+ (2.0*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,1,1)] = (-2.0*y1*z1)*gradients[2] + (4.0*y1*z1 + 8.0*y2*z2)*gradients[1]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
		+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
		+ (2.0*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,2,0)] = (-2.0*x1*x2 - 2.0*y1*y2)*gradients[2] + (4.0*z1*z2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
		+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
		+ (2.0*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,2,1)] = (2.0*x1_2 + 2.0*y1_2)*gradients[2] + (4.0*z1_2 + 8.0*z2_2 + 4.0*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*gradients[0]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
		+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
		+ (2.0*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	//DBG__Matrix_Dump(&DPDF[0],6,6);
}

void ImplicitHyperElasticFEMSolver::computeFirstPiolaStress(Eigen::MatrixXd &F, std::vector<double> &gradients, Eigen::MatrixXd &P)
{
	double x1 = F(0,0); double x2 = F(0,1);
	double y1 = F(1,0); double y2 = F(1,1);
	double z1 = F(2,0); double z2 = F(2,1);

	double x1_2 = x1*x1; double x2_2 = x2*x2;
	double y1_2 = y1*y1; double y2_2 = y2*y2;
	double z1_2 = z1*z1; double z2_2 = z2*z2;

	P = Eigen::MatrixXd(3,2);

	P(0,0) = (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*gradients[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*gradients[1] +  2.0*x1*gradients[0];
	P(0,1) = (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*gradients[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*x2*gradients[0];

	P(1,0) = (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*gradients[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*gradients[1] + 2.0*y1*gradients[0];
	P(1,1) = (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*gradients[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*y2*gradients[0];

	P(2,0) = (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*gradients[2] +  (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*gradients[1] + 2.0*z1*gradients[0];
	P(2,1) = (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*gradients[2] +  (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*z2*gradients[0];
}


//  This is the routine inspired from VEGA's implementation which does complex tensor manupulations.  
//  DGDF is obtained from DPDF by multiplying individual blocks of the DPDF tensor (each of dimension 3x2) by the bm matrix (2x2) 
//  to produce a (3x2) block. The exact types of block multiplied and the final result produced thereof is as follows:
//  | dP00_dF00 dP01_dF00 |                              | dG00_dF00 dG01_dF00 |
//  | dP10_dF00 dP11_dF00 |    *    | bm00 bm01 |     =  | dG10_dF00 dG11_dF00 |
//  | dP20_dF00 dP21_dF00 |         | bm10 bm11 |        | dG20_dF00 dG21)dF00 |
//  Basically repeat this process for F01, F10, F11, F20, F21 (Each 3x2 blocks will produce the corresponding blocks of dG_dF
//  Finally the final structure of the matrix is as follows (This is a ROW-MAJOR MATRIX) 
//  | dG00_dF00 dG00_dF01 dG00_dF10 dG00_dF11 dG00_dF20 dG00_dF21 |
//  | dG10_dF00 dG10_dF01 dG10_dF10 dG10_dF11 dG10_dF20 dG10_dF21 |
//  | dG20_dF00 dG20_dF01 dG20_dF10 dG20_dF11 dG20_dF20 dG20_dF21 | 
//  | dG01_dF00 dG01_dF01 dG01_dF10 dG01_dF11 dG01_dF20 dG01_dF21 |
//  | dG11_dF00 dG11_dF01 dG11_dF10 dG11_dF11 dG11_dF20 dG11_dF21 |
//  | dG21_dF00 dG21_dF01 dG21_dF10 dG21_dF11 dG21_dF20 dG21_dF21 |
void ImplicitHyperElasticFEMSolver::computeDGDF(std::vector<double> &DPDF, Eigen::Matrix2d bVec, std::vector<double> &DGDF)
{
	DGDF.resize(36,0);
	for (int abc = 0; abc < 2; abc++)
		for (int i = 0; i < 3; i++)
			for (int column = 0; column < 6; column++)
				for (int k = 0; k < 2; k++)
				{ 
					//std::cout << "{" << 18 * abc + 6 * i + column << "} [" << (2 * i + k) * 6 + column << "] ";
					//std::cout << "abc:" << abc << ", k:" << k << " ,i:" << i << " ,col:" << column << "\n";
					DGDF[18 * abc + 6 * i + column] += DPDF[(2 * i + k) * 6 + column] * (bVec(abc,k));
				}
}

//Finally we compute the product dG_dF * dF_dU to get dG_dU 
//  | dG00_dF00 dG00_dF01 dG00_dF10 dG00_dF11 dG00_dF20 dG00_dF21 |
//  | dG10_dF00 dG10_dF01 dG10_dF10 dG10_dF11 dG10_dF20 dG10_dF21 |
//  | dG20_dF00 dG20_dF01 dG20_dF10 dG20_dF11 dG20_dF20 dG20_dF21 | 
//  | dG01_dF00 dG01_dF01 dG01_dF10 dG01_dF11 dG01_dF20 dG01_dF21 |
//  | dG11_dF00 dG11_dF01 dG11_dF10 dG11_dF11 dG11_dF20 dG11_dF21 |
//  | dG21_dF00 dG21_dF01 dG21_dF10 dG21_dF11 dG21_dF20 dG21_dF21 |
//                         *
//   | dF_00\du_00 dF_00\du_01 dF_00\du_02 dF_00\du_10 dF_00\du_11 dF_00\du_12 dF_00\du_20 dF_00\du_21 dF_00\du_22 |
//   | dF_01\du_00 dF_01\du_01 dF_01\du_02 dF_01\du_10 dF_01\du_11 dF_01\du_12 dF_01\du_20 dF_01\du_21 dF_01\du_22 |
//   | dF_10\du_00 dF_10\du_01 dF_10\du_02 dF_10\du_10 dF_10\du_11 dF_10\du_12 dF_10\du_20 dF_10\du_21 dF_10\du_22 |
//   | dF_11\du_00 dF_11\du_01 dF_11\du_02 dF_11\du_10 dF_11\du_11 dF_11\du_12 dF_11\du_20 dF_11\du_21 dF_11\du_22 |
//   | dF_20\du_00 dF_20\du_01 dF_20\du_02 dF_20\du_10 dF_20\du_11 dF_20\du_12 dF_20\du_20 dF_20\du_21 dF_20\du_22 |
//   | dF_21\du_00 dF_21\du_01 dF_21\du_02 dF_21\du_10 dF_21\du_11 dF_21\du_12 dF_21\du_20 dF_21\du_21 dF_21\du_22 |
//                         =
//   | dG_00\du_00 dG_00\du_01 dG_00\du_02 dG_00\du_10 dG_00\du_11 dG_00\du_12 dG_00\du_20 dG_00\du_21 dG_00\du_22 |
//   | dG_10\du_00 dG_10\du_01 dG_10\du_02 dG_10\du_10 dG_10\du_11 dG_10\du_12 dG_10\du_20 dG_10\du_21 dG_10\du_22 |
//   | dG_20\du_00 dG_20\du_01 dG_20\du_02 dG_20\du_10 dG_20\du_11 dG_20\du_12 dG_20\du_20 dG_20\du_21 dG_20\du_22 |
//   | dG_01\du_00 dG_01\du_01 dG_01\du_02 dG_01\du_10 dG_01\du_11 dG_01\du_12 dG_01\du_20 dG_01\du_21 dG_01\du_22 |
//   | dG_11\du_00 dG_11\du_01 dG_11\du_02 dG_11\du_10 dG_11\du_11 dG_11\du_12 dG_11\du_20 dG_11\du_21 dG_11\du_22 |
//   | dG_21\du_00 dG_21\du_01 dG_21\du_02 dG_21\du_10 dG_21\du_11 dG_21\du_12 dG_21\du_20 dG_21\du_21 dG_21\du_22 |
//   Mind you these are the first 6 rows the elements of the last row are obtained as follows
//   dG_02\du_00 = -(dG_00\du_00 + dG_01\du_00) and so on 
void ImplicitHyperElasticFEMSolver::computeElementStiffnessMatrix(int el,std::vector<double> &DGDF,std::vector<double> &KELEM)
{
	double * dFdU = &dFdUs[54 * el];
	//DBG__Matrix_Dump(dFdU,6,9);
	// K is stored column-major (however, it doesn't matter because K is symmetric)
	// Now its time to initalize element K which is a 9x9 matrix  arranged as a double vector 
	KELEM.resize(81,0);
	for (int row = 0; row < 6; row++) {
		for (int column = 0; column < 9; column++) {
			double result = 0;
			for (int inner = 0; inner < 6; inner++) {
				//dGdF is 6x6, and dFdU is 6x9
				result += DGDF[6 * row + inner] * dFdU[9 * inner + column];
			}
			KELEM[9 * column + row] = result;
		}
	}

	for (int row = 0; row < 9; row++) {
		//17th column
		KELEM[9 * row + 6] = -KELEM[9 * row + 0] - KELEM[9 * row + 3] ;
		//8th column
		KELEM[9 * row + 7] = -KELEM[9 * row + 1] - KELEM[9 * row + 4] ;
		//9th column
		KELEM[9 * row + 8] = -KELEM[9 * row + 2] - KELEM[9 * row + 5] ;
	}

	/*for(int i=0;i<9;i++)
	{
	for(int j=0;j<9;j++)
	{
	std::cout << std::setw(4) << KELEM[9*i+j] << " ";
	}
	std::cout << "\n";
	}
	std::cout << "------------\n";*/
}

void ImplicitHyperElasticFEMSolver::addElementStiffnessMatrixToGlobalStiffnessMatrix(Cloth_Data* cloth,int el,std::vector<double> &KELEM)
{
	Triangles tri = cloth->getMesh()->get_triangle(el);
	std::vector<int> vertices(3);
	vertices[0] = tri.a;
	vertices[1] = tri.b;
	vertices[2] = tri.c;

	for (int vtxIndexA = 0; vtxIndexA < 3; vtxIndexA++)
	{
		for (int vtxIndexB = 0; vtxIndexB < 3; vtxIndexB++) {
			int vtxA = vertices[vtxIndexA];
			int vtxB = vertices[vtxIndexB];
			for (int i = 0; i < 3; i++) 
			{
				for (int j = 0; j < 3; j++) 
				{
					int row = 3 * vtxA + i;
					double * value = &KELEM[ELT(9, 3 * vtxIndexA + i, 3 * vtxIndexB + j)];
					int columnIndex = tangentStiffnessMatrix_->GetInverseIndex(row, 3 * vtxB + j);
					tangentStiffnessMatrix_->AddEntry(row, columnIndex, *value);
				}
			}
		}
	}
}


void ImplicitHyperElasticFEMSolver::addShearComponents(Cloth_Data *cloth)
{
	int bad=0;

	float young_modulus = cloth->get_property_obj()->get_youngs_modulus();
	float poisson = cloth->get_property_obj()->get_poisson_ratio();

	float lame_lambda = (young_modulus*poisson)/((1.0+poisson)*(1.0-(2.0*poisson)));
	float lame_mu = young_modulus/(2.0*(1.0+poisson));

	float h = cloth->get_property_obj()->get_timestep();

	int num_vertices = cloth->getMesh()->get_number_vertices();
	int num_triangles = cloth->getMesh()->get_number_triangles();

	for(int t=0;t<num_triangles;t++)
	{
		Triangles tri = cloth->getMesh()->get_triangle(t);

		Eigen::Vector3d pa = cloth->getMesh()->get_point_data( tri.a,lastFrameId_);
		Eigen::Vector3d pb = cloth->getMesh()->get_point_data( tri.b,lastFrameId_);
		Eigen::Vector3d pc = cloth->getMesh()->get_point_data( tri.c,lastFrameId_);

		Eigen::Matrix2d Bm = cloth->getDmInv(t);

		Eigen::Vector3d diff1 = pc-pa;
		Eigen::Vector3d diff2 = pc-pb;

		Eigen::MatrixXd tmp(3,2);
		tmp(0,0) = diff1.x(); tmp(0,1) = diff2.x();
		tmp(1,0) = diff1.y(); tmp(1,1) = diff2.y();
		tmp(2,0) = diff1.z(); tmp(2,1) = diff2.z();

		Eigen::MatrixXd F = tmp*Bm;

		//Assume here that we are using NeoHookean Materials 
		std::vector<double> gradients(3);
		std::vector<double> hessians(6);

		double IIIC = (F(0,1)*F(0,1)*F(1,0)*F(1,0)) - (2.0*F(0,0)*F(0,1)*F(1,0)*F(1,1)) + (F(0,0)*F(0,0)*F(1,1)*F(1,1))+ (F(0,1)*F(0,1)*F(2,0)*F(2,0)) + (F(1,1)*F(1,1)*F(2,0)*F(2,0)) - (2.0*F(0,0)*F(0,1)*F(2,0)*F(2,1)) - (2.0*F(1,0)*F(1,1)*F(2,0)*F(2,1)) + (F(0,0)*F(0,0)*F(2,1)*F(2,1)) + (F(1,0)*F(1,0)*F(2,1)*F(2,1));
		double IC = F(0,0)*F(0,0) + F(0,1)*F(0,1) + F(1,0)*F(1,0) + F(1,1)*F(1,1) + F(2,0)*F(2,0) + F(2,1)*F(2,1);

#ifdef NEOHOOKEAN_MATERIAL 
		gradients[0] = 0.5 * lame_mu;
		gradients[1] = 0.0;
		gradients[2] = (-0.5 * lame_mu + 0.25 * lame_lambda * log(IIIC)) / IIIC;

		hessians[0] = 0.0;
		hessians[1] = 0.0;
		hessians[2] = 0.0;
		hessians[3] = 0.0;
		hessians[4] = 0.0;
		hessians[5] = (0.25 * lame_lambda + 0.5 * lame_mu - 0.25 * lame_lambda * log(IIIC)) / (IIIC * IIIC);

#else
		gradients[0] = 0.25 * lame_lambda * (IC - 3.0) - 0.5 * lame_mu;
		gradients[1] = 0.25 * lame_mu;
		gradients[2] = 0.0;

		hessians[0] = 0.25 * lame_lambda;
		hessians[1] = 0.0;
		hessians[2] = 0.0;
		hessians[3] = 0.0;
		hessians[4] = 0.0;
		hessians[5] = 0.0;
#endif		


		//Add compression resistance gradient if needed 
#ifdef ADDCOMPRESSIONRESISTENCE
		double J = sqrt(IIIC);
		if (J < 1.0)
		{
			double compressionResistanceFactor =1.0;
			gradients[2] += -compressionResistanceFactor * (J - 1.0) * (J - 1.0) / (1728.0 * J);
			hessians[5] += compressionResistanceFactor * (1.0 - J) * (1.0 + J) / (3456.0 * J * J * J);
		}
#endif

		//Directly compute the piola kirchoff stress tensor 
		Eigen::MatrixXd P_1;
		computeFirstPiolaStress(F,gradients,P_1);

		Eigen::Vector2d ewan0 = cloth->getEdgeWeightedTriangleNormals(t,0);
		Eigen::Vector2d ewan1 = cloth->getEdgeWeightedTriangleNormals(t,1);
		Eigen::Vector2d ewan2 = cloth->getEdgeWeightedTriangleNormals(t,2);

		Eigen::MatrixXd internalForce_0 = P_1 * ewan0;
		Eigen::MatrixXd internalForce_1 = P_1 * ewan1;
		Eigen::MatrixXd internalForce_2 = P_1 * ewan2;

		internalForce_[3*tri.a + 0] += internalForce_0(0,0);
		internalForce_[3*tri.a + 1] += internalForce_0(1,0);
		internalForce_[3*tri.a + 2] += internalForce_0(2,0);

		internalForce_[3*tri.b + 0] += internalForce_1(0,0);
		internalForce_[3*tri.b + 1] += internalForce_1(1,0);
		internalForce_[3*tri.b + 2] += internalForce_1(2,0);

		internalForce_[3*tri.c + 0] += internalForce_2(0,0);
		internalForce_[3*tri.c + 1] += internalForce_2(1,0);
		internalForce_[3*tri.c + 2] += internalForce_2(2,0);

		//Now we will compute the stiffness matrices one by one 

		//First we need to compute the value of dPdF_atFHat which is a 4x4 tensor
		std::vector<double> dPdF;
		std::vector<double> dGdF;
		std::vector<double> KElem;

		//Next we compute teh dPdF which is a 6x6 tensor
		computeDPDF(F,gradients,hessians,dPdF);

		//Next we compute the dGdF which is a 6x6 tensor
		Eigen::Matrix2d bVec;
		bVec(0,0) = ewan0(0);
		bVec(0,1) = ewan0(1);
		bVec(1,0) = ewan1(0);
		bVec(1,1) = ewan1(1);
		computeDGDF(dPdF,bVec,dGdF);

		//Finally using DGDF and precomputede dFdUs we can compute the elementary stiffness matrix 
		computeElementStiffnessMatrix(t,dGdF,KElem);
		//Insert it in the correct place in the global matrix
		addElementStiffnessMatrixToGlobalStiffnessMatrix(cloth,t,KElem);

		//Test the tangent stiffness matrix 
#ifdef TESTELEMENTSTIFFNESSEIGENVALUE
		Eigen::MatrixXd KElemE(9,9);
		for(int r=0;r<9;r++)
		{
			for(int c=0;c<9;c++)
			{
				KElemE(r,c) = KElem[ELT(9,r,c)];
			}
		}
		Eigen::EigenSolver<Eigen::MatrixXd> eigSolver(KElemE);
		Eigen::VectorXd eigVals = eigSolver.eigenvalues().real();
		std::vector<double> evSorted(9);
		memcpy(&evSorted[0],eigVals.data(),sizeof(double)*9);
		std::sort(evSorted.begin(),evSorted.end());
		//std::cout << "[" << t << "] ";
		bool badTri = false;
		for(int i=0;i<6;i++)
		{
			if(fabs(evSorted[i])<1.0)
				continue;
			badTri = true;
		}
		if(badTri)
		{
			std::cout << "[" << t << "] ";
			for(int i=0;i<9;i++)
			{
				std::cout << evSorted[i] << ",";
			}
			std::cout << "\n";
			bad++;
		}
#endif

		/*if(t==1000)
		{
			for(int i=0;i<36;i++)
			{
				std::cout << "[" << i << "] " << std::setw(4) << dGdF[i] << "," << dGdF3[i] << "\n";
			}
		}*/


	}

#ifdef TESTELEMENTSTIFFNESSEIGENVALUE
	std::cout << "Bad:" << bad << " Total:" << num_triangles << "\n";
#endif
}

void ImplicitHyperElasticFEMSolver::addGravityComponents(Cloth_Data *cloth) {
	int num_vertices = cloth->getMesh()->get_number_vertices();
	//OMP_FOR
	for(int p=0;p<num_vertices;p++)
	{
		float gravForce = cloth->get_vertex_mass(p) * 9.8f; 
		externalForce_[3*p+1] -= gravForce;
	}

}

#ifdef IMPLICITNEWMARKSOLVER
void ImplicitHyperElasticFEMSolver::prepareNumericalEquations(Cloth_Data* cloth)
{

	int numVertices = cloth->getMesh()->get_number_vertices();
	float timeStep = cloth->get_property_obj()->get_timestep();

	//First prepare the numerical systems
	double dampingStiffness = cloth->get_property_obj()->getRayleighDampingCoeff();
	double dampingMass = cloth->get_property_obj()->getRayleighMassCoeff();

	tangentStiffnessMatrix_->ScalarMultiply(dampingStiffness, rayleighDampingMatrix_);
	rayleighDampingMatrix_->AddSubMatrix(dampingMass,*massMatrix_);
	rayleighDampingMatrix_->ScalarMultiplyAdd(alpha4_, tangentStiffnessMatrix_);
	tangentStiffnessMatrix_->AddSubMatrix(alpha1_, *massMatrix_);

	std::vector<double> residual(numVertices*3,0.0);

	massMatrix_->MultiplyVector(&qAccel_[0], &residual[0]);
  rayleighDampingMatrix_->MultiplyVectorAdd(qVel_.data(), &residual[0]);

	for(int i=0; i<numVertices*3; i++)
	{
		residual[i] += internalForce_[i] - externalForce_[i];	
		residual[i] *= (-1.0);
		residuals_[i] = residual[i];

	}

}
#else
void ImplicitHyperElasticFEMSolver::prepareNumericalEquations(Cloth_Data* cloth)
{

	int numVertices = cloth->getMesh()->get_number_vertices();
	float timeStep = cloth->get_property_obj()->get_timestep();

	//First prepare the numerical systems
	double dampingStiffness = cloth->get_property_obj()->getRayleighDampingCoeff();
	double dampingMass = cloth->get_property_obj()->getRayleighMassCoeff();

	tangentStiffnessMatrix_->ScalarMultiply(dampingStiffness, rayleighDampingMatrix_);
	rayleighDampingMatrix_->AddSubMatrix(dampingMass,*massMatrix_);
	std::vector<double> buffer(numVertices*3,0);
	std::vector<double> residual(numVertices*3,0.0);

	if(numIters_ != 0) 
	{
		// add K * (q_1 - q) to qresidual (will multiply by -h later)
		for(int i=0; i<numVertices*3; i++)
			buffer[i] = (q1_[i] - q_[i]);
		tangentStiffnessMatrix_->MultiplyVectorAdd(&buffer[0],&residual[0]);
	}

	*tangentStiffnessMatrix_ *= timeStep;
	*tangentStiffnessMatrix_ += *rayleighDampingMatrix_;
	tangentStiffnessMatrix_->MultiplyVectorAdd(qVel_.data(),&residual[0]);

	*tangentStiffnessMatrix_ *= timeStep;
	tangentStiffnessMatrix_->AddSubMatrix(1.0, *massMatrix_);

	for(int i=0; i<numVertices*3; i++)
	{
		residual[i] += internalForce_[i] - externalForce_[i];	
		residual[i] *= (-timeStep);
	}


	if (numIters_ != 0) // can skip on first iteration (zero contribution)
	{
		// add M * (qvel_1 - qvel) to qresidual
		for(int i=0; i<numVertices*3; i++)
			buffer[i] = qVel1_[i] - qVel_[i];
		massMatrix_->MultiplyVectorAdd(&buffer[0], &residual[0]);
	}

	for(int i=0; i<numVertices*3; i++)
		residuals_[i] = residual[i];
}
#endif

bool ImplicitHyperElasticFEMSolver::checkConvergence(Cloth_Data* cloth)
{
	bool converged = false;
	int n3 = (cloth->getMesh()->get_number_vertices())*3;
	
	//Subsize the matrix to form the constrained equations
	int numConstrainedDOFs = n3 - (3*numConstrainedVerts_);
	RemoveRows(n3, &constrainedResiduals_[0], &residuals_[0], numConstrainedVerts_*3, constrainedVerts_);

	/*for(int i=0;i<numConstrainedDOFs;i++)
	{
		std::cout << constrainedResiduals_[i] << residuals_[i] << "\n";
	}*/

	double error = 0.0;
  for(int i=0; i<numConstrainedDOFs; i++)
    error += constrainedResiduals_[i] * constrainedResiduals_[i];

	//std::cout << "I:" << numIters_ << " Error:" << error << "\n";

	//getchar();

	if (numIters_ == 0) 
	{
		error_ = error;
		errorQuotient_ = 1.0;
	}
	else
	{
		// error divided by the initial error, before performing this iteration
		errorQuotient_ = error / error_; 
	}

	std::cout << "I:" << numIters_ << " Error:" << errorQuotient_ << "\n";


	const double epsilon = 1.0e-7;

	if (errorQuotient_ < epsilon)
	{
		converged = true;
	}

	return converged;
}

void ImplicitHyperElasticFEMSolver::solverStep(Cloth_Data* cloth)
{
	float timeStep = cloth->get_property_obj()->get_timestep();
	int numVertices = cloth->getMesh()->get_number_vertices();
	int n3 = (cloth->getMesh()->get_number_vertices())*3;

	//First remove the rows from the tangent stiffness matrix 
	int numConstrainedDOFs = n3 - (3*numConstrainedVerts_);
	SparseMatrix *tempSpMatCopy = new SparseMatrix(*tangentStiffnessMatrix_);
	tempSpMatCopy->RemoveRowsColumns(numConstrainedVerts_*3,constrainedVerts_);

	double *resultConstrained = new double[numConstrainedDOFs];

#if defined USE_PARDISO_SOLVER
	PardisoSolver solver(tempSpMatCopy,7,0,0,0);
	solver.ComputeCholeskyDecomposition(tempSpMatCopy);
	int retVal = solver.SolveLinearSystem(resultConstrained,&constrainedResiduals_[0]);
#elif defined USE_CG_SOLVER
	CGSolver cgsolver(tempSpMatCopy);
	cgsolver.SolveLinearSystemWithJacobiPreconditioner(resultConstrained, &constrainedResiduals_[0], 1e-6, 10000);
#else
	int numGaussSeidelIterations = 8000;
	for(int iter=0; iter<numGaussSeidelIterations; iter++)
		tempSpMatCopy->DoOneGaussSeidelIteration(resultConstrained, &constrainedResiduals_[0]);
#endif

	//Expand the rows now to fill it 
	double *delV = new double[n3];
	InsertRows(n3,resultConstrained,delV,numConstrainedVerts_*3,constrainedVerts_);
	//Free the stuff we allocated
	delete(tempSpMatCopy);
	delete[] resultConstrained;

	for(int p=0;p<numVertices;p++)
	{
#ifdef IMPLICITNEWMARKSOLVER
		Eigen::Vector3d oldPos = cloth->getMesh()->get_point_data(p,0);
		q_[3*p+0] += delV[3*p+0];
		q_[3*p+1] += delV[3*p+1];
		q_[3*p+2] += delV[3*p+2];
		Eigen::Vector3d elemDelQ(q_[3*p+0],q_[3*p+1],q_[3*p+2]);
		Eigen::Vector3d newPos = oldPos + elemDelQ;
		cloth->setCurrentStepPositions(newPos,p);

		qAccel_[3*p+0] = alpha1_ * (q_[3*p+0] - q1_[3*p+0]) - alpha2_ * qVel1_[3*p+0] - alpha3_ * qAccel1_[3*p+0];
		qAccel_[3*p+1] = alpha1_ * (q_[3*p+1] - q1_[3*p+1]) - alpha2_ * qVel1_[3*p+1] - alpha3_ * qAccel1_[3*p+1];
		qAccel_[3*p+2] = alpha1_ * (q_[3*p+2] - q1_[3*p+2]) - alpha2_ * qVel1_[3*p+2] - alpha3_ * qAccel1_[3*p+2];

    qVel_[3*p+0] = alpha4_ * (q_[3*p+0] - q1_[3*p+0]) + alpha5_ * qVel1_[3*p+0] + alpha6_ * qAccel1_[3*p+0];
    qVel_[3*p+1] = alpha4_ * (q_[3*p+1] - q1_[3*p+1]) + alpha5_ * qVel1_[3*p+1] + alpha6_ * qAccel1_[3*p+1];
    qVel_[3*p+2] = alpha4_ * (q_[3*p+2] - q1_[3*p+2]) + alpha5_ * qVel1_[3*p+2] + alpha6_ * qAccel1_[3*p+2];

		Eigen::Vector3d elemVel(qVel_[3*p+0],qVel_[3*p+1],qVel_[3*p+2]);
		cloth->setCurrentStepVelocity(elemVel,p);
#else
		Eigen::Vector3d oldPos = cloth->getMesh()->get_point_data(p,0);

		qVel_[3*p+0] += delV[3*p+0];
		qVel_[3*p+1] += delV[3*p+1];
		qVel_[3*p+2] += delV[3*p+2];

		Eigen::Vector3d elemVel(qVel_[3*p+0],qVel_[3*p+1],qVel_[3*p+2]);
		cloth->setCurrentStepVelocity(elemVel,p);

		q_[3*p+0] += q1_[3*p+0] - q_[3*p+0] + timeStep * qVel_[3*p+0];
		q_[3*p+1] += q1_[3*p+1] - q_[3*p+1] + timeStep * qVel_[3*p+1];
		q_[3*p+2] += q1_[3*p+2] - q_[3*p+2] + timeStep * qVel_[3*p+2];


		Eigen::Vector3d elemDelQ(q_[3*p+0],q_[3*p+1],q_[3*p+2]);
		Eigen::Vector3d newPos = oldPos + elemDelQ;
		cloth->setCurrentStepPositions(newPos,p);
#endif
	}
}

void ImplicitHyperElasticFEMSolver::finalizeStep(Cloth_Data* cloth)
{
	int numVertices = cloth->getMesh()->get_number_vertices();
	for(int p=0;p<numVertices;p++)
	{
		Eigen::Vector3d currentPos(oldPosition_[3*p+0],oldPosition_[3*p+1],oldPosition_[3*p+2]);
		cloth->setCurrentStepPositions(currentPos,p);

		Eigen::Vector3d oldPos = cloth->getMesh()->get_point_data(p,0);

		Eigen::Vector3d elemVel(qVel_[3*p+0],qVel_[3*p+1],qVel_[3*p+2]);
		cloth->set_next_step_velocity(elemVel,p);

		Eigen::Vector3d elemDelQ(q_[3*p+0],q_[3*p+1],q_[3*p+2]);
		Eigen::Vector3d newPos = oldPos + elemDelQ;
		cloth->set_next_step_pos(newPos,p);
	}
}

#ifdef IMPLICITNEWMARKSOLVER
void ImplicitHyperElasticFEMSolver::initializeAlphas(Cloth_Data* cloth)
{
	double newMarkBeta = 0.25 ;
	double newMarkGamma = 0.5;
	double timeStep = cloth->get_property_obj()->get_timestep();

	alpha1_ = 1.0 / (newMarkBeta * timeStep * timeStep);
  alpha2_ = 1.0 / (newMarkBeta * timeStep);
  alpha3_ = (1.0 - 2.0 * newMarkBeta) / (2.0 * newMarkBeta);
  alpha4_ = newMarkGamma / (newMarkBeta * timeStep);
  alpha5_ = 1 - newMarkGamma/newMarkBeta;
  alpha6_ = (1.0 - newMarkGamma / (2.0 * newMarkBeta)) * timeStep;
}
#endif

void ImplicitHyperElasticFEMSolver::resetParameters() {
}

void ImplicitHyperElasticFEMSolver::generateConstantEnergyHessian(Cloth_Data *cloth) {
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
	//temp = kd*(*constantEnergyHessian_);
	//mKdQ_ = new SparseMatrix(temp);
	//temp = ke*(-h)*h*(*constantEnergyHessian_);
	//mhSquaredKeQ_ = new SparseMatrix(temp);
	//temp = ke*h*(*constantEnergyHessian_);
	//mhKeQ_ = new SparseMatrix(temp);
	//temp = kd*(-h)*(*constantEnergyHessian_);
	//mhKdQ_ = new SparseMatrix(temp);

	//Delete the residuals
	delete (spMatOutline);
	delete[] fakeIdentity;
}

void ImplicitHyperElasticFEMSolver::addQuadraticBending(Cloth_Data* cloth)
{
	// Assemble the force and velocity 
	int numVertices = cloth->getMesh()->get_number_vertices();
	double *x = new double[numVertices*3];
	std::vector<Eigen::Vector3d> posVec = cloth->getMesh()->getPositionVector(lastFrameId_);
	memcpy(x,&posVec[0],sizeof(double)*3*numVertices);

	//Find the elastic force and add it to LHS 
	double *fe = new double[numVertices*3];
	mKeQ_->MultiplyVector(x,fe);
	Eigen::Map<Eigen::VectorXd> fe_Map(fe,numVertices*3);
	internalForce_ -= fe_Map;

	//Add the jacobian to the stiffness matrix
	tangentStiffnessMatrix_->AddSubMatrix(-1.0,*mKeQ_,1);

	delete[] x;
	delete[] fe;

}