#include "ImplicitFEMSolver.h"
#include "PardisoSolver.h"

#define USECONJUGATEGRADIENTSOLVER 0

ImplicitFEMSolver::ImplicitFEMSolver(void)
	{
	}


ImplicitFEMSolver::~ImplicitFEMSolver(void)
	{
		 delete[] massMatrixDiagonal_;
		 delete LHSMatrix_;
		 delete[] constrainedVerts_;
	}

void ImplicitFEMSolver::initializeSparseMatrixFromOutline(Cloth_Data *cloth) {
	int numVertices = cloth->getMesh()->get_number_vertices();
	int numTriangles = cloth->getMesh()->get_number_triangles();
	SparseMatrixOutline *spMatOutline = new SparseMatrixOutline(numVertices*3);

	double *zero3x3 = new double[9];
	memset(zero3x3,0,sizeof(double)*9);

	//Initialize the diagonal 3x3 blocks
	for(int p=0; p<numVertices;p++) {
		spMatOutline->AddBlock3x3Entry(p,p,zero3x3);
	}
	//Initialize the offdiagonal 3x3 blocks
	for(int t=0;t<numTriangles;t++)
	{
		Triangles tri = cloth->getMesh()->get_triangle(t);
		int a = tri.a;
		int b = tri.b;
		int c = tri.c;

		//ab ac ba bc ca cb 
		spMatOutline->AddBlock3x3Entry(a,b,zero3x3);
		spMatOutline->AddBlock3x3Entry(a,c,zero3x3);
		spMatOutline->AddBlock3x3Entry(b,a,zero3x3);
		spMatOutline->AddBlock3x3Entry(b,c,zero3x3);
		spMatOutline->AddBlock3x3Entry(c,a,zero3x3);
		spMatOutline->AddBlock3x3Entry(c,b,zero3x3);
	}
	//Initialize the actual sparse matrix 
	LHSMatrix_= new SparseMatrix(spMatOutline);

	//Delete the outline 
	delete (spMatOutline);
	delete[] zero3x3;

	//Create the diagonal indices
	LHSMatrix_->BuildDiagonalIndices();

}


void ImplicitFEMSolver::advance_time_step(Cloth_Data *cloth) {
	//Set the last frame information
	lastFrameId_ = (cloth->getMesh()->get_number_frames()) - 1;
	
	//Prepare the LHS matrix for usage
	LHSMatrix_->ResetToZero();
	LHSMatrix_->AddDiagonalMatrix(massMatrixDiagonal_);
	
	//Prepare the RHS Vector for usage
	RHSVector_.setZero();

	//Add the physics components
	addShearComponents(cloth);
	//addShearComponentsCorotational(cloth);
	addGravityComponents(cloth);

	//Solve and report
	finalizeAndSolve(cloth);
}

void ImplicitFEMSolver::initialize(Cloth_Data *cloth) {
	initializeSparseMatrixFromOutline(cloth);
	
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

}

void ImplicitFEMSolver::addShearComponents(Cloth_Data *cloth) {
	float young_modulus = cloth->get_property_obj()->get_youngs_modulus();
	float poisson = cloth->get_property_obj()->get_poisson_ratio();

	float kd = 0.01;//cloth->get_property_obj()->get_damping_param();
	float h = cloth->get_property_obj()->get_timestep();

	int num_vertices = cloth->getMesh()->get_number_vertices();
	int num_triangles = cloth->getMesh()->get_number_triangles();

	//std::cout << num_vertices << " " << num_triangles << "\n";

	//@TODO: make this more generic by replacing it with a generic function 
	Eigen::Matrix3d E = isotropic_linear_elasticity(young_modulus, poisson);
	Eigen::Matrix3d E_Prime = -kd*(Eigen::Matrix3d::Identity());
	
	//OMP_FOR
	for(int t=0;t<num_triangles;t++)
	{
		Triangles tri = cloth->getMesh()->get_triangle(t);

		Eigen::Vector3d pa = cloth->getMesh()->get_point_data( tri.a,lastFrameId_);
		Eigen::Vector3d pb = cloth->getMesh()->get_point_data( tri.b,lastFrameId_);
		Eigen::Vector3d pc = cloth->getMesh()->get_point_data( tri.c,lastFrameId_);

		Eigen::Vector3d va = cloth->get_velocity(tri.a);
		Eigen::Vector3d vb = cloth->get_velocity(tri.b);
		Eigen::Vector3d vc = cloth->get_velocity(tri.c);

		float r_ua = (cloth->get_vertex_distribution(t))(0);
		float r_va = (cloth->get_vertex_distribution(t))(1);
		float r_ub = (cloth->get_vertex_distribution(t))(2);
		float r_vb = (cloth->get_vertex_distribution(t))(3);
		float r_uc = (cloth->get_vertex_distribution(t))(4);
		float r_vc = (cloth->get_vertex_distribution(t))(5);
		float d = (cloth->get_vertex_distribution(t))(6);

		Eigen::Vector3d U = (r_ua * pa) + (r_ub * pb) + (r_uc * pc);
		Eigen::Vector3d V = (r_va * pa) + (r_vb * pb) + (r_vc * pc);

		Eigen::Vector3d U_Prime = (r_ua * va) + (r_ub * vb) + (r_uc * vc);
		Eigen::Vector3d V_Prime = (r_va * va) + (r_vb * vb) + (r_vc * vc);

		Eigen::Vector3d strain_tensor_prime;
		strain_tensor_prime[0] = U.dot(U_Prime);
		strain_tensor_prime[1] = V.dot(V_Prime);
		strain_tensor_prime[2] = U_Prime.dot(V) + U.dot(V_Prime);

		Eigen::Vector3d stress_tensor_prime = E_Prime * strain_tensor_prime;

		Eigen::Vector3d strain_tensor;
		strain_tensor[0] = 0.5f*(U.dot(U)-1.0f);
		strain_tensor[1] = 0.5f*(V.dot(V)-1.0f);
		strain_tensor[2] = U.dot(V);

		//std::cout << strain_tensor(0) << "," << strain_tensor(2) << "," << strain_tensor(1) << "\n";

		Eigen::Vector3d stress_tensor = E * strain_tensor;
		stress_tensor += stress_tensor_prime;

		if(stress_tensor[0]<=0.0)
			stress_tensor[0] = 0.0;

		if(stress_tensor[1]<=0.0)
			stress_tensor[1] = 0.0;

		Eigen::Vector3d fA = (stress_tensor[0] * r_ua * U) + (stress_tensor[1] * r_va * V) + (stress_tensor[2] * r_ua * V ) + (stress_tensor[2] * r_va * U );
		fA *= (-fabs(d)*0.5f);
		Eigen::Vector3d fB = (stress_tensor[0] * r_ub * U) + (stress_tensor[1] * r_vb * V) + (stress_tensor[2] * r_ub * V ) + (stress_tensor[2] * r_vb * U );
		fB *= (-fabs(d)*0.5f);
		Eigen::Vector3d fC = (stress_tensor[0] * r_uc * U) + (stress_tensor[1] * r_vc * V) + (stress_tensor[2] * r_uc * V ) + (stress_tensor[2] * r_vc * U );
		fC *= (-fabs(d)*0.5f);

		RHSVector_[3*tri.a + 0] += fA[0];
		RHSVector_[3*tri.a + 1] += fA[1];
		RHSVector_[3*tri.a + 2] += fA[2];

		RHSVector_[3*tri.b + 0] += fB[0];
		RHSVector_[3*tri.b + 1] += fB[1];
		RHSVector_[3*tri.b + 2] += fB[2];	

		RHSVector_[3*tri.c + 0] += fC[0];
		RHSVector_[3*tri.c + 1] += fC[1];
		RHSVector_[3*tri.c + 2] += fC[2];

		 //Lets do the Jacobians one by one 
		Eigen::Matrix3d J_aa = E*(r_ua*U)*((r_ua*U).transpose()) + E*(r_va*V)*((r_va*V).transpose()) + E*(r_ua*V + r_va*U)*((r_ua*V + r_va*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_ua*r_ua + stress_tensor[1]*r_va*r_va + stress_tensor[2]*(r_ua*r_va + r_va*r_ua));
		J_aa *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_aa_prime = E_Prime*(r_ua*U)*((r_ua*U).transpose()) + E_Prime*(r_va*V)*((r_va*V).transpose()) + E_Prime*(r_ua*V + r_va*U)*((r_ua*V + r_va*U).transpose()) ;
		J_aa_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.a,3*tri.a,J_aa.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*tri.a,3*tri.a,J_aa_prime.data(),(-h));

		Eigen::Matrix3d J_ab = E*(r_ua*U)*((r_ub*U).transpose()) + E*(r_va*V)*((r_vb*V).transpose()) + E*(r_ua*V + r_va*U)*((r_ub*V + r_vb*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_ua*r_ub + stress_tensor[1]*r_va*r_vb + stress_tensor[2]*(r_ua*r_vb + r_va*r_ub));
		J_ab *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_ab_prime = E_Prime*(r_ua*U)*((r_ub*U).transpose()) + E_Prime*(r_va*V)*((r_vb*V).transpose()) + E_Prime*(r_ua*V + r_va*U)*((r_ub*V + r_vb*U).transpose()) ;
		J_ab_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.a,3*tri.b,J_ab.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*tri.a,3*tri.b,J_ab_prime.data(),(-h));

		Eigen::Matrix3d J_ac = E*(r_ua*U)*((r_uc*U).transpose()) + E*(r_va*V)*((r_vc*V).transpose()) + E*(r_ua*V + r_va*U)*((r_uc*V + r_vc*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_ua*r_uc + stress_tensor[1]*r_va*r_vc + stress_tensor[2]*(r_ua*r_vc + r_va*r_uc));
		J_ac *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_ac_prime = E_Prime*(r_ua*U)*((r_uc*U).transpose()) + E_Prime*(r_va*V)*((r_vc*V).transpose()) + E_Prime*(r_ua*V + r_va*U)*((r_uc*V + r_vc*U).transpose()) ;
		J_ac_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.a,3*tri.c,J_ac.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*tri.a,3*tri.c,J_ac_prime.data(),(-h));
		
		Eigen::Matrix3d J_ba = E*(r_ub*U)*((r_ua*U).transpose()) + E*(r_vb*V)*((r_va*V).transpose()) + E*(r_ub*V + r_vb*U)*((r_ua*V + r_va*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_ub*r_ua + stress_tensor[1]*r_vb*r_va + stress_tensor[2]*(r_ub*r_va + r_vb*r_ua));
		J_ba *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_ba_prime = E_Prime*(r_ub*U)*((r_ua*U).transpose()) + E_Prime*(r_vb*V)*((r_va*V).transpose()) + E_Prime*(r_ub*V + r_vb*U)*((r_ua*V + r_va*U).transpose()) ;
		J_ba_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.b,3*tri.a,J_ba.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*tri.b,3*tri.a,J_ba_prime.data(),(-h));

		Eigen::Matrix3d J_bb = E*(r_ub*U)*((r_ub*U).transpose()) + E*(r_vb*V)*((r_vb*V).transpose()) + E*(r_ub*V + r_vb*U)*((r_ub*V + r_vb*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_ub*r_ub + stress_tensor[1]*r_vb*r_vb + stress_tensor[2]*(r_ub*r_vb + r_vb*r_ub));
		J_bb *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_bb_prime = E_Prime*(r_ub*U)*((r_ub*U).transpose()) + E_Prime*(r_vb*V)*((r_vb*V).transpose()) + E_Prime*(r_ub*V + r_vb*U)*((r_ub*V + r_vb*U).transpose()) ;
		J_bb_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.b,3*tri.b,J_bb.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*tri.b,3*tri.b,J_bb_prime.data(),(-h));

		Eigen::Matrix3d J_bc = E*(r_ub*U)*((r_uc*U).transpose()) + E*(r_vb*V)*((r_vc*V).transpose()) + E*(r_ub*V + r_vb*U)*((r_uc*V + r_vc*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_ub*r_uc + stress_tensor[1]*r_vb*r_vc + stress_tensor[2]*(r_ub*r_vc + r_vb*r_uc));
		J_bc *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_bc_prime = E_Prime*(r_ub*U)*((r_uc*U).transpose()) + E_Prime*(r_vb*V)*((r_vc*V).transpose()) + E_Prime*(r_ub*V + r_vb*U)*((r_uc*V + r_vc*U).transpose()) ;
		J_bc_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.b,3*tri.c,J_bc.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*tri.b,3*tri.c,J_bc_prime.data(),(-h));

		Eigen::Matrix3d J_ca = E*(r_uc*U)*((r_ua*U).transpose()) + E*(r_vc*V)*((r_va*V).transpose()) + E*(r_uc*V + r_vc*U)*((r_ua*V + r_va*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_uc*r_ua + stress_tensor[1]*r_vc*r_va + stress_tensor[2]*(r_uc*r_va + r_vc*r_ua));
		J_ca *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_ca_prime = E_Prime*(r_uc*U)*((r_ua*U).transpose()) + E_Prime*(r_vc*V)*((r_va*V).transpose()) + E_Prime*(r_uc*V + r_vc*U)*((r_ua*V + r_va*U).transpose()) ;
		J_ca_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.c,3*tri.a,J_ca.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*tri.c,3*tri.a,J_ca_prime.data(),(-h));

		Eigen::Matrix3d J_cb = E*(r_uc*U)*((r_ub*U).transpose()) + E*(r_vc*V)*((r_vb*V).transpose()) + E*(r_uc*V + r_vc*U)*((r_ub*V + r_vb*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_uc*r_ub + stress_tensor[1]*r_vc*r_vb + stress_tensor[2]*(r_uc*r_vb + r_vc*r_ub));
		J_cb *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_cb_prime = E_Prime*(r_uc*U)*((r_ub*U).transpose()) + E_Prime*(r_vc*V)*((r_vb*V).transpose()) + E_Prime*(r_uc*V + r_vc*U)*((r_ub*V + r_vb*U).transpose()) ;
		J_cb_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.c,3*tri.b,J_cb.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*tri.c,3*tri.b,J_cb_prime.data(),(-h));

		Eigen::Matrix3d J_cc = E*(r_uc*U)*((r_uc*U).transpose()) + E*(r_vc*V)*((r_vc*V).transpose()) + E*(r_uc*V + r_vc*U)*((r_uc*V + r_vc*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_uc*r_uc + stress_tensor[1]*r_vc*r_vc + stress_tensor[2]*(r_uc*r_vc + r_vc*r_uc));
		J_cc *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_cc_prime = E_Prime*(r_uc*U)*((r_uc*U).transpose()) + E_Prime*(r_vc*V)*((r_vc*V).transpose()) + E_Prime*(r_uc*V + r_vc*U)*((r_uc*V + r_vc*U).transpose()) ;
		J_cc_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.c,3*tri.c,J_cc.data(),(-h*h));
		LHSMatrix_->Add3x3Entry(3*tri.c,3*tri.c,J_cc_prime.data(),(-h));

		//Eigen::Matrix3d J_ji = E*(r_uj*U)*((r_ui*U).transpose()) + E*(r_vj*V)*((r_vi*V).transpose()) + E*(r_uj*V + r_vj*U)*((r_ui*V + r_vi*U).transpose()) 
		//	                     + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_uj*r_ui + stress_tensor[1]*r_vj*r_vi + stress_tensor[2]*(r_uj*r_vi + r_vj*r_ui));

		//Now we will need to add the velocities in their correct places also 
		Eigen::Vector3d df_dx_a = h*(J_aa * va + J_ab * vb + J_ac * vc);
		Eigen::Vector3d df_dx_b = h*(J_ba * va + J_bb * vb + J_bc * vc);
		Eigen::Vector3d df_dx_c = h*(J_ca * va + J_cb * vb + J_cc * vc);

		RHSVector_[3*tri.a+0] += df_dx_a[0];
		RHSVector_[3*tri.a+1] += df_dx_a[1];
		RHSVector_[3*tri.a+2] += df_dx_a[2];

		RHSVector_[3*tri.b+0] += df_dx_b[0];
		RHSVector_[3*tri.b+1] += df_dx_b[1];
		RHSVector_[3*tri.b+2] += df_dx_b[2];

		RHSVector_[3*tri.c+0] += df_dx_c[0];
		RHSVector_[3*tri.c+1] += df_dx_c[1];
		RHSVector_[3*tri.c+2] += df_dx_c[2];

	}
} 

void ImplicitFEMSolver::addShearComponentsCorotational(Cloth_Data *cloth)
{
	float young_modulus = cloth->get_property_obj()->get_youngs_modulus();
	float poisson = cloth->get_property_obj()->get_poisson_ratio();

	float lame_lambda = (young_modulus*poisson)/((1.0+poisson)*(1.0-(2.0*poisson)));
	float lame_mu = young_modulus/(2.0*(1.0+poisson));

	float kd = 0.01;//cloth->get_property_obj()->get_damping_param();
	float h = cloth->get_property_obj()->get_timestep();

	int num_vertices = cloth->getMesh()->get_number_vertices();
	int num_triangles = cloth->getMesh()->get_number_triangles();

	//std::cout << num_vertices << " " << num_triangles << "\n";

	//@TODO: make this more generic by replacing it with a generic function 
	Eigen::Matrix3d E = isotropic_linear_elasticity(young_modulus, poisson);
	Eigen::Matrix3d E_Prime = -kd*(Eigen::Matrix3d::Identity());
	
	//OMP_FOR
	for(int t=0;t<num_triangles;t++)
	{
		Triangles tri = cloth->getMesh()->get_triangle(t);

		Eigen::Vector3d pa = cloth->getMesh()->get_point_data( tri.a,lastFrameId_);
		Eigen::Vector3d pb = cloth->getMesh()->get_point_data( tri.b,lastFrameId_);
		Eigen::Vector3d pc = cloth->getMesh()->get_point_data( tri.c,lastFrameId_);

		Eigen::Vector3d va = cloth->get_velocity(tri.a);
		Eigen::Vector3d vb = cloth->get_velocity(tri.b);
		Eigen::Vector3d vc = cloth->get_velocity(tri.c);

		float r_ua = (cloth->get_vertex_distribution(t))(0);
		float r_va = (cloth->get_vertex_distribution(t))(1);
		float r_ub = (cloth->get_vertex_distribution(t))(2);
		float r_vb = (cloth->get_vertex_distribution(t))(3);
		float r_uc = (cloth->get_vertex_distribution(t))(4);
		float r_vc = (cloth->get_vertex_distribution(t))(5);
		float d = (cloth->get_vertex_distribution(t))(6);

		Eigen::Vector3d U = (r_ua * pa) + (r_ub * pb) + (r_uc * pc);
		Eigen::Vector3d V = (r_va * pa) + (r_vb * pb) + (r_vc * pc);

		Eigen::MatrixXd F(3,2);
		F(0,0) = U(0); F(0,1) = V(0);
		F(1,0) = U(1); F(1,1) = V(1);
		F(2,0) = U(2); F(2,1) = V(2);

		Eigen::MatrixXd R_Mat = genericPolarDecomposition(F);

		//Calculate the first piola kirchoff stress tensor
		Eigen::Matrix2d RTFMI = R_Mat.transpose()*F - Eigen::Matrix2d::Identity();

		Eigen::MatrixXd P_F_1 = 2.0*lame_mu*(F-R_Mat) + lame_lambda*RTFMI.trace()*R_Mat;
		Eigen::MatrixXd P_F_2 = F.transpose()*P_F_1;

		Eigen::Vector3d U_Prime = (r_ua * va) + (r_ub * vb) + (r_uc * vc);
		Eigen::Vector3d V_Prime = (r_va * va) + (r_vb * vb) + (r_vc * vc);

		Eigen::Vector3d strain_tensor_prime;
		strain_tensor_prime[0] = U.dot(U_Prime);
		strain_tensor_prime[1] = V.dot(V_Prime);
		strain_tensor_prime[2] = U_Prime.dot(V) + U.dot(V_Prime);

		Eigen::Vector3d stress_tensor_prime = E_Prime * strain_tensor_prime;

		Eigen::Vector3d strain_tensor;
		strain_tensor[0] = 0.5f*(U.dot(U)-1.0f);
		strain_tensor[1] = 0.5f*(V.dot(V)-1.0f);
		strain_tensor[2] = U.dot(V);

		//std::cout << strain_tensor(0) << "," << strain_tensor(2) << "," << strain_tensor(1) << "\n";

		//Eigen::Vector3d stress_tensor = E * strain_tensor;
		Eigen::Vector3d stress_tensor;
		stress_tensor(0) = P_F_2(0,0);
		stress_tensor(1) = P_F_2(1,1);
		stress_tensor(2) = P_F_2(1,0);

		stress_tensor += stress_tensor_prime;

		Eigen::Vector3d fA = (stress_tensor[0] * r_ua * U) + (stress_tensor[1] * r_va * V) + (stress_tensor[2] * r_ua * V ) + (stress_tensor[2] * r_va * U );
		fA *= (-fabs(d)*0.5f);
		Eigen::Vector3d fB = (stress_tensor[0] * r_ub * U) + (stress_tensor[1] * r_vb * V) + (stress_tensor[2] * r_ub * V ) + (stress_tensor[2] * r_vb * U );
		fB *= (-fabs(d)*0.5f);
		Eigen::Vector3d fC = (stress_tensor[0] * r_uc * U) + (stress_tensor[1] * r_vc * V) + (stress_tensor[2] * r_uc * V ) + (stress_tensor[2] * r_vc * U );
		fC *= (-fabs(d)*0.5f);

		RHSVector_[3*tri.a + 0] += fA[0];
		RHSVector_[3*tri.a + 1] += fA[1];
		RHSVector_[3*tri.a + 2] += fA[2];

		RHSVector_[3*tri.b + 0] += fB[0];
		RHSVector_[3*tri.b + 1] += fB[1];
		RHSVector_[3*tri.b + 2] += fB[2];	

		RHSVector_[3*tri.c + 0] += fC[0];
		RHSVector_[3*tri.c + 1] += fC[1];
		RHSVector_[3*tri.c + 2] += fC[2];

		 //Lets do the Jacobians one by one 
		Eigen::Matrix3d J_aa = E*(r_ua*U)*((r_ua*U).transpose()) + E*(r_va*V)*((r_va*V).transpose()) + E*(r_ua*V + r_va*U)*((r_ua*V + r_va*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_ua*r_ua + stress_tensor[1]*r_va*r_va + stress_tensor[2]*(r_ua*r_va + r_va*r_ua));
		J_aa *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_aa_prime = E_Prime*(r_ua*U)*((r_ua*U).transpose()) + E_Prime*(r_va*V)*((r_va*V).transpose()) + E_Prime*(r_ua*V + r_va*U)*((r_ua*V + r_va*U).transpose()) ;
		J_aa_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.a,3*tri.a,J_aa.data(),(-h*h));
		//LHSMatrix_->Add3x3Entry(3*tri.a,3*tri.a,J_aa_prime.data(),(-h));

		Eigen::Matrix3d J_ab = E*(r_ua*U)*((r_ub*U).transpose()) + E*(r_va*V)*((r_vb*V).transpose()) + E*(r_ua*V + r_va*U)*((r_ub*V + r_vb*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_ua*r_ub + stress_tensor[1]*r_va*r_vb + stress_tensor[2]*(r_ua*r_vb + r_va*r_ub));
		J_ab *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_ab_prime = E_Prime*(r_ua*U)*((r_ub*U).transpose()) + E_Prime*(r_va*V)*((r_vb*V).transpose()) + E_Prime*(r_ua*V + r_va*U)*((r_ub*V + r_vb*U).transpose()) ;
		J_ab_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.a,3*tri.b,J_ab.data(),(-h*h));
		//LHSMatrix_->Add3x3Entry(3*tri.a,3*tri.b,J_ab_prime.data(),(-h));

		Eigen::Matrix3d J_ac = E*(r_ua*U)*((r_uc*U).transpose()) + E*(r_va*V)*((r_vc*V).transpose()) + E*(r_ua*V + r_va*U)*((r_uc*V + r_vc*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_ua*r_uc + stress_tensor[1]*r_va*r_vc + stress_tensor[2]*(r_ua*r_vc + r_va*r_uc));
		J_ac *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_ac_prime = E_Prime*(r_ua*U)*((r_uc*U).transpose()) + E_Prime*(r_va*V)*((r_vc*V).transpose()) + E_Prime*(r_ua*V + r_va*U)*((r_uc*V + r_vc*U).transpose()) ;
		J_ac_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.a,3*tri.c,J_ac.data(),(-h*h));
		//LHSMatrix_->Add3x3Entry(3*tri.a,3*tri.c,J_ac_prime.data(),(-h));
		
		Eigen::Matrix3d J_ba = E*(r_ub*U)*((r_ua*U).transpose()) + E*(r_vb*V)*((r_va*V).transpose()) + E*(r_ub*V + r_vb*U)*((r_ua*V + r_va*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_ub*r_ua + stress_tensor[1]*r_vb*r_va + stress_tensor[2]*(r_ub*r_va + r_vb*r_ua));
		J_ba *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_ba_prime = E_Prime*(r_ub*U)*((r_ua*U).transpose()) + E_Prime*(r_vb*V)*((r_va*V).transpose()) + E_Prime*(r_ub*V + r_vb*U)*((r_ua*V + r_va*U).transpose()) ;
		J_ba_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.b,3*tri.a,J_ba.data(),(-h*h));
		//LHSMatrix_->Add3x3Entry(3*tri.b,3*tri.a,J_ba_prime.data(),(-h));

		Eigen::Matrix3d J_bb = E*(r_ub*U)*((r_ub*U).transpose()) + E*(r_vb*V)*((r_vb*V).transpose()) + E*(r_ub*V + r_vb*U)*((r_ub*V + r_vb*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_ub*r_ub + stress_tensor[1]*r_vb*r_vb + stress_tensor[2]*(r_ub*r_vb + r_vb*r_ub));
		J_bb *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_bb_prime = E_Prime*(r_ub*U)*((r_ub*U).transpose()) + E_Prime*(r_vb*V)*((r_vb*V).transpose()) + E_Prime*(r_ub*V + r_vb*U)*((r_ub*V + r_vb*U).transpose()) ;
		J_bb_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.b,3*tri.b,J_bb.data(),(-h*h));
		//LHSMatrix_->Add3x3Entry(3*tri.b,3*tri.b,J_bb_prime.data(),(-h));

		Eigen::Matrix3d J_bc = E*(r_ub*U)*((r_uc*U).transpose()) + E*(r_vb*V)*((r_vc*V).transpose()) + E*(r_ub*V + r_vb*U)*((r_uc*V + r_vc*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_ub*r_uc + stress_tensor[1]*r_vb*r_vc + stress_tensor[2]*(r_ub*r_vc + r_vb*r_uc));
		J_bc *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_bc_prime = E_Prime*(r_ub*U)*((r_uc*U).transpose()) + E_Prime*(r_vb*V)*((r_vc*V).transpose()) + E_Prime*(r_ub*V + r_vb*U)*((r_uc*V + r_vc*U).transpose()) ;
		J_bc_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.b,3*tri.c,J_bc.data(),(-h*h));
		//LHSMatrix_->Add3x3Entry(3*tri.b,3*tri.c,J_bc_prime.data(),(-h));

		Eigen::Matrix3d J_ca = E*(r_uc*U)*((r_ua*U).transpose()) + E*(r_vc*V)*((r_va*V).transpose()) + E*(r_uc*V + r_vc*U)*((r_ua*V + r_va*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_uc*r_ua + stress_tensor[1]*r_vc*r_va + stress_tensor[2]*(r_uc*r_va + r_vc*r_ua));
		J_ca *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_ca_prime = E_Prime*(r_uc*U)*((r_ua*U).transpose()) + E_Prime*(r_vc*V)*((r_va*V).transpose()) + E_Prime*(r_uc*V + r_vc*U)*((r_ua*V + r_va*U).transpose()) ;
		J_ca_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.c,3*tri.a,J_ca.data(),(-h*h));
		//LHSMatrix_->Add3x3Entry(3*tri.c,3*tri.a,J_ca_prime.data(),(-h));

		Eigen::Matrix3d J_cb = E*(r_uc*U)*((r_ub*U).transpose()) + E*(r_vc*V)*((r_vb*V).transpose()) + E*(r_uc*V + r_vc*U)*((r_ub*V + r_vb*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_uc*r_ub + stress_tensor[1]*r_vc*r_vb + stress_tensor[2]*(r_uc*r_vb + r_vc*r_ub));
		J_cb *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_cb_prime = E_Prime*(r_uc*U)*((r_ub*U).transpose()) + E_Prime*(r_vc*V)*((r_vb*V).transpose()) + E_Prime*(r_uc*V + r_vc*U)*((r_ub*V + r_vb*U).transpose()) ;
		J_cb_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.c,3*tri.b,J_cb.data(),(-h*h));
		//LHSMatrix_->Add3x3Entry(3*tri.c,3*tri.b,J_cb_prime.data(),(-h));

		Eigen::Matrix3d J_cc = E*(r_uc*U)*((r_uc*U).transpose()) + E*(r_vc*V)*((r_vc*V).transpose()) + E*(r_uc*V + r_vc*U)*((r_uc*V + r_vc*U).transpose()) 
													 + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_uc*r_uc + stress_tensor[1]*r_vc*r_vc + stress_tensor[2]*(r_uc*r_vc + r_vc*r_uc));
		J_cc *= (-fabs(d)*0.5f);
		Eigen::Matrix3d J_cc_prime = E_Prime*(r_uc*U)*((r_uc*U).transpose()) + E_Prime*(r_vc*V)*((r_vc*V).transpose()) + E_Prime*(r_uc*V + r_vc*U)*((r_uc*V + r_vc*U).transpose()) ;
		J_cc_prime *= (-fabs(d)*0.5f);
		LHSMatrix_->Add3x3Entry(3*tri.c,3*tri.c,J_cc.data(),(-h*h));
		//LHSMatrix_->Add3x3Entry(3*tri.c,3*tri.c,J_cc_prime.data(),(-h));

		//Eigen::Matrix3d J_ji = E*(r_uj*U)*((r_ui*U).transpose()) + E*(r_vj*V)*((r_vi*V).transpose()) + E*(r_uj*V + r_vj*U)*((r_ui*V + r_vi*U).transpose()) 
		//	                     + Eigen::Matrix3d::Identity()*(stress_tensor[0]*r_uj*r_ui + stress_tensor[1]*r_vj*r_vi + stress_tensor[2]*(r_uj*r_vi + r_vj*r_ui));

		//Now we will need to add the velocities in their correct places also 
		Eigen::Vector3d df_dx_a = h*(J_aa * va + J_ab * vb + J_ac * vc);
		Eigen::Vector3d df_dx_b = h*(J_ba * va + J_bb * vb + J_bc * vc);
		Eigen::Vector3d df_dx_c = h*(J_ca * va + J_cb * vb + J_cc * vc);

		RHSVector_[3*tri.a+0] += df_dx_a[0];
		RHSVector_[3*tri.a+1] += df_dx_a[1];
		RHSVector_[3*tri.a+2] += df_dx_a[2];

		RHSVector_[3*tri.b+0] += df_dx_b[0];
		RHSVector_[3*tri.b+1] += df_dx_b[1];
		RHSVector_[3*tri.b+2] += df_dx_b[2];

		RHSVector_[3*tri.c+0] += df_dx_c[0];
		RHSVector_[3*tri.c+1] += df_dx_c[1];
		RHSVector_[3*tri.c+2] += df_dx_c[2];

	}
}

void ImplicitFEMSolver::addGravityComponents(Cloth_Data *cloth) {
	int num_vertices = cloth->getMesh()->get_number_vertices();
	OMP_FOR
	for(int p=0;p<num_vertices;p++)
	{
		float gravForce = cloth->get_vertex_mass(p) * 9.8f;
		RHSVector_[3*p+1] -= gravForce; 
	}

}

void ImplicitFEMSolver::finalizeAndSolve(Cloth_Data *cloth) {
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
	//Solve this using the Gauss Seidel Iterations
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
		//tempSpMatCopy->CheckLinearSystemSolution(resultConstrained,rhsConstrained); 
	}
	//Expand the rows now to fill it 
	double *delV = new double[num_vertices*3];
	InsertRows(num_vertices*3,resultConstrained,delV,numConstrainedVerts_*3,constrainedVerts_);
	//Free the stuff we allocated
	delete(tempSpMatCopy);
	delete[] rhsConstrained;
	delete[] resultConstrained;

	//std::cout << "#Iterations:" << code << "\n";

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

void ImplicitFEMSolver::resetParameters() {
	RHSVector_.setZero();
	LHSMatrix_->ResetToZero();
}