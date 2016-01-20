#pragma once

#define OPEN_MP
#ifdef OPEN_MP
#  include <omp.h>
#  ifndef _MSC_VER
#    define OMP_FOR _Pragma("omp parallel for")
#  else
#    define OMP_FOR __pragma(omp parallel for)
#  endif // vs compiler
#else // no omp
#  define OMP_FOR
#endif // #ifndef OPEN_MP

#include "dynamics_solver_interface.h"

#include <Eigen/Dense>

#include "sparseMatrix.h"
#include "cloth_data.h"
#include "geom_funcs.h"
#include "insertRows.h"
#include "CGSolver.h"
#include "PardisoSolver.h"


#define ADDCOMPRESSIONRESISTENCE 
#define NEOHOOKEAN_MATERIAL
#define USE_PARDISO_SOLVER
//#define USE_CG_SOLVER
#define ROWMAJOR6x6TENSOR
#define IMPLICITNEWMARKSOLVER
#define ADDBENDINGFORCES

#define CORNERV1 50
#define CORNERV2 2600

class ImplicitHyperElasticFEMSolver :
	public Dynamics_Solver_Interface
{
public:
	ImplicitHyperElasticFEMSolver(void);
	void computedFdU(Cloth_Data* cloth);
	~ImplicitHyperElasticFEMSolver(void);

	void initialize(Cloth_Data *cloth);
	void loadQAndQVels(Cloth_Data* cloth);
	void updateQAndQ1s(Cloth_Data* cloth);
	void initializeSparseMatrixFromOutline(Cloth_Data *cloth);
	void initializeMassMatrixFromOutline(Cloth_Data* cloth);
	void resetParameters();
	void generateConstantEnergyHessian(Cloth_Data* cloth);
	void addQuadraticBending(Cloth_Data* cloth);
	virtual void advance_time_step(Cloth_Data* cloth);

	int compute6x9TensorIndex(int i, int j, int m, int n);
	int compute6x6TensorIndex(int i, int j, int m, int n);

	void computeDPDF(Eigen::MatrixXd &F,std::vector<double> &gradients, std::vector<double> &hessians,std::vector<double> &DPDF);
	void computeFirstPiolaStress(Eigen::MatrixXd &F, std::vector<double> &gradients, Eigen::MatrixXd &P);
	void computeDGDF(std::vector<double> &DPDF, Eigen::Matrix2d bVec, std::vector<double> &DGDF);
	void computeElementStiffnessMatrix(int el, std::vector<double>& DGDF, std::vector<double>& KELEM);
	void addElementStiffnessMatrixToGlobalStiffnessMatrix(Cloth_Data* cloth, int el, std::vector<double>& KELEM);

	void addShearComponents(Cloth_Data* cloth);
	void addGravityComponents(Cloth_Data* cloth);

	void prepareNumericalEquations(Cloth_Data* cloth);
	bool checkConvergence(Cloth_Data* cloth);
	void solverStep(Cloth_Data* cloth);
	void finalizeStep(Cloth_Data* cloth);

#ifdef IMPLICITNEWMARKSOLVER
	void initializeAlphas(Cloth_Data *cloth);
	void computeInitialAccelerations(Cloth_Data* cloth);
#endif

private:

	SparseMatrix* tangentStiffnessMatrix_;

	int lastFrameId_;
	double *massMatrixDiagonal_;

	//Constrain information 
	//@TODO Later on shift to cloth_data 
	int numConstrainedVerts_;
	int *constrainedVerts_;

	std::vector<double> dDSdU;
	std::vector<double> dFdUs;

	SparseMatrix *rayleighDampingMatrix_;
	SparseMatrix *massMatrix_;

	Eigen::VectorXd internalForce_;
	Eigen::VectorXd externalForce_;

	std::vector<double> oldPosition_;
	
	Eigen::VectorXd q_;
	Eigen::VectorXd q1_;
	Eigen::VectorXd qVel_;
	Eigen::VectorXd qVel1_;
#ifdef IMPLICITNEWMARKSOLVER
	std::vector<double> qAccel_;
	std::vector<double> qAccel1_;
#endif


	std::vector<double> residuals_;
	std::vector<double> constrainedResiduals_;
	Eigen::VectorXd velocities_;

	int numIters_ ;
	double error_;
	double errorQuotient_;

	SparseMatrix *constantEnergyHessian_;
	SparseMatrix *mKeQ_;
	SparseMatrix *mKdQ_;

#ifdef IMPLICITNEWMARKSOLVER
	double alpha1_,alpha2_,alpha3_,alpha4_,alpha5_,alpha6_;
#endif

};

