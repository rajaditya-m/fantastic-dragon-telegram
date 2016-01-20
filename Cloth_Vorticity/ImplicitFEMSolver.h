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

#include <Eigen/Sparse>
#include "sparseMatrix.h"

#include "dynamics_solver_interface.h"
#include "cloth_data.h"
#include "CGSolver.h"
#include "insertRows.h"
#include "EigenHelpers.cpp"
//#include "debug_funcs.h"

class ImplicitFEMSolver :
	public Dynamics_Solver_Interface
	{
	public:
		ImplicitFEMSolver(void);
		~ImplicitFEMSolver(void);

		void initialize(Cloth_Data *cloth);
		void initializeSparseMatrixFromOutline(Cloth_Data *cloth);
		
		void resetParameters();

		virtual void advance_time_step(Cloth_Data* cloth);

		void addShearComponents(Cloth_Data* cloth);
	void addShearComponentsCorotational(Cloth_Data* cloth);
	void addGravityComponents(Cloth_Data* cloth);
		void addDampingComponents(Cloth_Data* cloth);
		void addBendingComponents(Cloth_Data* cloth);
		void finalizeAndSolve(Cloth_Data* cloth);


	private:
		Eigen::VectorXd RHSVector_;
		SparseMatrix *LHSMatrix_;
		int lastFrameId_;
		double *massMatrixDiagonal_;

		//Constrain information 
		//@TODO Later on shift to cloth_data 
		int numConstrainedVerts_;
		int *constrainedVerts_;



	};

