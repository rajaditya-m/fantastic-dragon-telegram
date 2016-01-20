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

#include <vector>
#include <algorithm>
#include <Eigen/Dense>

#include "debug_funcs.h"
#include "dynamics_solver_interface.h"
#include "cloth_data.h"	

class FEM_Solver :
	public Dynamics_Solver_Interface
{
public:
	FEM_Solver(void);
	~FEM_Solver(void);
	
	//The main funtion called by the external class
	virtual void advance_time_step(Cloth_Data* cloth);

	//Reset Function 
	void resetParameters();

	//Initialize vectors
	void initialize(Cloth_Data* data);

	//Force computation functions
	void compute_shear_forces(Cloth_Data* cloth);
	void compute_bending_forces(Cloth_Data* cloth);
	void compute_gravity_forces(Cloth_Data* cloth);
	void computeDampingForce(Cloth_Data* cloth);

	//functions to get the whole vectors 
	std::vector<Eigen::Vector3d>& getShearForce()				{ return shear_forces;	}
	std::vector<Eigen::Vector3d>& getBendingForce()				{ return bend_forces;	}
	std::vector<Eigen::Vector3d>& getAcceleration()				{ return acceleration_;	}
	std::vector<Eigen::Vector3d>& getDampingForce()				{ return dampingForce_;	}

	//Calculation functions 
	void get_next_velocity_positions(Cloth_Data* cloth);
private:
	std::vector<Eigen::Vector3d> shear_forces;
	std::vector<Eigen::Vector3d> bend_forces;
	std::vector<Eigen::Vector3d> gravity_forces;
	std::vector<Eigen::Vector3d> acceleration_;
	std::vector<Eigen::Vector3d> dampingForce_;
	int last_frame_id;
};

