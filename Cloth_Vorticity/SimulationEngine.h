#pragma once

#include "cloth_data.h"
#include "ImplicitFEMSolver.h"
#include "ImplicitHyperElasticFEMSolver.h"
#include "ImplicitMassSpringSolver.h"
#include "dynamics_solver_interface.h"
#include "body_data.h"
#include "CollisionEngine.h"


class SimulationEngine
{
public:
	SimulationEngine(Cloth_Data* cloth_data,ImplicitMassSpringSolver* solver,Body_Data* body_data,CollisionEngine* collEngine);
	SimulationEngine(Cloth_Data* cloth_data,ImplicitFEMSolver* solver,Body_Data* body_data,CollisionEngine* collEngine);
	SimulationEngine(Cloth_Data* cloth_data,ImplicitHyperElasticFEMSolver* solver,Body_Data* body_data,CollisionEngine* collEngine);
	~SimulationEngine(void);

	void generate_next_frame();
	void populatePerVertexBuffer();

public:
	

	Cloth_Data* clothData_;
	Dynamics_Solver_Interface* solver_;
	Body_Data* bodyData_;
	CollisionEngine* collEngine_;
};

