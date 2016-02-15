#include "SimulationEngine.h"


SimulationEngine::SimulationEngine(Cloth_Data* cloth_data,ImplicitFEMSolver* solver,Body_Data* body_data,CollisionEngine* collEngine)
{
	clothData_ = cloth_data;
	solver_ = solver;
	solver_->initialize(clothData_);
	bodyData_ = body_data;
	collEngine_ = collEngine;
}

SimulationEngine::SimulationEngine(Cloth_Data* cloth_data,ImplicitMassSpringSolver* solver,Body_Data* body_data,CollisionEngine* collEngine)
{
	clothData_ = cloth_data;
	solver_ = solver;
	solver_->initialize(clothData_);
	bodyData_ = body_data;
	collEngine_ = collEngine;
}

SimulationEngine::SimulationEngine(Cloth_Data* cloth_data,ImplicitHyperElasticFEMSolver* solver,Body_Data* body_data,CollisionEngine* collEngine)
{
	clothData_ = cloth_data;
	solver_ = solver;
	solver_->initialize(clothData_);
	bodyData_ = body_data;
	collEngine_ = collEngine;
}


SimulationEngine::~SimulationEngine(void)
{
}

void SimulationEngine::generate_next_frame()
{
	//Advance the time step to generate the next set of collision as per internal stuff
	solver_->advance_time_step(clothData_);
	
	//Then apply the collision engine
	collEngine_->resolveClothBodyCollision(clothData_,bodyData_);
	
	//Then resolve self collisions 
	//collEngine_->resolve_cloth_cloth_collisions(cloth_data);
	
	//Finalize the velocity and positions
	clothData_->finalize_velocity_position();

	//Control the flow of parameter to the cloth engine depending upon the type of rendering 
	populatePerVertexBuffer();
}

//We will implement this function later on as per the need. For now this does nothing. 
void SimulationEngine::populatePerVertexBuffer()
{
	/*switch(clothData_->getRenderMode())
	{
		case(HMAP_ACCLRN) : clothData_->setPerVertexVectorBuffer(solver_->getAcceleration()); break;
		case(HMAP_BENDING_FORCE) : clothData_->setPerVertexVectorBuffer(solver_->getBendingForce()); break;
		case(HMAP_SHEAR_FORCE) : clothData_->setPerVertexVectorBuffer(solver_->getShearForce()); break;
		case(HMAP_DAMPING_FORCE) : clothData_->setPerVertexVectorBuffer(solver_->getDampingForce()); break;
	}*/
}
