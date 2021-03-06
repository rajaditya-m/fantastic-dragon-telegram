#pragma once

#include "body_data.h"
#include "triangles.h"
#include "geom_funcs.h"
#include "constitutive_models.h"
#include "edge.h"
#include "triMesh.h"
#include "cloth_properties.h"
#include "RenderObject.h"
#include "global_typedefs.h"
#include "fem_solver.h"

#include <set>
#include <vector>
#include <list>

#include <Eigen\Dense>

#include <qdebug.h>

class Cloth_Data : public RenderObject
{
public:
	Cloth_Data( const char* name,
				const char* file_location_vertex,
				const char* file_location_mesh,
				const char* file_location_texture,
				const char* property_xml,
				const char* material_xml);
	Cloth_Data( const char* name,
				const char* file_location,
				const char* property_xml,
				const char* material_xml);
	~Cloth_Data(void);

	//Accessor methods
	ClothProperties* get_property_obj()								{ return cloth_prop_obj;				}
	Vector7d& get_vertex_distribution(int idx)						{ return vertex_distribution[idx];		}
	float get_vertex_mass(int idx)									{ return vertex_mass[idx];				}
	Eigen::Vector3d get_velocity(int idx)							{ return velocity[idx];					}
	std::vector<Eigen::Vector3d>& get_velocity_vector()				{ return velocity;						}
	Eigen::Vector2d getEdgeWeightedTriangleNormals(int t,int i) { return edgeWeightedTriangleNormals_[3*t+i];}

	Eigen::Matrix2d& getDmInv(int tidx) { return Dm_[tidx];}
	double getW(int tidx) { return W_[tidx];}

	TriMesh* getMesh()          { return mesh_;}

	//Mutator methods
	void set_next_step_pos(Eigen::Vector3d& pos, int idx)			{ next_step_position[idx] = pos ;		}
	void set_next_step_velocity(Eigen::Vector3d& vel, int idx)		{ next_step_velocity[idx] = vel;		}
	void setCurrentStepPositions(Eigen::Vector3d &pos, int idx);
	void setCurrentStepVelocity(Eigen::Vector3d &vel,int idx);
	void setPerVertexVectorBuffer(std::vector<Eigen::Vector3d> &v)	{ perVertexVectorBuffer_ = v;			}
	void finalize_velocity_position();

	void resetParameters();

	//Rendering functions 
	virtual void render(int frameId);

	//addition functions added by Guowei
	void setPositions(float* p){
		for (int i = 0; i < next_step_position.size(); i++){
			next_step_position[i][0] = p[3*i+0];
			next_step_position[i][1] = p[3*i+1];
			next_step_position[i][2] = p[3*i+2];
		}
	}
	void setVelocities(float* v){
		for (int i = 0; i < next_step_velocity.size(); i++){
			next_step_velocity[i][0] = v[3*i+0];
			next_step_velocity[i][1] = v[3*i+1];
			next_step_velocity[i][2] = v[3*i+2];
		}
	}

	float* getPositions(){
		p_.resize(next_step_position.size()*3);
		for (int i = 0; i < next_step_position.size(); i++){
			p_[3*i+0] = next_step_position[i][0];
			p_[3*i+1] = next_step_position[i][1];
			p_[3*i+2] = next_step_position[i][2];
		}

		return &p_[0];
	}
	float* getVelocities(){
		v_.resize(next_step_velocity.size()*3);
		for (int i = 0; i < next_step_velocity.size(); i++){
			v_[3*i+0] = next_step_velocity[i][0];
			v_[3*i+1] = next_step_velocity[i][1];
			v_[3*i+2] = next_step_velocity[i][2];
		}
		return &v_[0];
	}

private:
	//Contains the external fiel read params to the simulator
	ClothProperties* cloth_prop_obj;
	
	//Shear Distribution Precomputation
	std::vector<Vector7d> vertex_distribution;
	
	//Velocity and mass parameters
	std::vector<Eigen::Vector3d> velocity;
	std::vector<float> vertex_mass;

	//this contains the buffers that contain the next positions and velocities 
	std::vector<Eigen::Vector3d> next_step_velocity;
	std::vector<Eigen::Vector3d> next_step_position;

	std::vector<Eigen::Matrix2d> Dm_;
	std::vector<double> W_;
	std::vector<Eigen::Vector2d> edgeWeightedTriangleNormals_;

	//Additional Rendering Buffer (will be populated by SIM ENGINE)
	std::vector<Eigen::Vector3d> perVertexVectorBuffer_;


	std::vector<float> p_;
	std::vector<float> v_;
};

