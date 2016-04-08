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
#include <set>
#include <fstream>
#include <iostream>

#include "Eigen\Dense"

#include <QGLWidget>

#include "VolumeInformation.h"
#include "triangles.h"
#include "edge.h"
#include "global_typedefs.h"
#include "struct_decls.h"
#include "mesh_material.h"
#include "grid.h"
#include "math_func.h"


class TriMesh
{

public:
	//Constructor for Clothing (with UV)
	TriMesh(const char* meshName,
			const char* file_location_vertex,
			const char* file_location_mesh,
			const char* file_location_texture,
			const char* material_xml,
			float pdfs);
	TriMesh(const char* meshName,
			const char* file_location_obj,
			const char* material_xml,
			float pdfs);
	//Constructor for generic static mesh (without UV)
	TriMesh(const char* meshName,
			const char* file_location_vertex,
			const char* file_location_mesh,
			const char* material_xml,
			float pdfs);
	//Constructor for dynamic sequence of clothing 
	TriMesh(const char* meshName,
			const char* file_location_gen,
			const char* material_xml,
			int num_frames,
			float pdfs);
	TriMesh();
	TriMesh(float pdfs);
	~TriMesh(void);

	//This is used for non constructor initialization of the mesh 
	void initializeDynamicMesh(std::vector<Eigen::Vector3d> &pvec, std::vector<int> triangles, const char* meshName, const char* material_xml);

	//More important functions 
	int get_number_vertices() const 									{ return point_data[0].size() ;			}
	int get_number_triangles() const									{ return mesh_data.size();				}
	int get_number_edges() const										{ return edge_list.size();				}
	float get_average_edge_length()	const								{ return average_edge_length_;			}
	float get_longest_edge_length()	const								{ return longest_edge_length_;			}
	Pair_Vec3d& get_super_bounding_box()								{ return super_bounding_box;			}
	std::string getMeshName() const										{ return name_;							}
	std::vector<Eigen::Vector3d>& getPositionVector(int frameId) { return point_data[frameId];}

	//Accessor functions for this class 
	Triangles get_triangle(int tri_idx)	const							{ return mesh_data[tri_idx];			}
	bool is_part_of_mesh(int vert_idx)	const							{ return part_of_mesh[vert_idx];		}
	Eigen::Vector2d get_uv_data(int idx) const							{ return texture_data[idx];				}
	Eigen::Vector2d get_texture_data(int idx) const						{ return texture_data[idx];				}
	std::set<Edge>& get_edge_list()										{ return edge_list;						}
	std::vector<Edge>& getEdgeVector()                { return edgeVector_;}
	int getEdgeVectorSize()                           { return edgeVectorSize_;}
	Grid* get_grid_obj() const											{ return voxel_grid;					}
	bool isStatic() const												{ return static_;						}

	//Other accessor functions defined elsewhere
	int get_number_frames()	const ;										
	Eigen::Vector3d get_point_data(int idx, int frame = 0)	const ; 
	Eigen::Vector3d get_normal_data(int idx, int frame = 0)	const ;

	//Mutator functions for this class 
	void set_point_data(const Eigen::Vector3d p,int idx,int frame=0)	{ point_data[frame][idx] = p;			}
	void add_point_data_frame(std::vector<Eigen::Vector3d> &vec)		{ point_data.push_back(vec);			}
	void setStatic(bool st)												{ static_ = st;							}

	//Other functions
	void generate_frame_normals(int frame_idx) ;
	void add_frame_data(std::vector<Eigen::Vector3d> &vec);
	void resetMeshDefaults();
	int performMeshSelectionRayTest(double* rayStart, double* rayEnd, double *clickedWorldPos, double* selectedPos, int frameId);
	void saveFrameAsObj(const char* fileSuffix, int frameID);
	void saveAsOBJ(const char* fileSuffix);

	//Drawing Functions 
	void draw_mesh(int idx);
	void draw_wireframe(int idx);
	void draw_bounding_box(int idx);
	void draw_super_bounding_box();
	void draw_vector_valued_heatmap(std::vector<Eigen::Vector3d> vector_vertex_values,int frame_idx);
	
	void update_voxelization();					//This takes every frame data that is present and updates the voxelization
	void voxelize_frame(int frame_id);			//This takes one frame id and voxelizes it 


	//addition functions added by Guowei
	std::vector<int> getTriangles(){
		std::vector<int> t;
		t.clear();
		for (int i = 0; i < mesh_data.size(); i++){
			t.push_back(mesh_data[i].a);
			t.push_back(mesh_data[i].b);
			t.push_back(mesh_data[i].c);
		}
		return t;
	}
	std::vector<int> getEdges(){
		std::vector<int> e;
		e.clear();
		for (int i = 0; i < edgeVector_.size(); i++){
			e.push_back(edgeVector_[i].start);
			e.push_back(edgeVector_[i].end);
		}
		return e;
	}
	std::vector<int> getET(){
		std::vector<int> et;
		et.resize(2 * edgeVector_.size());
		for (int i = 0; i < edgeVector_.size(); i++){		
			et[2*i+0] = edgeVector_[i].adj_tri_1;
			et[2*i+0] = edgeVector_[i].adj_tri_2 ;
		}
		return et;
	}
	float* getPositions(){
		int currentFrame = point_data.size()-1;
		p_.resize(point_data[currentFrame].size()*3);
		for (int i = 0; i < point_data[currentFrame].size(); i++){
			p_[3*i+0] = (point_data[currentFrame][i][0]);
			p_[3*i+1] = (point_data[currentFrame][i][1]);
			p_[3*i+2] = (point_data[currentFrame][i][2]);
		}

		return &p_[0];
	}
	void* getOldPositions(){
		int previousFrame = point_data.size()-2;
		pOld_.resize(point_data[previousFrame].size()*3);
		for (int i = 0; i < point_data[previousFrame].size(); i++){
			pOld_[3*i+0] = (point_data[previousFrame][i][0]);
			pOld_[3*i+1] = (point_data[previousFrame][i][1]);
			pOld_[3*i+2] = (point_data[previousFrame][i][2]);
		}

		return &pOld_[0];	
	}
	void setPositions(float* p){
		int currentFrame = point_data.size()-1;
		for (int i = 0; i < point_data[currentFrame].size(); i++){
			point_data[currentFrame][i][0] = p[3*i+0];
			point_data[currentFrame][i][1] = p[3*i+1];
			point_data[currentFrame][i][2] = p[3*i+2];
		}
	}


private:
	//Mesh should have a name 
	std::string name_;

	//Data containign frame by frame stuff 
	Vector2d_Vec3d point_data;
	Vector2d_Vec3d normal_data;
	Vector_PairVec3d bounding_box;

	//Mesh specific Data
	std::vector<Triangles> mesh_data;
	std::vector<bool> part_of_mesh;
	std::set<Edge> edge_list;
	std::vector<Edge> edgeVector_;
	int edgeVectorSize_;
	Vector_Vec2d texture_data;

	const float CONST_SCALING;

	//Useful Statistics will be updated as and when needed
	float average_edge_length_;
	float longest_edge_length_;

	//Super bounding box measurements
	Pair_Vec3d super_bounding_box;
	
	//Other data
	bool has_uv_coods_;
	bool static_;

	//This should have some material properties 
	Mesh_Material* material;
	//@TODO : Support for texture rendering 

	//Voxelization wrt some grid boundaries
	bool voxelized_ ;
	Grid* voxel_grid;
	//Voxel to primitive DS
	std::vector< std::vector< std::vector <int> > > voxel_to_triangle;
	std::vector< std::vector< std::vector <int> > > voxel_to_point;
	std::vector< std::vector< std::vector <int> > > voxel_to_edge;
	//Primitive to voxel DS
	std::vector< std::vector< std::vector <Voxel> > > triangle_to_voxel;
	std::vector< std::vector< std::vector <Voxel> > > edge_to_voxel;
	std::vector< std::vector< Voxel > > point_to_voxel;

	std::vector<float> p_;
	std::vector<float> pOld_;

	std::vector<VolumeInformation> volInfo_;
};

