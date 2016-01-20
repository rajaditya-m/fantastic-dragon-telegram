#pragma once 

#include "global_typedefs.h"
#include <Eigen\Dense>
#include <vector>


class Grid
{
public:

	Grid()
	{
		grid_corner = Eigen::Vector3d::Zero();
		grid_corner_4f = Eigen::Vector4f::Zero();

		grid_resolution = inv_grid_resolution = 0.0;
	}

	//Mutator methods
	void set_grid_parameters(Eigen::Vector3d corner,float gr)
	{
		grid_corner = corner;

		grid_corner_4f = Eigen::Vector4f(grid_corner.x(),grid_corner.y(),grid_corner.z(),1.0f);

		grid_resolution  = gr;
		inv_grid_resolution = 1.0f/grid_resolution;
	}

	//This is the main function overloaded with voxels 
	Voxel get_voxel(Eigen::Vector3d point)
	{
		Eigen::Vector3d difference = point - grid_corner;
		float x_id = difference.x()*inv_grid_resolution;
		int x = (int)ceil(x_id);
		float y_id = difference.y()*inv_grid_resolution;
		int y = (int)ceil(y_id);
		float z_id = difference.z()*inv_grid_resolution;
		int z = (int)ceil(z_id);
		return Eigen::Vector3d(x,y,z);
	}

	Voxel get_voxel(Eigen::Vector4f point)
	{
		Eigen::Vector4f difference = point - grid_corner_4f;
		float x_id = difference.x()*inv_grid_resolution;
		int x = (int)ceil(x_id);
		float y_id = difference.y()*inv_grid_resolution;
		int y = (int)ceil(y_id);
		float z_id = difference.z()*inv_grid_resolution;
		int z = (int)ceil(z_id);
		return Eigen::Vector3d(x,y,z);
	}

	std::vector<Voxel> get_voxel(Eigen::Vector3d edge_vertex_1,Eigen::Vector3d edge_vertex_2)
	{
		std::vector<Voxel> edge_vector;
		Voxel edge_1 = this->get_voxel(edge_vertex_1);
		Voxel edge_2 = this->get_voxel(edge_vertex_2);
		//Now we have to classsify into cases 
		Voxel diff = edge_1 - edge_2;
		if(diff.x() > 1 || diff.x() < -1 || diff.y() > 1 || diff.y() < -1 ||diff.z() > 1 || diff.z() < -1)
		{
			std::cout << "[ERROR] Grid resolution is too coarse. Recommend reparametization.\n";
		}
		edge_vector.push_back(edge_1);
		edge_vector.push_back(edge_2);

		return edge_vector;
	}

	std::vector<Voxel> get_voxel(Eigen::Vector4f edge_vertex_1,Eigen::Vector4f edge_vertex_2)
	{
		std::vector<Voxel> edge_vector;
		Voxel edge_1 = this->get_voxel(edge_vertex_1);
		Voxel edge_2 = this->get_voxel(edge_vertex_2);
		//Now we have to classsify into cases 
		Voxel diff = edge_1 - edge_2;
		if(diff.x() > 1 || diff.x() < -1 || diff.y() > 1 || diff.y() < -1 ||diff.z() > 1 || diff.z() < -1)
		{
			std::cout << "[ERROR] Grid resolution is too coarse. Recommend reparametization.\n";
		}
		edge_vector.push_back(edge_1);
		edge_vector.push_back(edge_2);

		return edge_vector;
	}

	std::vector<Voxel> get_voxel(Eigen::Vector3d tri_vertex_1,Eigen::Vector3d tri_vertex_2,Eigen::Vector3d tri_vertex_3)
	{
		std::vector<Voxel> tri_vector;
		Voxel tri_1 = this->get_voxel(tri_vertex_1);
		Voxel tri_2 = this->get_voxel(tri_vertex_2);
		Voxel tri_3 = this->get_voxel(tri_vertex_3);

		Voxel diff_1 = tri_1 - tri_2;
		if(diff_1.x() > 1 || diff_1.x() < -1 || diff_1.y() > 1 || diff_1.y() < -1 ||diff_1.z() > 1 || diff_1.z() < -1)
		{
			std::cout << "[ERROR] Grid resolution is too coarse. Recommend reparametization.\n";
		}
		Voxel diff_2 = tri_2 - tri_3;
		if(diff_2.x() > 1 || diff_2.x() < -1 || diff_2.y() > 1 || diff_2.y() < -1 ||diff_2.z() > 1 || diff_2.z() < -1)
		{
			std::cout << "[ERROR] Grid resolution is too coarse. Recommend reparametization.\n";
		}
		Voxel diff_3 = tri_1 - tri_3;
		if(diff_3.x() > 1 || diff_3.x() < -1 || diff_3.y() > 1 || diff_3.y() < -1 ||diff_3.z() > 1 || diff_3.z() < -1)
		{
			std::cout << "[ERROR] Grid resolution is too coarse. Recommend reparametization.\n";
		}

		tri_vector.push_back(tri_1);
		tri_vector.push_back(tri_2);
		tri_vector.push_back(tri_3);

		return tri_vector;
	}
	
	std::vector<Voxel> get_voxel(Eigen::Vector4f tri_vertex_1,Eigen::Vector4f tri_vertex_2,Eigen::Vector4f tri_vertex_3)
	{
		std::vector<Voxel> tri_vector;
		Voxel tri_1 = this->get_voxel(tri_vertex_1);
		Voxel tri_2 = this->get_voxel(tri_vertex_2);
		Voxel tri_3 = this->get_voxel(tri_vertex_3);

		Voxel diff_1 = tri_1 - tri_2;
		if(diff_1.x() > 1 || diff_1.x() < -1 || diff_1.y() > 1 || diff_1.y() < -1 ||diff_1.z() > 1 || diff_1.z() < -1)
		{
			std::cout << "[ERROR] Grid resolution is too coarse. Recommend reparametization.\n";
		}
		Voxel diff_2 = tri_2 - tri_3;
		if(diff_2.x() > 1 || diff_2.x() < -1 || diff_2.y() > 1 || diff_2.y() < -1 ||diff_2.z() > 1 || diff_2.z() < -1)
		{
			std::cout << "[ERROR] Grid resolution is too coarse. Recommend reparametization.\n";
		}
		Voxel diff_3 = tri_1 - tri_3;
		if(diff_3.x() > 1 || diff_3.x() < -1 || diff_3.y() > 1 || diff_3.y() < -1 ||diff_3.z() > 1 || diff_3.z() < -1)
		{
			std::cout << "[ERROR] Grid resolution is too coarse. Recommend reparametization.\n";
		}

		tri_vector.push_back(tri_1);
		tri_vector.push_back(tri_2);
		tri_vector.push_back(tri_3);

		return tri_vector;
	}

private:
	Eigen::Vector3d grid_corner;
	Eigen::Vector4f grid_corner_4f;
	float grid_resolution;
	float inv_grid_resolution;
};