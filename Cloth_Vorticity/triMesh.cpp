#include "triMesh.h"
#include "geom_funcs.h"
#include "debug_funcs.h"

//Default constructor
TriMesh::TriMesh()
	:CONST_SCALING(1.0)
{
	
}

TriMesh::TriMesh(float pdfs)
	:CONST_SCALING(pdfs)
{
	
}

//This is the version that reads a single frame with a single filename with texture coordinates 
TriMesh::TriMesh(const char* meshName,
				 const char* file_location_vertex,
				 const char* file_location_mesh,
				 const char* file_location_texture,
				 const char* material_xml,
				 float pdfs)
	: CONST_SCALING(pdfs)
{
	name_ = meshName ;
	has_uv_coods_ = true;
	static_ = false;
	voxelized_ = false;
	voxel_grid = new Grid();

	std::string line;
	std::string first_iden = "";

	//Lets read the vertex binaries first
	Point_Struct temp;
	std::ifstream in_file(file_location_vertex,std::ios::in|std::ios::binary);
	in_file.read(reinterpret_cast<char* >(&temp),sizeof(Point_Struct));
	std::vector<Eigen::Vector3d> pointvec_f0;
	Eigen::Vector3d min_bb(99999.0f,99999.0f,99999.0f);
	Eigen::Vector3d max_bb(-99999.0f,-99999.0f,-99999.0f);
	while(in_file && !in_file.eof())
	{
		Eigen::Vector3d temp_point(temp.x*CONST_SCALING,temp.y*CONST_SCALING,temp.z*CONST_SCALING);
		pointvec_f0.push_back(temp_point);
		in_file.read(reinterpret_cast<char* >(&temp),sizeof(Point_Struct));
		if(temp_point.x()>max_bb.x())
			max_bb(0) = temp_point.x();
		if(temp_point.y()>max_bb.y())
			max_bb(1) = temp_point.y();
		if(temp_point.z()>max_bb.z())
			max_bb(2) = temp_point.z();
		if(temp_point.x()<min_bb.x())
			min_bb(0) = temp_point.x();
		if(temp_point.y()<min_bb.y())
			min_bb(1) = temp_point.y();
		if(temp_point.z()<min_bb.z())
			min_bb(2) = temp_point.z();
	}
	in_file.close();
	int num_vertices = pointvec_f0.size();
	part_of_mesh.resize(num_vertices,false);

	//Added the point data read so far
	point_data.push_back(pointvec_f0);
	
	//Add the bounding box and super bounding box (for a single frame this is the bounding box)
	bounding_box.push_back(std::make_pair(min_bb,max_bb));
	super_bounding_box = std::make_pair(min_bb,max_bb);

	//Now lets read the texture binaries 
	Texture_Struct temp_tx;
	in_file.open(file_location_texture,std::ios::in|std::ios::binary);
	in_file.read(reinterpret_cast<char* >(&temp_tx),sizeof(Texture_Struct));
	while(in_file && !in_file.eof())
	{
		Eigen::Vector2d temp_read_texture(temp_tx.u,temp_tx.v);
		texture_data.push_back(temp_read_texture);
		in_file.read(reinterpret_cast<char* >(&temp_tx),sizeof(Texture_Struct));
	}
	in_file.close();

	//Now lets read the triangle binaries
	average_edge_length_ = 0.0f;
	int num_counted_edges = 0;
	longest_edge_length_ = 0.0f;
	Triangle_Texture_Struct temp_t;
	in_file.open(file_location_mesh,std::ios::in|std::ios::binary);
	in_file.read(reinterpret_cast<char* >(&temp_t),sizeof(Triangle_Texture_Struct));

	std::vector<Eigen::Vector3d> normalvec_f0(num_vertices,Eigen::Vector3d::Zero());
	int restriction_counter = 0;
	while(in_file && !in_file.eof())
	{
		int vp1 = temp_t.a;
		int vp2 = temp_t.b;
		int vp3 = temp_t.c;

		int t1 = temp_t.at;
		int t2 = temp_t.bt;
		int t3 = temp_t.ct;


		Eigen::Vector3d v_ab = pointvec_f0[vp2-1] - pointvec_f0[vp1-1];
		Eigen::Vector3d v_ac = pointvec_f0[vp3-1] - pointvec_f0[vp1-1];
		Eigen::Vector3d v_bc = pointvec_f0[vp3-1] - pointvec_f0[vp2-1];

		float ab = v_ab.norm();
		float bc = v_bc.norm();
		float ac = v_ac.norm();

		if(ab>longest_edge_length_)
			longest_edge_length_ = ab;

		if(bc>longest_edge_length_)
			longest_edge_length_ = bc;
		
		if(ac>longest_edge_length_)
			longest_edge_length_ = ac;

		average_edge_length_ += (ab+bc+ac);
		num_counted_edges += 3;

		Eigen::Vector3d tri_norm_abc = Eigen::Vector3d::Zero();
		tri_norm_abc = v_ab.cross(v_ac);
		tri_norm_abc.normalize();

		normalvec_f0[vp1-1] += tri_norm_abc;
		normalvec_f0[vp2-1] += tri_norm_abc;
		normalvec_f0[vp3-1] += tri_norm_abc;

		Triangles temp_tri;
		temp_tri.add_textures(t1-1,t2-1,t3-1);
		temp_tri.add(vp1-1,vp2-1,vp3-1);
		temp_tri.addTriangleNormal(tri_norm_abc);
		mesh_data.push_back(temp_tri);

		part_of_mesh[vp1-1] = true;
		part_of_mesh[vp2-1] = true;
		part_of_mesh[vp3-1] = true;

		in_file.read(reinterpret_cast<char* >(&temp_t),sizeof(Triangle_Texture_Struct));
	}
	for(int p=0;p<pointvec_f0.size();p++)
	{
		if(!part_of_mesh[p])
			continue;
		normalvec_f0[p].normalize();
	}
	normal_data.push_back(normalvec_f0);

	//Calculate the values of the statistics
	average_edge_length_ /= num_counted_edges;

	//Populate the edge data-structure
	for(int t=0;t<mesh_data.size();t++)
	{
		int a = mesh_data[t].a;
		int b = mesh_data[t].b;
		int c = mesh_data[t].c;

		Edge e1;
		e1.add_edges(a,b);
		Edge e2;
		e2.add_edges(b,c);
		Edge e3;
		e3.add_edges(c,a);

		std::set<Edge>::iterator it = edge_list.find(e1);
		if(it==edge_list.end())
		{
			e1.add_triangle(t);
			edge_list.insert(e1);
		}
		else
		{
			const Edge &re = *it;
			Edge &e = const_cast<Edge&>(re);
			bool status = e.add_triangle(t);
			if(!status)
				std::cout <<  "[WARNING] Non-Manifold mesh detected.\n";
		}

		it = edge_list.find(e2);
		if(it==edge_list.end())
		{
			e2.add_triangle(t);
			edge_list.insert(e2);
		}
		else
		{
			const Edge &re = *it;
			Edge &e = const_cast<Edge&>(re);
			bool status = e.add_triangle(t);
			if(!status)
				std::cout << "[WARNING] Non-Manifold mesh detected.\n";
		}

		it = edge_list.find(e3);
		if(it==edge_list.end())
		{
			e3.add_triangle(t);
			edge_list.insert(e3);
		}
		else
		{
			const Edge &re = *it;
			Edge &e = const_cast<Edge&>(re);
			bool status = e.add_triangle(t);
			if(!status)
				std::cout << "[WARNING] Non-Manifold mesh detected.\n";
		}
	}

	//Now assign an edge index this is needed for the voxelization 
	std::set<Edge>::iterator it;
	int edge_counter = 0;
	for(it = edge_list.begin(); it != edge_list.end() ; it++)
	{
		const Edge &ce = *it;
		Edge &e = const_cast<Edge&>(ce);
		e.set_index(edge_counter);
		edge_counter++;
	}


	//This is the material property 
	material = new Mesh_Material(material_xml);


}

//This is the version that reads a single frame with a single filename with texture coordinates 
TriMesh::TriMesh(const char* meshName,
								 const char* file_location_obj,
								 const char* material_xml,
								 float pdfs):
	CONST_SCALING(pdfs){
			name_ = meshName ;
			has_uv_coods_ = true;
			static_ = false;
			voxelized_ = false;
			voxel_grid = new Grid();
			int num_vertices;
			std::vector<Eigen::Vector3d> pointvec_f0;
			std::vector<Eigen::Vector3d> normalvec_f0;
			Eigen::Vector3d min_bb(99999.0f,99999.0f,99999.0f);
			Eigen::Vector3d max_bb(-99999.0f,-99999.0f,-99999.0f);
			std::ifstream fileIn(file_location_obj);
			std::string line;
			std::stringstream ss;
			std::string firstIden;
			bool onceFlag = true;
			double x,y,z;
			int a,b,c;

			static_ = false;

			average_edge_length_ = 0.0f;
			int num_counted_edges = 0;
			longest_edge_length_ = 0.0f;
			while(std::getline(fileIn,line))
			{
			ss.clear();
			ss.str("");
			ss << line;
			ss >> firstIden;
			if(!firstIden.compare("v"))
			{
				ss >> x >> y >> z;
				Eigen::Vector3d temp_point(x*CONST_SCALING,y*CONST_SCALING,z*CONST_SCALING);
				pointvec_f0.push_back(temp_point);
				if(temp_point.x()>max_bb.x())
					max_bb(0) = temp_point.x();
				if(temp_point.y()>max_bb.y())
					max_bb(1) = temp_point.y();
				if(temp_point.z()>max_bb.z())
					max_bb(2) = temp_point.z();
				if(temp_point.x()<min_bb.x())
					min_bb(0) = temp_point.x();
				if(temp_point.y()<min_bb.y())
					min_bb(1) = temp_point.y();
				if(temp_point.z()<min_bb.z())
					min_bb(2) = temp_point.z();

				//Keep on updating this specific number everytme 
				num_vertices = pointvec_f0.size();
			}
			else if(!firstIden.compare("f"))
			{
				if(onceFlag) {
					normalvec_f0.resize(num_vertices,Eigen::Vector3d::Zero());
					part_of_mesh.resize(num_vertices,false);
					onceFlag = false;
				}
				std::string sa,sb,sc;
				ss >> sa >> sb >> sc;
				//proces the three numbers 
				std::size_t firstSlash = sa.find('/');
				std::string firstN = sa.substr(0,firstSlash);
				std::string rem = sa.substr(firstSlash+1);
				std::size_t secondSlash =rem.find('/');
				std::string secondN = rem.substr(0,secondSlash);
				std::string thirdN = rem.substr(secondSlash+1);
				
				int vp1 = atoi(firstN.c_str());
				int t1 = atoi(secondN.c_str());

				firstSlash = sb.find('/');
				firstN = sb.substr(0,firstSlash);
				rem = sb.substr(firstSlash+1);
				secondSlash =rem.find('/');
				secondN = rem.substr(0,secondSlash);
				thirdN = rem.substr(secondSlash+1);

				int vp2 = atoi(firstN.c_str());
				int t2 = atoi(secondN.c_str());

				firstSlash = sc.find('/');
				firstN = sc.substr(0,firstSlash);
				rem = sc.substr(firstSlash+1);
				secondSlash =rem.find('/');
				secondN = rem.substr(0,secondSlash);
				thirdN = rem.substr(secondSlash+1);

				int vp3 = atoi(firstN.c_str());
				int t3 = atoi(secondN.c_str());

				Eigen::Vector3d v_ab = pointvec_f0[vp2-1] - pointvec_f0[vp1-1];
				Eigen::Vector3d v_ac = pointvec_f0[vp3-1] - pointvec_f0[vp1-1];
				Eigen::Vector3d v_bc = pointvec_f0[vp3-1] - pointvec_f0[vp2-1];

				float ab = v_ab.norm();
				float bc = v_bc.norm();
				float ac = v_ac.norm();

				if(ab>longest_edge_length_)
					longest_edge_length_ = ab;

				if(bc>longest_edge_length_)
					longest_edge_length_ = bc;
		
				if(ac>longest_edge_length_)
					longest_edge_length_ = ac;

				average_edge_length_ += (ab+bc+ac);
				num_counted_edges += 3;

				Eigen::Vector3d tri_norm_abc = Eigen::Vector3d::Zero();
				tri_norm_abc = v_ab.cross(v_ac);
				tri_norm_abc.normalize();

				normalvec_f0[vp1-1] += tri_norm_abc;
				normalvec_f0[vp2-1] += tri_norm_abc;
				normalvec_f0[vp3-1] += tri_norm_abc;

				Triangles temp_tri;
				temp_tri.add_textures(t1-1,t2-1,t3-1);
				temp_tri.add(vp1-1,vp2-1,vp3-1);
				temp_tri.addTriangleNormal(tri_norm_abc);
				mesh_data.push_back(temp_tri);

				part_of_mesh[vp1-1] = true;
				part_of_mesh[vp2-1] = true;
				part_of_mesh[vp3-1] = true;

			}
			else if(!firstIden.compare("vt")) {
				ss >> x >> y;
				Eigen::Vector2d temp_read_texture(x,y);
				texture_data.push_back(temp_read_texture);
			}
		}

		point_data.push_back(pointvec_f0);
		volInfo_.push_back(VolumeInformation(point_data[0]));

		//Add the bounding box and super bounding box (for a single frame this is the bounding box)
		bounding_box.push_back(std::make_pair(min_bb,max_bb));
		super_bounding_box = std::make_pair(min_bb,max_bb);

		for(int p=0;p<pointvec_f0.size();p++)
		{
			if(!part_of_mesh[p])
				continue;
			normalvec_f0[p].normalize();
		}
		normal_data.push_back(normalvec_f0);

		//Calculate the values of the statistics
		average_edge_length_ /= num_counted_edges;

		//Populate the edge data-structure
		for(int t=0;t<mesh_data.size();t++)
		{
			int a = mesh_data[t].a;
			int b = mesh_data[t].b;
			int c = mesh_data[t].c;

			Edge e1;
			e1.add_edges(a,b);
			Edge e2;
			e2.add_edges(b,c);
			Edge e3;
			e3.add_edges(c,a);

			std::set<Edge>::iterator it = edge_list.find(e1);
			if(it==edge_list.end())
			{
				e1.add_triangle(t);
				edge_list.insert(e1);
			}
			else
			{
				const Edge &re = *it;
				Edge &e = const_cast<Edge&>(re);
				bool status = e.add_triangle(t);
				if(!status)
					std::cout <<  "[WARNING] Non-Manifold mesh detected.\n";
			}

			it = edge_list.find(e2);
			if(it==edge_list.end())
			{
				e2.add_triangle(t);
				edge_list.insert(e2);
			}
			else
			{
				const Edge &re = *it;
				Edge &e = const_cast<Edge&>(re);
				bool status = e.add_triangle(t);
				if(!status)
					std::cout << "[WARNING] Non-Manifold mesh detected.\n";
			}

			it = edge_list.find(e3);
			if(it==edge_list.end())
			{
				e3.add_triangle(t);
				edge_list.insert(e3);
			}
			else
			{
				const Edge &re = *it;
				Edge &e = const_cast<Edge&>(re);
				bool status = e.add_triangle(t);
				if(!status)
					std::cout << "[WARNING] Non-Manifold mesh detected.\n";
			}
		}

		//Now assign an edge index this is needed for the voxelization 
		std::set<Edge>::iterator it;
		int edge_counter = 0;
		for(it = edge_list.begin(); it != edge_list.end() ; it++)
		{
			const Edge &ce = *it;
			Edge &e = const_cast<Edge&>(ce);
			e.set_index(edge_counter);
			edge_counter++;
		}
		//This is the material property 
		material = new Mesh_Material(material_xml);

		//Last modification 
		edgeVector_ = std::vector<Edge>(edge_list.begin(),edge_list.end());
		edgeVectorSize_ = edgeVector_.size();
}

//Same version without the texture coordinates 
TriMesh::TriMesh(const char* meshName,
				 const char* file_location_vertex,
				 const char* file_location_mesh,
				 const char* material_xml,
				 float pdfs)
		: CONST_SCALING(pdfs)
{
	name_ = meshName ;
	has_uv_coods_ = false;
	voxelized_ = false;
	static_ = false;
	voxel_grid = new Grid();

	std::string line;
	std::string first_iden = "";

	//Lets read the vertex binaries first
	Point_Struct temp;
	std::ifstream in_file(file_location_vertex,std::ios::in|std::ios::binary);
	in_file.read(reinterpret_cast<char* >(&temp),sizeof(Point_Struct));
	std::vector<Eigen::Vector3d> pointvec_f0;
	Eigen::Vector3d min_bb(99999.0f,99999.0f,99999.0f);
	Eigen::Vector3d max_bb(-99999.0f,-99999.0f,-99999.0f);
	while(in_file && !in_file.eof())
	{
		Eigen::Vector3d temp_point(temp.x*CONST_SCALING,temp.y*CONST_SCALING,temp.z*CONST_SCALING);
		pointvec_f0.push_back(temp_point);
		in_file.read(reinterpret_cast<char* >(&temp),sizeof(Point_Struct));
		if(temp_point.x()>max_bb.x())
			max_bb(0) = temp_point.x();
		if(temp_point.y()>max_bb.y())
			max_bb(1) = temp_point.y();
		if(temp_point.z()>max_bb.z())
			max_bb(2) = temp_point.z();
		if(temp_point.x()<min_bb.x())
			min_bb(0) = temp_point.x();
		if(temp_point.y()<min_bb.y())
			min_bb(1) = temp_point.y();
		if(temp_point.z()<min_bb.z())
			min_bb(2) = temp_point.z();
	}
	in_file.close();
	int num_vertices = pointvec_f0.size();
	part_of_mesh.resize(num_vertices,false);

	//Added the point data read so far
	point_data.push_back(pointvec_f0);
	
	//Add the bounding box and super bounding box (for a single frame this is the bounding box)
	bounding_box.push_back(std::make_pair(min_bb,max_bb));
	super_bounding_box = std::make_pair(min_bb,max_bb);

	//Now lets read the triangle binaries
	average_edge_length_ = 0.0f;
	int num_counted_edges = 0;
	longest_edge_length_ = 0.0f;
	std::vector<Eigen::Vector3d> normalvec_f0(num_vertices,Eigen::Vector3d::Zero());
	in_file.open(file_location_mesh,std::ios::in|std::ios::binary);
	Triangle_Struct temp_t;
	in_file.read(reinterpret_cast<char* >(&temp_t),sizeof(Triangle_Struct));
	while(in_file && !in_file.eof())
	{
		int vp1 = temp_t.a;
		int vp2 = temp_t.b;
		int vp3 = temp_t.c;

		Eigen::Vector3d v_ab = pointvec_f0[vp2-1] - pointvec_f0[vp1-1];
		Eigen::Vector3d v_ac = pointvec_f0[vp3-1] - pointvec_f0[vp1-1];
		Eigen::Vector3d v_bc = pointvec_f0[vp3-1] - pointvec_f0[vp2-1];

		float ab = v_ab.norm();
		float bc = v_bc.norm();
		float ac = v_ac.norm();

		if(ab>longest_edge_length_)
			longest_edge_length_ = ab;

		if(bc>longest_edge_length_)
			longest_edge_length_ = bc;
		
		if(ac>longest_edge_length_)
			longest_edge_length_ = ac;

		average_edge_length_ += (ab+bc+ac);
		num_counted_edges += 3;

		Eigen::Vector3d tri_norm_abc = Eigen::Vector3d::Zero();
		tri_norm_abc = v_ab.cross(v_ac);
		tri_norm_abc.normalize();

		normalvec_f0[vp1-1] += tri_norm_abc;
		normalvec_f0[vp2-1] += tri_norm_abc;
		normalvec_f0[vp3-1] += tri_norm_abc;

		Triangles temp_tri;
		temp_tri.add(vp1-1,vp2-1,vp3-1);
		temp_tri.addTriangleNormal(tri_norm_abc);
		mesh_data.push_back(temp_tri);

		part_of_mesh[vp1-1] = true;
		part_of_mesh[vp2-1] = true;
		part_of_mesh[vp3-1] = true;

		in_file.read(reinterpret_cast<char* >(&temp_t),sizeof(Triangle_Struct));
	}

	for(int p=0;p<pointvec_f0.size();p++)
	{
		if(!part_of_mesh[p])
			continue;
		normalvec_f0[p].normalize();
	}

	normal_data.push_back(normalvec_f0);

	//Calculate the values of the statistics
	average_edge_length_ /= num_counted_edges;

	//Populate the edge data-structure
	for(int t=0;t<mesh_data.size();t++)
	{
		int a = mesh_data[t].a;
		int b = mesh_data[t].b;
		int c = mesh_data[t].c;

		Edge e1;
		e1.add_edges(a,b);
		Edge e2;
		e2.add_edges(b,c);
		Edge e3;
		e3.add_edges(c,a);

		std::set<Edge>::iterator it = edge_list.find(e1);
		if(it==edge_list.end())
		{
			e1.add_triangle(t);
			edge_list.insert(e1);
		}
		else
		{
			const Edge &re = *it;
			Edge &e = const_cast<Edge&>(re);
			bool status = e.add_triangle(t);
			if(!status)
				std::cout <<  "[WARNING] Non-Manifold mesh detected.\n";
		}

		it = edge_list.find(e2);
		if(it==edge_list.end())
		{
			e2.add_triangle(t);
			edge_list.insert(e2);
		}
		else
		{
			const Edge &re = *it;
			Edge &e = const_cast<Edge&>(re);
			bool status = e.add_triangle(t);
			if(!status)
				std::cout << "[WARNING] Non-Manifold mesh detected.\n";
		}

		it = edge_list.find(e3);
		if(it==edge_list.end())
		{
			e3.add_triangle(t);
			edge_list.insert(e3);
		}
		else
		{
			const Edge &re = *it;
			Edge &e = const_cast<Edge&>(re);
			bool status = e.add_triangle(t);
			if(!status)
				std::cout << "[WARNING] Non-Manifold mesh detected.\n";
		}
	}

	//Now assign an edge index this is needed for the voxelization 
	std::set<Edge>::iterator it;
	int edge_counter = 0;
	for(it = edge_list.begin(); it != edge_list.end() ; it++)
	{
		const Edge &ce = *it;
		Edge &e = const_cast<Edge&>(ce);
		e.set_index(edge_counter);
		edge_counter++;
	}

	//This is the material property 
	material = new Mesh_Material(material_xml);

	//Last modification 
	edgeVector_ = std::vector<Edge>(edge_list.begin(),edge_list.end());
	edgeVectorSize_ = edgeVector_.size();
}

//This is the version that reads k frames with the given generated filename
TriMesh::TriMesh(const char* meshName,
				 const char* file_location_gen,
				 const char* material_xml,
				 int num_frames,
				 float pdfs)
		: CONST_SCALING(pdfs)
{
	name_ = meshName ;
	has_uv_coods_ = false;
	voxelized_ = false;
	static_ = false;
	voxel_grid = new Grid();

	Eigen::Vector3d super_min_bb(99999.0f,99999.0f,99999.0f);
	Eigen::Vector3d super_max_bb(-99999.0f,-99999.0f,-99999.0f);

	std::string line;
	std::string first_iden = "";

	int per_file_vertex_count;
	int num_vertices;

	for(int file_idx = 0; file_idx < num_frames; file_idx++)
	{
		Eigen::Vector3d min_bb(99999.0f,99999.0f,99999.0f);
		Eigen::Vector3d max_bb(-99999.0f,-99999.0f,-99999.0f);

		std::stringstream svb;
		svb << file_location_gen <<  "_vertex_" << file_idx << ".bin" ;
		std::string file_location_vertex = svb.str();

		per_file_vertex_count = 0;

		Point_Struct temp;
		std::ifstream in_file(file_location_vertex.c_str(),std::ios::in|std::ios::binary);
		in_file.read(reinterpret_cast<char* >(&temp),sizeof(Point_Struct));
		std::vector<Eigen::Vector3d> pointvec_fi;

		while(in_file && !in_file.eof())
		{
			Eigen::Vector3d temp_point(temp.x*CONST_SCALING,temp.y*CONST_SCALING,temp.z*CONST_SCALING);
			pointvec_fi.push_back(temp_point);
			in_file.read(reinterpret_cast<char* >(&temp),sizeof(Point_Struct));
			if(temp_point.x()>max_bb.x())
				max_bb(0) = temp_point.x();
			if(temp_point.y()>max_bb.y())
				max_bb(1) = temp_point.y();
			if(temp_point.z()>max_bb.z())
				max_bb(2) = temp_point.z();
			if(temp_point.x()<min_bb.x())
				min_bb(0) = temp_point.x();
			if(temp_point.y()<min_bb.y())
				min_bb(1) = temp_point.y();
			if(temp_point.z()<min_bb.z())
				min_bb(2) = temp_point.z();
		}
		in_file.close();

		//Some measurememnts that are to be read just for the first time. 
		if(file_idx==0)
		{
			num_vertices = pointvec_fi.size();
			part_of_mesh.resize(num_vertices,false);
		}
		else
		{
			per_file_vertex_count = pointvec_fi.size();
			if(per_file_vertex_count != num_vertices)
			{
				std::cout << "[WARNING] Vertex count mismatch between initial frame and frame number# "<< file_idx << "\n";
			}
		}

		//Add the point data read so far
		point_data.push_back(pointvec_fi);
		int num_vertices = pointvec_fi.size();
		
		//Add the bounding box measurements
		bounding_box.push_back(std::make_pair(min_bb,max_bb));

		//Do the super bounding box measurements
		if(max_bb.x()>super_max_bb.x())
			super_max_bb(0) = max_bb.x();
		if(max_bb.y()>super_max_bb.y())
			super_max_bb(1) = max_bb.y();
		if(max_bb.z()>super_max_bb.z())
			super_max_bb(2) = max_bb.z();
		if(min_bb.x()<super_min_bb.x())
			super_min_bb(0) = min_bb.x();
		if(min_bb.y()<super_min_bb.y())
			super_min_bb(1) = min_bb.y();
		if(min_bb.z()<super_min_bb.z())
			super_min_bb(2) = min_bb.z();

		//Now read the triangle binaries
		average_edge_length_ = 0.0f;
		int num_counted_edges = 0;
		longest_edge_length_ = 0.0f;
		std::vector<Eigen::Vector3d> normalvec_fi(num_vertices,Eigen::Vector3d::Zero());
		Triangle_Struct temp_t;
		std::stringstream stm;
		stm << file_location_gen << "_mesh_" << file_idx << ".bin";
		std::string gen_file_name_mesh = stm.str();
		in_file.open(gen_file_name_mesh.c_str(),std::ios::in|std::ios::binary);
		in_file.read(reinterpret_cast<char* >(&temp_t),sizeof(Triangle_Struct));

		while(in_file && !in_file.eof())
		{
			int vp1 = temp_t.a;
			int vp2 = temp_t.b;
			int vp3 = temp_t.c;

			Eigen::Vector3d v_ab = pointvec_fi[vp2-1] - pointvec_fi[vp1-1];
			Eigen::Vector3d v_ac = pointvec_fi[vp3-1] - pointvec_fi[vp1-1];
			Eigen::Vector3d v_bc = pointvec_fi[vp3-1] - pointvec_fi[vp2-1];
		
			Eigen::Vector3d tri_norm_abc = Eigen::Vector3d::Zero();
			tri_norm_abc = v_ab.cross(v_ac);
			tri_norm_abc.normalize();

			normalvec_fi[vp1-1] += tri_norm_abc;
			normalvec_fi[vp2-1] += tri_norm_abc;
			normalvec_fi[vp3-1] += tri_norm_abc;

			if(file_idx==0)
			{
				Triangles temp_tri;
				temp_tri.add(vp1-1,vp2-1,vp3-1);
				temp_tri.addTriangleNormal(tri_norm_abc);
				mesh_data.push_back(temp_tri);

				part_of_mesh[vp1-1] = true;
				part_of_mesh[vp2-1] = true;
				part_of_mesh[vp3-1] = true;

				float ab = v_ab.norm();
				float bc = v_bc.norm();
				float ac = v_ac.norm();

				if(ab>longest_edge_length_)
					longest_edge_length_ = ab;

				if(bc>longest_edge_length_)
					longest_edge_length_ = bc;
		
				if(ac>longest_edge_length_)
					longest_edge_length_ = ac;

				average_edge_length_ += (ab+bc+ac);
				num_counted_edges += 3;
			}

			in_file.read(reinterpret_cast<char* >(&temp_t),sizeof(Triangle_Struct));
		}

		for(int p=0;p<pointvec_fi.size();p++)
		{
			if(!part_of_mesh[p])
				continue;
			normalvec_fi[p].normalize();
		}
		normal_data.push_back(normalvec_fi);

		if(file_idx==0)
		{
			//Calculate the values of the statistics
			average_edge_length_ /= num_counted_edges;

			//Populate the edge data-structure
			for(int t=0;t<mesh_data.size();t++)
			{
				int a = mesh_data[t].a;
				int b = mesh_data[t].b;
				int c = mesh_data[t].c;

				Edge e1;
				e1.add_edges(a,b);
				Edge e2;
				e2.add_edges(b,c);
				Edge e3;
				e3.add_edges(c,a);

				std::set<Edge>::iterator it = edge_list.find(e1);
				if(it==edge_list.end())
				{
					e1.add_triangle(t);
					edge_list.insert(e1);
				}
				else
				{
					const Edge &re = *it;
					Edge &e = const_cast<Edge&>(re);
					bool status = e.add_triangle(t);
					if(!status)
					std::cout <<  "[WARNING] Non-Manifold mesh detected.\n";
				}

				it = edge_list.find(e2);
				if(it==edge_list.end())
				{
					e2.add_triangle(t);
					edge_list.insert(e2);
				}
				else
				{
					const Edge &re = *it;
					Edge &e = const_cast<Edge&>(re);
					bool status = e.add_triangle(t);
					if(!status)
						std::cout << "[WARNING] Non-Manifold mesh detected.\n";
				}

				it = edge_list.find(e3);
				if(it==edge_list.end())
				{
					e3.add_triangle(t);
					edge_list.insert(e3);
				}
				else
				{
					const Edge &re = *it;
					Edge &e = const_cast<Edge&>(re);
					bool status = e.add_triangle(t);
					if(!status)
						std::cout << "[WARNING] Non-Manifold mesh detected.\n";
				}
			}
		}
	}

	//Now assign an edge index this is needed for the voxelization 
	std::set<Edge>::iterator it;
	int edge_counter = 0;
	for(it = edge_list.begin(); it != edge_list.end() ; it++)
	{
		const Edge &ce = *it;
		Edge &e = const_cast<Edge&>(ce);
		e.set_index(edge_counter);
		edge_counter++;
	}

	//Now seal off the super bounding box 
	super_bounding_box = std::make_pair(super_min_bb,super_max_bb);

	//This is the material property 
	material = new Mesh_Material(material_xml);

	//Last modification 
	edgeVector_ = std::vector<Edge>(edge_list.begin(),edge_list.end());
	edgeVectorSize_ = edgeVector_.size();
}

TriMesh::~TriMesh(void)
{
}

void TriMesh::initializeDynamicMesh(std::vector<Eigen::Vector3d>& pvec, std::vector<int> triangles, const char* meshName, const char* material_xml)
	{
	name_ = meshName ;
	has_uv_coods_ = false;
	static_ = false;
	voxelized_ = false;
	voxel_grid = new Grid();
	int num_vertices;
	int a,b,c;
	//std::vector<Eigen::Vector3d> pointvec_f0;
	std::vector<Eigen::Vector3d> normalvec_f0;
	Eigen::Vector3d min_bb(99999.0f,99999.0f,99999.0f);
	Eigen::Vector3d max_bb(-99999.0f,-99999.0f,-99999.0f);
	
	average_edge_length_ = 0.0f;
	int num_counted_edges = 0;
	longest_edge_length_ = 0.0f;
	static_ = false;

	std::cout <<"DDD";

	num_vertices = pvec.size();
	normalvec_f0.resize(num_vertices,Eigen::Vector3d::Zero());
	part_of_mesh.resize(num_vertices,false);

	//Find the bounding box here
	for(int i=0;i<num_vertices;i++)
	{
		Eigen::Vector3d temp_point = pvec[i];
		if(temp_point.x()>max_bb.x())
			max_bb(0) = temp_point.x();
		if(temp_point.y()>max_bb.y())
			max_bb(1) = temp_point.y();
		if(temp_point.z()>max_bb.z())
			max_bb(2) = temp_point.z();
		if(temp_point.x()<min_bb.x())
			min_bb(0) = temp_point.x();
		if(temp_point.y()<min_bb.y())
			min_bb(1) = temp_point.y();
		if(temp_point.z()<min_bb.z())
			min_bb(2) = temp_point.z();
	}

	//Do the triangle related operations here 
	for(int t=0;t<(triangles.size()/3);t++)
	{
		a = triangles[3*t+0];
		b = triangles[3*t+1];
		c = triangles[3*t+2];

		Eigen::Vector3d v_ab = pvec[b] - pvec[a];
		Eigen::Vector3d v_ac = pvec[c] - pvec[a];
		Eigen::Vector3d v_bc = pvec[c] - pvec[b];

		float ab = v_ab.norm();
		float bc = v_bc.norm();
		float ac = v_ac.norm();

		if(ab>longest_edge_length_)
			longest_edge_length_ = ab;

		if(bc>longest_edge_length_)
			longest_edge_length_ = bc;

		if(ac>longest_edge_length_)
			longest_edge_length_ = ac;

		average_edge_length_ += (ab+bc+ac);
		num_counted_edges += 3;

		Eigen::Vector3d tri_norm_abc = Eigen::Vector3d::Zero();
		tri_norm_abc = v_ab.cross(v_ac);
		tri_norm_abc.normalize();

		normalvec_f0[a] += tri_norm_abc;
		normalvec_f0[b] += tri_norm_abc;
		normalvec_f0[c] += tri_norm_abc;

		Triangles temp_tri;
		temp_tri.add(a,b,c);
		temp_tri.addTriangleNormal(tri_norm_abc);
		mesh_data.push_back(temp_tri);

		part_of_mesh[a] = true;
		part_of_mesh[b] = true;
		part_of_mesh[c] = true;

	}

	point_data.push_back(pvec);
	volInfo_.push_back(VolumeInformation(point_data[0]));

	//Add the bounding box and super bounding box (for a single frame this is the bounding box)
	bounding_box.push_back(std::make_pair(min_bb,max_bb));
	super_bounding_box = std::make_pair(min_bb,max_bb);

	for(int p=0;p<pvec.size();p++)
		{
		if(!part_of_mesh[p])
			continue;
		normalvec_f0[p].normalize();
		}
	normal_data.push_back(normalvec_f0);

	//Calculate the values of the statistics
	average_edge_length_ /= num_counted_edges;

	//Populate the edge data-structure
	for(int t=0;t<mesh_data.size();t++)
		{
		int a = mesh_data[t].a;
		int b = mesh_data[t].b;
		int c = mesh_data[t].c;

		Edge e1;
		e1.add_edges(a,b);
		Edge e2;
		e2.add_edges(b,c);
		Edge e3;
		e3.add_edges(c,a);

		std::set<Edge>::iterator it = edge_list.find(e1);
		if(it==edge_list.end())
			{
			e1.add_triangle(t);
			edge_list.insert(e1);
			}
		else
			{
			const Edge &re = *it;
			Edge &e = const_cast<Edge&>(re);
			bool status = e.add_triangle(t);
			if(!status)
				std::cout <<  "[WARNING] Non-Manifold mesh detected.\n";
			}

		it = edge_list.find(e2);
		if(it==edge_list.end())
			{
			e2.add_triangle(t);
			edge_list.insert(e2);
			}
		else
			{
			const Edge &re = *it;
			Edge &e = const_cast<Edge&>(re);
			bool status = e.add_triangle(t);
			if(!status)
				std::cout << "[WARNING] Non-Manifold mesh detected.\n";
			}

		it = edge_list.find(e3);
		if(it==edge_list.end())
			{
			e3.add_triangle(t);
			edge_list.insert(e3);
			}
		else
			{
			const Edge &re = *it;
			Edge &e = const_cast<Edge&>(re);
			bool status = e.add_triangle(t);
			if(!status)
				std::cout << "[WARNING] Non-Manifold mesh detected.\n";
			}
		}

	//Now assign an edge index this is needed for the voxelization 
	std::set<Edge>::iterator it;
	int edge_counter = 0;
	for(it = edge_list.begin(); it != edge_list.end() ; it++)
		{
		const Edge &ce = *it;
		Edge &e = const_cast<Edge&>(ce);
		e.set_index(edge_counter);
		edge_counter++;
		}
	//This is the material property 
	material = new Mesh_Material(material_xml);

	//Last modification 
	edgeVector_ = std::vector<Edge>(edge_list.begin(),edge_list.end());
	edgeVectorSize_ = edgeVector_.size();
	}


// Some of the accessor functions modified for use with static meshes
Eigen::Vector3d TriMesh::get_point_data(int idx, int frame)	const
{
	if(static_)
		frame = 0;
	return point_data[frame][idx];
}

Eigen::Vector3d TriMesh::get_normal_data(int idx, int frame)	const
{
	if(static_)
		frame = 0;
	return normal_data[frame][idx];
}

int TriMesh::get_number_frames() const
{
	if(static_)
		return 1;
	else 
		return point_data.size();
}

// Other body functions 
void TriMesh::generate_frame_normals(int frame_idx)
{
	int num_vertices = this->get_number_vertices();
	int num_triangles = this->get_number_triangles();
	int num_frames = this->get_number_frames();

	if(frame_idx >= num_frames)
	{
		std::cout << "[ERROR] Frame point data is not available for given frame to generate normals.\n";
		return ;
	}

	std::vector<Eigen::Vector3d> new_normal_vector(num_vertices,Eigen::Vector3d::Zero());

	OMP_FOR
	for(int t=0;t<num_triangles;t++)
	{
		Triangles tri = this->get_triangle(t);

		Eigen::Vector3d v_a = this->get_point_data(tri.a,frame_idx);
		Eigen::Vector3d v_b = this->get_point_data(tri.b,frame_idx);
		Eigen::Vector3d v_c = this->get_point_data(tri.c,frame_idx);

		Eigen::Vector3d v_ab = v_b - v_a ;
		Eigen::Vector3d v_ac = v_c - v_a ;
		Eigen::Vector3d v_bc = v_c - v_b ;

		Eigen::Vector3d tri_norm_abc = Eigen::Vector3d::Zero();
		tri_norm_abc = v_ab.cross(v_ac);
		tri_norm_abc.normalize();

		new_normal_vector[tri.a] += tri_norm_abc;
		new_normal_vector[tri.b] += tri_norm_abc;
		new_normal_vector[tri.c] += tri_norm_abc;

	}

	OMP_FOR
	for(int p=0;p<num_vertices;p++)
	{
		if(!(this->is_part_of_mesh(p)))
			continue;

		new_normal_vector[p].normalize();
	}

	//Now we need to be careful while addign this 
	int size_of_normal = normal_data.size();
	if(frame_idx==size_of_normal)
		normal_data.push_back(new_normal_vector);

	//@TODO: Generate this for arbitrary cases. Though I dont think we will need that 

}

void TriMesh::update_voxelization()
{
	//First check how many frames are present
	int num_frames_present = this->get_number_frames();
	//Next check how many frames have been voxelized 
	int num_frames_voxelized = voxel_to_triangle.size();
	//Now how many we have to voxelize 
	int num_frames_to_voxelize = num_frames_present-num_frames_voxelized;
	for(int i = 0; i <num_frames_to_voxelize ; i++)
	{
		int cur_frame_id = num_frames_voxelized+i;
		std::cout << "[INFO] Voxelizing frame # " <<  cur_frame_id << " for mesh :" << name_ << ".\n";

		//Compute the voxel location per frame
		std::vector<Voxel> point_voxel_curframe;
		point_voxel_curframe.resize(this->get_number_vertices());
		for(int p=0;p<this->get_number_vertices();p++)
		{
			Eigen::Vector3d pnt = this->get_point_data(p,cur_frame_id);
			Voxel pnt_voxel = voxel_grid->get_voxel(pnt);
			point_voxel_curframe[p] = pnt_voxel;
		}
		point_to_voxel.push_back(point_voxel_curframe);

		//Compute the voxel location per triangle
		std::vector< std::vector <Voxel> > triangle_voxel_curframe;
		triangle_voxel_curframe.resize(this->get_number_triangles());
		for(int t=0;t<this->get_number_triangles();t++)
		{
			Voxel vox_a = point_voxel_curframe[mesh_data[t].a];
			Voxel vox_b = point_voxel_curframe[mesh_data[t].b];
			Voxel vox_c = point_voxel_curframe[mesh_data[t].c];

			std::vector<Voxel> triangle_voxel;
			triangle_voxel.resize(3);
			triangle_voxel[0] = vox_a;
			triangle_voxel[1] = vox_b;
			triangle_voxel[2] = vox_c;

			triangle_voxel_curframe[t]=triangle_voxel;

		}
		triangle_to_voxel.push_back(triangle_voxel_curframe);

		//Compute the voxel location per edge
		std::vector < std::vector <Voxel> > edge_voxel_curframe;
		edge_voxel_curframe.resize(this->get_number_edges());
		std::set<Edge>::iterator it;
		for(it = edge_list.begin(); it != edge_list.end(); it++)
		{
			const Edge e  = *it;
			Voxel vox_a = point_voxel_curframe[e.get_start_vertex()];
			Voxel vox_b = point_voxel_curframe[e.get_end_vertex()];

			std::vector<Voxel> edge_voxel;
			edge_voxel.resize(2);
			edge_voxel[0] = vox_a;
			edge_voxel[1] = vox_b;

			edge_voxel_curframe[e.get_index()] = edge_voxel;
		}
		edge_to_voxel.push_back(edge_voxel_curframe);
	}

}

void TriMesh::voxelize_frame(int frame_id)
{
	int frame_till_voxelized = point_to_voxel.size() - 1;

	if(frame_till_voxelized<frame_id)
	{
		//Resize the voxel data structures
		point_to_voxel.resize(frame_id+1);
		triangle_to_voxel.resize(frame_id+1);
		edge_to_voxel.resize(frame_id+1);

		//@TODO:Resize the inverse data strcutures too
	}

	//Compute the voxel location per frame
	std::vector<Voxel> point_voxel_curframe;
	point_voxel_curframe.resize(this->get_number_vertices());
	for(int p=0;p<this->get_number_vertices();p++)
	{
		Eigen::Vector3d pnt = this->get_point_data(p,frame_id);
		Voxel pnt_voxel = voxel_grid->get_voxel(pnt);
		point_voxel_curframe[p] = pnt_voxel;
	}
	point_to_voxel[frame_id] = point_voxel_curframe;

	//Compute the voxel location per triangle
	std::vector< std::vector <Voxel> > triangle_voxel_curframe;
	triangle_voxel_curframe.resize(this->get_number_triangles());
	for(int t=0;t<this->get_number_triangles();t++)
	{
		Voxel vox_a = point_voxel_curframe[mesh_data[t].a];
		Voxel vox_b = point_voxel_curframe[mesh_data[t].b];
		Voxel vox_c = point_voxel_curframe[mesh_data[t].c];

		std::vector<Voxel> triangle_voxel;
		triangle_voxel.resize(3);
		triangle_voxel[0] = vox_a;
		triangle_voxel[1] = vox_b;
		triangle_voxel[2] = vox_c;

		triangle_voxel_curframe[t]=triangle_voxel;

	}
	triangle_to_voxel[frame_id] = triangle_voxel_curframe;

	//Compute the voxel location per edge
	std::vector < std::vector <Voxel> > edge_voxel_curframe;
	edge_voxel_curframe.resize(this->get_number_edges());
	std::set<Edge>::iterator it;
	for(it = edge_list.begin(); it != edge_list.end(); it++)
	{
		const Edge e  = *it;
		Voxel vox_a = point_voxel_curframe[e.get_start_vertex()];
		Voxel vox_b = point_voxel_curframe[e.get_end_vertex()];

		std::vector<Voxel> edge_voxel;
		edge_voxel.resize(2);
		edge_voxel[0] = vox_a;
		edge_voxel[1] = vox_b;

		edge_voxel_curframe[e.get_index()] = edge_voxel;
	}
	edge_to_voxel[frame_id] = edge_voxel_curframe;
}

//Dev notes: You will need to add 3 things when you add a new frame
// The frame data
// The normals 
// The bounding box calculations 
void TriMesh::add_frame_data(std::vector<Eigen::Vector3d> &vec)
{
	//Add this to the vector 
	point_data.push_back(vec);

	//Calculate the normals 
	int num_vertices = vec.size();
	int num_triangles = this->get_number_triangles();
	std::vector<Eigen::Vector3d> new_normal_vector(num_vertices,Eigen::Vector3d::Zero());
	for(int t=0;t<num_triangles;t++)
	{
		Triangles tri = this->get_triangle(t);

		Eigen::Vector3d v_a = vec[tri.a];
		Eigen::Vector3d v_b = vec[tri.b];
		Eigen::Vector3d v_c = vec[tri.c];

		Eigen::Vector3d v_ab = v_b - v_a ;
		Eigen::Vector3d v_ac = v_c - v_a ;
		Eigen::Vector3d v_bc = v_c - v_b ;

		Eigen::Vector3d tri_norm_abc = Eigen::Vector3d::Zero();
		tri_norm_abc = v_ab.cross(v_ac);
		tri_norm_abc.normalize();

		new_normal_vector[tri.a] += tri_norm_abc;
		new_normal_vector[tri.b] += tri_norm_abc;
		new_normal_vector[tri.c] += tri_norm_abc;

	}
	//Bounding box calculation 
	Eigen::Vector3d min_bb(99999.0f,99999.0f,99999.0f);
	Eigen::Vector3d max_bb(-99999.0f,-99999.0f,-99999.0f);
	for(int p=0;p<num_vertices;p++)
	{
		//Calculate the bounding box coordinates
		Eigen::Vector3d temp_point = vec[p];
		if(temp_point.x()>max_bb.x())
			max_bb(0) = temp_point.x();
		if(temp_point.y()>max_bb.y())
			max_bb(1) = temp_point.y();
		if(temp_point.z()>max_bb.z())
			max_bb(2) = temp_point.z();
		if(temp_point.x()<min_bb.x())
			min_bb(0) = temp_point.x();
		if(temp_point.y()<min_bb.y())
			min_bb(1) = temp_point.y();
		if(temp_point.z()<min_bb.z())
			min_bb(2) = temp_point.z();

		if(!(this->is_part_of_mesh(p)))
			continue;
		new_normal_vector[p].normalize();
	}
	normal_data.push_back(new_normal_vector);

	//Bounding box calculation 
	bounding_box.push_back(std::make_pair(min_bb,max_bb));
	Eigen::Vector3d super_min_bb = super_bounding_box.first;
	Eigen::Vector3d super_max_bb = super_bounding_box.second;
	if(max_bb.x()>super_max_bb.x())
		super_max_bb(0) = max_bb.x();
	if(max_bb.y()>super_max_bb.y())
		super_max_bb(1) = max_bb.y();
	if(max_bb.z()>super_max_bb.z())
		super_max_bb(2) = max_bb.z();
	if(min_bb.x()<super_min_bb.x())
		super_min_bb(0) = min_bb.x();
	if(min_bb.y()<super_min_bb.y())
		super_min_bb(1) = min_bb.y();
	if(min_bb.z()<super_min_bb.z())
		super_min_bb(2) = min_bb.z();
	super_bounding_box = std::make_pair(super_min_bb,super_max_bb);
}

//Drawing functions 
void TriMesh::draw_mesh(int frame_idx)
{
	if(static_)
		frame_idx = 0;

	material->set_material();
	//glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	glEnable (GL_LIGHTING);
	glEnable(GL_NORMALIZE);
	glDisable(GL_CULL_FACE);
	//This function draws the current render frame
	glBegin(GL_TRIANGLES);
	for(int t=0;t<this->get_number_triangles();t++)
	{
		Eigen::Vector3d p1 = point_data[frame_idx][mesh_data[t].a];
		Eigen::Vector3d p2 = point_data[frame_idx][mesh_data[t].b];
		Eigen::Vector3d p3 = point_data[frame_idx][mesh_data[t].c];


		Eigen::Vector3d n1 = normal_data[frame_idx][mesh_data[t].a];
		Eigen::Vector3d n2 = normal_data[frame_idx][mesh_data[t].b];
		Eigen::Vector3d n3 = normal_data[frame_idx][mesh_data[t].c];
		
		glNormal3f(n1.x(),n1.y(),n1.z());
		glVertex3f(p1.x(),p1.y(),p1.z());
		glNormal3f(n2.x(),n2.y(),n2.z());
		glVertex3f(p2.x(),p2.y(),p2.z());
		glNormal3f(n3.x(),n3.y(),n3.z());
		glVertex3f(p3.x(),p3.y(),p3.z());
	}
	glEnd();
	glDisable(GL_LIGHTING);
}

void TriMesh::draw_wireframe(int frame_idx)
{
	if(static_)
		frame_idx = 0;

	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	glDisable(GL_LIGHTING);
	glColor3f(1.0f,1.0f,1.0f);
	glBegin(GL_TRIANGLES);

	//This function draws the current render frame
	for(int t=0;t<this->get_number_triangles();t++)
	{
		Eigen::Vector3d p1 = point_data[frame_idx][mesh_data[t].a];
		Eigen::Vector3d p2 = point_data[frame_idx][mesh_data[t].b];
		Eigen::Vector3d p3 = point_data[frame_idx][mesh_data[t].c];

		


		Eigen::Vector3d n1 = normal_data[frame_idx][mesh_data[t].a];
		Eigen::Vector3d n2 = normal_data[frame_idx][mesh_data[t].b];
		Eigen::Vector3d n3 = normal_data[frame_idx][mesh_data[t].c];

		glNormal3f(n1.x(),n1.y(),n1.z());
		glVertex3f(p1.x(),p1.y(),p1.z());
		glNormal3f(n2.x(),n2.y(),n2.z());
		glVertex3f(p2.x(),p2.y(),p2.z());
		glNormal3f(n3.x(),n3.y(),n3.z());
		glVertex3f(p3.x(),p3.y(),p3.z());
	}
	glEnd();
	glEnable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
}

void TriMesh::draw_bounding_box(int frame_idx)
{
	if(static_)
		frame_idx = 0;

	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	glDisable(GL_LIGHTING);
	glColor3f(1.0f,1.0f,1.0f);

	Eigen::Vector3d p1 = bounding_box[frame_idx].first;
	Eigen::Vector3d p2 = bounding_box[frame_idx].second;
	float x1 = p1.x(); float y1 = p1.y(); float z1 = p1.z();
	float x2 = p2.x(); float y2 = p2.y(); float z2 = p2.z();

	glBegin(GL_QUADS);

	glVertex3f(x1,y1,z1);
	glVertex3f(x2,y1,z1);
	glVertex3f(x2,y2,z1);
	glVertex3f(x1,y2,z1);

	
	glVertex3f(x1,y1,z1);
	glVertex3f(x1,y2,z1);
	glVertex3f(x1,y2,z2);
	glVertex3f(x1,y1,z2);

	
	glVertex3f(x1,y2,z1);
	glVertex3f(x2,y2,z1);
	glVertex3f(x2,y2,z2);
	glVertex3f(x1,y2,z2);
	
	glVertex3f(x2,y1,z1);
	glVertex3f(x2,y2,z1);
	glVertex3f(x2,y2,z2);
	glVertex3f(x2,y1,z2);
	
	glVertex3f(x1,y1,z1);
	glVertex3f(x2,y1,z1);
	glVertex3f(x2,y1,z2);
	glVertex3f(x1,y1,z2);

	glVertex3f(x1,y1,z2);
	glVertex3f(x2,y1,z2);
	glVertex3f(x2,y2,z2);
	glVertex3f(x1,y2,z2);
	
	glEnd();
}

void TriMesh::draw_super_bounding_box()
{

	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	glDisable(GL_LIGHTING);
	glColor3f(1.0f,1.0f,1.0f);

	Eigen::Vector3d p1 = super_bounding_box.first;
	Eigen::Vector3d p2 = super_bounding_box.second;
	float x1 = p1.x(); float y1 = p1.y(); float z1 = p1.z();
	float x2 = p2.x(); float y2 = p2.y(); float z2 = p2.z();

	glBegin(GL_QUADS);

	glVertex3f(x1,y1,z1);
	glVertex3f(x2,y1,z1);
	glVertex3f(x2,y2,z1);
	glVertex3f(x1,y2,z1);

	
	glVertex3f(x1,y1,z1);
	glVertex3f(x1,y2,z1);
	glVertex3f(x1,y2,z2);
	glVertex3f(x1,y1,z2);

	
	glVertex3f(x1,y2,z1);
	glVertex3f(x2,y2,z1);
	glVertex3f(x2,y2,z2);
	glVertex3f(x1,y2,z2);
	
	glVertex3f(x2,y1,z1);
	glVertex3f(x2,y2,z1);
	glVertex3f(x2,y2,z2);
	glVertex3f(x2,y1,z2);
	
	glVertex3f(x1,y1,z1);
	glVertex3f(x2,y1,z1);
	glVertex3f(x2,y1,z2);
	glVertex3f(x1,y1,z2);

	glVertex3f(x1,y1,z2);
	glVertex3f(x2,y1,z2);
	glVertex3f(x2,y2,z2);
	glVertex3f(x1,y2,z2);
	
	glEnd();
}

void TriMesh::draw_vector_valued_heatmap(std::vector<Eigen::Vector3d> vector_vertex_values,int frame_idx)
{
	if(static_)
		frame_idx = 0;

	//For color map visualization
	int r_c_idx[5] = {0.0,1.0,0.0,0.0,1.0};
	int g_c_idx[5] = {0.0,1.0,1.0,1.0,0.0};
	int b_c_idx[5] = {1.0,0.0,0.0,1.0,0.0};

	//Convert the vector to scalar
	std::vector<float> scalar_values;
	scalar_values.resize(vector_vertex_values.size());
	float min = 99999.9f;
	float max = -99999.9f;
	std::vector<Eigen::Vector3d>::iterator it;
	unsigned int counter = 0;
	for(it = vector_vertex_values.begin(); it != vector_vertex_values.end() ; ++it)
	{
		Eigen::Vector3d cur_vec = *it;
		float scalar_value = cur_vec.norm();
		scalar_values[counter] = scalar_value;
		counter++;
		if(scalar_value>max)
			max = scalar_value;
		if(scalar_value<min)
			min = scalar_value;
	}

	//Get the color values to the scalars 
	std::vector<Eigen::Vector3d> vertex_colors;
	vertex_colors.resize(vector_vertex_values.size());
	std::vector<float>::iterator it_sc;
	counter = 0;
	for(it_sc = scalar_values.begin(); it_sc != scalar_values.end() ; it_sc++)
	{
		int col_idx = 0;
		double attrition_fac = (*it_sc)/(max-min);
		if(attrition_fac<0.25)
			col_idx = 0;
		else if(attrition_fac<0.5)
			col_idx = 1;
		else if(attrition_fac<0.75)
			col_idx = 2;
		else
			col_idx = 3;
		float r_c = r_c_idx[col_idx]*(1-attrition_fac)+r_c_idx[col_idx+1]*attrition_fac;
		float g_c = g_c_idx[col_idx]*(1-attrition_fac)+g_c_idx[col_idx+1]*attrition_fac;
		float b_c = b_c_idx[col_idx]*(1-attrition_fac)+b_c_idx[col_idx+1]*attrition_fac;
		vertex_colors[counter] = Eigen::Vector3d(r_c,g_c,b_c);
		counter++;
	}

	//Render the triangles
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	glDisable(GL_LIGHTING);
	glBegin(GL_TRIANGLES);
	std::vector<Triangles>::iterator it_tri;
	for(it_tri = mesh_data.begin(); it_tri != mesh_data.end() ; ++it_tri)
	{
		Eigen::Vector3d p1 = point_data[frame_idx][it_tri->a];
		Eigen::Vector3d p2 = point_data[frame_idx][it_tri->b];
		Eigen::Vector3d p3 = point_data[frame_idx][it_tri->c];

		Eigen::Vector3d c1 = vertex_colors[it_tri->a];
		Eigen::Vector3d c2 = vertex_colors[it_tri->b];
		Eigen::Vector3d c3 = vertex_colors[it_tri->c];
	
		glColor3f(c1.x(),c1.y(),c1.z());
		glVertex3f(p1.x(),p1.y(),p1.z());
		glColor3f(c2.x(),c2.y(),c2.z());
		glVertex3f(p2.x(),p2.y(),p2.z());
		glColor3f(c3.x(),c3.y(),c3.z());
		glVertex3f(p3.x(),p3.y(),p3.z());

	}
}

void TriMesh::resetMeshDefaults()
{
	//Points
	std::vector<Eigen::Vector3d> firstFramePoints = point_data[0];
	point_data.clear();
	point_data.push_back(firstFramePoints);

	//Normals
	std::vector<Eigen::Vector3d> firstFrameNormals = normal_data[0];
	normal_data.clear();
	normal_data.push_back(firstFrameNormals);

	//@TODO : Bounding Boxes Calculation

	//@TODO : Voxelization too possibly
}

int TriMesh::performMeshSelectionRayTest(double* rayStart, double* rayEnd, double *clickedWorldPos, double* selectedPos, int frameId) {

	if(static_)
		frameId = 0;

	int numTriangles = this->get_number_triangles();
	std::vector<std::pair<double,int> > intersectedTriangles;
	for(int i = 0; i<numTriangles; i++) {
		Eigen::Vector3d v0 = point_data[frameId][mesh_data[i].a];
		Eigen::Vector3d v1 = point_data[frameId][mesh_data[i].b];
		Eigen::Vector3d v2 = point_data[frameId][mesh_data[i].c];
		double u,v,t;
		bool intersects = segmentTriangleIntersection(rayStart,rayEnd,v0.data(),v1.data(),v2.data(),0.001,u,v,t);
		if(intersects) {
			intersectedTriangles.push_back(std::pair<double,int>(t,i));
		}
	}
	if (intersectedTriangles.size() > 0) {
		std::sort(intersectedTriangles.begin(), intersectedTriangles.end());
		int triIdx =  intersectedTriangles.front().second;
		Eigen::Vector3d v0 = point_data[frameId][mesh_data[triIdx].a];
		Eigen::Vector3d v1 = point_data[frameId][mesh_data[triIdx].b];
		Eigen::Vector3d v2 = point_data[frameId][mesh_data[triIdx].c];
		double minDist[3] = {pointLineDistance(rayStart, rayEnd, v0.data()),pointLineDistance(rayStart, rayEnd, v1.data()),pointLineDistance(rayStart, rayEnd, v2.data())};
		int selectedVertex = mesh_data[triIdx].a;
		int minIdx = 0;
		if (minDist[1] < minDist[0]) {
			selectedVertex = mesh_data[triIdx].b;
			minIdx = 1;
		}
		if (minDist[2] < minDist[minIdx]) {
			minIdx =2;
			selectedVertex = mesh_data[triIdx].c;
		}
		selectedPos[0] = point_data[frameId][selectedVertex].x();
		selectedPos[1] = point_data[frameId][selectedVertex].y();
		selectedPos[2] = point_data[frameId][selectedVertex].z();
		return selectedVertex; 
	} else {
		return -1;
	}


}

void TriMesh::saveFrameAsObj(const char* filePrefix, int frameID) {
	std::stringstream ss;
	ss << filePrefix << "\\mesh_" << frameID << ".obj";
	std::string s = ss.str();
	std::cout << "Writing OBJ FILE to " << s.c_str() << "\n";
	std::ofstream fileOut(s.c_str());
	for(int p=0;p<point_data[frameID].size();p++) {
		fileOut << "v " << point_data[frameID][p][0] << " " << point_data[frameID][p][1] << " " << point_data[frameID][p][2] << "\n";
	}
	for(int t=0;t<mesh_data.size();t++) {
		fileOut << "f " << (mesh_data[t].a)+1 << " " << (mesh_data[t].b)+1 << " " << (mesh_data[t].c)+1 << "\n";
	}
	fileOut.close();

}

void TriMesh::saveAsOBJ(const char* filePrefix) {
	if(static_) {
		saveFrameAsObj(filePrefix,0);
	}
	else {
		for(int i=0;i<point_data.size();i++) {
			saveFrameAsObj(filePrefix,i);
		}
	}
}



//added by Guowei

