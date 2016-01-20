#include "fem_solver.h"


FEM_Solver::FEM_Solver(void)
{
}


FEM_Solver::~FEM_Solver(void)
{
}

void FEM_Solver::advance_time_step(Cloth_Data* cloth)
{
	//First get the farme id we will process
	last_frame_id = (cloth->getMesh()->get_number_frames()) - 1;

	//Compute in-plane deformation forces
	this->compute_shear_forces(cloth);
	//Compute bending forces 
	this->compute_bending_forces(cloth);
	//Compute graavity forces
	this->compute_gravity_forces(cloth);
	//Compute damping forces
	this->computeDampingForce(cloth);
	//Generate the velocities
	this->get_next_velocity_positions(cloth);

}

void FEM_Solver::compute_shear_forces(Cloth_Data* cloth)
{
	float young_modulus = cloth->get_property_obj()->get_youngs_modulus();
	float poisson = cloth->get_property_obj()->get_poisson_ratio();

	int num_vertices = cloth->getMesh()->get_number_vertices();
	int num_triangles = cloth->getMesh()->get_number_triangles();

	//std::cout << num_vertices << " " << num_triangles << "\n";

	//@TODO: make this more generic by replacing it with a generic function 
	Eigen::Matrix3d E = isotropic_linear_elasticity(young_modulus, poisson);

	//Check size if needed resize then fill 
	if(shear_forces.size()<num_vertices)
	{
		shear_forces.resize(num_vertices);
	}
	std::fill(shear_forces.begin(),shear_forces.end(),Eigen::Vector3d::Zero());
	
	OMP_FOR
	for(int t=0;t<num_triangles;t++)
	{
		Triangles tri = cloth->getMesh()->get_triangle(t);

		Eigen::Vector3d x0_4 = cloth->getMesh()->get_point_data( tri.a,last_frame_id);
		Eigen::Vector3d x1_4 = cloth->getMesh()->get_point_data( tri.b,last_frame_id);
		Eigen::Vector3d x2_4 = cloth->getMesh()->get_point_data( tri.c,last_frame_id);

		Eigen::Vector3d pa(x0_4.x(),x0_4.y(),x0_4.z());
		Eigen::Vector3d pb(x1_4.x(),x1_4.y(),x1_4.z());
		Eigen::Vector3d pc(x2_4.x(),x2_4.y(),x2_4.z());

		float r_ua = (cloth->get_vertex_distribution(t))(0);
		float r_va = (cloth->get_vertex_distribution(t))(1);
		float r_ub = (cloth->get_vertex_distribution(t))(2);
		float r_vb = (cloth->get_vertex_distribution(t))(3);
		float r_uc = (cloth->get_vertex_distribution(t))(4);
		float r_vc = (cloth->get_vertex_distribution(t))(5);
		float d = (cloth->get_vertex_distribution(t))(6);

		//if(t==34)
	//		std::cout << "Triangle# " << t << r_ua << " " << r_va << " " << r_ub << " " << r_vb << " " << r_uc << " " << r_vc << "\n"; 

		Eigen::Vector3d U = (r_ua * pa) + (r_ub * pb) + (r_uc * pc);
		Eigen::Vector3d V = (r_va * pa) + (r_vb * pb) + (r_vc * pc);

		Eigen::Vector3d strain_tensor;
		strain_tensor[0] = 0.5f*(U.dot(U)-1.0f);
		strain_tensor[1] = 0.5f*(V.dot(V)-1.0f);
		strain_tensor[2] = U.dot(V);

		Eigen::Vector3d stress_tensor = E * strain_tensor;
		//if(t%21==0)
		//	DBG__Vector_Dump(strain_tensor,t,"Strain tensor");

		Eigen::Vector3d fA = (stress_tensor[0] * r_ua * U) + (stress_tensor[1] * r_va * V) + (stress_tensor[2] * r_ua * V ) + (stress_tensor[2] * r_va * U );
		fA *= (-fabs(d)*0.5f);
		Eigen::Vector3d fB = (stress_tensor[0] * r_ub * U) + (stress_tensor[1] * r_vb * V) + (stress_tensor[2] * r_ub * V ) + (stress_tensor[2] * r_vb * U );
		fB *= (-fabs(d)*0.5f);
		Eigen::Vector3d fC = (stress_tensor[0] * r_uc * U) + (stress_tensor[1] * r_vc * V) + (stress_tensor[2] * r_uc * V ) + (stress_tensor[2] * r_vc * U );
		fC *= (-fabs(d)*0.5f);

		//if(t%251==0)
		//{
			//DBG__Vector_Dump(strain_tensor,t,"STRAIN");
			//DBG__Vector_Dump(U,t,"U");
			//DBG__Vector_Dump(fA,t,"F_a");
			//DBG__Vector_Dump(fB,t,"F_b");
			//DBG__Vector_Dump(fC,t,"F_c");
		//}

		shear_forces[tri.a] += fA;
		shear_forces[tri.b] += fB;
		shear_forces[tri.c] += fC;
	}
	
}

void FEM_Solver::compute_bending_forces(Cloth_Data* cloth)
{
	int num_triangles = cloth->getMesh()->get_number_triangles();
	int num_vertices = cloth->getMesh()->get_number_vertices();

	float k_bend = cloth->get_property_obj()->get_k_bend();
	//std::cout << k_bend << "\n";
	//Check size if needed resize then fill 
	if(bend_forces.size()<num_vertices)
	{
		bend_forces.resize(num_vertices);
	}
	std::fill(bend_forces.begin(),bend_forces.end(),Eigen::Vector3d::Zero());

	std::set<Edge> edge_list = cloth->getMesh()->get_edge_list();
	std::set<Edge>::iterator it;

	std::vector<Edge> edgeVector = cloth->getMesh()->getEdgeVector();
	int sizeOfEdgeVector = cloth->getMesh()->getEdgeVectorSize();

	int counter = 0;
	OMP_FOR
	for(int i=0;i<sizeOfEdgeVector;i++)
	{
		int t1 = edgeVector[i].get_tri_1();
		int t2 = edgeVector[i].get_tri_2();
		int v3 = edgeVector[i].get_start_vertex();
		int v4 = edgeVector[i].get_end_vertex();

		if(t2==-1)
			continue;
		counter++;
		Triangles tri1 = cloth->getMesh()->get_triangle(t1);
		Triangles tri2 = cloth->getMesh()->get_triangle(t2);

		int v1,v2;

		//Resolve the point dependencies to get x1 and x2 
		int t1_v0 = tri1.a;
		int t1_v1 = tri1.b;
		int t1_v2 = tri1.c;

		if((t1_v0 == v3 && t1_v1 == v4)||(t1_v0 == v4 && t1_v1 == v3))
			v1 = t1_v2;
		else if((t1_v1 == v3 && t1_v2 == v4)||(t1_v1 == v4 && t1_v2 == v3))
			v1 = t1_v0;
		else 
			v1 = t1_v1;

		int t2_v0 = tri2.a;
		int t2_v1 = tri2.b;
		int t2_v2 = tri2.c;

		if((t2_v0 == v3 && t2_v1 == v4)||(t2_v0 == v4 && t2_v1 == v3))
			v2 = t2_v2;
		else if((t2_v1 == v3 && t2_v2 == v4)||(t2_v1 == v4 && t2_v2 == v3))
			v2 = t2_v0;
		else 
			v2 = t2_v1;

		//std::cout << v1 << " " << v2 << " " << v3 << " " << v4 << "\n";

		Eigen::Vector3d x1 = cloth->getMesh()->get_point_data(v1,last_frame_id);
		Eigen::Vector3d x2 = cloth->getMesh()->get_point_data(v2,last_frame_id);
		Eigen::Vector3d x3 = cloth->getMesh()->get_point_data(v3,last_frame_id);
		Eigen::Vector3d x4 = cloth->getMesh()->get_point_data(v4,last_frame_id);

		Eigen::Vector3d e_v = (x4-x3);
		float e = e_v.norm();

		Eigen::Vector3d a_n1_v = (x1-x3).cross(x1-x4);
		float a_n1 = a_n1_v.norm();

		Eigen::Vector3d a_n2_v = (x2-x4).cross(x2-x3);
		float a_n2 = a_n2_v.norm();

		Eigen::Vector3d n1 = a_n1_v.normalized();
		Eigen::Vector3d n2 = a_n2_v.normalized();
		Eigen::Vector3d ee = e_v.normalized();
		/*
		if(counter==12)
		{
			DBG__Vector_Dump(n1,counter,"N1");
			DBG__Vector_Dump(n2,counter,"N2");
			DBG__Vector_Dump(ee,counter,"EE");

		}
		*/
		float sin_theta = (n1.cross(n2)).dot(ee);
		bool sintheta_positive = (sin_theta<0.0f)?false:true;
	
		float term_below_sqr = 0.5f*(1.0f-(n1.dot(n2)));
		//std::cout << n1.dot(n2) << "\n";
		if(term_below_sqr<0.0)
			term_below_sqr = 1.0e-8;
		float sin_theta_half = sqrt(term_below_sqr);
		sin_theta_half = (sintheta_positive ? sin_theta_half : -sin_theta_half);

		float magnitude_shear = k_bend * sin_theta_half * ((e*e)/(a_n1+a_n2));

		//Now we calculate the values of the u1,u2,u2,u4
		Eigen::Vector3d u1 = (e/(a_n1*a_n1))*a_n1_v;
		Eigen::Vector3d f1 = magnitude_shear * u1;

		Eigen::Vector3d u2 = (e/(a_n2*a_n2))*a_n2_v;
		Eigen::Vector3d f2 = magnitude_shear * u2;

		Eigen::Vector3d u3 = (((x1-x4).dot(e_v))/(e*a_n1*a_n1))*a_n1_v + (((x2-x4).dot(e_v))/(e*a_n2*a_n2))*a_n2_v;
		Eigen::Vector3d f3 = magnitude_shear * u3;

		Eigen::Vector3d u4 = -(((x1-x3).dot(e_v))/(e*a_n1*a_n1))*a_n1_v - (((x2-x3).dot(e_v))/(e*a_n2*a_n2))*a_n2_v;
		Eigen::Vector3d f4 = magnitude_shear * u4;

		///if(counter==12)
		//{
		//	DBG__Vector_Dump(u1,counter,"U1");
		//	DBG__Vector_Dump(u2,counter,"U2");
		//	DBG__Vector_Dump(u3,counter,"U3");
			

		//}

		bend_forces[v1] += f1;
		bend_forces[v2] += f2;
		bend_forces[v3] += f3;
		bend_forces[v4] += f4;

		//DBG__Vector_Dump(bend_forces[v1],t1,"BEND");
		//DBG__Vector_Dump(bend_forces[v2],t1,"BEND");
		//DBG__Vector_Dump(bend_forces[v3],t1,"BEND");
		//DBG__Vector_Dump(bend_forces[v4],t1,"BEND");
	}
}

void FEM_Solver::compute_gravity_forces(Cloth_Data* cloth)
{
	int num_vertices = cloth->getMesh()->get_number_vertices();

	//Check size if needed resize then fill 
	if(gravity_forces.size()<num_vertices)
	{
		gravity_forces.resize(num_vertices);
	}
	
	OMP_FOR
	for(int p=0;p<num_vertices;p++)
	{
		float grav_force = cloth->get_vertex_mass(p) * 9.8f;
		gravity_forces[p] = Eigen::Vector3d(0.0,-grav_force,0.0);
	}

}

void FEM_Solver::computeDampingForce(Cloth_Data* cloth)
{
	int numTriangles = cloth->getMesh()->get_number_triangles();
	int numVertices = cloth->getMesh()->get_number_vertices();

	for(int p=0;p<numVertices;p++)
	{
		dampingForce_[p] = Eigen::Vector3d::Zero();
	}

	OMP_FOR
	for(int t=0; t<numTriangles; t++)
	{
		Triangles tri = cloth->getMesh()->get_triangle(t);

		Eigen::Vector3d velocityV0 = cloth->get_velocity(tri.a);
		Eigen::Vector3d velocityV1 = cloth->get_velocity(tri.b);
		Eigen::Vector3d velocityV2 = cloth->get_velocity(tri.c);

		//Implement a simple damping model 
		float bendC = cloth->get_property_obj()->get_damping_param();
		Eigen::Vector3d d12 = bendC*(velocityV1-velocityV2);
		Eigen::Vector3d d20 = bendC*(velocityV2-velocityV0);
		Eigen::Vector3d d01 = bendC*(velocityV0-velocityV1);

		dampingForce_[tri.b] += d12;
		dampingForce_[tri.c] -= d12;

		dampingForce_[tri.c] += d20;
		dampingForce_[tri.a] -= d20;

		dampingForce_[tri.a] += d01;
		dampingForce_[tri.b] -= d01;
	}
}

void FEM_Solver::get_next_velocity_positions(Cloth_Data* cloth)
{
	float time_step = cloth->get_property_obj()->get_timestep();
	int num_vertices = cloth->getMesh()->get_number_vertices();
	if(acceleration_.size()<num_vertices)
	{
		acceleration_.resize(num_vertices);
	}
	int lastVClicked = cloth->getLastClickedVertex();
	Eigen::Vector3d uiForce = cloth->getCurrentUIForce();
	OMP_FOR
	for(int p=0;p<num_vertices;p++)
	{
	  //DBG__Vector_Dump(shear_forces[p],p,"Shear");
		Eigen::Vector3d tot_force = shear_forces[p]+gravity_forces[p]+bend_forces[p]+dampingForce_[p];
		if(p==lastVClicked) {
			tot_force += uiForce;
		}
		float vertex_mass = cloth->get_vertex_mass(p);
		acceleration_[p] = tot_force * (1.0/vertex_mass);
		Eigen::Vector3d old_pos = cloth->getMesh()->get_point_data(p,last_frame_id);
		Eigen::Vector3d old_vel = cloth->get_velocity(p);
		Eigen::Vector3d new_vel = old_vel + (time_step * acceleration_[p]);
		Eigen::Vector3d new_pos = old_pos + (time_step * new_vel);
		if(p==2600 || p==50) {
			new_pos = old_pos;
			new_vel = Eigen::Vector3d::Zero();
		}
		cloth->set_next_step_pos(new_pos,p);
		cloth->set_next_step_velocity(new_vel,p);
	}
	//std::cout << "___________________________________________________\n";
	///DBG__Vector_Dump(cloth->get_velocity(cloth->getMesh()->get_triangle(0).a),0,"TR-A");
	//DBG__Vector_Dump(cloth->get_velocity(cloth->getMesh()->get_triangle(0).b),0,"TR-B");
	//DBG__Vector_Dump(cloth->get_velocity(cloth->getMesh()->get_triangle(0).c),0,"TR-C");
	//std::cout << "___________________________________________________\n";

}

void FEM_Solver::resetParameters()
{
	std::fill(shear_forces.begin(),shear_forces.end(),Eigen::Vector3d::Zero());
	std::fill(bend_forces.begin(),bend_forces.end(),Eigen::Vector3d::Zero());
	std::fill(gravity_forces.begin(),gravity_forces.end(),Eigen::Vector3d::Zero());
	std::fill(acceleration_.begin(),acceleration_.end(),Eigen::Vector3d::Zero());
	std::fill(dampingForce_.begin(),dampingForce_.end(),Eigen::Vector3d::Zero());
}

void FEM_Solver::initialize(Cloth_Data* cloth)
{
	int numVertices = cloth->getMesh()->get_number_vertices();
	shear_forces.resize(numVertices,Eigen::Vector3d::Zero());
	bend_forces.resize(numVertices,Eigen::Vector3d::Zero());
	gravity_forces.resize(numVertices,Eigen::Vector3d::Zero());
	acceleration_.resize(numVertices,Eigen::Vector3d::Zero());
	dampingForce_.resize(numVertices,Eigen::Vector3d::Zero());
}