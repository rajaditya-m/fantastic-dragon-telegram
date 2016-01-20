#pragma once

#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>




class ClothProperties
{
public:
	//The constructor to read the XML file and parse the result 
	ClothProperties(const char* property_xml);
	~ClothProperties(void);

	//Accessor methods 
	float get_youngs_modulus() const			{ return youngs_modulus;		}
	float get_poisson_ratio() const				{ return poisson_ratio;			}
	float get_timestep() const					{ return time_step;				}
	float get_thickness() const					{ return thickness;				}
	float get_max_solver_iteration() const		{ return max_solver_iterations;	}
	float get_friction() const					{ return friction_coeff;		}
	float get_damping_param() const				{ return damping_param;			}
	float get_density() const					{ return density;				}
	float get_k_bend() const					{ return k_bend;				}
	float getKStiffness() const { return kStiffness_;}
	float getQuadBendStiffness() const { return qBStiffness_;}
	float getQuadBendDamping()  const { return qBDamping_;}
	std::string getDataDumpLoc() const { return dataDumpLocation_;}
	float getRayleighDampingCoeff() const { return rayDampCoeff_;}
	float getRayleighMassCoeff() const { return rayMassCoeff_;}

private:
	float youngs_modulus;
	float poisson_ratio;
	float time_step;
	float thickness;
	float max_solver_iterations;
	float friction_coeff;
	float damping_param;
	float k_bend;
	float density;
	float kStiffness_;
	std::string dataDumpLocation_;
	float qBStiffness_;
	float qBDamping_;
	float rayDampCoeff_;
	float rayMassCoeff_;

};