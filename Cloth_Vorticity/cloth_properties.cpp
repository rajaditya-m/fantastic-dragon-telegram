#include "cloth_properties.h"


ClothProperties::ClothProperties(const char* property_xml)
{
	boost::property_tree::ptree pt;
	boost::property_tree::read_xml(property_xml,pt);

	//Now read everything one by one 
	youngs_modulus = pt.get<float>("property.youngs_modulus");
	poisson_ratio = pt.get<float>("property.poisson_ratio");
	time_step = pt.get<float>("property.time_step");
	thickness = pt.get<float>("property.thickness");
	max_solver_iterations = pt.get<float>("property.max_iter");
	friction_coeff = pt.get<float>("property.friction_coeff");
	damping_param = pt.get<float>("property.damping_param");
	density = pt.get<float>("property.density");
	k_bend = pt.get<float>("property.k_bend");
	kStiffness_ = pt.get<float>("property.k_stiffness");
	qBDamping_ = pt.get<float>("property.quad_damping");
	qBStiffness_ = pt.get<float>("property.quad_stiffness");
	dataDumpLocation_ = pt.get<std::string>("property.data_dump_location");
	rayDampCoeff_ = pt.get<float>("property.rayleigh_damp_coeff");
	rayMassCoeff_ = pt.get<float>("property.rayleigh_mass_coeff");
}


ClothProperties::~ClothProperties(void)
{
}
