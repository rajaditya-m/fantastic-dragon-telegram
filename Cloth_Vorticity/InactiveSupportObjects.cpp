#include "InactiveSupportObjects.h"


InactiveSupportObjects::InactiveSupportObjects(const char* name, const char* file_location, const char* material_xml)
	:RenderObject(UNDEFINED,false){
		mesh_ = new TriMesh(name,file_location, material_xml,1.0);
		mesh_->setStatic(true);
}


InactiveSupportObjects::~InactiveSupportObjects(void)
	{
	}
