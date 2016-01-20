#include "body_data.h"
#include "triMesh.h"
#include "RenderObject.h"


Body_Data::Body_Data(const char* name,const char* file_location,const char* material_xml,int seq_limit,int start_frame)
	:RenderObject(COLL_OBJ,false)
{
	mesh_ = new TriMesh(name,file_location,material_xml,seq_limit,1.0);
}

Body_Data::Body_Data(const char* name, const char* file_location_vertex, const char* file_location_mesh, const char* file_location_texture, const char* material_xml)
	:RenderObject(COLL_OBJ,false)
{
	mesh_ = new TriMesh(name,file_location_vertex,file_location_mesh,file_location_texture,material_xml,1.0);
	mesh_->setStatic(true);
}

Body_Data::Body_Data(const char* name, const char* file_location, const char* material_xml) 
	:RenderObject(COLL_OBJ,false){
		mesh_ = new TriMesh(name,file_location, material_xml,1.0);
		mesh_->setStatic(true);
}

Body_Data::~Body_Data(void)
{
}

