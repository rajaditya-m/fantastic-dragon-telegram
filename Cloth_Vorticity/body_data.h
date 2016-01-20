#pragma once

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

#include <QGLWidget>

#include <Eigen\Dense>

#include "triangles.h"
#include "struct_decls.h"
#include "geom_funcs.h"
#include "triMesh.h"
#include "RenderObject.h"

#define PREDEF_SCALING 15.0

class Body_Data : public RenderObject
{
public:
	Body_Data(const char* name,const char* file_location,const char* material_xml,int seq_limit,int start_frame=0);
	Body_Data(const char* name, const char* file_location_vertex, const char* file_location_mesh, const char* file_location_texture, const char* material_xml);
	Body_Data(const char* name, const char* file_location, const char* material_xml);
	~Body_Data(void);

};

