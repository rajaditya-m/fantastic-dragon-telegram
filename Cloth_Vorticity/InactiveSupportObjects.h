#pragma once

#include "RenderObject.h"

#define PREDEF_SCALING 15.0

class InactiveSupportObjects : public RenderObject
	{
	public:
		InactiveSupportObjects(const char* name, const char* file_location, const char* material_xml);
		~InactiveSupportObjects(void);
	};

