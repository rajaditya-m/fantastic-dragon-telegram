#pragma once

#include "cloth_data.h"
#include "body_data.h"

class CollisionEngine
{
public:
	CollisionEngine(void);
	~CollisionEngine(void);

	//Collision Resolution between body and cloth 
	void resolveClothBodyCollision(Cloth_Data* cloth,Body_Data* body);
	void resolveClothClothCollision(Cloth_Data* cloth);

};

