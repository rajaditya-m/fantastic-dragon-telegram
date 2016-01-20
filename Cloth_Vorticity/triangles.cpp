#include "triangles.h"


Triangles::Triangles(void)
{
	a = b = c = -1;
}

Triangles::Triangles(int _a, int _b, int _c)
{
	a = _a;
	b = _b;
	c = _c;
}

void Triangles::add(int _a,int _b, int _c)
{
	a = _a;
	b = _b;
	c = _c;
}

void Triangles::add_textures(int _at,int _bt,int _ct)
{
	at = _at;
	bt = _bt;
	ct = _ct;
}

void Triangles::add_lmt_value(int x,int y)
{
	lmt_x = x;
	lmt_y = y;
}

Triangles::~Triangles(void)
{
}
