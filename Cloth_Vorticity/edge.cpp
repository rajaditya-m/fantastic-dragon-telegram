#include "edge.h"



Edge::Edge(void)
{
	start = end = -1;
	adj_tri_1 = adj_tri_2 = -1;
	index = -1;
}

void Edge::add_edges(int a,int b)
{
	if(a<b) 
	{
		start = a;
		end = b;
	}
	else
	{
		start = b;
		end = a;
	}
}

bool Edge::add_triangle(int tri_idx)
{
	if(adj_tri_1==-1)
	{
		adj_tri_1 = tri_idx;
		return true;
	}
	else if(adj_tri_2==-1)
	{
		adj_tri_2 = tri_idx;
		return true;
	}
	else
	{
		return false;
	}
}

bool Edge::operator< (const Edge e) const
{
	if(this->start < e.start)
		return true;
	else if(this->start == e.start)
	{
		if(this->end < e.end)
			return true;
		else
			return false;
	}
	else return false;
}

void Edge::set_index(int idx)
{
	index = idx;
}

int Edge::get_tri_1()const 
{
	return adj_tri_1;
}

int Edge::get_tri_2()const 
{
	return adj_tri_2;
}

int Edge::get_start_vertex() const
{
	return start;
}

int Edge::get_end_vertex() const
{
	return end;
}

Edge::~Edge(void)
{
}
