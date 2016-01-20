#pragma once

class Edge
{
public:
	Edge(void);
	~Edge(void);
	void add_edges(int a,int b);
	bool operator< (const Edge e) const;
	bool add_triangle(int tri_idx);
	void set_index(int idx);
	int get_tri_1() const;
	int get_tri_2() const;
	int get_start_vertex() const;
	int get_end_vertex() const;
	int get_index() const					{ return index;			}
public:
	int start;
	int end;
	int adj_tri_1;
	int adj_tri_2;
	int index;
};