#pragma once

#include <Eigen\Dense>

inline float triangle_area(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c)
{
	Eigen::Vector3d vBA = b-a;
	Eigen::Vector3d vCA = c-a;
	Eigen::Vector3d crossPdk = vBA.cross(vCA);
	float area = crossPdk.norm()*0.5f;
	return area;
}

inline Eigen::Vector2f compute_outward_normal(Eigen::Vector2f a, Eigen::Vector2f b, Eigen::Vector2f c)
{
	Eigen::Vector2f ab = b-a;
	Eigen::Vector2f ac = c-a;
	Eigen::Vector2f N(ab.y(),-ab.x());

	float d = N.dot(ac);

	if(d>0.0)
		N = -N;

	N.normalized();
	
	return N;
}

inline Eigen::Vector3d computeTriangleNormals(Eigen::Vector3d a,Eigen::Vector3d b,Eigen::Vector3d c)
{
	Eigen::Vector3d ab = b-a;
	Eigen::Vector3d ac = c-a;
	Eigen::Vector3d normal = ab.cross(ac);
	normal = normal.normalized();
	return normal;
}


inline Eigen::Matrix3d compute_rotation_matrix_a_to_b(Eigen::Vector3d norm_in_a,Eigen::Vector3d norm_in_b)
{
	float cos_theta = norm_in_a.dot(norm_in_b);
	
	Eigen::Vector3d unit_axis = norm_in_a.cross(norm_in_b);
	unit_axis.normalized();

	float x = unit_axis.x();
	float y = unit_axis.y();
	float z = unit_axis.z();

	float c = cos_theta;
	float term_under_sqrt = 1.0f-c*c;
	if(term_under_sqrt<0.0)
		std::cout << "[ERROR] Error while computing the rotation matrix\n";
	float s = sqrt(term_under_sqrt);
	float C = 1-c;

	Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Zero();

	rot_mat(0,0) = x*x*C+c;
	rot_mat(0,1) = x*y*C-z*s;
	rot_mat(0,2) = x*z*C+y*s;

	rot_mat(1,0) = y*x*C+z*s;
	rot_mat(1,1) = y*y*C+c;
	rot_mat(1,2) = y*z*C-x*s;

	rot_mat(2,0) = z*x*C-y*s;
	rot_mat(2,1) = z*y*C+x*s;
	rot_mat(2,2) = z*z*C+c;

	return rot_mat;
}

inline bool compare_point_leq(Eigen::Vector3d a, Eigen::Vector3d b)
{
	for(int i=0;i<3;i++)
	{
		if(a(i)>b(i))
			return false;
	}
	return true;
}

inline bool compare_point_geq(Eigen::Vector3d a, Eigen::Vector3d b)
{
	for(int i=0;i<3;i++)
	{
		if(a(i)<b(i))
			return false;
	}
	return true;
}

inline bool segmentTriangleIntersection(double* segment_start, double* segment_end,
																 double* tri_v0, double* tri_v1, double* tri_v2,
																 double threshold, double& u, double& v, double& t) {
	const double kEpsilon = 1e-10;
	Eigen::Vector3d d = Eigen::Map<Eigen::Vector3d>(segment_end) - Eigen::Map<Eigen::Vector3d>(segment_start);
	Eigen::Vector3d e1,e2,pvec;
	e1 = Eigen::Map<Eigen::Vector3d>(tri_v1) - Eigen::Map<Eigen::Vector3d>(tri_v0);
	e2 = Eigen::Map<Eigen::Vector3d>(tri_v2) - Eigen::Map<Eigen::Vector3d>(tri_v0);
	pvec = d.cross(e2);
	double det = e1.dot(pvec);
	if(std::abs(det) <= kEpsilon) 
		return false;

	double invDet = 1.0/det;

	Eigen::Vector3d tvec = Eigen::Map<Eigen::Vector3d>(segment_start) - Eigen::Map<Eigen::Vector3d>(tri_v0);

	u = tvec.dot(pvec) * invDet;

	if(u <-threshold || 1.0 + threshold < u) 
		return false;

	Eigen::Vector3d qvec = tvec.cross(e1);
	v = d.dot(qvec) * invDet;
	if (-threshold > v || 1.0 + threshold < u + v )
		return false;
	t = e2.dot(qvec) * invDet;
	if (t < -threshold || t >= 1.0 + threshold)
		return false;

	return true;
}

inline double pointLineDistance(double *p0, double *p1, double *p, double* length = NULL) {
	Eigen::Vector3d vec1 = Eigen::Map<Eigen::Vector3d>(p1) - Eigen::Map<Eigen::Vector3d>(p0);
	vec1.normalize();
	Eigen::Vector3d vec2 = Eigen::Map<Eigen::Vector3d>(p) - Eigen::Map<Eigen::Vector3d>(p0);
	double projectionLength = vec1.dot(vec2);
	if (length) {
		*length = projectionLength;
	}
	double vec2Norm = vec2.norm();
	double distance = 0;
	distance = std::sqrt((vec2Norm*vec2Norm) - (projectionLength*projectionLength));
	return distance;
}

inline Eigen::Vector3d findOrthonormalVector(Eigen::Vector3d &vec)
{
  // find smallest abs component of v
  int smallestIndex = 0;
  for(int dim=1; dim<3; dim++)
    if (fabs(vec[dim]) < fabs(vec[smallestIndex]))
      smallestIndex = dim;

  Eigen::Vector3d axis(0.0, 0.0, 0.0);
  axis[smallestIndex] = 1.0;

  // this cross-product will be non-zero (as long as v is not zero)
  Eigen::Vector3d result = vec.cross(axis);
	result.normalize();
  return result;
}

inline void barycentricCoods(Eigen::Vector3d p, Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c, float &u, float &v, float &w)
{
    Eigen::Vector3d v0 = b - a, v1 = c - a, v2 = p - a;
    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0f - v - w;
}