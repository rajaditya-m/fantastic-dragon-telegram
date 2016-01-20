#pragma once 

#include <cmath> 

inline bool is_equal(float x, float y)
{
	const float epsilon = 1.0e-10;
	return std::abs(x - y) <= epsilon * std::abs(x);
}

inline bool is_equal(double x, double y)
{
  const double epsilon = 1.0e-10;
  return std::abs(x - y) <= epsilon * std::abs(x);
}

inline double inverseSinc(double y) {
	int numIters = 100;
	double x = M_PI;
	for(int i=0;i<numIters;i++) {
		double fx = sin(x) - (x*y);
		double d_fx = (cos(x) - y);
		double xnew = x - (fx/d_fx);
		x = xnew;
	}
	return x;
}

inline double derivativeInverseSinc(double x) {
	double sincX = inverseSinc(x);
	double res = sincX/(cos(sincX)-x);
	return res;
}

inline double derivativeInverseSinc(double x, double sincX) {
	double res = sincX/(cos(sincX)-x);
	return res;
}

inline double derivativeSinc(double x) {
	double res = (cos(x)/x) - (sin(x)/(x*x));
	return res;
}

inline double cotTheta(const Eigen::Vector3d v, const Eigen::Vector3d w)
{
  const double cosTheta = v.dot(w);
  const double sinTheta = v.cross(w).norm();
  return (cosTheta / sinTheta);
}
 