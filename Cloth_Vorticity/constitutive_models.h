#pragma once 

#include <Eigen\Dense>
#include <Eigen\SVD>

inline Eigen::Matrix2f corotational_elasticity(Eigen::Matrix2f &m,float lame_param_lambda, float lame_param_mu)
{
	//Compute polar decomposition of the matrix
	Eigen::JacobiSVD<Eigen::Matrix2f> svd(m);
	Eigen::Matrix2f u = svd.matrixU();
	Eigen::Matrix2f v = svd.matrixV();
	Eigen::Vector2f s = svd.singularValues();
	
	Eigen::Matrix2f I = Eigen::Matrix2f::Identity();
	Eigen::Matrix2f s_mat = Eigen::Matrix2f::Zero();
	s_mat(0,0) = s(0);
	s_mat(1,1) = s(1);

	Eigen::Matrix2f r = v * s_mat * v.transpose();
	

	Eigen::Matrix2f p_f = (2.0*lame_param_mu*(m-r)) + (lame_param_lambda*(r.transpose()*m - I).trace()*r);

	return p_f;

}


inline Eigen::Matrix2f linear_elasticity(Eigen::Matrix2f &m, float lame_param_lambda,float lame_param_mu)
{
	Eigen::Matrix2f I = Eigen::Matrix2f::Identity();

	Eigen::Matrix2f p_f = (lame_param_mu*(m + m.transpose() - 2.0*I)) + (lame_param_lambda * ((m-I).trace()) * I); 

	return p_f;
}

inline Eigen::Matrix3d isotropic_linear_elasticity(float ym,float poisson)
{
	Eigen::Matrix3d E = Eigen::Matrix3d::Zero();
	E(0,0) = 1.0f; E(0,1) = poisson; 
	E(1,0) = poisson ; E(1,1) = 1.0f;
	E(2,2) = (1.0f-poisson)*0.5f;
	E *= (ym/(1.0f-(poisson*poisson))) ;
	//std::cout << E << "\n";
	return E;

}

inline Eigen::MatrixXd genericPolarDecomposition(Eigen::MatrixXd J)
{
	Eigen::MatrixXd U_Squared = J.transpose()*J;
	Eigen::EigenSolver<Eigen::MatrixXd> es(U_Squared);

	double lambda1 = es.eigenvalues()[0].real();
	double lambda2 = es.eigenvalues()[1].real();
	
	Eigen::VectorXd v1 = es.eigenvectors().col(0).real();
	Eigen::VectorXd v2 = es.eigenvectors().col(1).real();

	Eigen::MatrixXd U = std::sqrt(lambda1)*v1*v1.transpose() + std::sqrt(lambda2)*v2*v2.transpose();
	Eigen::MatrixXd U_Inv = U.inverse();
	Eigen::MatrixXd R  =J * U_Inv;

	return R;
}