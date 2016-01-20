#include <Eigen\Dense>
#include <Eigen\Sparse>

inline void insert3x3MatrixBlock(Eigen::SparseMatrix<double,Eigen::RowMajor> &spMat, Eigen::Matrix3d &block, int rowStart, int colStart) {
	spMat.coeffRef(rowStart,colStart) += block(0,0);
	spMat.coeffRef(rowStart,colStart+1) += block(0,1);
	spMat.coeffRef(rowStart,colStart+2) += block(0,2);
	
	spMat.coeffRef(rowStart+1,colStart) += block(1,0);
	spMat.coeffRef(rowStart+1,colStart+1) += block(1,1);
	spMat.coeffRef(rowStart+1,colStart+2) += block(1,2);

	spMat.coeffRef(rowStart+2,colStart) += block(2,0);
	spMat.coeffRef(rowStart+2,colStart+1) += block(2,1);
	spMat.coeffRef(rowStart+2,colStart+2) += block(2,2);

}

inline void insert3x3MatrixBlock(Eigen::SparseMatrix<double,Eigen::RowMajor> &spMat, Eigen::Matrix3d &block, int rowStart, int colStart, double constMult) {
	spMat.coeffRef(rowStart,colStart) += (constMult*block(0,0));
	spMat.coeffRef(rowStart,colStart+1) += (constMult*block(0,1));
	spMat.coeffRef(rowStart,colStart+2) += (constMult*block(0,2));
	
	spMat.coeffRef(rowStart+1,colStart) += (constMult*block(1,0));
	spMat.coeffRef(rowStart+1,colStart+1) += (constMult*block(1,1));
	spMat.coeffRef(rowStart+1,colStart+2) += (constMult*block(1,2));

	spMat.coeffRef(rowStart+2,colStart) += (constMult*block(2,0));
	spMat.coeffRef(rowStart+2,colStart+1) += (constMult*block(2,1));
	spMat.coeffRef(rowStart+2,colStart+2) += (constMult*block(2,2));

}