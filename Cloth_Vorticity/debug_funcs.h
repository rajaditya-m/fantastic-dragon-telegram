#include <iostream>
#include <fstream>
#include <Eigen\Dense>

inline void DBG__Vector_Dump(Eigen::Vector3d vec)
{
	std::cout << "[" << vec[0] << "," << vec[1] << "," << vec[2] << "]\n";
}

inline void DBG__Vector_Dump(Eigen::Vector3d vec,int seq_num,const char* header="")
{
	std::cout << header << " Seq#" << seq_num << " [" << vec[0] << "," << vec[1] << "," << vec[2] << "]\n";
}

//Assume column major storage
inline void DBG__Matrix_Dump(double* rawSpace, int rows, int cols)
{
	for(int r=0;r<rows;r++)
	{
	  std::cout << "[";
		for(int c=0;c<cols;c++)
		{
			std::cout << rawSpace[c*rows+r] << " ";
		}
		std::cout << "]\n";
	}
	std::cout << "\n";
}

inline void DBG__Save_To_File(int numElems, double *vec, const char *fileName) {
	std::ofstream fileOut(fileName);
	for(int i=0;i<numElems;i++) {
		fileOut << vec[i] << "\n";
	}
	fileOut.close();
}
