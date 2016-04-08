#pragma once

#include <Eigen\Dense>
#include <vector>

class VolumeInformation
	{
	public:
		VolumeInformation(void);
		//The most important constructor
		VolumeInformation(std::vector<Eigen::Vector3d> pointData);
		~VolumeInformation(void);

		void createMitsubaVOLFile(const char* prefix, int id);


	private:
		//first we need the information for bounding box 
		float xmin_;
		float xmax_;
		float ymin_;
		float ymax_;
		float zmin_;
		float zmax_;
	};

