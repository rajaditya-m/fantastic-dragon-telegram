#include "VolumeInformation.h"
#include <fstream>
#define GRIDRESOLUTION 1000

VolumeInformation::VolumeInformation(void)
	{
	}


VolumeInformation::~VolumeInformation(void)
	{
	}

VolumeInformation::VolumeInformation(std::vector<Eigen::Vector3d> pointData)
{
  // First we find the bounding box of our common grid elements 
  float xmn=9999.99,ymn=9999.9,zmn=9999.9;
	float xmx=-9999.99,ymx=-9999.99,zmx=-9999.99;
	for(int i=0;i<pointData.size();i++)
	{
		Eigen::Vector3d pdata = pointData[i];
		if(pdata.x()<xmn)
			xmn = pdata.x();
		if(pdata.y()<ymn)
			ymn = pdata.y();
		if(pdata.z()<zmn)
			zmn = pdata.z();
		if(pdata.x()>xmx)
			xmx = pdata.x();
		if(pdata.y()>ymx)
			ymx = pdata.y();
		if(pdata.z()>zmx)
			zmx = pdata.z();
	}
	xmin_ = xmn;
	ymin_ = ymn;
	zmin_ = zmn;
	xmax_ = xmx;
	ymax_ = ymx;
	zmax_ = zmx;

}

void VolumeInformation::createMitsubaVOLFile(const char* prefix, int id)
{
  //First create the filename and the filepointer
  char fileName[512];
	sprintf(fileName,"%z%z",prefix,id);
	std::ofstream outfile(fileName,std::ios::binary|std::ios::out);

	//The ascii characters ;VOL'
	char header[] = {'V','O','L'};
	outfile.write(reinterpret_cast<char*>(&header),sizeof(&header));

	//Next we create the voxelization as needed 
	float xDelta = (xmax_ - xmin_)/GRIDRESOLUTION;
	float yDelta = (ymax_ - ymin_)/GRIDRESOLUTION;
	float zDelta = (zmax_ - zmin_)/GRIDRESOLUTION;

	for(int xCount = 0;xCount < GRIDRESOLUTION;xCount++)
	{
		for(int yCount = 0;yCount < GRIDRESOLUTION;yCount++)
		{
			for(int zCount = 0;zCount < GRIDRESOLUTION;zCount++)
			{
				
			}
		}
	}


	//Dotn forget to close the file pointer 
	outfile.close();
}
