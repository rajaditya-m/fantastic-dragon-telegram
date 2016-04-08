#include "YarnOverlay.h"


YarnOverlay::YarnOverlay(void)
	:RenderObject(YARN_OVERLAY,false)
	{
	}

YarnOverlay::YarnOverlay(int numWeftYarns, int numWarpYarns,Cloth_Data* clothData):
	RenderObject(YARN_OVERLAY,false)
{
	numCrossingNodes_ = numWarpYarns*numWeftYarns;
	std::vector<Eigen::Vector3d> pd = clothData->getMesh()->getPositionVector(0);
	//For now just hack this by using the bounding box parameters
	Pair_Vec3d bbData =  clothData->getMesh()->get_super_bounding_box();

	Eigen::Vector3d minBB = bbData.first;
	Eigen::Vector3d maxBB = bbData.second;

	//Assume that everythning is in the XZ plane 
	double minX = minBB.x();
	double minZ = minBB.z();
	double minY = minBB.y();

	double maxX = maxBB.x();
	double maxZ = maxBB.z();

	//Now we do the main assignment of the yarns 
	//Assume that Warp yarns run in the X direction 
	// Weft yarsn run in the Z direction
	// Everythin in column order 

	/*
	(minX,minZ) ----WARP------(maxX,minZ)
	|                                 |
	|                                 |
	|                                 |
	WEFT     -------WARP--------   WEFT
	|                                 |
	|                                 |
	|                                 |
	(minX,maxZ) ----WARP0-----(maxX,maxZ)
	*/
	
	double warpDelta = (maxX-minX)/(numWarpYarns-1);
	double weftDelta = (maxZ-minZ)/(numWeftYarns-1);

	std::vector<Eigen::Vector3d> pvec;
	for(int i=0;i<numWarpYarns;i++)
	{
		for(int j=0;j<numWeftYarns;j++)
		{
			Eigen::Vector3d pntVec(minX+i*warpDelta,minY,minZ+j*weftDelta);
			pvec.push_back(pntVec);
		}
	}

	//Now time to figure out the rendering triangles 
	std::vector<int> tri;
	for(int i=0;i<numWarpYarns-1;i++)
	{
		for(int j=0;j<numWeftYarns-1;j++)
		{
			int a = i*(numWarpYarns)+j;
			int b = i*(numWarpYarns)+(j+1);
			int c = (i+1)*(numWarpYarns)+(j+1);
			int d = (i+1)*(numWarpYarns)+j;
			tri.push_back(a);
			tri.push_back(b);
			tri.push_back(c);
			tri.push_back(a);
			tri.push_back(c);
			tri.push_back(d);
		}
	}
	mesh_ = new TriMesh();
	mesh_->initializeDynamicMesh(pvec,tri,"Yarn","body_material.xml");

	//Now is the part where we generate the barycentric coordinates for the triangles 
	int numClothTriangles = clothData->getMesh()->get_number_triangles();
	for(int i=0;i<numCrossingNodes_;i++)
	{
		Eigen::Vector3d p = pvec[i];
		Eigen::Vector3d bary;
		int baryTri;
		float u,v,w;
		for(int t=0;t<numClothTriangles;t++)
		{
			Triangles tC = clothData->getMesh()->get_triangle(t);
			Eigen::Vector3d a = clothData->getMesh()->get_point_data(tC.a);
			Eigen::Vector3d b = clothData->getMesh()->get_point_data(tC.b);
			Eigen::Vector3d c = clothData->getMesh()->get_point_data(tC.c);
			barycentricCoods(p,a,b,c,u,v,w);
			if(u>=0.0 && v>=0.0 && w>=0.0)
			{
				bary  = Eigen::Vector3d(u,v,w);
				baryTri = t;
				break;
			}
		}
		barycentricInformation_.push_back(std::pair<int,Eigen::Vector3d>(baryTri,bary));

	}
	//std::vector<Triangles> clothTri = clothData->getMesh()->get_triangle();

}


YarnOverlay::~YarnOverlay(void)
	{
	}

void YarnOverlay::addFrameInformation(Cloth_Data* clothData)
{
	int numFrames = clothData->getMesh()->get_number_frames();
	//Compute all the barycentric points 
	std::vector<Eigen::Vector3d> pointVec;
	for(int i=0;i<numCrossingNodes_;i++)
	{
		Eigen::Vector3d baryCoods = barycentricInformation_[i].second;
		int tri = barycentricInformation_[i].first;
		
		Triangles tCloth = clothData->getMesh()->get_triangle(tri);
		Eigen::Vector3d a = clothData->getMesh()->get_point_data(tCloth.a,numFrames-1);
		Eigen::Vector3d b = clothData->getMesh()->get_point_data(tCloth.b,numFrames-1);
		Eigen::Vector3d c = clothData->getMesh()->get_point_data(tCloth.c,numFrames-1);
		Eigen::Vector3d baryT = baryCoods.x()*a + baryCoods.y()*b + baryCoods.z()*c;
		DBG__Vector_Dump(baryT);
		pointVec.push_back(baryT);
	}
	mesh_->add_frame_data(pointVec);
}
