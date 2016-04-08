#pragma once

#include <vector>
#include <utility>
#include <Eigen\Dense>
#include "cloth_data.h"
#include "RenderObject.h"

class YarnOverlay : public RenderObject
	{

	public:
		YarnOverlay(void);
		//This constructor is too basic we need something better to handle any cloth orientation 
		YarnOverlay(int numWeftYarns, int numWarpYarns, Cloth_Data* clothData);
		~YarnOverlay(void);

		//The update function for yarn overlay 
		void addFrameInformation(Cloth_Data* clothData);

	private:
		int numCrossingNodes_;
		std::vector<std::pair<int,Eigen::Vector3d> > barycentricInformation_;
	};

