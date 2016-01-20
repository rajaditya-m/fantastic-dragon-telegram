#include "RenderObject.h"

RenderObject::RenderObject(SceneObjectId id,bool rotInv)
{
	sceneId_ = id;
	rotationInvariant_ = rotInv;
	lastClickedVertex_ = -1;
	lastClickedPosition_ = Eigen::Vector3d::Zero();
}

RenderObject::~RenderObject()
{
}

void RenderObject::render(int frameId)
{
	if(render_)
	{
	 
		switch(renderMode_)
		{
		case WIREFRAME : mesh_->draw_wireframe(frameId); break;
		case SHADING : mesh_->draw_mesh(frameId);  break; 
		default : std::cout << "[ERROR] The selected mode of rendering is not supported for this object.\n"; break;
		}
	}
}

int RenderObject::performObjectSelectionRayTest(double* rayStart, double* rayEnd, double *clickedWorldPos, double* selectedPos, int frameId) {
	return mesh_->performMeshSelectionRayTest(rayStart, rayEnd, clickedWorldPos, selectedPos,frameId);
}

void RenderObject::updateUIForce(Eigen::Vector3d currPosition, int vertexId, int frameId) {
	Eigen::Vector3d origin = mesh_->get_point_data(vertexId,frameId);
	Eigen::Vector3d diff = currPosition-origin;
	double compliance = 1.0;
	diff *= compliance;
	UIForce_ = diff;
	std::cout << "Vertex Clicked: " << vertexId << " [" << diff[0] << "," << diff[1] << "," << diff[2] << "]\n";
}

void RenderObject::resetUIForce() {
	UIForce_ = Eigen::Vector3d::Zero();
	lastClickedPosition_ = Eigen::Vector3d::Zero();
	lastClickedVertex_ = -1;
}