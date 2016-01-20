#pragma once

#include "triMesh.h"
#include "mesh_material.h"
#include "global_typedefs.h"

//A renderable object has the following things in it :
// A mesh and a material xml file  


class RenderObject 
{
public:
	RenderObject(SceneObjectId id,bool rotInv);
	~RenderObject();

	virtual void render(int frameId);

	void setMesh(TriMesh* mesh)					{ mesh_ = mesh;					}
	TriMesh* getMesh() const					{ return mesh_;					}

	void setRenderable(bool render)				{ render_ = render;				}
	bool isRenderable() const					{ return render_;				}

	void setLastClickedVertex(int p)   { lastClickedVertex_ = p;}
	int getLastClickedVertex()         { return lastClickedVertex_;}

	void setLastClickedPosition(Eigen::Vector3d pos) { lastClickedPosition_ = pos;}
	Eigen::Vector3d getLastClickedPosition()      { return lastClickedPosition_ ;}

	void setRenderMode(RenderMode rendermode)	{ renderMode_ = rendermode;		}
	RenderMode getRenderMode() const			{ return renderMode_;			}

	SceneObjectId getSceneObjectId()			{ return sceneId_;				}

	bool isRotationInvariant() const			{ return rotationInvariant_;	}

	int performObjectSelectionRayTest(double* rayStart, double* rayEnd, double *clickedWorldPos, double* selectedPos, int frameId);

	void updateUIForce(Eigen::Vector3d currPosition, int vertexId, int frameId);

	Eigen::Vector3d getCurrentUIForce()   { return UIForce_;} 

	void resetUIForce();



protected:
	TriMesh* mesh_;
	bool render_;
	SceneObjectId sceneId_;
	RenderMode renderMode_;
	bool rotationInvariant_;
	int lastClickedVertex_;
	Eigen::Vector3d lastClickedPosition_;
	Eigen::Vector3d UIForce_;
};