#pragma once

#define RADPERDEG 0.0174533
#define TEX_ROWS 256
#define TEX_COLS 256

#include <QGLWidget>
#include <QtOpenGL>
#include <gl/glu.h>
#include <math.h>
#include <vector>
#include <chrono>
#include <boost\filesystem.hpp>
#include <boost\algorithm\string.hpp>
#include "RenderObject.h"


//This will contain all the data required to render the scene object along with the requisite 
// Information 
//@TODO : get onmi present lighting information 
//@TODO : find better ways to input the eye position right now we will work with the defaults 
//@TODO : support multiple lighting position

class Scene :public QGLWidget
{
public:
	int getXRot() const						{ return xRot_;				}
	int getYRot() const						{ return yRot_;				}
	int getZRot() const						{ return zRot_;				}
	float getScaling() const				{ return scaling_;			}
	QPoint getLastPosition() const			{ return lastPos_;			}
	int getCurrentRenderFrame() const		{ return renderFrame_;		} 
	bool renderGround() const				{ return renderGround_;		}
	bool renderAxes() const					{ return renderAxes_;		}
	int getLastVertexClicked() const { return lastVertexClicked_;}

	
	//void setLastWorldPosition(Eigen::Vector3d v)  { lastWorldPos_ = v;}
	void setXRot(int x)						{ xRot_ = x;				}
	void setYRot(int y)						{ yRot_ = y;				}
	void setZRot(int z)						{ zRot_ = z;				}
	void setScaling(float s)				{ scaling_ = s;				}
	void setLastPosition(QPoint p)			{ lastPos_ = p;				}
	void setCurrentRenderFrame(int f)		{ renderFrame_ = f;			}
	void setRenderGround(bool rg)			{ renderGround_ = rg;		}
	void setRenderAxes(bool ra)				{ renderAxes_ = ra;			}
	void setAllRotations(int x, int y, int z);
	void setLastVertexClicked(int x)   { lastVertexClicked_ = x;}


	void addRenderObject(RenderObject* obj);

	int performSelectionRayTest(double* rayStart, double* rayEnd, double *clickedWorldPos, double* selectedPos); 
	void updateDraggingForce(double x, double y, double z);
	void resetDraggingForce();

	//Save as functions 
	void saveObjs(int code);


	//update the prefixes
	void updateOBJPrefixes(std::string &prefix);


	//Rendering Functions 
	void renderScene();
	void setupScene();
	void resizeScene(int width, int height);

	//Helper Objects in the screen 
	void drawAxes(float length);
	void drawGradientBackGround() ;
	void drawGround();

	Scene(void);
	~Scene(void);

private:
	GLubyte* readAsciiPPMImage(char* fName,int* height,int* width);
	void Arrow(float x1,float y1,float z1,float x2,float y2,float z2,float D);

private:
	//General Mouse Manupulation Options 
	int xRot_;
	int yRot_;
	int zRot_;
	float scaling_;
	QPoint lastPos_;
	int lastVertexClicked_;
	int clickedVertexHandle_;
	Eigen::Vector3d arrowDest_;

	//Vector of Rendering Objects 
	std::vector<RenderObject*> renderObj_;

	//Lighting Information (At this point single lighting is supported) 
	Eigen::Vector3d lightPos1_;

	//Eye Position Information 
	Eigen::Vector3d eyePos_;
	Eigen::Vector3d centerPos_;
	Eigen::Vector3d upDirection_;

	//Perspective Projection Settings 
	Eigen::Vector3d perspective_;

	//Floor texture data 
	GLuint groundTex_;
	GLubyte* textureData_;
	//GLubyte textureData_[TEX_ROWS][TEX_COLS][4];

	//Current render frame
	int renderFrame_;

	//Some generic boolean flipping 
	bool renderAxes_;
	bool renderGround_;

	//Some more saving information 
	std::string clothPrefix_;
	std::string bodyPrefix_;
	bool dirsCreated_;
};

