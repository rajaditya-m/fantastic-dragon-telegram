#pragma once

#include <QGLWidget>
#include <QtOpenGL>
#include <gl/glu.h>
#include <math.h>

#include "global_typedefs.h"
#include "body_data.h"
#include "cloth_data.h"
#include "InactiveSupportObjects.h"
#include "ImplicitFEMSolver.h"
#include "ImplicitHyperElasticFEMSolver.h"
//#include "ImplicitMassSpringSolver.h"
#include "SimulationEngine.h"
#include "Scene.h"
#include "CollisionEngine.h"
#include "YarnOverlay.h"


class GLWidget : public QGLWidget
{
	Q_OBJECT

public:
	GLWidget(QWidget *parent=0);
	~GLWidget();

	QSize sizeHint() const;
	QSize minimumSizeHint() const;

protected:
	void initializeGL();
	void paintGL();
	void resizeGL(int width,int height);

	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);

public slots:
	void setRenderingFrameNumber(double frame_number);
	void setRenderingFrameNumber(int frame_number);
	void startAnimation();
	void pauseAnimation();
	void resetAnimation();
	void nextFrame();
	
	void bodyDisplayToggled(bool checked);
	void clothDisplayToggled(bool checked);
	void renderGroundToggled(bool checked);
	void renderAxesToggled(bool checked);
	void setBodyRenderInShading(bool checked);
	void setBodyRenderInWireframe(bool checked);
	void setClothRenderInShading(bool checked);
	void setClothRenderInWireframe(bool checked);
	void setClothRenderInHeatmapVelocity(bool checked);
	void setClothRenderInHeatmapAcceleration(bool checked);
	void setClothRenderInHeatmapBendingForce(bool checked);
	void setClothRenderInHeatmapDampingForce(bool checked);
	void setClothRenderInHeatmapShearForce(bool checked);
	void setSaveAsOBJCloth() ;
	void setSaveAsOBJBody();
	void setSaveAsOBJAll();
signals:
	void frameNumberUpdated(double frame_nos);
	void frameNumberUpdated(int frame_nos);
	void animFrameNumberUpdated(int frame_nos);


private:
	Body_Data *body_information_;
	Cloth_Data *cloth_information_;
	InactiveSupportObjects* support_;
	//ImplicitFEMSolver* fem_solver_;
	ImplicitHyperElasticFEMSolver* fem_solver_;
	ImplicitMassSpringSolver* mass_spring_solver_;
	SimulationEngine* sim_engine_;
	CollisionEngine* collisionEngine_;
	YarnOverlay *yarnOverlay_;


	//Animation timer
	QTimer* animation_timer;

	//Scence Data 
	Scene* scene_;
};
