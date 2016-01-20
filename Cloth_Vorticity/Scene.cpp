#include "Scene.h"


Scene::Scene(void)
{
	//Scene data
	xRot_ = yRot_ = zRot_ = 0;
	scaling_ = 1.0f;

	//Eye Data Default 
	eyePos_ = Eigen::Vector3d(0.0f,13.0f,8.0f);
	centerPos_ = Eigen::Vector3d(0.0f,0.0f,0.0f);
	upDirection_ = Eigen::Vector3d(0.0f,1.0f,0.0f);

	//Lighting Data Default 
	lightPos1_ = Eigen::Vector3d(0.0f,10.0f,3.0f);

	//Current RenderFrame
	renderFrame_ = 0;

	//Default rendering data 
	renderGround_ = true;
	renderAxes_ = false;

	textureData_ = NULL;

	lastVertexClicked_ = -1;
	clickedVertexHandle_ = 0;

	//Some more information
	const char* suffix = "F:\\Cloth_Sim_Data";
	//Get the current date and time also and create a folder 
	std::chrono::time_point<std::chrono::system_clock> timeNow;
	timeNow = std::chrono::system_clock::now();
	std::time_t timeNowT = std::chrono::system_clock::to_time_t(timeNow);
	const char* timeStr = std::ctime(&timeNowT);
	std::string timeS(timeStr);
	boost::replace_all(timeS," ","_");
	boost::replace_all(timeS,":","_");
	boost::trim(timeS);
	std::stringstream ssc;
	ssc << suffix << "\\" << timeS << "\\Cloth";
	clothPrefix_ = ssc.str();
	std::stringstream ssb;
	ssb << suffix << "\\" << timeS << "\\Body";
	bodyPrefix_ = ssb.str();
	dirsCreated_ = false;
}


Scene::~Scene(void)
{
}

void Scene::setAllRotations(int x,int y, int z)
{
	xRot_ = x;
	yRot_ = y;
	zRot_ = z;
}

void Scene::renderScene()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
	glLoadIdentity();
		glTranslatef(0.0, 0.0, -8.0);
		glRotatef(xRot_/16.0, 1.0, 0.0, 0.0);
		glRotatef(yRot_/16.0, 0.0, 1.0, 0.0);
	gluLookAt(eyePos_.x(),eyePos_.y(),eyePos_.z(),centerPos_.x(),centerPos_.y(),centerPos_.z(),upDirection_.x(),upDirection_.y(),upDirection_.z());
	glScalef(scaling_,scaling_,scaling_);

	drawGradientBackGround();
	if(renderAxes_)
		drawAxes(16.0);
	if(renderGround_)
		drawGround();
	

	std::vector<RenderObject*>::iterator it;
	for(it = renderObj_.begin(); it!= renderObj_.end(); ++it)
	{
		(*it)->render(renderFrame_);
	}

	if(lastVertexClicked_!=-1) {
		Eigen::Vector3d originArrow = renderObj_[clickedVertexHandle_]->getMesh()->get_point_data(lastVertexClicked_,renderFrame_);
		//glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
		//glDisable(GL_LIGHTING);
		glBegin(GL_LINES);
		glColor3f(0.0,0.0,1.0);
		glVertex3d(arrowDest_[0],arrowDest_[1],arrowDest_[2]);
		glVertex3d(originArrow[0],originArrow[1],originArrow[2]);
		glEnd();
	}
	
}

void Scene::addRenderObject(RenderObject* obj)
{
	renderObj_.push_back(obj);
}

void Scene::setupScene()
{
	//Set up lighting and perspective data
	qglClearColor(Qt::black);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_CULL_FACE);
		glShadeModel(GL_SMOOTH);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		static GLfloat lightPosition[4] = { lightPos1_.x(),lightPos1_.y(),lightPos1_.z(),1.0};
		glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

	//Setup Ground Textures
	int imgHt,imgWd;
	textureData_ = readAsciiPPMImage("textures/woodfloortex.ppm",&imgHt,&imgWd);
	glGenTextures(1, &groundTex_);
	glBindTexture(GL_TEXTURE_2D, groundTex_);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,TEX_COLS,TEX_ROWS,0,GL_RGB,GL_UNSIGNED_BYTE,textureData_);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}

void Scene::resizeScene(int width, int height)
{
	glViewport(0,0,width,height);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
	gluPerspective(60, 1, .1, 1000);
		glMatrixMode(GL_MODELVIEW);
}

//@TODO: Custom Axis Labellign Colors
void Scene::drawAxes(float length)
{
		glPushMatrix();
	glColor3f(1.0,0.0,0.0);
		Arrow(0,0,0, length,0,0, 0.2);
		glPopMatrix();

		glPushMatrix();
	glColor3f(0.0,1.0,0.0);
		Arrow(0,0,0, 0,length,0, 0.2);
		glPopMatrix();

		glPushMatrix();
	glColor3f(0.0,0.0,1.0);
		Arrow(0,0,0, 0,0,length, 0.2);
		glPopMatrix();
	
}

//@TODO: Custom Gradient Background
void Scene::drawGradientBackGround() 
{
	// Meshlab background colors (given by Xiaofeng )
	float lower_color[3] = {115 / 255.0f, 115 / 255.0f, 230 / 255.0f};
	float upper_color[3] = {0 / 255.0f, 0 / 255.0f, 0 / 255.0f};

	GLboolean lighting_enabled, depth_test_enabled;
	glGetBooleanv(GL_DEPTH_TEST, &depth_test_enabled);
	glGetBooleanv(GL_LIGHTING, &lighting_enabled);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);

	glBegin(GL_QUADS);
	glColor3fv(lower_color);
	glVertex3f(-1.0,-1.0, 1);
	glVertex3f(1.0,-1.0, 1);
	glColor3fv(upper_color);
	glVertex3f(1.0, 1.0, 1);
	glVertex3f(-1.0, 1.0, 1);
	glEnd();

	if (lighting_enabled)
	{ 
		glEnable(GL_LIGHTING);
	}
	if (depth_test_enabled)
	{ 
		glEnable(GL_DEPTH_TEST);
	}

	glPopMatrix(); // Pop model view matrix
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}

//@TODO : Custom Texture Function Selection
void Scene::drawGround()
{
	glPushMatrix();

	glEnable(GL_TEXTURE_2D);

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glBindTexture(GL_TEXTURE_2D,groundTex_);

	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	glBegin(GL_QUADS);

	glColor3f(0.5,0.5,0.5);
	glTexCoord2f(0, 0);
	glVertex3f(50.0,0.0,-50.0);
	glTexCoord2f(1, 0);
	glVertex3f(-50.0,0.0,-50.0);
	glTexCoord2f(1, 1);
	glVertex3f(-50.0,0.0,50.0);
	glTexCoord2f(0,1);
	glVertex3f(50.0,0.0,50.0);

	glEnd();

	glDisable(GL_TEXTURE_2D);

	glPopMatrix();
}

//Private Functions for internal usage only

void Scene::Arrow(float x1,float y1,float z1,float x2,float y2,float z2,float D)
{
	float x=x2-x1;
	float y=y2-y1;
	float z=z2-z1;
	float L=sqrt(x*x+y*y+z*z);

	GLUquadricObj *quadObj;

	glPushMatrix ();

	glTranslated(x1,y1,z1);

	if((x!=0.)||(y!=0.))
	{
		glRotated(atan2(y,x)/RADPERDEG,0.,0.,1.);
		glRotated(atan2(sqrt(x*x+y*y),z)/RADPERDEG,0.,1.,0.);
	} 
	else if (z<0)
	{
		glRotated(180,1.,0.,0.);
	}
	
	glTranslatef(0,0,L-4*D);

	//Sketches the end cone
	quadObj = gluNewQuadric ();
	gluQuadricDrawStyle (quadObj, GLU_FILL);
	gluQuadricNormals (quadObj, GLU_SMOOTH);
	gluCylinder(quadObj, 0.65*D, 0.0, 4*D, 32, 1);
	gluDeleteQuadric(quadObj);
	
	glTranslatef(0,0,-L+4*D);

	//Sketches the rod
	quadObj = gluNewQuadric ();
	gluQuadricDrawStyle (quadObj, GLU_FILL);
	gluQuadricNormals (quadObj, GLU_SMOOTH);
	gluCylinder(quadObj, 0.25*D, 0.25*D, L-4*D, 32, 1);
	gluDeleteQuadric(quadObj);
		
	//Sketches the end sealing disc
	quadObj = gluNewQuadric ();
	gluQuadricDrawStyle (quadObj, GLU_FILL);
	gluQuadricNormals (quadObj, GLU_SMOOTH);
	gluDisk(quadObj, 0.0, 0.25*D, 32, 1);
	gluDeleteQuadric(quadObj);
	
		glPopMatrix ();
}

GLubyte* Scene::readAsciiPPMImage(char* fName,int* height,int* width) 
{
		
	FILE* in = fopen(fName, "r"); 

	int tht,twt;

			int  ccv; 
			char header[100]; 
			fscanf(in, "%s %d %d %d", header, &twt, &tht, &ccv); 
			int r, g, b; 

			*height = tht;
			*width = twt;

			GLubyte* texImg = new GLubyte[tht*twt*3];
		GLubyte readImage[500][500][3]; 

			for (int i=tht-1; i>=0; i--)
			{
				for (int j=0; j<twt; j++)
			{
						fscanf(in, "%d %d %d", &r, &g, &b); 
						readImage[i][j][0] = (GLubyte)r; 
						readImage[i][j][1] = (GLubyte)g; 
						readImage[i][j][2] = (GLubyte)b; 
				}
			}
				
			for (int i=0; i<tht; i++)
			{
				for ( int j=0; j<twt; j++)
				{
					texImg[i*twt*3+j*3+0] = readImage[i][j][0];
				texImg[i*twt*3+j*3+1] = readImage[i][j][1];
				texImg[i*twt*3+j*3+2] = readImage[i][j][2];
				}
			}	  
			fclose(in);
			return texImg; 
}

int Scene::performSelectionRayTest(double* rayStart, double* rayEnd, double *clickedWorldPos, double* selectedPos) {
	//Interate over the renderable object 
	std::vector<RenderObject*>::iterator it;
	int counter = 0;
	for(it = renderObj_.begin(); it!= renderObj_.end(); ++it,counter++)
	{
		bool isVisible = (*it)->isRenderable();
		SceneObjectId scId = (*it)->getSceneObjectId();
		if(isVisible && scId== CLOTH) {
			int vertexSelected = (*it)->performObjectSelectionRayTest(rayStart, rayEnd, clickedWorldPos, selectedPos,renderFrame_);
			if(vertexSelected != -1) {
				lastVertexClicked_ = vertexSelected;
				Eigen::Vector3d clickedWorldPosEig(clickedWorldPos[0],clickedWorldPos[1],clickedWorldPos[2]);
				arrowDest_ = clickedWorldPosEig;
				(*it)->setLastClickedVertex(lastVertexClicked_);
				(*it)->setLastClickedPosition(clickedWorldPosEig);
				clickedVertexHandle_ = counter;
				return vertexSelected;
			}
		}
	}
	return -1;
}

void Scene::updateDraggingForce(double x, double y, double z) {
	Eigen::Vector3d pos(x,y,z);
	arrowDest_ = pos;
	renderObj_[clickedVertexHandle_]->updateUIForce(pos,lastVertexClicked_,renderFrame_);
}

void Scene::resetDraggingForce() {
	renderObj_[clickedVertexHandle_]->resetUIForce();
	lastVertexClicked_ = -1;
	clickedVertexHandle_ = -1;
}

void Scene::saveObjs(int code) {
	//These codes are internally implemented 
	// 0 - all 1 - cloth 2 - body
	if(!dirsCreated_) {
		boost::filesystem::create_directories(clothPrefix_.c_str());
		boost::filesystem::create_directories(bodyPrefix_.c_str());
		dirsCreated_ = true;
	}

	std::vector<RenderObject*>::iterator it;
	for(it = renderObj_.begin(); it!= renderObj_.end(); ++it) {
		SceneObjectId scId = (*it)->getSceneObjectId();
		if(code==0) {
			if(scId==COLL_OBJ)
				(*it)->getMesh()->saveAsOBJ(bodyPrefix_.c_str());
			else
				(*it)->getMesh()->saveAsOBJ(clothPrefix_.c_str());
		} else if(code==1 && scId==CLOTH) {
			(*it)->getMesh()->saveAsOBJ(clothPrefix_.c_str());
		} else if(code==2 && scId==COLL_OBJ) {
			(*it)->getMesh()->saveAsOBJ(bodyPrefix_.c_str());
		}
	}
}

void Scene::updateOBJPrefixes(std::string &prefix)
{
	boost::trim(prefix);
	const char* suffix = prefix.c_str();
	//Get the current date and time also and create a folder 
	std::chrono::time_point<std::chrono::system_clock> timeNow;
	timeNow = std::chrono::system_clock::now();
	std::time_t timeNowT = std::chrono::system_clock::to_time_t(timeNow);
	const char* timeStr = std::ctime(&timeNowT);
	std::string timeS(timeStr);
	boost::replace_all(timeS," ","_");
	boost::replace_all(timeS,":","_");
	boost::trim(timeS);
	std::stringstream ssc;
	ssc << suffix << "\\" << timeS << "\\Cloth";
	clothPrefix_ = ssc.str();
	std::stringstream ssb;
	ssb << suffix << "\\" << timeS << "\\Body";
	bodyPrefix_ = ssb.str();
	dirsCreated_ = false;

	std::cout << clothPrefix_ << "\n";
	std::cout << bodyPrefix_ << "\n";

}
