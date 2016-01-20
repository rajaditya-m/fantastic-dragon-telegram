#include <QGLWidget>
#include <QtOpenGL>
#include <gl/glu.h>

inline double getPixelDepth(int pixel_position_x, int pixel_position_y) {
	GLint view_port[4];
	glGetIntegerv(GL_VIEWPORT, view_port);
	float depth;
	glReadPixels(pixel_position_x, view_port[3] - pixel_position_y, 1, 1,
							 GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
	return (double)depth;
}

inline void getPixelWorldPosition(int pixel_position_x, int pixel_position_y, double* world_pos) {

	double depth = getPixelDepth(pixel_position_x,pixel_position_y);
	GLint view_port[4]; // viewport dimensions+pos
	GLdouble projectioin_matrix[16];
	GLdouble model_view_matrix[16];
	glGetDoublev (GL_MODELVIEW_MATRIX, model_view_matrix);
	glGetDoublev (GL_PROJECTION_MATRIX, projectioin_matrix);
	glGetIntegerv(GL_VIEWPORT, view_port);
	//view[3]-cursorY = conversion from upper left (0,0) to lower left (0,0)
	//Unproject 2D Screen coordinates into wonderful world coordinates
	gluUnProject((GLdouble) pixel_position_x, (GLdouble) view_port[3] - pixel_position_y, depth,
							 model_view_matrix, projectioin_matrix, view_port,
							 world_pos + 0, world_pos + 1, world_pos + 2);
}

inline void getSelectionRay(int pixel_position_x, int pixel_position_y, double *starting_point, double *ending_point) {
	double objX, objY, objZ; // holder for world coordinates
	GLint view_port[4]; // viewport dimensions+pos
	GLdouble projectioin_matrix[16];
	GLdouble model_view_matrix[16];
	glGetDoublev (GL_MODELVIEW_MATRIX, model_view_matrix);
	glGetDoublev (GL_PROJECTION_MATRIX, projectioin_matrix);
	glGetIntegerv(GL_VIEWPORT, view_port);
	//view[3]-cursorY = conversion from upper left (0,0) to lower left (0,0)
	//Unproject 2D Screen coordinates into wonderful world coordinates
	gluUnProject((GLdouble) pixel_position_x, (GLdouble) view_port[3] - pixel_position_y, 1,
							 model_view_matrix, projectioin_matrix, view_port, &objX, &objY, &objZ);
	ending_point[0] = objX;
	ending_point[1] = objY;
	ending_point[2] = objZ;
	gluUnProject((GLdouble) pixel_position_x, (GLdouble) view_port[3] - pixel_position_y, 0,
							 model_view_matrix, projectioin_matrix, view_port, &objX, &objY, &objZ);
	starting_point[0] = objX;
	starting_point[1] = objY;
	starting_point[2] = objZ;
}
