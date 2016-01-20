#include "cloth_vorticity.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	Cloth_Vorticity w;
	w.show();
	return a.exec();
}
