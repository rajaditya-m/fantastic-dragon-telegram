#include "cloth_vorticity.h"

Cloth_Vorticity::Cloth_Vorticity(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	ui_custom.setUpUi(this,ui);

}

Cloth_Vorticity::~Cloth_Vorticity()
{

}
