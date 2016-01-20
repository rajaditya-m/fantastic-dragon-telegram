#ifndef CLOTH_VORTICITY_H
#define CLOTH_VORTICITY_H

#include <QtWidgets/QMainWindow>
#include "ui_cloth_vorticity.h"
#include "ui_custom_mods.h"
#include "glwidget.h"

class Cloth_Vorticity : public QMainWindow
{
	Q_OBJECT

public:
	Cloth_Vorticity(QWidget *parent = 0);
	~Cloth_Vorticity();

private:
	Ui::Cloth_VorticityClass ui;
	Ui::ui_custom_mods ui_custom;

};

#endif // CLOTH_VORTICITY_H
