#pragma once

#include "ui_cloth_vorticity.h"
#include <QtWidgets/QActionGroup>

namespace Ui {

class ui_custom_mods : public Ui_Cloth_VorticityClass
{
public:
	QActionGroup* body_render_action_group;
	QActionGroup* cloth_render_action_group;
	void setUpUi(QMainWindow *CVC,Ui::Cloth_VorticityClass parent)
	{
		//For showing the body
		parent.actionShow_Body->setCheckable(true);
		parent.actionShow_Body->setChecked(true);
		QObject::connect(parent.actionShow_Body, SIGNAL(toggled(bool)), parent.widget, SLOT(bodyDisplayToggled(bool)));

		//For showing the cloth 
		parent.actionShow_Clothing->setCheckable(true);
		parent.actionShow_Clothing->setChecked(true); 
		QObject::connect(parent.actionShow_Clothing, SIGNAL(toggled(bool)), parent.widget, SLOT(clothDisplayToggled(bool)));

		//For showing the grid 
		parent.actionShow_Co_ord_Axis->setCheckable(true);
		QObject::connect(parent.actionShow_Co_ord_Axis, SIGNAL(toggled(bool)), parent.widget, SLOT(renderAxesToggled(bool)));

		//For showing the ground 
		parent.actionShow_Ground->setCheckable(true);
		parent.actionShow_Ground->setChecked(true);
		QObject::connect(parent.actionShow_Ground, SIGNAL(toggled(bool)), parent.widget, SLOT(renderGroundToggled(bool)));

		//This is the render group for the cloth shader
		cloth_render_action_group = new QActionGroup(CVC);
		parent.actionShaded_Cloth->setCheckable(true);
		parent.actionWireframe_Cloth->setCheckable(true);
		parent.actionHeatmap_Velocity->setCheckable(true);
		parent.actionHeatmap_Acceleration->setCheckable(true);
		parent.actionHeatmap_Bending_Force->setCheckable(true);
		parent.actionHeatmap_Damping_Force->setCheckable(true);
		parent.actionHeatmap_Shear_Force->setCheckable(true);
		cloth_render_action_group->addAction(parent.actionShaded_Cloth);
		cloth_render_action_group->addAction(parent.actionWireframe_Cloth);
		cloth_render_action_group->addAction(parent.actionHeatmap_Velocity);
		cloth_render_action_group->addAction(parent.actionHeatmap_Acceleration);
		cloth_render_action_group->addAction(parent.actionHeatmap_Bending_Force);
		cloth_render_action_group->addAction(parent.actionHeatmap_Damping_Force);
		cloth_render_action_group->addAction(parent.actionHeatmap_Shear_Force);
		cloth_render_action_group->setExclusive(true);

		parent.actionShaded_Cloth->setChecked(true);	
		QObject::connect(parent.actionShaded_Cloth, SIGNAL(toggled(bool)), parent.widget, SLOT(setClothRenderInShading(bool)));
		QObject::connect(parent.actionWireframe_Cloth, SIGNAL(toggled(bool)), parent.widget, SLOT(setClothRenderInWireframe(bool)));
		QObject::connect(parent.actionHeatmap_Velocity, SIGNAL(toggled(bool)), parent.widget, SLOT(setClothRenderInHeatmapVelocity(bool)));
		QObject::connect(parent.actionHeatmap_Acceleration, SIGNAL(toggled(bool)), parent.widget, SLOT(setClothRenderInHeatmapAcceleration(bool)));
		QObject::connect(parent.actionHeatmap_Bending_Force, SIGNAL(toggled(bool)), parent.widget, SLOT(setClothRenderInHeatmapBendingForce(bool)));
		QObject::connect(parent.actionHeatmap_Damping_Force, SIGNAL(toggled(bool)), parent.widget, SLOT(setClothRenderInHeatmapDampingForce(bool)));
		QObject::connect(parent.actionHeatmap_Shear_Force, SIGNAL(toggled(bool)), parent.widget, SLOT(setClothRenderInHeatmapShearForce(bool)));
		
		//This is the group for body shader
		body_render_action_group = new QActionGroup(CVC);
		parent.actionShaded_Body->setCheckable(true);
		parent.actionWireframe_Body->setCheckable(true);
		body_render_action_group->addAction(parent.actionShaded_Body);
		body_render_action_group->addAction(parent.actionWireframe_Body);
		body_render_action_group->setExclusive(true);
		parent.actionShaded_Body->setChecked(true);
		QObject::connect(parent.actionShaded_Body, SIGNAL(toggled(bool)), parent.widget, SLOT(setBodyRenderInShading(bool)));
		QObject::connect(parent.actionWireframe_Body, SIGNAL(toggled(bool)), parent.widget, SLOT(setBodyRenderInWireframe(bool)));

		//This is for the play option 
		QObject::connect(parent.actionPlay, SIGNAL(triggered()), parent.widget, SLOT(startAnimation()));

		//This is for the pauase animation 
		QObject::connect(parent.actionPause, SIGNAL(triggered()), parent.widget, SLOT(pauseAnimation()));
	
		//This is for reset
		QObject::connect(parent.actionReset, SIGNAL(triggered()), parent.widget, SLOT(resetAnimation()));

		//This is for the save option (cloth)
		QObject::connect(parent.actionSave_Cloth, SIGNAL(triggered()), parent.widget, SLOT(setSaveAsOBJCloth()));

		//This is for the save option (body)
		QObject::connect(parent.actionBody, SIGNAL(triggered()), parent.widget, SLOT(setSaveAsOBJBody()));

		//This is for the save option (all)
		QObject::connect(parent.actionAll, SIGNAL(triggered()), parent.widget, SLOT(setSaveAsOBJAll()));

	}
};

}

