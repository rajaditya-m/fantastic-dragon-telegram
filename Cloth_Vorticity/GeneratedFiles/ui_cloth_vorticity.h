/********************************************************************************
** Form generated from reading UI file 'cloth_vorticity.ui'
**
** Created by: Qt User Interface Compiler version 5.4.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CLOTH_VORTICITY_H
#define UI_CLOTH_VORTICITY_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QSlider>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>
#include "glwidget.h"

QT_BEGIN_NAMESPACE

class Ui_Cloth_VorticityClass
{
public:
    QAction *actionShow_Body;
    QAction *actionShow_Clothing;
    QAction *actionShaded_Cloth;
    QAction *actionWireframe_Cloth;
    QAction *actionShaded_Body;
    QAction *actionWireframe_Body;
    QAction *actionShow_Co_ord_Axis;
    QAction *actionShow_Ground;
    QAction *actionVelocity;
    QAction *actionTotal_Force;
    QAction *actionStretch_Force;
    QAction *actionDamping_Force;
    QAction *actionBending_Force;
    QAction *actionGravity_Forrce;
    QAction *actionAcceleration;
    QAction *actionPlay;
    QAction *actionPause;
    QAction *actionStop;
    QAction *actionHeatmap_Velocity;
    QAction *actionHeatmap_Acceleration;
    QAction *actionHeatmap_Shear_Force;
    QAction *actionHeatmap_Bending_Force;
    QAction *actionHeatmap_Damping_Force;
    QAction *actionReset;
    QAction *actionSave_Cloth;
    QAction *actionBody;
    QAction *actionAll;
    QWidget *centralWidget;
    GLWidget *widget;
    QSlider *horizontalSlider;
    QDoubleSpinBox *doubleSpinBox;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuView;
    QMenu *menuRender;
    QMenu *menuHeatMap_Cloth;
    QMenu *menuAnimation;
    QMenu *menuActions;
    QMenu *menuSave_As_OBJ;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *Cloth_VorticityClass)
    {
        if (Cloth_VorticityClass->objectName().isEmpty())
            Cloth_VorticityClass->setObjectName(QStringLiteral("Cloth_VorticityClass"));
        Cloth_VorticityClass->resize(1189, 778);
        actionShow_Body = new QAction(Cloth_VorticityClass);
        actionShow_Body->setObjectName(QStringLiteral("actionShow_Body"));
        actionShow_Body->setCheckable(false);
        actionShow_Clothing = new QAction(Cloth_VorticityClass);
        actionShow_Clothing->setObjectName(QStringLiteral("actionShow_Clothing"));
        actionShow_Clothing->setCheckable(false);
        actionShaded_Cloth = new QAction(Cloth_VorticityClass);
        actionShaded_Cloth->setObjectName(QStringLiteral("actionShaded_Cloth"));
        actionShaded_Cloth->setCheckable(false);
        actionWireframe_Cloth = new QAction(Cloth_VorticityClass);
        actionWireframe_Cloth->setObjectName(QStringLiteral("actionWireframe_Cloth"));
        actionWireframe_Cloth->setCheckable(false);
        actionShaded_Body = new QAction(Cloth_VorticityClass);
        actionShaded_Body->setObjectName(QStringLiteral("actionShaded_Body"));
        actionShaded_Body->setCheckable(false);
        actionWireframe_Body = new QAction(Cloth_VorticityClass);
        actionWireframe_Body->setObjectName(QStringLiteral("actionWireframe_Body"));
        actionWireframe_Body->setCheckable(false);
        actionShow_Co_ord_Axis = new QAction(Cloth_VorticityClass);
        actionShow_Co_ord_Axis->setObjectName(QStringLiteral("actionShow_Co_ord_Axis"));
        actionShow_Ground = new QAction(Cloth_VorticityClass);
        actionShow_Ground->setObjectName(QStringLiteral("actionShow_Ground"));
        actionVelocity = new QAction(Cloth_VorticityClass);
        actionVelocity->setObjectName(QStringLiteral("actionVelocity"));
        actionTotal_Force = new QAction(Cloth_VorticityClass);
        actionTotal_Force->setObjectName(QStringLiteral("actionTotal_Force"));
        actionStretch_Force = new QAction(Cloth_VorticityClass);
        actionStretch_Force->setObjectName(QStringLiteral("actionStretch_Force"));
        actionDamping_Force = new QAction(Cloth_VorticityClass);
        actionDamping_Force->setObjectName(QStringLiteral("actionDamping_Force"));
        actionBending_Force = new QAction(Cloth_VorticityClass);
        actionBending_Force->setObjectName(QStringLiteral("actionBending_Force"));
        actionGravity_Forrce = new QAction(Cloth_VorticityClass);
        actionGravity_Forrce->setObjectName(QStringLiteral("actionGravity_Forrce"));
        actionAcceleration = new QAction(Cloth_VorticityClass);
        actionAcceleration->setObjectName(QStringLiteral("actionAcceleration"));
        actionPlay = new QAction(Cloth_VorticityClass);
        actionPlay->setObjectName(QStringLiteral("actionPlay"));
        actionPause = new QAction(Cloth_VorticityClass);
        actionPause->setObjectName(QStringLiteral("actionPause"));
        actionStop = new QAction(Cloth_VorticityClass);
        actionStop->setObjectName(QStringLiteral("actionStop"));
        actionHeatmap_Velocity = new QAction(Cloth_VorticityClass);
        actionHeatmap_Velocity->setObjectName(QStringLiteral("actionHeatmap_Velocity"));
        actionHeatmap_Acceleration = new QAction(Cloth_VorticityClass);
        actionHeatmap_Acceleration->setObjectName(QStringLiteral("actionHeatmap_Acceleration"));
        actionHeatmap_Shear_Force = new QAction(Cloth_VorticityClass);
        actionHeatmap_Shear_Force->setObjectName(QStringLiteral("actionHeatmap_Shear_Force"));
        actionHeatmap_Bending_Force = new QAction(Cloth_VorticityClass);
        actionHeatmap_Bending_Force->setObjectName(QStringLiteral("actionHeatmap_Bending_Force"));
        actionHeatmap_Damping_Force = new QAction(Cloth_VorticityClass);
        actionHeatmap_Damping_Force->setObjectName(QStringLiteral("actionHeatmap_Damping_Force"));
        actionReset = new QAction(Cloth_VorticityClass);
        actionReset->setObjectName(QStringLiteral("actionReset"));
        actionSave_Cloth = new QAction(Cloth_VorticityClass);
        actionSave_Cloth->setObjectName(QStringLiteral("actionSave_Cloth"));
        actionBody = new QAction(Cloth_VorticityClass);
        actionBody->setObjectName(QStringLiteral("actionBody"));
        actionAll = new QAction(Cloth_VorticityClass);
        actionAll->setObjectName(QStringLiteral("actionAll"));
        centralWidget = new QWidget(Cloth_VorticityClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        widget = new GLWidget(centralWidget);
        widget->setObjectName(QStringLiteral("widget"));
        widget->setGeometry(QRect(9, 9, 1171, 651));
        horizontalSlider = new QSlider(centralWidget);
        horizontalSlider->setObjectName(QStringLiteral("horizontalSlider"));
        horizontalSlider->setGeometry(QRect(9, 680, 1101, 22));
        horizontalSlider->setMaximum(3000);
        horizontalSlider->setOrientation(Qt::Horizontal);
        horizontalSlider->setTickPosition(QSlider::TicksBelow);
        horizontalSlider->setTickInterval(5);
        doubleSpinBox = new QDoubleSpinBox(centralWidget);
        doubleSpinBox->setObjectName(QStringLiteral("doubleSpinBox"));
        doubleSpinBox->setGeometry(QRect(1120, 680, 62, 21));
        doubleSpinBox->setWrapping(true);
        doubleSpinBox->setDecimals(0);
        doubleSpinBox->setMaximum(3000);
        Cloth_VorticityClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(Cloth_VorticityClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1189, 21));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QStringLiteral("menuFile"));
        menuView = new QMenu(menuBar);
        menuView->setObjectName(QStringLiteral("menuView"));
        menuRender = new QMenu(menuBar);
        menuRender->setObjectName(QStringLiteral("menuRender"));
        menuHeatMap_Cloth = new QMenu(menuRender);
        menuHeatMap_Cloth->setObjectName(QStringLiteral("menuHeatMap_Cloth"));
        menuAnimation = new QMenu(menuBar);
        menuAnimation->setObjectName(QStringLiteral("menuAnimation"));
        menuActions = new QMenu(menuBar);
        menuActions->setObjectName(QStringLiteral("menuActions"));
        menuSave_As_OBJ = new QMenu(menuActions);
        menuSave_As_OBJ->setObjectName(QStringLiteral("menuSave_As_OBJ"));
        Cloth_VorticityClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(Cloth_VorticityClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        Cloth_VorticityClass->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(Cloth_VorticityClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        Cloth_VorticityClass->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuView->menuAction());
        menuBar->addAction(menuRender->menuAction());
        menuBar->addAction(menuAnimation->menuAction());
        menuBar->addAction(menuActions->menuAction());
        menuRender->addAction(actionShow_Body);
        menuRender->addAction(actionShow_Clothing);
        menuRender->addSeparator();
        menuRender->addAction(actionShow_Co_ord_Axis);
        menuRender->addAction(actionShow_Ground);
        menuRender->addSeparator();
        menuRender->addAction(menuHeatMap_Cloth->menuAction());
        menuRender->addSeparator();
        menuRender->addAction(actionShaded_Body);
        menuRender->addAction(actionWireframe_Body);
        menuHeatMap_Cloth->addAction(actionShaded_Cloth);
        menuHeatMap_Cloth->addAction(actionWireframe_Cloth);
        menuHeatMap_Cloth->addSeparator();
        menuHeatMap_Cloth->addAction(actionHeatmap_Velocity);
        menuHeatMap_Cloth->addAction(actionHeatmap_Acceleration);
        menuHeatMap_Cloth->addSeparator();
        menuHeatMap_Cloth->addAction(actionHeatmap_Shear_Force);
        menuHeatMap_Cloth->addAction(actionHeatmap_Bending_Force);
        menuHeatMap_Cloth->addAction(actionHeatmap_Damping_Force);
        menuAnimation->addAction(actionPlay);
        menuAnimation->addAction(actionPause);
        menuAnimation->addAction(actionReset);
        menuActions->addAction(menuSave_As_OBJ->menuAction());
        menuSave_As_OBJ->addAction(actionSave_Cloth);
        menuSave_As_OBJ->addAction(actionBody);
        menuSave_As_OBJ->addAction(actionAll);

        retranslateUi(Cloth_VorticityClass);
        QObject::connect(doubleSpinBox, SIGNAL(valueChanged(double)), widget, SLOT(setRenderingFrameNumber(double)));
        QObject::connect(horizontalSlider, SIGNAL(valueChanged(int)), widget, SLOT(setRenderingFrameNumber(int)));
        QObject::connect(widget, SIGNAL(frameNumberUpdated(double)), doubleSpinBox, SLOT(setValue(double)));
        QObject::connect(widget, SIGNAL(frameNumberUpdated(int)), horizontalSlider, SLOT(setValue(int)));
        QObject::connect(widget, SIGNAL(animFrameNumberUpdated(int)), horizontalSlider, SLOT(setValue(int)));

        QMetaObject::connectSlotsByName(Cloth_VorticityClass);
    } // setupUi

    void retranslateUi(QMainWindow *Cloth_VorticityClass)
    {
        Cloth_VorticityClass->setWindowTitle(QApplication::translate("Cloth_VorticityClass", "Cloth_Vorticity", 0));
        actionShow_Body->setText(QApplication::translate("Cloth_VorticityClass", "Show Body", 0));
        actionShow_Clothing->setText(QApplication::translate("Cloth_VorticityClass", "Show Clothing", 0));
        actionShaded_Cloth->setText(QApplication::translate("Cloth_VorticityClass", "Shaded Cloth", 0));
        actionWireframe_Cloth->setText(QApplication::translate("Cloth_VorticityClass", "Wireframe Cloth", 0));
        actionShaded_Body->setText(QApplication::translate("Cloth_VorticityClass", "Shaded Body", 0));
        actionWireframe_Body->setText(QApplication::translate("Cloth_VorticityClass", "Wireframe Body", 0));
        actionShow_Co_ord_Axis->setText(QApplication::translate("Cloth_VorticityClass", "Show Co-ord Axis", 0));
        actionShow_Ground->setText(QApplication::translate("Cloth_VorticityClass", "Show Ground", 0));
        actionVelocity->setText(QApplication::translate("Cloth_VorticityClass", "Velocity", 0));
        actionTotal_Force->setText(QApplication::translate("Cloth_VorticityClass", "Total Force", 0));
        actionStretch_Force->setText(QApplication::translate("Cloth_VorticityClass", "Stretch Force", 0));
        actionDamping_Force->setText(QApplication::translate("Cloth_VorticityClass", "Damping Force", 0));
        actionBending_Force->setText(QApplication::translate("Cloth_VorticityClass", "Bending Force", 0));
        actionGravity_Forrce->setText(QApplication::translate("Cloth_VorticityClass", "Gravity Force", 0));
        actionAcceleration->setText(QApplication::translate("Cloth_VorticityClass", "Acceleration", 0));
        actionPlay->setText(QApplication::translate("Cloth_VorticityClass", "Play", 0));
        actionPause->setText(QApplication::translate("Cloth_VorticityClass", "Pause", 0));
        actionStop->setText(QApplication::translate("Cloth_VorticityClass", "Stop", 0));
        actionHeatmap_Velocity->setText(QApplication::translate("Cloth_VorticityClass", "Heatmap - Velocity", 0));
        actionHeatmap_Acceleration->setText(QApplication::translate("Cloth_VorticityClass", "Heatmap - Acceleration", 0));
        actionHeatmap_Shear_Force->setText(QApplication::translate("Cloth_VorticityClass", "Heatmap - Shear Force", 0));
        actionHeatmap_Bending_Force->setText(QApplication::translate("Cloth_VorticityClass", "Heatmap - Bending Force", 0));
        actionHeatmap_Damping_Force->setText(QApplication::translate("Cloth_VorticityClass", "Heatmap - Damping Force", 0));
        actionReset->setText(QApplication::translate("Cloth_VorticityClass", "Reset ", 0));
        actionSave_Cloth->setText(QApplication::translate("Cloth_VorticityClass", "Cloth", 0));
        actionBody->setText(QApplication::translate("Cloth_VorticityClass", "Body", 0));
        actionAll->setText(QApplication::translate("Cloth_VorticityClass", "All", 0));
        menuFile->setTitle(QApplication::translate("Cloth_VorticityClass", "File", 0));
        menuView->setTitle(QApplication::translate("Cloth_VorticityClass", "View", 0));
        menuRender->setTitle(QApplication::translate("Cloth_VorticityClass", "Render", 0));
        menuHeatMap_Cloth->setTitle(QApplication::translate("Cloth_VorticityClass", "Cloth Rendering", 0));
        menuAnimation->setTitle(QApplication::translate("Cloth_VorticityClass", "Animation", 0));
        menuActions->setTitle(QApplication::translate("Cloth_VorticityClass", "Actions", 0));
        menuSave_As_OBJ->setTitle(QApplication::translate("Cloth_VorticityClass", "Save As OBJ", 0));
    } // retranslateUi

};

namespace Ui {
    class Cloth_VorticityClass: public Ui_Cloth_VorticityClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CLOTH_VORTICITY_H
