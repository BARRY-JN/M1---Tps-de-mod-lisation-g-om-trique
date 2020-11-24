/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.13.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include <meshviewerwidget.h>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QHBoxLayout *horizontalLayout;
    QWidget *widget_2;
    QVBoxLayout *verticalLayout;
    QPushButton *pushButton_chargement;
    QFrame *frame;
    QVBoxLayout *verticalLayout_2;
    QLabel *label;
    QWidget *widget_4;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_4;
    QPushButton *pushButton_vertexMoins;
    QLabel *labelVertex;
    QPushButton *pushButton_vertexPlus;
    QWidget *widget_3;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_5;
    QPushButton *pushButton_edgeMoins;
    QLabel *labelEdge;
    QPushButton *pushButton_edgePlus;
    QWidget *widget;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_6;
    QPushButton *pushButton_faceMoins;
    QLabel *labelFace;
    QPushButton *pushButton_facePlus;
    QPushButton *pushButton_voisinage;
    QFrame *frame_3;
    QHBoxLayout *horizontalLayout_6;
    QPushButton *pushButton_bordure;
    QFrame *frame_2;
    QVBoxLayout *verticalLayout_3;
    QWidget *widget_5;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_v1_chemin;
    QSpinBox *spinBox_v1_chemin;
    QLabel *label_v2_chemin;
    QSpinBox *spinBox_v2_chemin;
    QPushButton *pushButton_afficherChemin;
    QSpacerItem *verticalSpacer;
    MeshViewerWidget *displayWidget;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(632, 408);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        horizontalLayout = new QHBoxLayout(centralWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        widget_2 = new QWidget(centralWidget);
        widget_2->setObjectName(QString::fromUtf8("widget_2"));
        QSizePolicy sizePolicy(QSizePolicy::Maximum, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(150);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(widget_2->sizePolicy().hasHeightForWidth());
        widget_2->setSizePolicy(sizePolicy);
        widget_2->setMinimumSize(QSize(150, 0));
        verticalLayout = new QVBoxLayout(widget_2);
        verticalLayout->setSpacing(4);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(3, 3, 3, 3);
        pushButton_chargement = new QPushButton(widget_2);
        pushButton_chargement->setObjectName(QString::fromUtf8("pushButton_chargement"));
        pushButton_chargement->setMinimumSize(QSize(200, 0));

        verticalLayout->addWidget(pushButton_chargement);

        frame = new QFrame(widget_2);
        frame->setObjectName(QString::fromUtf8("frame"));
        frame->setFrameShape(QFrame::StyledPanel);
        frame->setFrameShadow(QFrame::Raised);
        verticalLayout_2 = new QVBoxLayout(frame);
        verticalLayout_2->setSpacing(1);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(2, 2, 2, 2);
        label = new QLabel(frame);
        label->setObjectName(QString::fromUtf8("label"));

        verticalLayout_2->addWidget(label);

        widget_4 = new QWidget(frame);
        widget_4->setObjectName(QString::fromUtf8("widget_4"));
        horizontalLayout_2 = new QHBoxLayout(widget_4);
        horizontalLayout_2->setSpacing(1);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(2, 2, 2, 2);
        label_4 = new QLabel(widget_4);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        QSizePolicy sizePolicy1(QSizePolicy::Maximum, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(label_4->sizePolicy().hasHeightForWidth());
        label_4->setSizePolicy(sizePolicy1);

        horizontalLayout_2->addWidget(label_4);

        pushButton_vertexMoins = new QPushButton(widget_4);
        pushButton_vertexMoins->setObjectName(QString::fromUtf8("pushButton_vertexMoins"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(pushButton_vertexMoins->sizePolicy().hasHeightForWidth());
        pushButton_vertexMoins->setSizePolicy(sizePolicy2);

        horizontalLayout_2->addWidget(pushButton_vertexMoins);

        labelVertex = new QLabel(widget_4);
        labelVertex->setObjectName(QString::fromUtf8("labelVertex"));
        QSizePolicy sizePolicy3(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(labelVertex->sizePolicy().hasHeightForWidth());
        labelVertex->setSizePolicy(sizePolicy3);
        labelVertex->setAlignment(Qt::AlignCenter);

        horizontalLayout_2->addWidget(labelVertex);

        pushButton_vertexPlus = new QPushButton(widget_4);
        pushButton_vertexPlus->setObjectName(QString::fromUtf8("pushButton_vertexPlus"));
        sizePolicy2.setHeightForWidth(pushButton_vertexPlus->sizePolicy().hasHeightForWidth());
        pushButton_vertexPlus->setSizePolicy(sizePolicy2);

        horizontalLayout_2->addWidget(pushButton_vertexPlus);


        verticalLayout_2->addWidget(widget_4);

        widget_3 = new QWidget(frame);
        widget_3->setObjectName(QString::fromUtf8("widget_3"));
        sizePolicy3.setHeightForWidth(widget_3->sizePolicy().hasHeightForWidth());
        widget_3->setSizePolicy(sizePolicy3);
        horizontalLayout_3 = new QHBoxLayout(widget_3);
        horizontalLayout_3->setSpacing(1);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(2, 2, 2, 2);
        label_5 = new QLabel(widget_3);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        sizePolicy1.setHeightForWidth(label_5->sizePolicy().hasHeightForWidth());
        label_5->setSizePolicy(sizePolicy1);

        horizontalLayout_3->addWidget(label_5);

        pushButton_edgeMoins = new QPushButton(widget_3);
        pushButton_edgeMoins->setObjectName(QString::fromUtf8("pushButton_edgeMoins"));
        sizePolicy2.setHeightForWidth(pushButton_edgeMoins->sizePolicy().hasHeightForWidth());
        pushButton_edgeMoins->setSizePolicy(sizePolicy2);

        horizontalLayout_3->addWidget(pushButton_edgeMoins);

        labelEdge = new QLabel(widget_3);
        labelEdge->setObjectName(QString::fromUtf8("labelEdge"));
        sizePolicy3.setHeightForWidth(labelEdge->sizePolicy().hasHeightForWidth());
        labelEdge->setSizePolicy(sizePolicy3);
        labelEdge->setAlignment(Qt::AlignCenter);

        horizontalLayout_3->addWidget(labelEdge);

        pushButton_edgePlus = new QPushButton(widget_3);
        pushButton_edgePlus->setObjectName(QString::fromUtf8("pushButton_edgePlus"));
        sizePolicy2.setHeightForWidth(pushButton_edgePlus->sizePolicy().hasHeightForWidth());
        pushButton_edgePlus->setSizePolicy(sizePolicy2);
        pushButton_edgePlus->setMinimumSize(QSize(10, 0));

        horizontalLayout_3->addWidget(pushButton_edgePlus);


        verticalLayout_2->addWidget(widget_3);

        widget = new QWidget(frame);
        widget->setObjectName(QString::fromUtf8("widget"));
        horizontalLayout_4 = new QHBoxLayout(widget);
        horizontalLayout_4->setSpacing(1);
        horizontalLayout_4->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(2, 2, 2, 2);
        label_6 = new QLabel(widget);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        sizePolicy1.setHeightForWidth(label_6->sizePolicy().hasHeightForWidth());
        label_6->setSizePolicy(sizePolicy1);

        horizontalLayout_4->addWidget(label_6);

        pushButton_faceMoins = new QPushButton(widget);
        pushButton_faceMoins->setObjectName(QString::fromUtf8("pushButton_faceMoins"));
        sizePolicy2.setHeightForWidth(pushButton_faceMoins->sizePolicy().hasHeightForWidth());
        pushButton_faceMoins->setSizePolicy(sizePolicy2);

        horizontalLayout_4->addWidget(pushButton_faceMoins);

        labelFace = new QLabel(widget);
        labelFace->setObjectName(QString::fromUtf8("labelFace"));
        sizePolicy3.setHeightForWidth(labelFace->sizePolicy().hasHeightForWidth());
        labelFace->setSizePolicy(sizePolicy3);
        labelFace->setAlignment(Qt::AlignCenter);

        horizontalLayout_4->addWidget(labelFace);

        pushButton_facePlus = new QPushButton(widget);
        pushButton_facePlus->setObjectName(QString::fromUtf8("pushButton_facePlus"));
        sizePolicy2.setHeightForWidth(pushButton_facePlus->sizePolicy().hasHeightForWidth());
        pushButton_facePlus->setSizePolicy(sizePolicy2);

        horizontalLayout_4->addWidget(pushButton_facePlus);


        verticalLayout_2->addWidget(widget);

        pushButton_voisinage = new QPushButton(frame);
        pushButton_voisinage->setObjectName(QString::fromUtf8("pushButton_voisinage"));

        verticalLayout_2->addWidget(pushButton_voisinage);


        verticalLayout->addWidget(frame);

        frame_3 = new QFrame(widget_2);
        frame_3->setObjectName(QString::fromUtf8("frame_3"));
        frame_3->setFrameShape(QFrame::StyledPanel);
        frame_3->setFrameShadow(QFrame::Raised);
        horizontalLayout_6 = new QHBoxLayout(frame_3);
        horizontalLayout_6->setSpacing(1);
        horizontalLayout_6->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        horizontalLayout_6->setContentsMargins(1, 1, 1, 1);
        pushButton_bordure = new QPushButton(frame_3);
        pushButton_bordure->setObjectName(QString::fromUtf8("pushButton_bordure"));

        horizontalLayout_6->addWidget(pushButton_bordure);


        verticalLayout->addWidget(frame_3);

        frame_2 = new QFrame(widget_2);
        frame_2->setObjectName(QString::fromUtf8("frame_2"));
        frame_2->setFrameShape(QFrame::StyledPanel);
        frame_2->setFrameShadow(QFrame::Raised);
        verticalLayout_3 = new QVBoxLayout(frame_2);
        verticalLayout_3->setSpacing(2);
        verticalLayout_3->setContentsMargins(11, 11, 11, 11);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        verticalLayout_3->setContentsMargins(1, 2, 2, 2);
        widget_5 = new QWidget(frame_2);
        widget_5->setObjectName(QString::fromUtf8("widget_5"));
        horizontalLayout_5 = new QHBoxLayout(widget_5);
        horizontalLayout_5->setSpacing(1);
        horizontalLayout_5->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        horizontalLayout_5->setContentsMargins(1, 1, 1, 1);
        label_v1_chemin = new QLabel(widget_5);
        label_v1_chemin->setObjectName(QString::fromUtf8("label_v1_chemin"));
        sizePolicy1.setHeightForWidth(label_v1_chemin->sizePolicy().hasHeightForWidth());
        label_v1_chemin->setSizePolicy(sizePolicy1);

        horizontalLayout_5->addWidget(label_v1_chemin);

        spinBox_v1_chemin = new QSpinBox(widget_5);
        spinBox_v1_chemin->setObjectName(QString::fromUtf8("spinBox_v1_chemin"));

        horizontalLayout_5->addWidget(spinBox_v1_chemin);

        label_v2_chemin = new QLabel(widget_5);
        label_v2_chemin->setObjectName(QString::fromUtf8("label_v2_chemin"));
        sizePolicy1.setHeightForWidth(label_v2_chemin->sizePolicy().hasHeightForWidth());
        label_v2_chemin->setSizePolicy(sizePolicy1);

        horizontalLayout_5->addWidget(label_v2_chemin);

        spinBox_v2_chemin = new QSpinBox(widget_5);
        spinBox_v2_chemin->setObjectName(QString::fromUtf8("spinBox_v2_chemin"));

        horizontalLayout_5->addWidget(spinBox_v2_chemin);


        verticalLayout_3->addWidget(widget_5);

        pushButton_afficherChemin = new QPushButton(frame_2);
        pushButton_afficherChemin->setObjectName(QString::fromUtf8("pushButton_afficherChemin"));

        verticalLayout_3->addWidget(pushButton_afficherChemin);


        verticalLayout->addWidget(frame_2);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);


        horizontalLayout->addWidget(widget_2);

        displayWidget = new MeshViewerWidget(centralWidget);
        displayWidget->setObjectName(QString::fromUtf8("displayWidget"));

        horizontalLayout->addWidget(displayWidget);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 632, 22));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "MainWindow", nullptr));
        pushButton_chargement->setText(QCoreApplication::translate("MainWindow", "Charger OBJ", nullptr));
        label->setText(QCoreApplication::translate("MainWindow", "S\303\251lections courantes :", nullptr));
        label_4->setText(QCoreApplication::translate("MainWindow", "v : ", nullptr));
        pushButton_vertexMoins->setText(QCoreApplication::translate("MainWindow", "-", nullptr));
        labelVertex->setText(QCoreApplication::translate("MainWindow", "-1", nullptr));
        pushButton_vertexPlus->setText(QCoreApplication::translate("MainWindow", "+", nullptr));
        label_5->setText(QCoreApplication::translate("MainWindow", "e : ", nullptr));
        pushButton_edgeMoins->setText(QCoreApplication::translate("MainWindow", "-", nullptr));
        labelEdge->setText(QCoreApplication::translate("MainWindow", "-1", nullptr));
        pushButton_edgePlus->setText(QCoreApplication::translate("MainWindow", "+", nullptr));
        label_6->setText(QCoreApplication::translate("MainWindow", "f : ", nullptr));
        pushButton_faceMoins->setText(QCoreApplication::translate("MainWindow", "-", nullptr));
        labelFace->setText(QCoreApplication::translate("MainWindow", "-1", nullptr));
        pushButton_facePlus->setText(QCoreApplication::translate("MainWindow", "+", nullptr));
        pushButton_voisinage->setText(QCoreApplication::translate("MainWindow", "Passer en mode voisinage", nullptr));
        pushButton_bordure->setText(QCoreApplication::translate("MainWindow", "Afficher la bordure", nullptr));
        label_v1_chemin->setText(QCoreApplication::translate("MainWindow", "  v1 : ", nullptr));
        label_v2_chemin->setText(QCoreApplication::translate("MainWindow", "  v2 : ", nullptr));
        pushButton_afficherChemin->setText(QCoreApplication::translate("MainWindow", "Afficher le chemin", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
