#include "mainwindow.h"
#include "ui_mainwindow.h"


//Calcul de l'aire de la face faceID
float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    FaceHandle fh = _mesh->face_handle(faceID);

    HalfedgeHandle heh = _mesh->halfedge_handle(fh);
    VertexHandle vh_1 = _mesh->to_vertex_handle(heh);
    heh = _mesh->next_halfedge_handle(heh);
    VertexHandle vh_2 = _mesh->to_vertex_handle(heh);
    heh = _mesh->next_halfedge_handle(heh);
    VertexHandle vh_3 = _mesh->to_vertex_handle(heh);
    //(Base x Hauteur)/2
    return ( (_mesh->point(vh_2) - _mesh->point(vh_1)) % (_mesh->point(vh_3) - _mesh->point(vh_1)) ).norm() / 2;
}

//l'angle signé entre les normales des faces voisines à l'arête allant de vertID0 à vertID1
float MainWindow::angleFF(MyMesh* _mesh, int faceID0,  int faceID1, int vertID0, int vertID1)
{
    FaceHandle fh0 = _mesh->face_handle(faceID0);
    FaceHandle fh1 = _mesh->face_handle(faceID1);

    VertexHandle vh0 = _mesh->vertex_handle(vertID0);
    VertexHandle vh1 = _mesh->vertex_handle(vertID1);

    if(((_mesh->normal(fh0) % _mesh->normal(fh1)) | (_mesh->point(vh1) - _mesh->point(vh0))) < 0){
        return (std::acos(_mesh->normal(fh0) | _mesh->normal(fh1)))*-1;
    }else{
        return std::acos(_mesh->normal(fh0) | _mesh->normal(fh1));
    }
}

//l'angle absolu entre les deux arêtes d"une face depuis un de ses sommets
float MainWindow::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
{
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    FaceHandle fh = _mesh->face_handle(faceID);

    HalfedgeHandle heh = _mesh->halfedge_handle(fh);

    VertexHandle vh1;
    VertexHandle vh2;
    VertexHandle vh3;

    for(int i = 0; i < 3; i++){
        vh1 = _mesh->to_vertex_handle(heh);
        if(vh1 == vh){
            vh2 = _mesh->from_vertex_handle(heh);
            heh = _mesh->next_halfedge_handle(heh);
            vh3 = _mesh->to_vertex_handle(heh);
            break;
        }
        heh = _mesh->next_halfedge_handle(heh);
    }

    return std::acos((_mesh->point(vh2) - _mesh->point(vh1)).normalize()|(_mesh->point(vh3) - _mesh->point(vh1)).normalize());
}

//Aire barycentrique d'un sommet
float MainWindow::Bary(MyMesh* _mesh, int vertexID){
    VertexHandle vh = _mesh->vertex_handle(vertexID);

    float Bary = 0;
    FaceHandle fh;
    for (MyMesh::VertexFaceIter vf_it=_mesh->vf_iter(vh); vf_it.is_valid(); ++vf_it)
    {
        fh = *vf_it;
        Bary += faceArea(_mesh, fh.idx());
    }

    return Bary/3;
}

//Somme des angleFF autour d'un sommet
float MainWindow::angleH_Curv(MyMesh* _mesh, int vertexID){
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    float angle = 0;
    HalfedgeHandle heh;
    FaceHandle fh0;
    FaceHandle fh1;
    VertexHandle vh0;
    VertexHandle vh1;
    for (MyMesh::VertexOHalfedgeIter voh_it=_mesh->voh_iter(vh); voh_it.is_valid(); voh_it++)
    {
        heh = *voh_it;

        fh0 = _mesh->face_handle(heh);
        if(!fh0.is_valid()) continue;

        fh1 = _mesh->opposite_face_handle(heh);
        if(!fh1.is_valid()) continue;

        vh0 = _mesh->from_vertex_handle(heh);
        vh1 = _mesh->to_vertex_handle(heh);

        angle += angleFF(_mesh, fh0.idx(), fh1.idx(), vh0.idx(), vh1.idx()) * ( _mesh->point(vh1) - _mesh->point(vh0) ).norm();
    }
    return angle;
}

//Coloration de tout le maillage à l'aide de la formule de la Courbure moyenne en un certain point
void MainWindow::H_Curv(MyMesh* _mesh)
{
    float aireBary;
	float angle;
	VertexHandle vh_cur;
	for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert!=_mesh->vertices_end(); ++curVert)
	{
		vh_cur = *curVert;
		aireBary = Bary(_mesh, vh_cur.idx());
		angle = angleH_Curv(_mesh, vh_cur.idx());

		_mesh->data(vh_cur).value = (1 / (4*aireBary) ) * angle;
	}
}

//Somme des angleEE
float MainWindow::angleK_Curv(MyMesh* _mesh, int vertexID){
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    float angle = 0;
    FaceHandle fh;
    for (MyMesh::VertexFaceIter vf_it=_mesh->vf_iter(vh); vf_it.is_valid(); vf_it++)
    {
        fh = *vf_it;
        angle += angleEE(_mesh, vh.idx(), fh.idx());
    }
    return angle;
}


//Coloration de tout le maillage à l'aide de la formule de la Courbure gaussienne
void MainWindow::K_Curv(MyMesh* _mesh)
{
    float Bary2;
    float angle;
    VertexHandle vh_cur;
    for (MyMesh::VertexIter curVert=_mesh->vertices_begin(); curVert!=_mesh->vertices_end(); curVert++)
    {
        vh_cur = *curVert;
        Bary2 = Bary(_mesh, vh_cur.idx());
        angle = angleK_Curv(_mesh, vh_cur.idx());

        _mesh->data(vh_cur).value = (1/Bary2)*(2*M_PI - angle);
    }
}
/* **** fin de la partie à compléter **** */



/* **** début de la partie boutons et IHM **** */
void MainWindow::on_pushButton_H_clicked()
{
    H_Curv(&mesh);
    displayMesh(&mesh, DisplayMode::TemperatureMap); // permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

void MainWindow::on_pushButton_K_clicked()
{
    K_Curv(&mesh);
    displayMesh(&mesh, DisplayMode::TemperatureMap); // permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

/*  Cette fonction est à utiliser UNIQUEMENT avec le fichier testAngleArea.obj
    Elle est appelée par le bouton "Test angles/aires"

    Elle permet de vérifier les fonctions faceArea, angleFF et angleEE.
    Elle doit afficher :

    Aire de la face 0 : 2
    Aire de la face 1 : 2
    Angle entre les faces 0 et 1 : 1.5708
    Angle entre les faces 1 et 0 : -1.5708
    Angle au sommet 1 sur la face 0 : 0.785398 */
void MainWindow::on_pushButton_angleArea_clicked()
{
    qDebug() << "Aire de la face 0 :" << faceArea(&mesh, 0);
    qDebug() << "Aire de la face 1 :" << faceArea(&mesh, 1);

    qDebug() << "Angle entre les faces 0 et 1 :" << angleFF(&mesh, 0, 1, 1, 2);
    qDebug() << "Angle entre les faces 1 et 0 :" << angleFF(&mesh, 1, 0, 1, 2);

    qDebug() << "Angle au sommet 1 sur la face 0 :" << angleEE(&mesh, 1, 0);
    qDebug() << "Angle au sommet 3 sur la face 1 :" << angleEE(&mesh, 3, 1);
}

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}
/* **** fin de la partie boutons et IHM **** */

/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, DisplayMode mode)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(mode == DisplayMode::TemperatureMap)
    {
        QVector<float> values;
        for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
            values.append(fabs(_mesh->data(*curVert).value));
        qSort(values);

        float range = values.at(values.size()*0.8);

        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if(mode == DisplayMode::Normal)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if(mode == DisplayMode::ColorShading)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}
