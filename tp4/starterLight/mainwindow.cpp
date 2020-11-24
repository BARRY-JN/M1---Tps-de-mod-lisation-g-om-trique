#include "mainwindow.h"
#include "ui_mainwindow.h"

/* **** début de la partie à compléter **** */

//Colorie l'arête sélectionné
void MainWindow::showEdgeSelection(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);

    /* **** à compléter ! (Partie 1) ****
     * cette fonction utilise la variables de sélection edgeSelection qui est l'ID de
     * l'arête sélectionnée et qui est égale à -1 si la sélection est vide
     */

    if(vertexSelection>-1)
    {
        _mesh->set_color(_mesh->vertex_handle(vertexSelection), MyMesh::Color(0, 0, 255));
        _mesh->data(_mesh->vertex_handle(vertexSelection)).thickness = 12;
    }

    if(edgeSelection>-1)
    {
        _mesh->set_color(_mesh->edge_handle(edgeSelection), MyMesh::Color(254, 0, 0));
        _mesh->data(_mesh->edge_handle(edgeSelection)).thickness = 12;
    }

    if(faceSelection>-1)
        _mesh->set_color(_mesh->face_handle(faceSelection), MyMesh::Color(0, 0, 255));
        _mesh->data(_mesh->edge_handle(faceSelection)).thickness = 12;

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

//Marque les arêtes à supprimer
void MainWindow::collapseEdge(MyMesh* _mesh, int edgeID)
{
    /* **** à compléter ! (Partie 1) ****
     * cette fonction utilise l'opérateur collapse pour supprimer l'arête d'index edgeID
     * Attention à ne pas oublier garbage_collection() !
     */
    EdgeHandle eh = _mesh->edge_handle(edgeID);
    HalfedgeHandle heh = _mesh->halfedge_handle(eh, 0);

    MyMesh::VertexHandle v1 = _mesh->to_vertex_handle(heh);
    MyMesh::VertexHandle v0 = _mesh->from_vertex_handle(heh);
    MyMesh::Point point0 = _mesh->point(v0);
    MyMesh::Point point1 = _mesh->point(v1);

    MyMesh::Point point = (point1 + point0)/2;
    _mesh->point(v1) = point;
    mesh.collapse(heh);

    // permet de nettoyer le maillage et de garder la cohérence des indices après un collapse
    _mesh->garbage_collection();
}

// fonction pratique pour faire des tirages aléatoires
int randInt(int low, int high){return qrand() % ((high + 1) - low) + low;}

void MainWindow::decimation(MyMesh* _mesh, int percent, QString method)
{
    /* **** à compléter ! (Partie 2 et 3) ****
     * Cette fonction supprime des arêtes jusqu'à atteindre un pourcentage d'arêtes restantes, selon un critère donné
     * percent : pourcentage de l'objet à garder
     * method  : la méthode à utiliser parmis : "Aléatoire", "Par taille", "Par angle", "Par planéité"
     */

	int start_n_edge = _mesh->n_edges();
    if(method == "Aléatoire")
    {
        while(_mesh->n_edges() > start_n_edge * static_cast<float>(percent) / 100){
           int rand ;
           while(true){
               rand = randInt(0,_mesh->n_edges());
               if(_mesh->is_collapse_ok(_mesh->halfedge_handle(_mesh->edge_handle(static_cast <unsigned int>(rand)),0))){
                    break;
               }

           }
           collapseEdge(_mesh,rand);
       }

    }
    else if(method == "Par taille")
    {
        QMultiMap <float, int> ed;
        for (MyMesh::EdgeIter  cur = _mesh->edges_begin();cur != _mesh->edges_end();cur++)
        {


            MyMesh::HalfedgeHandle heh = _mesh->halfedge_handle(*cur,0);
            MyMesh::VertexHandle v1 = _mesh->to_vertex_handle(heh);
            MyMesh::VertexHandle v0 = _mesh->from_vertex_handle(heh);
            MyMesh::Point point0 = _mesh->point(v0);
            MyMesh::Point point1 = _mesh->point(v1);

            float size = (point1 + point0).norm();
            ed.insert(size,(*cur).idx());
        }

        int i = 0;
        while(_mesh->n_edges() > start_n_edge * static_cast<float>(percent) / 100 && i < ed.size()){
           if(_mesh->edge_handle(static_cast <unsigned int>(ed.first())).is_valid()){
                if(_mesh->is_collapse_ok(_mesh->halfedge_handle(_mesh->edge_handle((static_cast <unsigned int>(ed.first()))),0))){
                    collapseEdge(_mesh,ed.first());
                }
            }
            i++;
            ed.remove(ed.key(ed.first()));
        }
    }
    else if(method == "Par angle")
    {
        QMultiMap <float, int> ed;
        for(MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
        {
            FaceHandle fh = *curFace;
            for (MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(fh); curVertex.is_valid(); curVertex ++){

                VertexHandle vh = *curVertex;

                QVector<VertexHandle> listePoints;
                QVector<VertexHandle> listePointsOnFace;
                QVector<float> vectors;

                //Les points voisins de notre vertex
                for (MyMesh::VertexVertexIter curVertex = _mesh->vv_iter(vh); curVertex.is_valid(); curVertex ++){
                    VertexHandle v = *curVertex;
                    listePoints.append(v);
                }

                //les points voisins qui appartiennent aussi a la meme face
                for(MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(fh); curVertex.is_valid(); curVertex ++){
                    VertexHandle v = *curVertex;
                    if(listePoints.contains(v) && v.idx() !=vh.idx()){
                        listePointsOnFace.append(v);
                    }
                }

                //Nouveaux vecteurs avec les nouveaux points
                for(int i=0; i<listePointsOnFace.size();i++)
                {
                    vectors.append((_mesh->point(listePointsOnFace[i])[0])-(_mesh->point(vh)[0]));
                    vectors.append((_mesh->point(listePointsOnFace[i])[1])-(_mesh->point(vh)[1]));
                    vectors.append((_mesh->point(listePointsOnFace[i])[2])-(_mesh->point(vh)[2]));

                }

                //Normalisation des vecteurs
                QVector<float> norme;
                float tmp = vectors[0]*vectors[0]+vectors[1]*vectors[1]+vectors[2]*vectors[2];

                norme.append(sqrt(tmp));

                tmp = vectors[3]*vectors[3]+vectors[4]*vectors[4]+vectors[5]*vectors[5];
                norme.append(sqrt(tmp));

                for(int i=0; i<listePointsOnFace.size(); i++){
                    vectors[i*3] = vectors[i*3]/norme[i];
                    vectors[i*3+1] = vectors[i*3+1]/norme[i];
                    vectors[i*3+2] = vectors[i*3+2]/norme[i];
                }

                float a = vectors[0]*vectors[0+3];
                float b = vectors[1]*vectors[1+3];
                float c = vectors[2]*vectors[2+3];
                float prodScal = a+b+c;
                abs(prodScal);

                float angle = acos(prodScal);
                ed.insert(angle, vh.idx());
             }
        }

        int i = 0;
        while(_mesh->n_edges()>start_n_edge*static_cast<float>(percent)/100 && i<ed.size()){
            if(_mesh->edge_handle(static_cast <unsigned int>(ed.first())).is_valid()){
                if(_mesh->is_collapse_ok(_mesh->halfedge_handle(_mesh->edge_handle((static_cast <unsigned int>(ed.first()))),0))){
                    collapseEdge(_mesh,ed.first());
                }
            }
            i++;
            ed.remove(ed.key(ed.first()));
        }

    }
    else if(method == "Par planéité"){
        unsigned nb_aretes_del = static_cast<unsigned>(_mesh->n_edges() * static_cast<unsigned>(percent) / 100);
        std::vector<float> planeite = std::vector<float>(_mesh->n_edges(), 0.0f);

        float angle_dihedral;
        VertexHandle vh_to;
        VertexHandle vh_from;
        HalfedgeHandle heh;
        float cpt;
        for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++){
            heh = _mesh->halfedge_handle(*curEdge, 0);
            vh_to = _mesh->to_vertex_handle(heh);
            vh_from = _mesh->from_vertex_handle(heh);
            cpt = 0;
            angle_dihedral = 0;
            for (MyMesh::VertexEdgeIter curVEdge = _mesh->ve_iter(vh_to); curVEdge.is_valid(); curVEdge ++){
                if(!_mesh->is_boundary(*curVEdge)){
                    if(_mesh->calc_dihedral_angle(*curVEdge) == 0.0f) angle_dihedral += static_cast<float>(M_PI);
                    else angle_dihedral += _mesh->calc_dihedral_angle(*curVEdge);
                    cpt++;
                }
            }
            for (MyMesh::VertexEdgeIter curVEdge = _mesh->ve_iter(vh_from); curVEdge.is_valid(); curVEdge ++){
                if(!_mesh->is_boundary(*curVEdge)){
                    if(_mesh->calc_dihedral_angle(*curVEdge) == 0.0f){
                        angle_dihedral += static_cast<float>(M_PI);
                    }else angle_dihedral += _mesh->calc_dihedral_angle(*curVEdge);
                    cpt++;
                }
            }
            planeite.at(static_cast<unsigned long>(curEdge->idx())) = std::abs(angle_dihedral)/cpt;
        }

        unsigned maxElementIndex;
        EdgeHandle eh;
        for(unsigned i = 0; i < nb_aretes_del;){
            maxElementIndex = static_cast<unsigned>(std::max_element(planeite.begin(), planeite.end()) - planeite.begin());

            eh = _mesh->edge_handle(maxElementIndex);
            planeite.at(static_cast<unsigned long>(eh.idx())) = -1.0f;

            heh = _mesh->halfedge_handle(eh, 0);
            if(!_mesh->is_collapse_ok(heh)) continue;


            vh_to = _mesh->to_vertex_handle(heh);
            vh_from = _mesh->from_vertex_handle(heh);
            if(!_mesh->is_boundary(eh)){
                if(_mesh->is_boundary(vh_from) || _mesh->is_boundary(vh_to)) continue;
            }
            collapseEdge(_mesh, eh.idx());

            if(!_mesh->is_boundary(eh)){
                i += 3;
            }
            else i += 2;
        }
        _mesh->garbage_collection();
        qDebug() << "nb aretes finale " << _mesh->n_edges();
    }
    else
    {
        qDebug() << "Méthode inconnue !!!";
    }

}

/* **** début de la partie boutons et IHM **** */
void MainWindow::updateEdgeSelectionIHM()
{
    /* **** à compléter ! (Partie 3) ****
     * Cette fonction met à jour l'interface, les critères pourrons être affichés dans la zone de texte pour les vérifier
     */

    QString infos = "";
    infos = infos + "Surface : " + QString::number(0) + "\n";
    infos = infos + "C1 : " + QString::number(0) + "\n";
    infos = infos + "C2 : " + QString::number(0) + "\n";
    infos = infos + "C3 : " + QString::number(0) + "\n";
    ui->infoEdgeSelection->setPlainText(infos);

    ui->labelEdge->setText(QString::number(edgeSelection));

    // on montre la nouvelle sélection
    showEdgeSelection(&mesh);
}
/* **** fin de la partie à compléter **** */

void MainWindow::on_pushButton_edgeMoins_clicked()
{
    // mise à jour de l'interface
    edgeSelection = edgeSelection - 1;
    updateEdgeSelectionIHM();
}

void MainWindow::on_pushButton_edgePlus_clicked()
{
    // mise à jour de l'interface
    edgeSelection = edgeSelection + 1;
    updateEdgeSelectionIHM();
}

void MainWindow::on_pushButton_delSelEdge_clicked()
{
    // on supprime l'arête d'indice edgeSelection
    collapseEdge(&mesh, edgeSelection);

    // on actualise la sélection
    showEdgeSelection(&mesh);
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

void MainWindow::on_pushButton_decimate_clicked()
{
    decimation(&mesh, ui->horizontalSlider->value(), ui->comboBox->currentText());
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
    edgeSelection = -1;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

