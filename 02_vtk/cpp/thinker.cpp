#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;

class CalcNode
{
friend class CalcMesh;
protected:
    double x, y, z;          // current coordinates
    double smth;              // scalar field
    double vx, vy, vz;        // velocity
public:
    CalcNode() : x(0), y(0), z(0), smth(0), vx(0), vy(0), vz(0) {}
    CalcNode(double x, double y, double z) : x(x), y(y), z(z), smth(0), vx(0), vy(0), vz(0) {}
};


class Element
{
friend class CalcMesh;
protected:
    unsigned long nodesIds[4];
};

class CalcMesh
{
protected:
    vector<CalcNode> nodes;
    vector<Element> elements;
    vector<double> initialX, initialY, initialZ;   // initial positions
    double pivotX, pivotY, pivotZ;                  // centroid for stretching

public:
    CalcMesh(const vector<double>& nodesCoords, const vector<size_t>& tetrsPoints)
    {
        size_t numNodes = nodesCoords.size() / 3;
        nodes.resize(numNodes);
        initialX.resize(numNodes);
        initialY.resize(numNodes);
        initialZ.resize(numNodes);

        double sumX = 0, sumY = 0, sumZ = 0;
        for (size_t i = 0; i < numNodes; ++i) {
            double x = nodesCoords[i*3];
            double y = nodesCoords[i*3 + 1];
            double z = nodesCoords[i*3 + 2];
            nodes[i] = CalcNode(x, y, z);
            initialX[i] = x;
            initialY[i] = y;
            initialZ[i] = z;
            sumX += x;
            sumY += y;
            sumZ += z;
        }
        pivotX = sumX / numNodes;
        pivotY = sumY / numNodes;
        pivotZ = sumZ / numNodes;

        size_t numTetrs = tetrsPoints.size() / 4;
        elements.resize(numTetrs);
        for (size_t i = 0; i < numTetrs; ++i) {
            elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
        }
    }

    void updatePositions(double time)
    {
        double A = 37.5;               // amplitude (jump height)
        double g = 50.0;                // gravity
        double T = sqrt(8.0 * A / g);   // period for a complete bounce

        double t_mod = fmod(time, T);
        double frac = t_mod / T;
        // Parabolic motion: z_offset = A * (1 - (2*frac - 1)^2)
        double z_offset = A * (1.0 - pow(2.0*frac - 1.0, 2.0));

        double stretchAmplitude = 0.1;
        double scale = 1.0 + stretchAmplitude * (z_offset / A);

        for (size_t i = 0; i < nodes.size(); ++i) {
            double stretchedZ = pivotZ + (initialZ[i] - pivotZ) * scale;
            nodes[i].x = initialX[i];
            nodes[i].y = initialY[i];
            nodes[i].z = stretchedZ + z_offset;
        }
    }

    void updateVelocityField(double time)
    {
        double A = 37.5;
        double g = 50.0;
        double T = sqrt(8.0 * A / g);
        double t_mod = fmod(time, T);
        double frac = t_mod / T;
        // vz = derivative of z_offset: d/dt [A*(1 - (2*frac-1)^2)] = - (4A/T) * (2*frac - 1)
        double vz = - (4.0 * A / T) * (2.0*frac - 1.0);

        for (size_t i = 0; i < nodes.size(); ++i) {
            nodes[i].vx = 0.0;
            nodes[i].vy = 0.0;
            nodes[i].vz = vz;
        }
        cout << "Time " << time << ", vz = " << vz << ", z_offset = " << (nodes[0].z - initialZ[0]) << endl;
    }

    void updateScalarField(double time)
    {
        double waveSpeed = 5.0;
        double waveFreq = 0.5;

        for (size_t i = 0; i < nodes.size(); ++i) {
            double dx = nodes[i].x - pivotX;
            double dy = nodes[i].y - pivotY;
            double dz = nodes[i].z - pivotZ;
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            double wave = sin(2 * M_PI * (waveFreq * time - dist / waveSpeed));
            nodes[i].smth = wave;
        }
    }

    void snapshot(unsigned int snap_number, double time)
    {
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("thoughtWave");
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        for (size_t i = 0; i < nodes.size(); ++i) {
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);
            double v[3] = {nodes[i].vx, nodes[i].vy, nodes[i].vz};
            vel->InsertNextTuple(v);
            smth->InsertNextValue(nodes[i].smth);
        }

        unstructuredGrid->SetPoints(dumpPoints);
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(smth);

        for (size_t i = 0; i < elements.size(); ++i) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId(0, elements[i].nodesIds[0]);
            tetra->GetPointIds()->SetId(1, elements[i].nodesIds[1]);
            tetra->GetPointIds()->SetId(2, elements[i].nodesIds[2]);
            tetra->GetPointIds()->SetId(3, elements[i].nodesIds[3]);
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        string fileName = "thinker-step-" + to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};

int main()
{
    const unsigned int GMSH_TETR_CODE = 4;
    double tau = 0.05;
    int numSteps = 100;

    gmsh::initialize();
    gmsh::model::add("thinker");

    try {
        gmsh::merge("lowest-poly-thinker.stl");
    } catch(...) {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return -1;
    }

    gmsh::model::mesh::removeDuplicateNodes();

    double angle = 60.0;
    bool includeBoundary = false;
    bool forceParametrizablePatches = true;
    double curveAngle = 180.0;

    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                        forceParametrizablePatches,
                                        curveAngle * M_PI / 180.);
    gmsh::model::mesh::createGeometry();

    std::vector<std::pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);
    if (surfaces.empty()) {
        gmsh::logger::write("No surfaces found – aborting");
        gmsh::finalize();
        return -2;
    }
    std::vector<int> surfaceTags;
    for (auto& s : surfaces) surfaceTags.push_back(s.second);
    int loop = gmsh::model::geo::addSurfaceLoop(surfaceTags);
    gmsh::model::geo::addVolume({loop});
    gmsh::model::geo::synchronize();

    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", "3");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    try {
        gmsh::model::mesh::generate(3);
    } catch(...) {
        gmsh::logger::write("3D meshing failed – check STL or parameters");
        gmsh::finalize();
        return -3;
    }

    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);

    std::vector<std::size_t>* tetras = nullptr;
    for (size_t i = 0; i < elementTypes.size(); ++i) {
        if (elementTypes[i] == GMSH_TETR_CODE) {
            tetras = &elementNodeTags[i];
            break;
        }
    }
    if (!tetras) {
        cerr << "No tetrahedral elements found!" << endl;
        gmsh::finalize();
        return -4;
    }

    cout << "Nodes: " << nodeTags.size() << ", Tetras: " << tetras->size()/4 << endl;

    CalcMesh mesh(nodesCoord, *tetras);

    gmsh::finalize();

    double time = 0.0;
    for (int step = 0; step < numSteps; ++step) {
        cout << "Step " << step << ", time = " << time << endl;

        mesh.updatePositions(time);
        mesh.updateVelocityField(time);
        mesh.updateScalarField(time);

        mesh.snapshot(step, time);

        time += tau;
    }

    return 0;
}