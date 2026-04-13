#ifndef COMMON_H
#define COMMON_H

#include <dolfin.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkTetra.h>
#include <vector>

using namespace dolfin;

void write_vtu(std::shared_ptr<const Mesh> mesh, const Function& u, const std::string& filename) {
    auto vtkGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto temperature = vtkSmartPointer<vtkDoubleArray>::New();
    temperature->SetName("temperature");

    std::vector<double> vertex_values;
    u.compute_vertex_values(vertex_values, *mesh);

    const std::vector<double>& coords = mesh->coordinates();
    for (std::size_t i = 0; i < mesh->num_vertices(); ++i) {
        points->InsertNextPoint(coords[3*i], coords[3*i+1], coords[3*i+2]);
        temperature->InsertNextValue(vertex_values[i]);
    }
    vtkGrid->SetPoints(points);
    vtkGrid->GetPointData()->AddArray(temperature);

    const std::vector<unsigned int>& cells = mesh->cells();
    for (std::size_t i = 0; i < mesh->num_cells(); ++i) {
        auto tetra = vtkSmartPointer<vtkTetra>::New();
        for (std::size_t j = 0; j < 4; ++j) tetra->GetPointIds()->SetId(j, cells[4*i + j]);
        vtkGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
    }

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(vtkGrid);
    writer->Write();
}

class CoolingRings : public SubDomain {
public:
    std::vector<double> centers;
    double width;
    CoolingRings(std::vector<double> c, double w = 0.3) : centers(c), width(w) {}
    bool inside(const Array<double>& x, bool on_boundary) const {
        if (!on_boundary) return false;
        double phi = std::atan2(x[1], x[0]);
        for (double c : centers) {
            double diff = std::abs(phi - c);
            while (diff > M_PI) diff = std::abs(diff - 2.0 * M_PI);
            if (diff < width) return true;
        }
        return false;
    }
};

class HeatSource : public Expression {
public:
    double xc, yc, zc, A, s;
    HeatSource() : xc(0), yc(0), zc(0), A(1500.0), s(0.15) {}
    void eval(Array<double>& values, const Array<double>& x) const {
        double r2 = std::pow(x[0]-xc, 2) + std::pow(x[1]-yc, 2) + std::pow(x[2]-zc, 2);
        values[0] = A * std::exp(-r2 / (s*s));
    }
};

#endif
