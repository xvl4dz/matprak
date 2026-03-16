#include <iostream>
#include <vector>
#include <cmath>
#include <gmsh.h>
#include <set>

class HollowTorus {
public:
    HollowTorus(double R, double r_outer, double r_inner)
        : m_R(R)
        , m_r_outer(r_outer)
        , m_r_inner(r_inner) 
    {}

    int createGeometry() const {
        int outer_torus = createSolidTorus(m_R, m_r_outer, 1);
        int inner_torus = createSolidTorus(m_R, m_r_inner, 2);

        std::vector<std::pair<int, int>> out;
        std::vector<std::pair<int, int>> object = {{3, outer_torus}};
        std::vector<std::pair<int, int>> tool   = {{3, inner_torus}};
        std::vector<std::vector<std::pair<int, int>>> outMap;
        gmsh::model::occ::cut(object, tool, out, outMap);

        gmsh::model::occ::synchronize();

        if (out.empty() || out[0].first != 3) {
            throw std::runtime_error("Boolean cut did not produce a volume.");
        }
        return out[0].second;
    }

private:
    double m_R;
    double m_r_outer;
    double m_r_inner;

    static int createSolidTorus(double R, double r, int tag = 0) {
        return gmsh::model::occ::addTorus(0, 0, 0, R, r, tag, 2*M_PI);
    }
};

double computeMeshSize(double thickness, int minElements = 4) {
    return thickness / minElements;
}

int main(int argc, char** argv) {
    gmsh::initialize();

    gmsh::model::add("hollow_torus");

    const double majorRadius   = 5.0;      // R
    const double outerRadius   = 1.5;      // outer minor radius
    const double innerRadius   = 1.2;      // inner minor radius
    const double wallThickness = outerRadius - innerRadius; // 0.3

    if (wallThickness <= 0) {
        std::cerr << "Error: outer radius must be larger than inner radius.\n";
        return 1;
    }

    // Compute mesh size
    double meshSize = computeMeshSize(wallThickness, 4);
    std::cout << "Wall thickness = " << wallThickness
              << ", using mesh size = " << meshSize << std::endl;

    gmsh::option::setNumber("Mesh.MeshSizeMax", meshSize);

    // Create the hollow torus geometry
    HollowTorus torus(majorRadius, outerRadius, innerRadius);
    int volumeTag;
    try {
        volumeTag = torus.createGeometry();
    } catch (const std::exception& e) {
        std::cerr << "Geometry creation failed: " << e.what() << std::endl;
        return 1;
    }

    // Add a physical group
    gmsh::model::addPhysicalGroup(3, {volumeTag}, 1);
    gmsh::model::setPhysicalName(3, 1, "HollowTorus");

    gmsh::model::mesh::generate(3);

    gmsh::write("hollow_torus.msh");

    std::set<std::string> args(argv, argv + argc);
    if (args.find("-nopopup") == args.end()) {
        gmsh::fltk::run();
    }

    gmsh::finalize();
    return 0;
}