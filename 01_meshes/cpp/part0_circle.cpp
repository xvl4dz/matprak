#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
    gmsh::initialize();
    gmsh::model::add("cylinder");

    double lc = 1e-2;
    double r = 0.3;
    double cx = 0.5, cy = 0.5;

    gmsh::model::geo::addPoint(cx, cy, 0, lc, 1);

    gmsh::model::geo::addPoint(cx + r, cy, 0, lc, 2);   // right
    gmsh::model::geo::addPoint(cx, cy + r, 0, lc, 3);   // top
    gmsh::model::geo::addPoint(cx - r, cy, 0, lc, 4);   // left
    gmsh::model::geo::addPoint(cx, cy - r, 0, lc, 5);   // bottom

    gmsh::model::geo::addCircleArc(2, 1, 3, 1);
    gmsh::model::geo::addCircleArc(3, 1, 4, 2);
    gmsh::model::geo::addCircleArc(4, 1, 5, 3);
    gmsh::model::geo::addCircleArc(5, 1, 2, 4);

    gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);

    gmsh::model::geo::synchronize();

    int ps = gmsh::model::addPhysicalGroup(2, {1});
    gmsh::model::setPhysicalName(2, ps, "My circle");

    gmsh::model::mesh::generate(2);
    gmsh::write("circle.msh");

    std::set<std::string> args(argv, argv + argc);
    if (!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}