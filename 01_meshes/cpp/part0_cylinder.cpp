#include <set>
#include <vector>
#include <gmsh.h>

int main(int argc, char **argv)
{
    gmsh::initialize();
    gmsh::model::add("cylinder_builtin");

    double lc = 1e-2;
    double r = 0.3;
    double cx = 0.5, cy = 0.5;
    double h = 1.0;

    // ---------- Bottom circle (z = 0) ----------
    int centerBottom = 1;
    gmsh::model::geo::addPoint(cx, cy, 0, lc, centerBottom);

    // Points on bottom circle (counter‑clockwise)
    int pBottomRight = 2;   gmsh::model::geo::addPoint(cx + r, cy, 0, lc, pBottomRight);
    int pBottomTop   = 3;   gmsh::model::geo::addPoint(cx, cy + r, 0, lc, pBottomTop);
    int pBottomLeft  = 4;   gmsh::model::geo::addPoint(cx - r, cy, 0, lc, pBottomLeft);
    int pBottomBottom= 5;   gmsh::model::geo::addPoint(cx, cy - r, 0, lc, pBottomBottom);

    // Arcs (counter‑clockwise)
    gmsh::model::geo::addCircleArc(pBottomRight, centerBottom, pBottomTop,    1);
    gmsh::model::geo::addCircleArc(pBottomTop,   centerBottom, pBottomLeft,   2);
    gmsh::model::geo::addCircleArc(pBottomLeft,  centerBottom, pBottomBottom, 3);
    gmsh::model::geo::addCircleArc(pBottomBottom, centerBottom, pBottomRight, 4);

    // Curve loop and bottom surface
    gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);

    // ---------- Top circle (z = h) ----------
    int centerTop = 6;
    gmsh::model::geo::addPoint(cx, cy, h, lc, centerTop);

    // Points on top circle (same order as bottom)
    int pTopRight = 7;   gmsh::model::geo::addPoint(cx + r, cy, h, lc, pTopRight);
    int pTopTop   = 8;   gmsh::model::geo::addPoint(cx, cy + r, h, lc, pTopTop);
    int pTopLeft  = 9;   gmsh::model::geo::addPoint(cx - r, cy, h, lc, pTopLeft);
    int pTopBottom= 10;  gmsh::model::geo::addPoint(cx, cy - r, h, lc, pTopBottom);

    // Arcs (counter‑clockwise)
    gmsh::model::geo::addCircleArc(pTopRight, centerTop, pTopTop,    5);
    gmsh::model::geo::addCircleArc(pTopTop,   centerTop, pTopLeft,   6);
    gmsh::model::geo::addCircleArc(pTopLeft,  centerTop, pTopBottom, 7);
    gmsh::model::geo::addCircleArc(pTopBottom, centerTop, pTopRight, 8);

    // Curve loop and top surface
    gmsh::model::geo::addCurveLoop({5, 6, 7, 8}, 2);
    gmsh::model::geo::addPlaneSurface({2}, 2);

    // ---------- Vertical lines connecting corresponding points ----------
    int lineRight = 9;  gmsh::model::geo::addLine(pBottomRight, pTopRight, lineRight);
    int lineTop   = 10; gmsh::model::geo::addLine(pBottomTop,   pTopTop,   lineTop);
    int lineLeft  = 11; gmsh::model::geo::addLine(pBottomLeft,  pTopLeft,  lineLeft);
    int lineBottom= 12; gmsh::model::geo::addLine(pBottomBottom, pTopBottom, lineBottom);

    // ---------- Four lateral ruled surfaces ----------
    gmsh::model::geo::addCurveLoop({1, 10, -5, -9}, 3);
    gmsh::model::geo::addSurfaceFilling({3}, 3);

    gmsh::model::geo::addCurveLoop({2, 11, -6, -10}, 4);
    gmsh::model::geo::addSurfaceFilling({4}, 4);

    gmsh::model::geo::addCurveLoop({3, 12, -7, -11}, 5);
    gmsh::model::geo::addSurfaceFilling({5}, 5);

    gmsh::model::geo::addCurveLoop({4, 9, -8, -12}, 6);
    gmsh::model::geo::addSurfaceFilling({6}, 6);

    gmsh::model::geo::synchronize();

    // Create a surface loop from all six surfaces
    std::vector<int> surfaceTags = {1, 2, 3, 4, 5, 6};
    int surfaceLoop = gmsh::model::geo::addSurfaceLoop(surfaceTags, 1);

    // Create volume from the surface loop
    int vol = gmsh::model::geo::addVolume({surfaceLoop});

    gmsh::model::mesh::generate(3);
    gmsh::write("cylinder.msh");

    std::set<std::string> args(argv, argv + argc);
    if (!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}