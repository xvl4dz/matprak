#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("cube");

  double lc = 2e-2;

  for (int x = 0; x <= 1; ++x)
    for (int y = 0; y <= 1; ++y)
      for (int z = 0; z <= 1; ++z) 
        gmsh::model::geo::addPoint(x, y, z, lc*(1 + 3*x)*(1 + 3*y)*(1 + 3*z), 100*x + 10*y + z);

  gmsh::model::geo::addLine(0, 1, 1);
  gmsh::model::geo::addLine(1, 11, 2);
  gmsh::model::geo::addLine(11, 10, 3);
  gmsh::model::geo::addLine(10, 0, 4);

  gmsh::model::geo::addLine(0, 100, 5);
  gmsh::model::geo::addLine(100, 110, 6);
  gmsh::model::geo::addLine(110, 10, 7);

  gmsh::model::geo::addLine(101, 1, 8);
  gmsh::model::geo::addLine(101, 100, 9);
  gmsh::model::geo::addLine(111, 101, 10);
  gmsh::model::geo::addLine(11, 111, 11);
  gmsh::model::geo::addLine(110, 111, 12);

  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
  gmsh::model::geo::addPlaneSurface({1}, 1);

  gmsh::model::geo::addCurveLoop({4, 5, 6, 7}, 2);
  gmsh::model::geo::addPlaneSurface({2}, 2);

  gmsh::model::geo::addCurveLoop({2, 11, 10, 8}, 3);
  gmsh::model::geo::addPlaneSurface({3}, 3);

  gmsh::model::geo::addCurveLoop({6, 12, 10, 9}, 4);
  gmsh::model::geo::addPlaneSurface({4}, 4);

  gmsh::model::geo::addCurveLoop({3, -7, 12, -11}, 5);
  gmsh::model::geo::addPlaneSurface({5}, 5);
  
  gmsh::model::geo::addCurveLoop({1, -8, 9, -5}, 6);
  gmsh::model::geo::addPlaneSurface({6}, 6);

  gmsh::model::geo::addSurfaceLoop({1, 2, 3, 4, 5, 6}, 1);
  gmsh::model::geo::addVolume({1});

  gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(2);

  gmsh::write("cube.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}

