#include <set>
#include <cmath>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("thinker");

  try {
    gmsh::merge("../lowest-poly-thinker.stl");
  } catch(...) {
    gmsh::logger::write("Could not load STL mesh: bye!");
    gmsh::finalize();
    return 0;
  }

  gmsh::model::mesh::removeDuplicateNodes();


  double angle = 60;
  bool includeBoundary = false;
  bool forceParametrizablePatches = true;
  double curveAngle = 180;

  gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                      forceParametrizablePatches,
                                      curveAngle * M_PI / 180.);

  gmsh::model::mesh::createGeometry();

  std::vector<std::pair<int, int> > s;
  gmsh::model::getEntities(s, 2);
  if (s.empty()) {
    gmsh::logger::write("No surfaces found – aborting volume creation");
    gmsh::finalize();
    return 1;
  }

  std::vector<int> sl;
  for(auto surf : s) sl.push_back(surf.second);
  int l = gmsh::model::geo::addSurfaceLoop(sl);
  gmsh::model::geo::addVolume({l});

  gmsh::model::geo::synchronize();

  bool funny = false;
  int f = gmsh::model::mesh::field::add("MathEval");
  gmsh::model::mesh::field::setString(f, "F", funny ? "2*sin((x+y)/5) + 3" : "4");
  gmsh::model::mesh::field::setAsBackgroundMesh(f);

  try {
    gmsh::model::mesh::generate(3);
  } catch(...) {
    gmsh::logger::write("3D meshing failed – generating 2D surface mesh instead");
    gmsh::model::mesh::generate(2);
  }

  gmsh::write("thinker.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();
  return 0;
}