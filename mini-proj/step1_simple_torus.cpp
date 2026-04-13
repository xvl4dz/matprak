#include "common.h"
#include "Heat.h"

int main()
{
  auto mesh = std::make_shared<Mesh>("torus.xml");
  auto V = std::make_shared<Heat::FunctionSpace>(mesh);
  Function u(V);
  u = Constant(0.0);
  write_vtu(mesh, u, "torus_shape.vtu");
  std::cout << "Mesh written to torus_shape.vtu" << std::endl;
  return 0;
}
