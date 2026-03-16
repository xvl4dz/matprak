#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

// MODIFIED: Two heat sources
class TwinSource : public Expression {
  void eval(Array<double>& values, const Array<double>& x) const {
    double s1 = exp(-(pow(x[0]-0.2, 2) + pow(x[1]-0.2, 2)) / 0.01);
    double s2 = exp(-(pow(x[0]-0.8, 2) + pow(x[1]-0.8, 2)) / 0.01);
    values[0] = 15 * (s1 + s2);
  }
};

// MODIFIED: Dirichlet on ALL FOUR edges (on_boundary)
class ClosedBox : public SubDomain {
  bool inside(const Array<double>& x, bool on_boundary) const {
    return on_boundary; 
  }
};

int main() {
  // MODIFIED: Doubled resolution (64x64)
  auto mesh = std::make_shared<UnitSquareMesh>(64, 64);
  auto V = std::make_shared<Poisson::FunctionSpace>(mesh);

  auto u0 = std::make_shared<Constant>(0.0);
  auto boundary = std::make_shared<ClosedBox>();
  DirichletBC bc(V, u0, boundary);

  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  L.f = std::make_shared<TwinSource>(); //
  L.g = std::make_shared<Constant>(0.0);

  Function u(V);
  solve(a == L, u, bc);

  File file("edited_result.pvd");
  file << u;
  return 0;
}