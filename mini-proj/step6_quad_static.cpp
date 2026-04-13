#include "common.h"
#include "Heat.h"
#include <cmath>

int main() {
    auto mesh = std::make_shared<Mesh>("torus.xml");
    auto V = std::make_shared<Heat::FunctionSpace>(mesh);

    auto f = std::make_shared<HeatSource>();
    f->xc = 1.5; f->yc = 0.0; f->zc = 0.0;
    f->A = 1500.0; f->s = 0.15;
    
    std::vector<double> angles = {M_PI/4, 3*M_PI/4, 5*M_PI/4, 7*M_PI/4};
    auto rings = std::make_shared<CoolingRings>(angles);
    DirichletBC bc(V, std::make_shared<Constant>(0.0), rings);

    auto a = std::make_shared<Heat::BilinearForm>(V, V);
    auto L = std::make_shared<Heat::LinearForm>(V);
    a->k = std::make_shared<Constant>(1.0);
    a->dt = std::make_shared<Constant>(1e6);
    L->dt = std::make_shared<Constant>(1e6);
    L->u0 = std::make_shared<Constant>(0.0);
    L->f = f;

    Function u(V);
    solve(*a == *L, u, bc);
    write_vtu(mesh, u, "step6_static_4rings.vtu");
    return 0;
}
