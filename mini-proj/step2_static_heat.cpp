#include "common.h"
#include "Heat.h"

int main() {
    auto mesh = std::make_shared<Mesh>("torus.xml");
    auto V = std::make_shared<Heat::FunctionSpace>(mesh);

    double source_angle = 0.0; 
    double ring_offset = M_PI; 

    auto boundary = std::make_shared<CoolingRings>(std::vector<double>{source_angle + ring_offset}, 0.4);
    auto bc = std::make_shared<DirichletBC>(V, std::make_shared<Constant>(0.0), boundary);

    const double R = 1.5;
    auto f = std::make_shared<HeatSource>();
    f->xc = R * std::cos(source_angle);
    f->yc = R * std::sin(source_angle);
    f->A = 1500.0; 

    auto a = std::make_shared<Heat::BilinearForm>(V, V);
    auto L = std::make_shared<Heat::LinearForm>(V);
    a->k = std::make_shared<Constant>(1.0);
    a->dt = std::make_shared<Constant>(1e8);
    L->dt = std::make_shared<Constant>(1e8);
    L->u0 = std::make_shared<Constant>(0.0);
    L->f = f;

    Function u(V);
    solve(*a == *L, u, *bc);
    write_vtu(mesh, u, "static_source.vtu");
    return 0;
}
