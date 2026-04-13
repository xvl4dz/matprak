#include "common.h"
#include "Heat.h"
#include <cmath>

int main() {
    auto mesh = std::make_shared<Mesh>("torus.xml");
    auto V = std::make_shared<Heat::FunctionSpace>(mesh);
    auto u0 = std::make_shared<Function>(V); *u0 = Constant(0.0);
    
    double dt = 0.05, T = 50.0, R = 1.5;
    auto a = std::make_shared<Heat::BilinearForm>(V, V);
    auto L = std::make_shared<Heat::LinearForm>(V);
    a->k = std::make_shared<Constant>(1.0);
    a->dt = std::make_shared<Constant>(dt);
    L->dt = std::make_shared<Constant>(dt);

    auto f = std::make_shared<HeatSource>();
    f->xc = R; f->yc = 0; f->s = 0.15;

    auto rings = std::make_shared<CoolingRings>(std::vector<double>{M_PI});
    DirichletBC bc(V, std::make_shared<Constant>(0.0), rings);

    for (double t = 0; t <= T + 1e-8; t += dt) {
        f->A = (std::fmod(t, 10.0) < 5.0) ? 1500.0 : 0.0;
        
        L->u0 = u0; L->f = f;
        Function u(V);
        solve(*a == *L, u, bc);
        *u0 = u;
        
        int step = static_cast<int>(std::round(t * 100));
        write_vtu(mesh, u, "step4_strobe_" + std::to_string(step) + ".vtu");
    }
    return 0;
}
