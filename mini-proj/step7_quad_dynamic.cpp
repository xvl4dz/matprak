#include "common.h"
#include "Heat.h"

int main() {
    auto mesh = std::make_shared<Mesh>("torus.xml");
    auto V = std::make_shared<Heat::FunctionSpace>(mesh);
    auto u0 = std::make_shared<Function>(V); *u0 = Constant(0.0);
    
    double dt = 0.1, T = 60.0, R = 1.5;
    double omega = (2.0 * M_PI) / 20.0; 

    auto a = std::make_shared<Heat::BilinearForm>(V, V);
    auto L = std::make_shared<Heat::LinearForm>(V);
    a->k = std::make_shared<Constant>(1.0);
    a->dt = std::make_shared<Constant>(dt);
    L->dt = std::make_shared<Constant>(dt);

    auto f = std::make_shared<HeatSource>();
    f->s = 0.15;

    std::vector<double> static_angles = {M_PI/4, 3*M_PI/4, 5*M_PI/4, 7*M_PI/4};
    auto rings = std::make_shared<CoolingRings>(static_angles, 0.25);
    DirichletBC bc(V, std::make_shared<Constant>(0.0), rings);

    for (double t = 0; t <= T + 1e-8; t += dt) {
        double current_phi = omega * t;
        f->xc = R * std::cos(current_phi);
        f->yc = R * std::sin(current_phi);

        double ramp = std::min(1.0, t / 1.0);
        f->A = 1500.0 * ramp;

        L->u0 = u0; L->f = f;
        Function u(V);
        solve(*a == *L, u, bc);
        *u0 = u;
        
        if (std::abs(std::fmod(t, 1.0)) < dt/2.0) {
            write_vtu(mesh, u, "step7_dynamic_" + std::to_string(int(t)) + ".vtu");
        }
    }
    return 0;
}
