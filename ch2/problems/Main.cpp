#include "World.hpp"
#include "Output.hpp"
#include "PotentialSolver.hpp"
#include "Species.hpp"
#include "Config.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
    std::string config_path = (argc > 1) ? argv[1] : "config.txt";
    Config cfg = Config::load(config_path);

    double3 origin = {cfg.origin_x, cfg.origin_y, cfg.origin_z};
    double3 xm = {origin[0] + cfg.xcount*cfg.dx,
                  origin[1] + cfg.ycount*cfg.dy,
                  origin[2] + cfg.zcount*cfg.dz};

    World world(cfg.xcount+1, cfg.ycount+1, cfg.zcount+1);
    world.setExtents(origin, xm);
    world.setTime(cfg.dt, cfg.num_ts);

    std::vector<Species> species;
    species.reserve(2);
    species.push_back(Species("O+", 16*Const::AMU, Const::QE, world));
    species.push_back(Species("e-", Const::ME, -1*Const::QE, world));

    int3 np_eles_grid = {cfg.xcount+1, cfg.ycount+1, cfg.zcount+1};
    int3 np_ions_grid = {2*cfg.xcount+1, 2*cfg.ycount+1, 2*cfg.zcount+1};

    species[0].loadParticlesBoxQS(world.getX0(), world.getXm(), cfg.ni, np_ions_grid);
    species[1].loadParticlesBoxQS(world.getX0(), world.getXc(), cfg.ne, np_eles_grid);

    for (Species &sp: species) {
        sp.advance();
        sp.computeNumberDensity();
    }

    world.computeRho(species);


    PotentialSolver solver(world, 10000, 1e-6);
    solver.solve();
    solver.computeEF();
    solver.logPotential();
    solver.logElectricField();

    while (world.advanceTime()) {
        for (Species &sp: species) {
            sp.advance();
            sp.computeNumberDensity();
        }

        world.computeRho(species);

        solver.solve();
        solver.computeEF();

        Output::screenOutput(world,species);
        Output::diagOutput(world,species);

		/*periodically write out results*/
        if (world.getTs()%100==0 || world.isLastTimeStep())
			Output::fields(world, species);
    }

	std::cout<<"Simulation took "<<world.getWallTime()<<" seconds";
}
