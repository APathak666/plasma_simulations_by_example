#include "World.hpp"
#include "Output.hpp"
#include "PotentialSolver.hpp"
#include "Species.hpp"

int main() {
    /* ================= problem configuration ================= */
    int xcount = 20, ycount = 20, zcount = 20;          // cell counts per axis (nodes = cells + 1)
    double dx = 0.01, dy = 0.01, dz = 0.01;             // grid spacing (cell size)
    double3 origin = {-0.1, -0.1, 0.0};                 // domain origin (x0)
    double n0 = 0;                                    // reference number density (1/m^3)
    double dt = 2e-10;                                   // time step
    int num_ts = 1;                                  // number of time steps
    /* ========================================================= */

    double3 xm = {origin[0] + xcount*dx,
                  origin[1] + ycount*dy,
                  origin[2] + zcount*dz};

    World world(xcount+1, ycount+1, zcount+1);
    world.setExtents(origin, xm);
    world.setTime(dt, num_ts);

    std::vector<Species> species;
    species.reserve(2);
    species.push_back(Species("O+", 16*Const::AMU, Const::QE, world));
    species.push_back(Species("e-", Const::ME, -1*Const::QE, world));

    int3 np_eles_grid = {xcount+1, ycount+1, zcount+1};
    int3 np_ions_grid = {2*xcount+1, 2*ycount+1, 2*zcount+1};

    species[0].loadParticlesBoxQS(world.getX0(), world.getXm(), n0, np_ions_grid);
    species[1].loadParticlesBoxQS(world.getX0(), world.getXc(), n0, np_eles_grid);

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
