#include "World.hpp"
#include "Output.hpp"
#include "PotentialSolver.hpp"
#include "Species.hpp"

int main(int argc, char* argv[]) {
    World world(21, 21, 41);
    world.setExtents({-0.1, -0.1, 0.0}, {0.1, 0.1, 0.4});
    world.setTime(2e-10, 10000);

    std::vector<Species> species;
    species.reserve(2);
    species.push_back(Species("O+", 16*Const::AMU, Const::QE, world));
    species.push_back(Species("e-", Const::ME, -1*Const::QE, world));

    double phi_sphere = -100;
    if (argc > 1) phi_sphere = atof(argv[1]);
    world.addSphere({0, 0, 0.15}, 0.05, phi_sphere);
    world.addInlet();

    // int3 np_eles_grid = {21, 21, 41};
    // int3 np_ions_grid = {41, 41, 81};

    // species[0].loadParticlesBoxQS(world.getX0(), world.getXm(), 1e11, np_ions_grid);
    // species[1].loadParticlesBoxQS(world.getX0(), world.getXc(), 1e11, np_eles_grid);

    for (Species &sp: species) {
        sp.advance();
        sp.computeNumberDensity();
    }

    world.computeRho(species);


    PotentialSolver solver(world, 10000, 1e-6);
    solver.solve();
    solver.computeEF();

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
