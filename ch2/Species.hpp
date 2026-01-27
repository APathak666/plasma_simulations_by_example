#ifndef _SPECIES_H
#define _SPECIES_H

#include "World.hpp"
#include "Field.hpp"
#include <string>
#include <vector>

struct Particle {
    double3 pos;
    double3 vel;
    double mpw;

    Particle(double3 x, double3 v, double mpw) : pos(x), vel(v), mpw(mpw) {}
};

class Species {
    public:
        Species(std::string name, double mass, double charge, World &world) :
            name(name), mass(mass), charge(charge), den(world.ni, world.nj, world.nk), world(world) {}
        size_t getNp() const { return particles.size(); }
        void advance();
        void computeNumberDensity();
        void addParticle(double3 pos, double3 vel, double mpw);
        void loadParticlesBoxQS(double3 x1, double3 x2, double number_den, int3 num_sim);
        double getRealCount();
        double3 getMomentum();
        double getKE();

        const std::string name;
        const double mass;
        const double charge;
        Field den;
        std::vector<Particle> particles;

    protected:
        World &world;
};

#endif
