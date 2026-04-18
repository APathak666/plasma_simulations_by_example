#include "Species.hpp"

void Species::addParticle(double3 pos, double3 vel, double mpw) {
    double dt = world.getDt();
    double3 lc = world.XtoL(pos);
    double3 ef_i = world.ef.gather(lc);
    vel -= (0.5*charge*dt/mass)*ef_i;

    particles.emplace_back(pos, vel, mpw);
}

void Species::computeNumberDensity() {
    den = 0;
    for (Particle &particle: particles) {
        double3 lc = world.XtoL(particle.pos);
        den.scatter(lc, particle.mpw);
    }

    den /= world.node_vol;
}

void Species::loadParticlesBoxQS(double3 x1, double3 x2, double number_den, int3 num_sim) {
    double3 dh;
    for (int i = 0; i < 3; i++) dh[i] = (x2[i]-x1[i])/(num_sim[i]-1);

    double mpw = number_den*(x2[0]-x1[0])*(x2[1]-x1[1])*(x2[2]-x1[2]);
    mpw /= (num_sim[0]-1)*(num_sim[1]-1)*(num_sim[2]-1);  /* this is to account for half-weighing of particles at the edge */

    for (int i = 0; i < num_sim[0]; i++)
        for (int j = 0; j < num_sim[1]; j++)
            for (int k = 0; k < num_sim[2]; k++) {
                double3 pos;
                pos[0] = x1[0] + i*dh[0];
                pos[1] = x1[1] + j*dh[1];
                pos[2] = x1[2] + k*dh[2];

                if (pos[0] == x2[0]) pos[0] -= (1e-4)*dh[0];
                if (pos[1] == x2[1]) pos[1] -= (1e-4)*dh[1];
                if (pos[2] == x2[2]) pos[2] -= (1e-4)*dh[2];

                double w = 1;
                if (i == 0 || i == num_sim[0]-1) w*= 0.5;
                if (j == 0 || j == num_sim[1]-1) w*= 0.5;
                if (k == 0 || k == num_sim[2]-1) w*= 0.5;

                double3 vel;
                vel = 0;
                addParticle(pos, vel, mpw*w);
            }
}

void Species::advance() {
    double dt = world.getDt();
    double3 x0 = world.getX0();
    double3 xm = world.getXm();

    for (Particle &particle: particles) {
        double3 lc = world.XtoL(particle.pos);
        double3 ef_i = world.ef.gather(lc);
        particle.vel += (charge*dt/mass)*ef_i;
        particle.pos += particle.vel*dt;

        for (int i = 0; i < 3; i++) {
            if (particle.pos[i] < x0[i]) {
                particle.pos[i] = 2*x0[i] - particle.pos[i];
                particle.vel[i] = -1.0*particle.vel[i];
            }
            else if (particle.pos[i] > xm[i]) {
                particle.pos[i] = 2*xm[i] - particle.pos[i];
                particle.vel[i] = -1.0*particle.vel[i];
            }
        }
    }
}

double Species::getRealCount() {
    double real_count = 0;
    for (Particle &particle: particles) real_count += particle.mpw;
    return real_count;
}

double3 Species::getMomentum() {
    double3 mom;
    mom = 0;

    for (Particle &particle: particles) mom += particle.mpw*particle.vel;
    return mom*mass;
}

double Species::getKE() {
    double ke = 0;
    for (Particle &particle: particles)
        ke += particle.mpw*(particle.vel[0]*particle.vel[0] + particle.vel[1]*particle.vel[1] + particle.vel[2]*particle.vel[2]);

    return 0.5*ke*mass;
}
