#include "World.hpp"
#include "Species.hpp"

World::World(int ni, int nj, int nk) : ni{ni}, nj{nj}, nk{nk}, nn{ni, nj, nk},
                                        phi(ni, nj, nk), rho(ni, nj, nk), node_vol(ni, nj, nk),
                                        ef(ni, nj, nk) {
    time_start =  std::chrono::high_resolution_clock::now();	//save starting time point
}

void World::setExtents(const double3 &_x0, const double3 _xm) {
    x0 = _x0;
    xm = _xm;
    xc = 0.5*(x0+xm);

    for (int i = 0; i < 3; i++) {
        dh[i] = (xm(i) - x0(i))/(nn(i) - 1);
    }

    computeNodeVolumes();
}

/*returns elapsed wall time in seconds*/
double World::getWallTime() {
  auto time_now = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_delta = time_now-time_start;
  return time_delta.count();
}

void World::computeNodeVolumes() {

    double vol = dh[0]*dh[1]*dh[2];

    for (int i = 0; i < ni; i++)
        for (int j = 0; j < nj; j++)
            for (int k = 0; k < nk; k++) {
                node_vol[i][j][k] = vol;
                if (i == 0 || i == ni-1) node_vol[i][j][k] /= 2;
                if (j == 0 || j == nj-1) node_vol[i][j][k] /= 2;
                if (k == 0 || k == nk-1) node_vol[i][j][k] /= 2;
            }
}

void World::computeRho(std::vector<Species> &species) {
    rho = 0;
    for (Species &sp: species) {
        if (!sp.charge) continue;
        rho += sp.charge*sp.den;
    }
}

double3 World::XtoL(double3 x) const {
    double3 lc;
    for (int i = 0; i < 3; i++) lc[i] = (x[i]-x0(i))/dh(i);

    return lc;
}

double World::getPE() {
    double pe = 0;

    for (int i = 0; i < ni; i++)
        for (int j = 0; j < nj; j++)
            for (int k = 0; k < nk; k++)
                pe += (ef[i][j][k][0]*ef[i][j][k][0] + ef[i][j][k][1]*ef[i][j][k][1] + ef[i][j][k][2]*ef[i][j][k][2])*node_vol[i][j][k];


    return 0.5*Const::EPS_0*pe;
}
