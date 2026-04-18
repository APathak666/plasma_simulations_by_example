#ifndef _WORLD_H
#define _WORLD_H

// #include "Species.hpp"
#include "Field.hpp"
#include <vector>
#include <chrono>

class Species;

namespace Const
{
	const double EPS_0 = 8.85418782e-12;  	// C/(V*m), vacuum permittivity
	const double QE = 1.602176565e-19;		// C, electron charge
	const double AMU = 1.660538921e-27;		// kg, atomic mass unit
	const double ME = 9.10938215e-31;		// kg, electron mass
	const double K = 1.380648e-23;			// J/K, Boltzmann constant
	const double PI = 3.141592653;			// pi
	const double EvToK = QE/K;				// 1eV in K ~ 11604
}

class World {
    public:
        World(int ni, int nj, int nk);
        void setExtents(const double3 &_x0, const double3 _xm);  /* set mesh span, recompute spacing */
        double3 getX0() const {return double3(x0);}
        double3 getXm() const {return double3(xm);}
        double3 getXc() const {return double3(xc);}
        double3 getDh() const {return double3(dh);}
        void computeNodeVolumes();
        void computeRho(std::vector<Species> &species);
        double3 XtoL(double3 x) const;

    	/*functions for accessing time information*/
        int getTs() const {return ts;}
        double getTime() const {return time;}
        double getWallTime();  /*returns wall time in seconds*/
        double getDt() const {return dt;}
        bool isLastTimeStep() const {return ts==num_ts-1;}
        /*sets time step and number of time steps*/
        void setTime(double dt, int num_ts) {this->dt=dt;this->num_ts=num_ts;}        
        /*advances to the next time step, returns true as long as more time steps remain*/
        bool advanceTime() {time+=dt;ts++;return ts<=num_ts;}
        bool inBounds(double3 pos) {
            for (int i=0;i<3;i++)
                if (pos[i]<x0[i] || pos[i]>=xm[i]) return false;
            return true;
        }
        double getPE();

        const int3 nn;        /* number of nodes */
        const int ni, nj, nk;   /* number of nodes in individual variables */
        Field phi;
        Field rho;
        Field node_vol;
        Field3 ef;

    protected:
        double3 x0;   /* origin vector */
        double3 xm;   /* diagonally opposite end vector of cube */
        double3 dh;   /* element spacing */
        double3 xc;   /* centroid vector */
        double dt = 0;		//time step length
        double time = 0;	//physical time
        int ts = -1;		//current time step
        int num_ts = 0;		//number of time steps
        std::chrono::time_point<std::chrono::high_resolution_clock> time_start;
};

#endif
