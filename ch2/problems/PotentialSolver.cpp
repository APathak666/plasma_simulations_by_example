#include "PotentialSolver.hpp"
#include <cmath>
#include <fstream>

bool PotentialSolver::solve() {
    Field &phi = world.phi;
    Field &rho = world.rho;

    double3 dh = world.getDh();
    double idx2 = 1/(dh[0]*dh[0]);
    double idy2 = 1/(dh[1]*dh[1]);
    double idz2 = 1/(dh[2]*dh[2]);

    /* Dirichlet BCs on z-faces */
    for (int i = 0; i < world.ni; i++)
        for (int j = 0; j < world.nj; j++) {
            phi[i][j][0]           = 0.0;   // z_min
            phi[i][j][world.nk-1]  = 10.0;  // z_max
        }

    bool converged = false;

    std::ofstream conv_out("convergence.csv");
    conv_out << "iter,residual\n";

    Field phi_new(world.ni, world.nj, world.nk);

    for (int it  = 0; it < max_solver_it; it++) {

        for (int i = 0; i < world.ni; i++)
            for (int j = 0; j < world.nj; j++)
                for (int k = 1; k < world.nk-1; k++) {
                    phi_new[i][j][k] = 0.0;

                    if (i == 0) phi_new[i][j][k] += idx2*2*phi[i+1][j][k];
                    else if (i == world.ni-1) phi_new[i][j][k] += idx2*2*phi[i-1][j][k];
                    else phi_new[i][j][k] += idx2*(phi[i-1][j][k]+phi[i+1][j][k]);

                    if (j == 0) phi_new[i][j][k] += idy2*2*phi[i][j+1][k];
                    else if (j == world.nj-1) phi_new[i][j][k] += idy2*2*phi[i][j-1][k];
                    else phi_new[i][j][k] += idy2*(phi[i][j-1][k]+phi[i][j+1][k]);

                    phi_new[i][j][k] += idz2*(phi[i][j][k-1]+phi[i][j][k+1]);
                    phi_new[i][j][k] += rho[i][j][k]/Const::EPS_0;
                    phi_new[i][j][k] /= (2*idx2+2*idy2+2*idz2);

                    /*SOR*/
                    // phi[i][j][k] = phi[i][j][k] + 1.4*(phi_new-phi[i][j][k]);
                }

        for (int i = 0; i < world.ni; i++)
            for (int j = 0; j < world.nj; j++)
                for (int k = 1; k < world.nk-1; k++)
                    phi[i][j][k] = phi_new[i][j][k];

        if (!(it%25)) {
            double R_sum = 0;

            for (int i = 0; i < world.ni; i++)
                for (int j = 0; j < world.nj; j++)
                    for (int k = 1; k < world.nk-1; k++) {
                        double R_ijk = -phi[i][j][k]*(2*idx2+2*idy2+2*idz2);

                        if (i == 0) R_ijk += idx2*2*phi[i+1][j][k];
                        else if (i == world.ni-1) R_ijk += idx2*2*phi[i-1][j][k];
                        else R_ijk += idx2*(phi[i-1][j][k]+phi[i+1][j][k]);

                        if (j == 0) R_ijk += idy2*2*phi[i][j+1][k];
                        else if (j == world.nj-1) R_ijk += idy2*2*phi[i][j-1][k];
                        else R_ijk += idy2*(phi[i][j-1][k]+phi[i][j+1][k]);

                        R_ijk += idz2*(phi[i][j][k-1]+phi[i][j][k+1]);
                        R_ijk += rho[i][j][k]/Const::EPS_0;

                        R_sum += R_ijk*R_ijk;
                    }

            R_sum /= world.ni*world.nj*world.nk;
            conv_out << it << "," << sqrt(R_sum) << "\n";

			if (sqrt(R_sum) < tolerance) {
				std::cout << "GS solver converged after " << it << " iterations\n";
                converged = true;
                break;
            }
        }
    }

    if (!converged) std::cerr << "GS solver failed to converge\n";
    return converged;
}

void PotentialSolver::computeEF() {
    Field &phi = world.phi;
    double3 dh = world.getDh();
    double idx = 1/(2*dh[0]);
    double idy = 1/(2*dh[1]);
    double idz = 1/(2*dh[2]);

    for (int i = 0; i < world.ni; i++)
        for (int j = 0; j < world.nj; j++)
            for (int k = 0; k < world.nk; k++) {
                double3 &ef = world.ef[i][j][k]; //reference to (i,j,k) ef vec3

				/*x component*/
				if (i==0)	/*forward*/
					// ef[0] = -(-3*phi[i][j][k]+4*phi[i+1][j][k]-phi[i+2][j][k])*idx;
                    ef[0] = 0;
				else if (i==world.ni-1)  /*backward*/
					// ef[0] = -(phi[i-2][j][k]-4*phi[i-1][j][k]+3*phi[i][j][k])*idx;
                    ef[0] = 0;
                else  /*central*/
					ef[0] = -(phi[i+1][j][k] - phi[i-1][j][k])*idx;

				/*y component*/
				if (j==0)
					// ef[1] = -(-3*phi[i][j][k] + 4*phi[i][j+1][k]-phi[i][j+2][k])*idy;
                    ef[1] = 0;
                else if (j==world.nj-1)
					// ef[1] = -(phi[i][j-2][k] - 4*phi[i][j-1][k] + 3*phi[i][j][k])*idy;
                    ef[1] = 0;
				else
					ef[1] = -(phi[i][j+1][k] - phi[i][j-1][k])*idy;

				/*z component*/
				if (k==0)
					ef[2] = -(-3*phi[i][j][k] + 4*phi[i][j][k+1]-phi[i][j][k+2])*idz;
                    // ef[2] = 0;
                else if (k==world.nk-1)
					ef[2] = -(phi[i][j][k-2] - 4*phi[i][j][k-1]+3*phi[i][j][k])*idz;
                    // ef[2] = 0;
				else
					ef[2] = -(phi[i][j][k+1] - phi[i][j][k-1])*idz;
            }
}

void PotentialSolver::logPotential() {
    std::ofstream out("potential.csv");
    if (!out) { std::cerr << "Could not open potential.csv\n"; return; }
    out << "x,y,z,phi,rho\n";

    double3 x0 = world.getX0();
    double3 dh = world.getDh();

    for (int i = 0; i < world.ni; i++)
        for (int j = 0; j < world.nj; j++)
            for (int k = 0; k < world.nk; k++) {
                double x = x0[0] + i*dh[0];
                double y = x0[1] + j*dh[1];
                double z = x0[2] + k*dh[2];
                out << x << "," << y << "," << z << ","
                    << world.phi[i][j][k] << "," << world.rho[i][j][k] << "\n";
            }
}

void PotentialSolver::logElectricField() {
    std::ofstream out("ef.csv");
    if (!out) { std::cerr << "Could not open ef.csv\n"; return; }
    out << "x,y,z,ef_x,ef_y,ef_z\n";

    double3 x0 = world.getX0();
    double3 dh = world.getDh();

    for (int i = 0; i < world.ni; i++)
        for (int j = 0; j < world.nj; j++)
            for (int k = 0; k < world.nk; k++) {
                double x = x0[0] + i*dh[0];
                double y = x0[1] + j*dh[1];
                double z = x0[2] + k*dh[2];
                out << x << "," << y << "," << z << ","
                    << world.ef[i][j][k][0] << "," << world.ef[i][j][k][1] << ","
                    << world.ef[i][j][k][2] << "\n";
            }
}
