#include "PotentialSolver.hpp"
#include <cmath>

bool PotentialSolver::solve() {
    Field &phi = world.phi;
    Field &rho = world.rho;

    double3 dh = world.getDh();
    double idx2 = 1/(dh[0]*dh[0]);
    double idy2 = 1/(dh[1]*dh[1]);
    double idz2 = 1/(dh[2]*dh[2]);

    bool converged = false;

    for (int it  = 0; it < max_solver_it; it++) {
        for (int i = 1; i < world.ni-1; i++)
            for (int j = 1; j < world.nj-1; j++)
                for (int k = 1; k < world.nk-1; k++) {
                    double phi_new = (rho[i][j][k]/Const::EPS_0 +
                                    idx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
                                    idy2*(phi[i][j-1][k]+phi[i][j+1][k]) +
                                    idz2*(phi[i][j][k-1]+phi[i][j][k+1]))/(2*idx2+2*idy2+2*idz2);

                    /*SOR*/
                    phi[i][j][k] = phi[i][j][k] + 1.4*(phi_new-phi[i][j][k]);
                }

        if (!(it%25)) {
            double R_sum = 0;

            for (int i = 1; i < world.ni-1; i++)
                for (int j = 1; j < world.nj-1; j++)
                    for (int k = 1; k < world.nk-1; k++) {
                        double R_ijk = -phi[i][j][k]*(2*idx2+2*idy2+2*idz2) +
									rho[i][j][k]/Const::EPS_0 +
									idx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
									idy2*(phi[i][j-1][k]+phi[i][j+1][k]) +
									idz2*(phi[i][j][k-1]+phi[i][j][k+1]);

                        R_sum += R_ijk*R_ijk;
                    }

            R_sum /= world.ni*world.nj*world.nk;

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
					ef[0] = -(-3*phi[i][j][k]+4*phi[i+1][j][k]-phi[i+2][j][k])*idx;	
				else if (i==world.ni-1)  /*backward*/
					ef[0] = -(phi[i-2][j][k]-4*phi[i-1][j][k]+3*phi[i][j][k])*idx;	
				else  /*central*/
					ef[0] = -(phi[i+1][j][k] - phi[i-1][j][k])*idx;	

				/*y component*/
				if (j==0)
					ef[1] = -(-3*phi[i][j][k] + 4*phi[i][j+1][k]-phi[i][j+2][k])*idy;
				else if (j==world.nj-1)
					ef[1] = -(phi[i][j-2][k] - 4*phi[i][j-1][k] + 3*phi[i][j][k])*idy;
				else
					ef[1] = -(phi[i][j+1][k] - phi[i][j-1][k])*idy;

				/*z component*/
				if (k==0)
					ef[2] = -(-3*phi[i][j][k] + 4*phi[i][j][k+1]-phi[i][j][k+2])*idz;
				else if (k==world.nk-1)
					ef[2] = -(phi[i][j][k-2] - 4*phi[i][j][k-1]+3*phi[i][j][k])*idz;
				else
					ef[2] = -(phi[i][j][k+1] - phi[i][j][k-1])*idz;
            }
}
