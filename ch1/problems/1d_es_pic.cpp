#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

namespace Const {
    const double Q_E = 1.602176565e-19; 	/* single positive charge */
    const double N_POSITIVE = 1e12;   		/* number density of +ve ions, m^(-3) */
	const double EPS_0 = 8.85418782e-12;	/* epsilon0 */
	const double SOR_W = 1.4;				/* weight of new value in SOR */
	const int N_I = 21;						/* number of nodes */
};

using namespace std;
using namespace Const;
using dvector = vector<double>;

void logPotential(double x0, double dx, const dvector &phi, const dvector &rho);
void logElectricField(double x0, double dx, const dvector &ef, const dvector &rho);
bool solvePotentialGS(double dx, dvector &phi, const dvector &rho, int max_it=10000);
void computeEF(double dx, dvector &ef, const dvector &phi);

int main() {
    const double x0 = 0;
    const double xm = 0.1;
    double dx = (xm-x0)/(N_I-1);

    dvector phi(N_I);
    dvector rho(N_I, Q_E*N_POSITIVE);   /* total charge density = number density*single charge = Q_E*n_i */
    dvector ef(N_I, 0);

	solvePotentialGS(dx, phi, rho);
	logPotential(x0, dx, phi, rho);
	computeEF(dx, ef, phi);
	logElectricField(x0, dx, ef, rho);

    return 0;
}

void logPotential(double x0, double dx, const dvector &phi, const dvector &rho) {
	ofstream out("potential.csv");
	if (!out) { cerr << "Could not open potential.csv\n"; return; }
	out << "x,phi,phi_a\n";
	for (int i = 0; i < N_I; i++) {
		double x = x0 + i*dx;
		out << x << "," << phi[i] << "," << rho[i]*(0.01 - (x*x))/(2*EPS_0) << "\n";
	}
}

void logElectricField(double x0, double dx, const dvector &ef, const dvector &rho) {
	ofstream out("ef.csv");
	if (!out) { cerr << "Could not open ef.csv\n"; return; }
	out << "x,ef,ef_a\n";
	for (int i = 0; i < N_I; i++) {
		double x = x0 + i*dx;
		out << x << "," << ef[i] << "," << rho[i]*x/EPS_0 << "\n";
	}
}

bool solvePotentialGS(double dx, dvector &phi, const dvector &rho, int max_it) {
	int ni = phi.size();
	double dx2 = dx*dx;
	double R_sum;

	dvector phi_new(ni);
	for (int i = 0; i < ni; i++) {
		phi_new[i] = 0;
		phi[i] = 0;
	}

	// phi[0] = 0;
	phi[ni-1] = 0;

	ofstream conv_out("convergence.csv");
	ofstream phi_out("phi_history.csv");
	conv_out << "iter,residual\n";
	phi_out << "iter";
	for (int i = 0; i < ni; i++) phi_out << ",x" << i;
	phi_out << "\n";

	for (int k = 0; k < max_it; k++) {

		// phi[0] = phi[1];
		for (int i = 1; i < ni - 1; i++) {
			phi_new[i] = (phi[i-1] + phi[i+1] + rho[i]*dx2/EPS_0)/2;	/* Gauss-Siedel scheme, 1.52 */
			// phi[i] = phi[i] + SOR_W*(phi_new[i] - phi[i]);					/* Successive over-relaxation */
		}
		// phi_new[0] = (4*phi[1] - phi[2])/3;
		phi_new[0] = phi[1] + (rho[0]*dx2)/(2*EPS_0);

		for (int i = 0; i < ni; i++) phi[i] = phi_new[i];

		/* compute residual every iteration for logging */
		R_sum = 0;
		for (int i = 1; i < ni - 1; i++) {
			double Ri = (phi[i-1] - 2*phi[i] + phi[i+1])/dx2 + (rho[i]/EPS_0);
			R_sum += Ri*Ri;
		}
		R_sum = sqrt(R_sum / ni);
		conv_out << k << "," << R_sum << "\n";

		/* log phi every 10 iterations */
		if (!(k % 10)) {
			phi_out << k;
			for (int i = 0; i < ni; i++) phi_out << "," << phi[i];
			phi_out << "\n";
		}

		if (R_sum < 1e-6) {
			cout << "GS solver converged after " << k << " iterations\n";
			return true;
		}
	}

	cout << "GS solver failed to converge\n";
	return false;
}

void computeEF(double dx, dvector &ef, const dvector &phi) {
	int ni = ef.size();

	for (int i = 1; i < ni - 1; i++)
		ef[i] = (phi[i-1] - phi[i+1])/(2*dx);

	/* Second-order difference on boundaries for better accuracy */
	// ef[0] = (3*phi[0] - 4*phi[1] + phi[2])/(2*dx);
	ef[0] = 0;
	ef[ni-1] = (-phi[ni-3] + 4*phi[ni-2] - 3*phi[ni-1])/(2*dx);
}
