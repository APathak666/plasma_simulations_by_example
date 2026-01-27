#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

namespace Const {
    const double Q_E = 1.602176565e-19; 	/* single positive charge */
    const double N_POSITIVE = 1e12;   		/* number density of +ve ions, m^(-3) */
	const double EPS_0 = 8.85418782e-12;	/* epsilon0 */
	const double M_ELEC = 9.10938215e-31;	/* mass of electron in kg */
	const double SOR_W = 1.4;				/* weight of new value in SOR */
	const int N_I = 21;						/* number of nodes */
	const int MAX_ITER = 4000;				/* max number of timesteps for sim */
};

using namespace std;
using namespace Const;
using dvector = vector<double>;

bool outputCSV(double x0, double dx, const dvector &phi, const dvector &rho, const dvector &ef);
bool solvePotentialGS(double dx, dvector &phi, const dvector &rho, int max_it=5000);
void computeEF(double dx, dvector &ef, const dvector &phi);
double XtoL(double x, double x0, double dx) { return (x-x0)/dx; }
double gather(double li, const dvector &field);					/* general function, used to gather any field (ef, phi, etc) */

int main() {
    const double x0 = 0;
    const double xm = 0.1;
    double dx = (xm-x0)/(N_I-1);

    dvector phi(N_I);
    dvector rho(N_I, Q_E*N_POSITIVE);   /* total charge density = number density*single charge = Q_E*n_i */
    dvector ef(N_I, 0);

	solvePotentialGS(dx, phi, rho);
	computeEF(dx, ef, phi);

	double m = M_ELEC;
	double q = -Q_E;
	double x = 4*dx;	/* initial position */
	double v = 0;

	double dt = 1e-10;

	double li = XtoL(x, x0, dx);
	double ef_i = gather(li, ef);
	v -= dt*q*ef_i/(2*m);			/* velocity rewind */

	double phi_max = phi[0];
	for (int i = 0; i < N_I; i++) if (phi[i] > phi_max) phi_max = phi[i];

	ofstream out("trace.csv");
	if (!out) { cerr << "Failed to open trace\n"; return -1; }
	out << "time,x,v,KE,PE\n";

	double phi_p, ke, pe;
	for (int ts = 0; ts < MAX_ITER; ts++) {
		double li = XtoL(x, x0, dx);
		double ef_i = gather(li, ef);
		v += dt*q*ef_i/m;
		x += v*dt;

		phi_p = gather(XtoL(x, x0, dx), phi);
		ke = 0.5*m*v*v/Q_E;
		pe = q*(phi_p - phi_max)/Q_E;
		out <<ts*dt << "," << x << "," << v << "," << ke << "," << pe << "\n";

		if (!(ts%1000)) cout << "ts: " << ts << "\tx: " << x << "\tv: " << v << "\tKE: " << ke << "\tPE: " << pe << "\n";
	}

	cout << "ts: " << MAX_ITER << "\tx: " << x << "\tv: " << v << "\tKE: " << ke << "\tPE: " << pe << "\n";

	outputCSV(x0, dx, phi, rho, ef);

    return 0;
}

bool outputCSV(double x0, double dx, const dvector &phi, const dvector &rho, const dvector &ef) {
	ofstream out("results.txt");	/* open file for writing */
	if (!out) {
		cerr<<"Could not open output file!"<<endl; 
		return false;
	}
	
	out<<"x\tphi\trho\t\tef\n";		/* write header */
	for (int i=0;i<phi.size();i++) {
		out<<x0+i*dx; /* write i-th position */
		out<<"\t"<<phi[i]<<"\t"<< rho[i]<<"\t"<<ef[i]; /* write values */
		out<<"\n";	/* new line, not using endl to avoid buffer flush */
	}
	
	/* file closed automatically when "out" variable is destroyed */
	return true;
}

bool solvePotentialGS(double dx, dvector &phi, const dvector &rho, int max_it) {
	int ni = phi.size();
	double dx2 = dx*dx;
	double R_sum;

	phi[0] = 0;
	phi[ni-1] = 0;

	for (int k = 0; k < max_it; k++) {
		for (int i = 1; i < ni - 1; i++) {
			double phi_i_new = (phi[i-1] + phi[i+1] + rho[i]*dx2/EPS_0)/2;	/* Gauss-Siedel scheme, 1.52 */
			phi[i] = phi[i] + SOR_W*(phi_i_new - phi[i]);					/* Successive over-relaxation */
		}

		/* convergence check every 50 iter (measure norm of error matrix: 1.54) */
		if (!(k%50)) {
			R_sum = 0;
			for (int i = 1; i < ni - 1; i++) {
				double Ri = (phi[i-1] - 2*phi[i] + phi[i+1])/dx2 + (rho[i]/EPS_0);
				R_sum += Ri*Ri;
			}

			R_sum /= ni;
			// cout << R_sum << endl;

			if (sqrt(R_sum) < 1e-6) {
				cout << "GS solver converged after " << k << " iterations\n";
				return true;
			}
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
	ef[0] = (3*phi[0] - 4*phi[1] + phi[2])/(2*dx);
	ef[ni-1] = (-phi[ni-3] + 4*phi[ni-2] - 3*phi[ni-1])/(2*dx);
}

double gather(double li, const dvector &field) {
	int i = (int) li;
	double di = li - i;
	return field[i]*(1-di) + field[i+1]*di;
}
