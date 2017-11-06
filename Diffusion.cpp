// One Dimensional Diffusion by Floyd Wiseman
// translated into C++ by R. Stauber

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;
const double PI = 2 * acos(0.0);

// ************************ DECLARE VARIABLES *********************************

double n_x = 0.1870133;
double n_y = 0.2093135;
double n_z = 0.3698417;
double s_x = -1.0;
double s_y = 1.0;
double s_z = 1.0;
double energy = 0;
double beta_b, beta_d;

double t, p, s, t_c, c_b, d, d_c, m_d, m_b, m_j, m_k, r_ave, lam, z_b;
double del_v, del_z, del_zd, x_run, r_run, fr, sum_a, speed, adj;

double r_b, r_d, v_b, v_d;
double n_b, n_d, n_sb, n_zb, n_sd, sum_bi, it, itm, n_bd, coll_2b, coll_bd, n_ref;
double n_steps, steps, coll_2d, o_o, o;
double time, delta_t, dx, dy, dz, r_check, dist, phi;
float x[95000], y[95000], z[95000], xc[95000], yc[95000], zc[95000];
float vx[95000], vy[95000], vz[95000], vxc[95000], vyc[95000], vzc[95000];
float  x_j[95000], y_j[95000], z_j[95000], x_k[95000], y_k[95000], z_k[95000];
int tag[9500], n_cb[1000], hit[9500], n_cd[1000];
int sum_ter, n_pair, k_b, o_b;


// ************************ FUNCTION DEFINITIONS ******************************
double getDouble(string message) {
	double myNumber = 0;
	string input = "";
	while (true) {
		cout << message << ":  ";
		getline(cin, input);
		// This code converts from string to number safely.
		stringstream myStream(input);
		if (myStream >> myNumber)
			break;
		cout << "Invalid, please try again" << endl;
	}
	// cout << "You entered: " << myNumber << endl << endl;
	return myNumber;
}

void GetInputs() {
	t = getDouble("Input the temperature in Kelvin");
	p = getDouble("Input the pressure in atm for the bath molecules");
	del_z = pow((0.132626*t / p), (1.0 / 3.0));				// cell size
	cout << "The cell size in nm is " << del_z << endl;
	n_sb = getDouble("Input the number of bath molecules along one of\n"
		"the non-diffusing (x and y) coordinates");
	s = del_z * n_sb;											// length of a side
	cout << "The xy area of the container is " << pow(s, 2) << " nm^2" << endl;
	n_zb = getDouble("Input the number of bath molecules along the diffusing coordinate");
	d = del_z * n_zb;											// bin length
	cout << "The bin length is " << d << " nm" << endl;
	n_b = n_sb * n_sb * n_zb;						// number of bath particles per bin
	cout << "The number of bath particles in each bin is " << n_b << endl;
	m_b = getDouble("Input the molar mass of the bath molecules in kg/mol");
	m_d = getDouble("Input the molar mass of the diffusing molecules in kg/mol");
	r_b = getDouble("Input the molecular radius of the bath molecules in nm");
	r_d = getDouble("Input the molecular radius of the diffusing molecules in nm");
	v_b = (4.994347e9*sqrt(t / m_b));
	v_d = (4.994347e9*sqrt(t / m_d));
	beta_b = (6.01361e-20)*m_b / t;
	beta_d = (6.01361e-20)*m_d / t;
	lam = 0.043373 * t / (p * sqrt(1. + m_d / m_b)*(r_b + r_d)*(r_b + r_d)); // mean free path
	cout << "The mean free path for the dilute diffusing molecules is " << lam << endl;
	cout << "The collision probability within a bin is " << 1. - exp(-d / (lam*cos(1.))) << endl;
	while (true) {
		n_sd = getDouble("Input the number of diffusing molecules along the x coordinate");
		if (n_sd < (s / (2.5*r_d))) break;
		cout << "The number of diffusing molecules is too large.  Try again" << endl;
	}
	del_zd = s / n_sd;
	n_d = n_sd * n_sd;
	cout << "The number of diffusing molecules is " << n_d << endl;
	cout << "They are spaced " << del_zd << " nm apart in the xy plane" << endl;
}

// ********************************  MAIN  *************************************
int main()
{
	cout << "One-Dimensional Diffusion Model\n";
	GetInputs();
	// *** Particle placements ***
	int q, k;
	q = k = 1;
	for (int i = 1; i <= n_sd; i++) {
		for (int j = 1; j <= n_sd; j++) {
			x[q] = (float(i) - 0.5) * del_zd + (0.125 * r_d / d) * cos(sqrt(float(j + q)));
			y[q] = (float(j) - 0.5) * del_zd + (0.125 * r_d / d) * sin(sqrt(float(j + i)));
			z[q] = (0.5 - float(k)) * del_z +  (0.25 *  r_d / d) * sin(sqrt(float(i + q)));
			if (k == 5) k = 0;
			k++;
			q++;
		}
	}
	q = 1 + n_d;
	for (int i = 1; i <= n_sb; i++) {
		for (int j = 1; j <= n_sb; j++) {
			for (int k = 1; k <= n_zb; k++) {
				x[q] = (float(i) - 0.5) * del_z + (0.125 * r_b / d) * sin(sqrt(float(j + q)));
				y[q] = (float(j) - 0.5) * del_z + (0.125 * r_b / d) * cos(sqrt(float(k + q)));
				z[q] = (0.5 - float(k)) * del_z + (0.25 *  r_b / d) * cos(sqrt(float(i + k)));
				q++;
			}
		}
	}
	// *** Velocity Vectors ***
	double del_v = 0.00005 / sqrt(beta_d);
	for (int i = 1; i <= n_d; i++) {
		n_x += 0.2730911;
		if (n_x >= 1.0) n_x -= 1.0;
		n_y += 0.2300787;
		if (n_y >= 1.0) n_y -= 1.0;
		n_z += 0.2911071;
		if (n_z >= 1.0) n_z -= 1.0;

		vx[i] = 0.5 * del_v;
		vy[i] = 0.5 * del_v;
		vz[i] = 0.5 * del_v;

		sum_a = 0.0;
		while (sum_a < n_x) {
			sum_a = sum_a + 2.0 * sqrt(beta_d / PI) * del_v * exp(-beta_d * vx[i] * vx[i]);
			vx[i] += del_v;
			s_x = -s_x;
		}
		sum_a = 0.0;
		while (sum_a < n_y) {
			sum_a = sum_a + 2.0 * sqrt(beta_d / PI) * del_v * exp(-beta_d * vy[i] * vy[i]);
			vy[i] += del_v;
			s_y = -s_y;
		}
		sum_a = 0.0;
		while (sum_a < n_z) {
			sum_a = sum_a + 2.0 * sqrt(beta_d / PI) * del_v * exp(-beta_d * vz[i] * vz[i]);
			vz[i] += del_v;
			s_z = -s_z;
		}
		vx[i] = s_x * vx[i];
		vy[i] = s_y * vy[i];
		vz[i] = s_z * vz[i];

		energy += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
	}
	t_c = m_d * energy / (2.49434e19 * n_d);
	adj = t / t_c;
	energy = 0.0;
	del_v = 0.00005 / sqrt(beta_b);
	for (int i = 1; i <= n_d; i++) {
		vx[i] = sqrt(adj) * vx[i];
		vy[i] = sqrt(adj) * vy[i];
		vz[i] = sqrt(adj) * vz[i];
	}
	for (int i = n_d + 1; i <= n_d + n_b; i++) {
		n_x += 0.2199431;
		if (n_x >= 1.0) n_x -= 1.0;
		n_y += 0.3170933;
		if (n_y >= 1.0) n_y -= 1.0;
		n_z += 0.2609717;
		if (n_z >= 1.0) n_z -= 1.0;

		vx[i] = 0.5 * del_v;
		vy[i] = 0.5 * del_v;
		vz[i] = 0.5 * del_v;

		sum_a = 0.0;
		while (sum_a < n_x) {
			sum_a = sum_a + 2.0 * sqrt(beta_b / PI) * del_v * exp(-beta_b * vx[i] * vx[i]);
			vx[i] += del_v;
			s_x = -s_x;
		}
		sum_a = 0.0;
		while (sum_a < n_y) {
			sum_a = sum_a + 2.0 * sqrt(beta_b / PI) * del_v * exp(-beta_b * vy[i] * vy[i]);
			vy[i] += del_v;
			s_y = -s_y;
		}
		sum_a = 0.0;
		while (sum_a < n_z) {
			sum_a = sum_a + 2.0 * sqrt(beta_b / PI) * del_v * exp(-beta_b * vz[i] * vz[i]);
			vz[i] += del_v;
			s_z = -s_z;
		}
		vx[i] = s_x * vx[i];
		vy[i] = s_y * vy[i];
		vz[i] = s_z * vz[i];

		energy += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
	}
	t_c = m_b * energy / (2.49434e19 * n_b);
	adj = t / t_c;
	for (int i = n_d + 1; i <= n_d + n_b; i++) {
		vx[i] = sqrt(adj) * vx[i];
		vy[i] = sqrt(adj) * vy[i];
		vz[i] = sqrt(adj) * vz[i];
	}

	r_ave = getDouble("Input the average distance travelled by a bath molecule"
		" in nm per time step");
	delta_t = r_ave / v_b;
	cout << "The time step in seconds is " << delta_t << endl;

	cout << "Press enter to exit";
	cin.ignore();
	return 0;
}