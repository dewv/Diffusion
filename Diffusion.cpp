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
int tag[95000], n_cb[1000], hit[95000], n_cd[1000];
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

// ********************************  MAIN  *************************************
int main()
{
	cout << "One-Dimensional Diffusion Model\n";
	// *** Get Inputs ***
	t = getDouble("Input the temperature in Kelvin");
	p = getDouble("Input the pressure in atm for the bath molecules");
	del_z = pow((0.132626 * t / p), (1.0 / 3.0));				// cell size
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
	lam = 0.043373 * t / (p * sqrt(1. + m_d / m_b) * (r_b + r_d) * (r_b + r_d)); // mean free path
	cout << "The mean free path for the dilute diffusing molecules is " << lam << endl;
	cout << "The collision probability within a bin is " << 1. - exp(-d / (lam*cos(1.))) << endl;
	while (true) {
		n_sd = getDouble("Input the number of diffusing molecules along the x coordinate");
		if (n_sd < (s / (2.5 * r_d))) break;
		cout << "The number of diffusing molecules is too large.  Try again" << endl;
	}
	del_zd = s / n_sd;
	n_d = n_sd * n_sd;
	cout << "The number of diffusing molecules is " << n_d << endl;
	cout << "They are spaced " << del_zd << " nm apart in the xy plane" << endl;
	// *** Particle placements ***
	int q, k;
	q = k = 0;
	for (int i = 1; i <= n_sd; i++) {
		for (int j = 1; j <= n_sd; j++) {
			x[q] = (float(i) - 0.5) * del_zd + (0.125 * r_d / d) * cos(sqrt(float(j + q)));
			cout << x[q] << endl;
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
	z_b = 60.08 * r_b *r_b * n_b * v_b * p / t;
	cout << "The b-b collision frequency in sec^-1 is " << z_b << endl;
	n_steps = getDouble("Input the number of time steps before diffusion is turned on");
	steps = getDouble("Input the number of time steps before outputting data");

	// *** Initialize trajectory sim variables ***
	it = 1;
	itm = 1;
	k_b = 1;
	coll_2b = 0;
	coll_bd = 0;
	coll_2d = 0;
	n_ref = 0;
	sum_bi = 0;
	sum_ter = 0;
	time = 0;
	x_run = 0;
	r_run = 0;
	fr = 0;
	for (int i = 0; i < 95000; i++) {
		hit[i] = 0;
		tag[i] = 0;
	}
	
	cout << "Beginning trajectory simulation" << endl;

	for (o_o = 0; o_o < 20; o_o++) {
		for (o = 0; o < steps; o++) {
			if (it == n_steps) {
				q = 0;
				for (int i = 0; i < n_sd; i++) {
					for (int j = 0; j < n_sd; j++) {
						x[q] = (float(i) - .5)*del_zd;
						y[q] = (float(j) - .5)*del_zd;
						z[q] = r_d;
						vz[q] = abs(vz[q]);
						q++;
					}
				}
			}
			// *** Translation equations ***
			for (int i = 0, i < n_d + k_b * n_b; i++) {
				x[i] = x[i] + vx[i] * delta_t;
				y[i] = y[i] + vy[i] * delta_t;
				z[i] = z[i] + vz[i] * delta_t;
			}
			// *** Wall recoil ***
			if (it < n_steps) {
				// Compartmentalized confinement in the z direction
				for (int i = 0; i < n_d; i++) {
					if ((z[i] < -0.5 * del_z && vz[i] < 0) ||
						(z[i] > -r_d         && vz[i] > 0)) {
						if (tag[i] == 0) vz[i] = -vz[i];
					}
				}
				for (int i = n_d; i < n_d + k_b * n_b) {
					if ((z[i] < r_b + r_d && vz[i] < 0) ||
						(z[i] > d - r_b    && vz[i] > 0)) {
						if (tag[i] == 0) vz[i] = -vz[i];
					}
				}
			}
			else {
				// Full container confinement in the z direction
				for (int i = 0; i < n_d; i++) {
					if ((z[i] < r_d					 && vz[i] < 0) ||
						(z[i] > float(k_b) * d - r_d && vz[i] > 0)) {
						if (tag[i] == 0) vz[i] = -vz[i];
					}
				}
				for (int i = n_d; i < n_d + k_b * n_b) {
					if ((z[i] < r_b					  && vz[i] < 0) ||
						(z[i] >= float(k_b) * d - r_b && vz[i] > 0)) {
						if (tag[i] == 0) vz[i] = -vz[i];
					}
				}
			}
			for (int i = 0; i < n_d + k_b * n_b) {
				// Confinement in the x and y directions
				if (i < n_d) r_check = r_d;
				else rcheck = r_b;
				if ((x[i] <= r_check     && vx[i] < 0) ||
					(x[i] >= s - r_check && vx[i] > 0)) {
					if (tag[i] == 0) vx[i] = -vx[i];
				}
				if (y[i] <= r_check) && vy[i] < 0) ||
				(y[i] >= s - r_check && vy[i] > 0)) {
				if (tag[i] == 0) vy[i] = -vy[i];
				}
			}
			// *** Molecular recoil module
			j = k = 0;
			while (j < n_d + k_b * n_b) {
				while (k < n_d + k_b * n_b) {
					dx = x[j] - x[k];
					dy = y[j] - y[k];
					dz = z[j] - z[k];
					dist = dx*dx + dy*dy + dz*dz;
					if (sqrt(dist) > 1.5 * (r_b + r_d)) continue; // goto 100 (i.e. to next k)
						if ((j < n_d && k < n_d)) {
							r_check = 2.0 * r_d;
							m_j = m_d;
							m_k = m_d;
							if (sqrt(dist) <= r_check && tag[j] == 0 && tag[k] == 0 && time == 0) {
								coll_2d++;
							}
						}
					// p. 7
					if (j <= n_d && k > n_d) {
						r_check = r_b + r_d;
						m_j = m_d;
						m_k = m_b;
						if (sqrt(dist) <= r_check && tag[j] == 0 && tag[k] == 0 && time > 0) {
							coll_bd++;
						}
					}
					if (j > n_d && k <= n_d) {
						r_check = r_b + r_d;
						m_j = m_b;
						m_k = m_d;
						if (sqrt(dist) <= r_check && tag[j] == 0 and tag[k] == 0 && time > 0) {
							collbd++;
						}
					}
					if (j > n_d && k > n_d) {
						r_check = 2 * r_b;
						m_j = m_b;
						m_k = m_d;
						if (sqrt(dist) <= r_check && tag[j] == 0 and tag[k] == 0 && time > 0) {
							coll2b++;
						}
					}
					if (sqrt(dist) <= r_check && tag[j] == 0 and tag[k] == 0) {
						// molecular recoil
						if (j <= n_d && k <= n_d && it >= n_steps) break;
						phi = 2 * (dx * (vx[j] - vx[k]) + dy * (vy[j] - vy[k]) + dz * (vz[j] - vz[k])) / (dist * (m_j + m_k));
						if (time > 0) { 
							// reflection criteria
							if (j > n_d && k <= n_d && vz[k] / (vz[k] + phi * m_j * dz) < 0) n_ref++;
							if (j <= n_d && k > n_d && vz[j] / (vz[j] + phi * m_k * dz) < 0) n_ref++;
						}
						vx[k] += phi * m_j * dx;
						vy[k] += phi * m_j * dy;
						vz[k] += phi * m_j * dz;
						vx[j] += phi * m_k * dx;
						vy[j] += phi * m_k * dy;
						vz[j] += phi * m_k * dz;
						hit[k]++;
						hit[j]++;
						sum_bi++;
					}
					// p.8
					if (sqrt(dist) < r_check) {
						x_j[j] = x[j] + vx[j] * delta_t;
						y_j[j] = y[j] + vy[j] * delta_t;
						z_j[j] = z[j] + vz[j] * delta_t;
						x_k[j] = x[k] + vx[k] * delta_t;
						y_k[j] = y[k] + vy[k] * delta_t;
						z_k[j] = z[k] + vz[k] * delta_t;
						dist = sqrt(pow((x_j[j] - x_k[k]), 2) + pow((y_j[j] - y_k[k]), 2) + pow((z_j[j] - z_k[k]), 2));
						if (dist <= r_check && tag[j] == 0 && tag[k] == 0) {
							tag[j] = 1;
							tag[k] = 1;
						}
						else if (dist > r_check && tag[j] == 1 && tag[k] == 1) {
							tag[j] = 0;
							tag[k] = 0;
						}
						if (it > n_steps) {
							for (int i = 0; i < n_d; i++) {
								x_run += abs(vz[i]);
								r_run += sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
								fr = x_run / r_run;
							}
						}
					}
					// *** 100 ***
				} // end of k while loop in molecular recoil
			} // end of j while loop in molecular recoil
			for (int i = 0; i < n_d + k_b *n_b) {
				if (hit[i] > 1) sum_ter++;
				hit[i] = 0;
			}
			if (n_steps == it) {
				j = 0;
				for (i = n_d; i < n_d + n_b) {
					j++;
					vxc[j] = vx[i];
					vyc[j] = vy[i];
					vzc[j] = vz[i];
					xc[j] = x[i];
					yc[j] = y[i];
					zc[j] = z[i];
				}
			}
			for (i = 0; i < n_d; i++) {
				if (z[i] > (float(k_b) * d - r_d)) {
					k = 0;
					for (j = n_d + k_b * n_b; j < n_d + (k_b + 1) * n_b; j++) {
						k++;
						vx[j] = vxc[k];
						vy[j] = vyc[k];
						vz[j] = vzc[k];
						x[j] = xc[k];
						y[j] = yc[k];
						z[j] = zc[k] + float(k_b) * d;
					}
					k_b++;
				}
			}
			if (it >= n_steps) {
				time += delta_t;
			}
			it++;
			itm += 1;
			if (itm == 100) {
				cout << "The time step number is " << it << endl;
				if (time == 0) cout << "Spin-up phase" << endl;
				else cout << "Diffusing phase" << endl;
				cout << "<distance> (nm) for the bath molecules is " << float(it) * delta_t * v_b << endl;
				cout << "<distance> (nm) for the diffusing molecules is " << float(it) * delta_t * v_d << endl;
				cout << "The cumulative number of bimolecular collisions is " << sum_bi << endl;
				cout << "The cumulative number of termolecular collisions is " << sum_ter << endl;
				if (time == 0) cout << "The cumulative number of d-d collisions is " << coll_2d << endl;
				cout << "The number of bins is " << k_b << endl;
				itm = 0;
			}

		} // end of o for loop
		if (time > 0) {
			cout << "The time is " << time << endl;
			cout << "The average collision frequency for the bath molecules is " << float(coll_2b) / time << endl;
			cout << "The average collision frequency between the bath and diffusing molecules is " << float(coll_bd) / time << endl;
			if (coll_bd > 0) cout << "The refection coefficient is " << float(n_ref) / float(coll_bd) << endl;
			cout << "The ratio between the distance traveled in one coordinate to the total distance travelled is " << fr << endl;
			cout << "     Bin #     [diffusing molecules]/nm^-3     [Total]" << endl;
			for (i = 0; i < k_b; i++) {
				n_cd[i] = 0;
				n_cb[i] = 0;
			}
			for (j = 0; j < n_d; j++) {
				for (i = 0; i < k_b; i++) {
					if (z[j] >= float(i - 1)*d && z[j] <= float(i) * d) {
						n_cd[i]++;
						exit;
					}
				}
			}
			for (j = n_d; j < n_d + k_b * n_b; j++) {
				for (i = 0; i < k_b; i++) {
					if (z[j] >= float(i - 1) * d && z[j] <= float(i) * d) {
						n_cb[i]++;
						exit;
					}
				}
			}
			for (i = 0; i < k_b; i++) {
				cout << i + 1 << "          " << float(n_cd[i]) / (s * s * d) << "          " << float(n_cd[i] + n_cb[i]) / (s * s * d) << endl;
			}
		} // end of time > 0 if
		n_pair = 0;
		for (i = 0; i < n_d + k_b * n_b; i++) {
			for (j = i; j < n_d + k_b * n_b) {
				if (i > n_d && j > n_d) r_check = 2 * r_b;
				else if (i <= n_d && j <= n_d) r_check = 0;
				else r_check = r_b + r_d;
				if (sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2) + pow(z[i] - z[j], 2)) <= r_check) n_pair++;
			}
		}
		cout << "The number of molecules currently in contact are " << 2 * float(n_pair) << endl;
		string input = "";
		while (true) {
			cout << "Do you want to continue the simulation? (y or n)" << ":  ";
			getline(cin, input);
			stringstream myStream(input);
			if (myStream == 'y' or myStream == 'n')
				break;
			cout << "Invalid, please try again" << endl;
		}
		cout << "You entered: " << myStream << endl << endl;
		if (myStream == 'n') break;
		else steps = getDouble("Input the new number of time steps");
	} // end of o_o for loop
	cout << "Press enter to exit";
	cin.ignore();
	return 0;
}