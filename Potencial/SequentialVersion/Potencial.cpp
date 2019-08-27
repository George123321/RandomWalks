// Potencial.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <windows.h>

constexpr long int N = 1000; // Number of particles;
#define ntrial 150000 // The number of total evolutions
#define ds 0.1
#define Nbar 100 // The number of samples
#define k 1.0
#define m 1.0
#define h 1.0
#define OutputDataCoord "output/output.csv"
#define OutputDataNE "output/outputNE.csv"

double left = -1;
double right = 1;
double D = h * h / (2 * m);
double dt = ds * ds / (2 * D);
int nave = ntrial / Nbar;

double V(double x) {
	return 0.5 * k * x * x;
}

void walk() {
	FILE* f;
	FILE* fNE;

	f = fopen(OutputDataCoord, "w");
	fprintf(f, "x\n");
	fclose(f);
	f = fopen(OutputDataCoord, "a");

	fNE = fopen(OutputDataNE, "w");
	fprintf(fNE, "N,E_r\n");
	fclose(fNE);
	fNE = fopen(OutputDataNE, "a");

	std::vector<double> xsite(N);
	double Vref = 0;

	std::mt19937 gen;
	gen.seed(time(0));
	std::uniform_real_distribution<> urd(0, 1);

	for (size_t i = 0; i < xsite.size(); i++) {
		xsite[i] = right - (right - left) * urd(gen);
		Vref += V(xsite[i]);
	}
	Vref = Vref / xsite.size();

	if (f) {
		for (size_t i = 0; i < xsite.size(); i++) {
			fprintf(f, "%lf\n", xsite[i]);
		}
	}
	if (fNE) {
		fprintf(fNE, "%d,%f\n", xsite.size(), Vref);
	}

	for (int trial = 0; trial < ntrial; trial++) {
		std::vector<double> xsiteNew;
		for (size_t i = 0; i < xsite.size(); i++) {
			double dV = (V(xsite[i]) - Vref) * dt;
			double r = urd(gen);
			if (dV < 0) {
				if (r < -dV) {
					xsiteNew.insert(xsiteNew.end(), xsite[i]);
				}
			}
			else {
				if (r < dV) {
					xsite.erase(xsite.begin() + i);
					i = i - 1;
				}
			}
			// Shift the pedestrian
			if (!(dV > 0 and r < dV)) { // If the pedestrian is alive
				if (urd(gen) < 0.5) {
					xsite[i] = xsite[i] + ds;
				}
				else {
					xsite[i] = xsite[i] - ds;
				}
			}
		}
		xsite.insert(xsite.end(), xsiteNew.begin(), xsiteNew.end());

		double V_ref_new = 0;
		long int dN = xsite.size() - N;
		long int N_current = xsite.size();
		std::vector<double> xsiteNew1;
		if (dN >= 0) {
			for (size_t i = 0; i < xsite.size(); i++) {
				double r = urd(gen);
				if (r < dN / (N_current + 0.0)) {
					xsite.erase(xsite.begin() + i);
					i = i - 1;
				}
				else {
					V_ref_new += V(xsite[i]);
				}
			}
		}
		else if (dN < 0) {
			for (size_t i = 0; i < xsite.size(); i++) {
				double r = urd(gen);
				if (r < -dN / (N_current + 0.0)) {
					xsiteNew1.insert(xsiteNew1.end(), xsite[i]);
					V_ref_new += 2 * V(xsite[i]);
				}
				else {
					V_ref_new += V(xsite[i]);
				}
			}
			xsite.insert(xsite.end(), xsiteNew1.begin(), xsiteNew1.end());
		}
		Vref = V_ref_new / xsite.size();

		if (f) {
			for (size_t i = 0; i < xsite.size(); i++) {
				fprintf(f, "%lf\n", xsite[i]);
			}
		}
		if (!f) {
			printf("Something is wrong! Not f\n");
		}
		if (fNE) {
			fprintf(fNE, "%d,%f\n", xsite.size(), Vref);
		}
		else { printf("Something is wrong! Not NE\n"); }

		if ((trial + 1) % nave == 0) {
			double perc = (double)(trial + 1) / (double)ntrial * 100;
			std::cout << perc << " percent of work completed" << std::endl;
		}
	}

}

int main() {
	walk();
}
