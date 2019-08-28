// Potencial.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//#include "pch.h"
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include "omp.h"


#define ntrial 350 // The number of total evolutions
#define ds 0.1
#define Nbar 100 // The number of samples
#define k 1.0
#define m 1.0
#define h 1.0
#define VeryLargeValue 1000

const long int N = 100000; // Number of particles;
double left = -1;
double right = 1;
double D = h * h / (2 * m);
double dt = ds * ds / (2 * D);
int nave = ntrial / Nbar;

double V(double x) {
	return 0.5 * k * x * x;
}

static bool sort_using_greater_than(double u, double v)
{
	return u > v;
}

void fill_xsite(std::vector<double>& xsite) {
#pragma omp parallel
	{
		std::mt19937 gen;
		gen.seed(time(0) * (omp_get_thread_num() + 2));
		std::uniform_real_distribution<> urd(0, 1);
		#pragma omp for schedule(dynamic, 1000)
		for (size_t i = 0; i < xsite.size(); i++) {
			xsite[i] = right - (right - left) * urd(gen);
		}
	}
}

double count_Vmean(std::vector<double>& xsite) {
	double Vmean = 0;
	#pragma omp parallel for reduction(+:Vmean)
	for (size_t i = 0; i < xsite.size(); i++) {
		Vmean += V(xsite[i]);
	}
	Vmean = Vmean / xsite.size();
	return Vmean;
}

std::vector<double> time_shift(std::vector<double>& xsite, double Vref, int trial) {
	std::vector<double> xsiteNew;
	std::vector<double> NumbersToDelete;
	#pragma omp parallel
	{	
		std::mt19937 gen;
		gen.seed(time(0) * (omp_get_thread_num()+2) * (trial+1));
		std::uniform_real_distribution<> urd(0, 1);
		#pragma omp for schedule(dynamic, 100)
		for (size_t i = 0; i < xsite.size(); i++) { 
			double dV = (V(xsite[i]) - Vref) * dt;
			double r = urd(gen);
			/*
#pragma omp critical(print)
			{
				std::cout << "ThNum: " << omp_get_thread_num() << " : " << r << std::endl;
			}
			*/
			if (dV < 0) {
				if (r < -dV) {
					#pragma omp critical(insert)
					{
						xsiteNew.insert(xsiteNew.end(), xsite[i]);
					}
				}
			}
			else {
				if (r < dV) {
					xsite[i] = VeryLargeValue;
				}
			}
		}
	}
	for (size_t i = 0; i < xsite.size(); i++) {
		if (xsite[i] == VeryLargeValue) {
			xsite.erase(xsite.begin() + i);
			i -= 1;
		}
	}
	return xsiteNew;
}

void space_shift(std::vector<double>& xsite, int trial) {
#pragma omp parallel
	{
		std::mt19937 gen;
		gen.seed(time(0) * (omp_get_thread_num() + 2) * (trial + 1));
		std::uniform_real_distribution<> urd(0, 1);
		#pragma omp for schedule(dynamic, 1000)
		for (size_t i = 0; i < xsite.size(); i++) {
			if (urd(gen) < 0.5) {
				xsite[i] = xsite[i] + ds;
			}
			else {
				xsite[i] = xsite[i] - ds;
			}
		}
	}
	
}

void fix_particle_number(std::vector<double>& xsite, int trial) {
	double V_ref_new = 0;
	long int dN = xsite.size() - N;
	long int N_current = xsite.size();
	std::vector<double> xsiteNew;
	std::vector<double> NumbersToDelete;
	if (dN >= 0) {
		#pragma omp parallel 
		{
			std::mt19937 gen;
			gen.seed(time(0) * (omp_get_thread_num() + 2) * (trial + 1));
			std::uniform_real_distribution<> urd(0, 1);
			#pragma omp for schedule(dynamic, 1000)
			for (size_t i = 0; i < xsite.size(); i++) {
				double r = urd(gen);
				if (r < dN / (N_current + 0.0)) {
					xsite[i] = VeryLargeValue;
				}
			}
		}
		for (size_t i = 0; i < xsite.size(); i++) {
			if (xsite[i] == VeryLargeValue) {
				xsite.erase(xsite.begin() + i);
				i -= 1;
			}
		}
	}
	else if (dN < 0) {
		#pragma omp parallel
		{
			std::mt19937 gen;
			gen.seed(time(0) * (omp_get_thread_num() + 2) * (trial + 1));
			std::uniform_real_distribution<> urd(0, 1);
			#pragma omp for schedule(dynamic, 1000)
			for (size_t i = 0; i < xsite.size(); i++) {
				double r = urd(gen);
				if (r < -dN / (N_current + 0.0)) {
					#pragma omp critical(toDublicate)
					{
						xsiteNew.insert(xsiteNew.end(), xsite[i]);
					}
				}
			}
		}
		xsite.insert(xsite.end(), xsiteNew.begin(), xsiteNew.end());
	}
}

void walk() {
	/*
	FILE *f;
	f = fopen(R"(C:\Users\George\Desktop\git_projects\C_plus_plus\SchrodingerEq\output\output.csv)", "w");
	fprintf(f, "x\n");
	fclose(f);
	*/
	std::vector<double> xsite(N);  // Vector with particle coordinates
	// Filling xsite
	fill_xsite(xsite);
	double Vref = count_Vmean(xsite);
	std::vector<double> Vref_arr(ntrial+1);
	Vref_arr[0] = Vref;

	for (int trial = 0; trial < ntrial; trial++) {
		auto xsiteNew = time_shift(xsite, Vref, trial);
		space_shift(xsite, trial);
		xsite.insert(xsite.end(), xsiteNew.begin(), xsiteNew.end());
		fix_particle_number(xsite, trial);
		
		Vref = count_Vmean(xsite);
		Vref_arr[trial + 1] = Vref;

		/*
		if ((trial + 1) % nave == 0) {
			double perc = (double)(trial + 1) / (double)ntrial * 100;
			std::cout << perc << " percent of work completed" << std::endl;
		std::cout << Vref << std::endl;
		}*/
	}
	/*
	double result = 0;
	for (int i = 0; i < 50; i++) {
		result += Vref_arr[ntrial - 50 + i];
	}
	result = result / 50;
	*/
	double result = Vref_arr[Vref_arr.size() - 1];
	std::cout << "Final result is " << result << std::endl;
	/*
	f = fopen(R"(C:\Users\George\Desktop\git_projects\C_plus_plus\SchrodingerEq\output\output.csv)", "a");
	if (f) {
		for (int i = 0; i < xsite.size(); i++) {
			fprintf(f, "%lf\n", xsite[i]);
		}
	}*/
}

int main() {
	walk();
}
