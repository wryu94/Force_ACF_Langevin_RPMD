//
//  main.cpp
//  Force_ACF_remastered
//
//  Created by Harry Ryu on 2017. 5. 11..
//  Copyright © 2017년 Harry Ryu. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;

int P = 6; // Number of beads
double h_bar = 1;
double mass = 1;
int N = 10000; // Number of timesteps taken
int N_samp = 1000; // Number of sampling from the canonical distribution
double beta = 8;
double beta_P = beta/P;
double dt = 0.01; // Timestep
double quad = 0.5;
double cub = 0.1;
double quart = 0.01;

double velocity_dist()
// Box Muller algorithm that returns Gaussian random number with 0 mean and nonunitary variance of 1/(beta_P * mass)
// Return velocity in unit of Å/fs
{
    double first, v1, v2, rsq, fac;
    int again = 1;
    while (again)
    {
        v1 = 2.0*rand()/RAND_MAX - 1.0;
        v2 = 2.0*rand()/RAND_MAX - 1.0;
        rsq = v1*v1 + v2*v2;
        if (rsq < 1.0 & rsq != 0.0) again = 0;
    }
    fac = sqrt(-2.0*log(rsq)/rsq);
    first = v1 * fac * sqrt(1/(beta_P * mass));
    return first;
}

double pot_deriv(double x)
{
    double value;
    value = 2 * quad * x + 3 * cub * x * x + 4 * quart * x * x * x;
    return value;
}

void force_calc(double *F, int time, double **matrix, double *x_normal_0, double *v_normal_0, double *freq)
{
    for (int k = 0 ; k < P ; ++k)
    {
        F[k] = 0.0;
        for (int l = 1 ; l < P ; ++l)
        {
            F[k] = F[k] + pot_deriv(matrix[k][l] * (x_normal_0[l] * cos(freq[l] * time * dt) + (v_normal_0[l] / freq[l]) * sin(freq[l] * time * dt)));
        }
    }
}

void corr_calc(double corr, double *force_time, double *force_init)
{
    corr = 0.0;
    for (int k = 0 ; k < P ; ++k)
    {
        for (int m = 0 ; m < P ; ++m)
        {
            corr = corr + (force_init[k] * (force_time[m] - force_init[m]));
        }
    }
}

void print(double *inp)
{
    for (int i = 0 ; i < P ; ++i)
    {
        cout << inp[i] << " ";
    }
    cout << endl;
}

int main()
{
    double** O = new double*[P]; // Coordinate transform matrix, from bead to normal mode
    double** O_T = new double*[P]; // Coordinate transform matrix, from normal mode to bead
    for (int i = 0 ; i < P ; ++i)
    {
        O[i] = new double[P];
        O_T[i] = new double[P];
    }
    
    for (int n = 0 ; n < P ; ++n)
    {
        O[0][n] = sqrt(1/static_cast<double>(P));
    }
    
    for (int m = 1 ; m < P/2 ; ++m)
    {
        for (int n = 0 ; n < P ; ++n)
        {
            O[m][n] = sqrt(2.0/static_cast<double>(P)) * cos(2 * M_PI * m * (n+1) / P);
        }
    }
    
    for (int n = 0 ; n < P ; ++n)
    {
        if (n % 2 == 0)
        {
            O[P/2][n] = -sqrt(1/static_cast<double>(P));
        }
        if (n % 2 == 1)
            O[P/2][n] = sqrt(1/static_cast<double>(P));
    }
    
    for (int m = P/2 + 1 ; m < P ; ++m)
    {
        for (int n = 0 ; n < P ; ++n)
        {
            O[m][n] = sqrt(2.0/static_cast<double>(P)) * sin( 2 * M_PI * m * (n+1) / P);
        }
    }
    
    for (int i = 0 ; i < P ; ++i)
    {
        for (int j = 0 ; j < P ; ++j)
        {
            O_T[i][j] = O[j][i];
        }
    }
    
    double* alpha = new double[P];
    for (int k = 0 ; k < P ; ++k)
    {
        alpha[k] = (2 / (beta_P * h_bar)) * sin(k * M_PI / P);
    }
    
    double* x_0 = new double[P]; // Initial normal mode position
    double* v_0 = new double[P]; // Initial normal mode velocity
    for (int i = 0 ; i < P ; ++i)
    // Initialize at position 1
    {
        x_0[i] = 1.0;
    }
    
    double** partial_force = new double*[N]; // partial_force[time][number_of_bead]
    for (int i = 0 ; i < N ; ++i)
    {
        partial_force[i] = new double[P];
    }
    
    double** correlation = new double*[N_samp];
    for (int i = 0 ; i < N_samp ; ++i)
    {
        correlation[i] = new double[N];
    }
    
    double* corr_tot = new double[N];
    
    // ------------------------------------------------------------------
    // ------------------------------------------------------------------
    for (int rep = 0 ; rep < N_samp ; ++rep)
    {
        for (int i = 0 ; i < P ; ++i)
        {
            v_0[i] = velocity_dist(); // Initialize velocity every run
        }
        
        for (int time = 0 ; time < N ; ++time)
        {
            force_calc(partial_force[time], time, O_T, x_0, v_0, alpha);
            
            correlation[rep][time] = 0.0;
            for (int k = 0 ; k < P ; ++k)
            {
                for (int m = 0 ; m < P ; ++m)
                {
                    correlation[rep][time] = correlation[rep][time] + (partial_force[0][k] * (partial_force[time][m] - partial_force[0][m]));
                }
            }
        }
        
        cout << "Run " << rep + 1 << " out of " << N_samp << endl;
    }
    
    for (int i = 0 ; i < N ; ++i)
    {
        corr_tot[i] = 0.0;
        for (int k = 0 ; k < N_samp ; ++k)
        {
            corr_tot[i] = corr_tot[i] + correlation[k][i] / N_samp;
        }
    }
    
    for (int i = 0 ; i < N ; ++i)
    {
        cout << corr_tot[i] << endl;
    }
    
    
    for (int i = 0 ; i < P ; ++i)
    {
        delete[] O[i];
        delete[] O_T[i];
    }
    delete[] O;
    delete[] O_T;
    
    for (int i = 0 ; i < N_samp ; ++i)
    {
        delete[] correlation[i];
    }
    delete[] correlation;
    
    for (int i = 0 ; i < N ; ++i)
    {
        delete[] partial_force[i];
    }
    delete[] partial_force;
    
    delete[] corr_tot;
    
    return 0;
}
