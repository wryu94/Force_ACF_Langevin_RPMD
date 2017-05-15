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

int P = 32; // Number of beads
double h_bar = 1;
double mass = 1;
int N = 50000; // Number of timesteps taken
double dt = 0.01; // Timestep
int corr_length = 10000;

int N_samp = 100; // Number of sampling from the canonical distribution
double beta = 8; // Temperature
double beta_P = beta/P; // Scaled temperature
double quad = 0.0;
double cub = 0.0;
double quart = 0.25;

double velocity_dist()
// Box Muller algorithm that samples velocity from the canonical distribution with variance P/(beta*mass)
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
// Derivative of the potential
{
    double value;
    value = 2 * quad * x + 3 * cub * x * x + 4 * quart * x * x * x;
    return value;
}

void cent_const_bead(double *bead, int time, double **matrix, double *x_normal_0, double *v_normal_0, double *freq)
// Position of the bead with centroid constrained at 0 using the harmonic motion of the noncentroid modes
{
    for (int k = 0 ; k < P ; ++k)
    {
        bead[k] = 0.0;
        for (int l = 1 ; l < P ; ++l)
        {
            bead[k] = bead[k] + matrix[k][l] * (x_normal_0[l] * cos(freq[l] * time * dt) + (v_normal_0[l] / freq[l]) * sin(freq[l] * time * dt));
        }
    }
}

void print(double *inp)
// Print function
{
    for (int i = 0 ; i < P ; ++i)
    {
        cout << inp[i] << " ";
    }
    cout << endl;
}

double *calc_corr(double *q)
// Calculate a correlation function of length corr_length for trajectory with length N
{
    double *corr = new double[corr_length];
    for (int n = 0 ; n < corr_length ; ++n)
    {
        corr[n] = 0.0;
        for (int i = 0 ; i < (N-n) ; ++i)
        {
            corr[n] = corr[n] + q[i] * q[i+n] / (N - n);
        }
    }
    return corr;
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
    
    double* alpha = new double[P]; // sqrt of the eigenvalues of the characteristic spring matrix
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
    
    double* restrained_bead = new double[P]; // centroid constrained bead position
    
    double** force_on_cent = new double*[N_samp]; // Force on centroid for each run and each timestep
    for (int i = 0 ; i < N_samp ; ++i)
    {
        force_on_cent[i] = new double[N];
    }
    
    double** correlation = new double*[N_samp]; // Correlation of the force
    for (int i = 0 ; i < N_samp ; ++i)
    {
        correlation[i] = new double[corr_length];
    }
    
    double* corr_tot = new double[N_samp]; // Average of the correlation
    
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
            cent_const_bead(restrained_bead, time, O_T, x_0, v_0, alpha); // Update the position of the restrained bead
            
            force_on_cent[rep][time] = 0.0; // Update the force on centroid
            for (int k = 0 ; k < P ; ++k)
            {
                force_on_cent[rep][time] = force_on_cent[rep][time] - pot_deriv(restrained_bead[k]);
            }
        }
        
        correlation[rep] = calc_corr(force_on_cent[rep]);
        
        cout << "Run " << rep + 1 << " out of " << N_samp << endl;
    }
    
    /*
    for (int rep = 0 ; rep < N_samp ; ++rep)
    {
        for (int i = 0 ; i < corr_length ; ++i)
        {
            cout << correlation[rep][i] << endl;
        }
        cout << endl;
    }
    */
    
    double initial = 0.0;
    for (int k = 0 ; k < N_samp ; ++k)
    {
        initial = initial + correlation[k][0] / N_samp;
    }
    
    for (int i = 0 ; i < corr_length ; ++i)
    {
        corr_tot[i] = 0.0;
        for (int k = 0 ; k < N_samp ; ++k)
        {
            corr_tot[i] = corr_tot[i] + correlation[k][i] / N_samp ;
        }
        corr_tot[i] = corr_tot[i] - initial;
    }
     
    /*
    for (int i = 0 ; i < corr_length ; ++i)
    {
        cout << corr_tot[i] << endl;
    }
    */
    
    // ------------------
    // Exporting force acf
    ofstream myfile;
    myfile.open ("force_acf.txt");
    for (int i = 0 ; i < corr_length ; ++i)
    {
        myfile << i * dt << "   " << corr_tot[i] << endl;
    }
    myfile.close();
    
    // ------------------
    for (int i = 0 ; i < P ; ++i)
    {
        delete[] O[i];
        delete[] O_T[i];
    }
    delete[] O;
    delete[] O_T;
    
    delete[] alpha;
    delete[] x_0;
    delete[] v_0;
    
    delete[] restrained_bead;
    
    for (int i = 0 ; i < N_samp ; ++i)
    {
        delete[] force_on_cent[i];
    }
    delete[] force_on_cent;
    
    for (int i = 0 ; i < N_samp ; ++i)
    {
        delete[] correlation[i];
    }
    delete[] correlation;
    
    delete[] corr_tot;
    
    return 0;
}
