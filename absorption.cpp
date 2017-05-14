//
//  main.cpp
//  Absorption_spectrum
//
//  Created by Harry Ryu on 2017. 4. 24..
//  Copyright © 2017년 Harry Ryu. All rights reserved.
//

//  This code calculates the absorption spectrum of RPMD system using generalized Lanagevin equation formalism

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
double beta = 8;
double beta_P = beta/P;
double dt = 0.01; // Timestep

int main()
{
    int N = 10000; // How long the correlation function is defined for
    
    int N_abs = N; // How long the spectrum will be constructed for
    
    double *correlation = new double[N];
    double *time = new double[N];
    
    double w_bar_squared = 1.0;
    // PMF frequency with harmonic approximation
    // force ACF initial value is always 0, so it's just the bare fequency
    
    FILE * pFile;
    
    pFile = fopen ("force_acf.txt", "r"); // Calling in force ACF
    for (int i = 0 ; i < N ; ++i)
    {
        fscanf(pFile, "%lf", &time[i]);
        fscanf(pFile, "%lf", &correlation[i]);
    }
    fclose (pFile);
    
    /*
    for (int i = 0 ; i < N ; ++i) // Checking whether we read it in correctly or not
    {
        cout << time[i] << "  " << correlation[i] << endl;
    }
    */
    
    double dw = (2 * M_PI)/time[N-1]; // Smallest unit of frequency
    double dt = time[1]; // Timestep
    
    double *ker = new double[N]; // Friction kernel
    for (int i = 0 ; i < N ; ++i)
    {
        ker[i] = (beta_P/mass) * correlation[i];
    }
    
    double *gamma_1 = new double[N_abs]; // Real part of half Fourier transform of friction kernel
    double *gamma_2 = new double[N_abs]; // Imaginary part of half Fourier transform of friction kernel
    
    for (int i = 0 ; i < N_abs; ++i)
    // FT of correlation function defined as P dimensional array
    // i is indexing of frequency
    {
        gamma_1[i] = 0;
        gamma_2[i] = 0;
        for (int k = 0 ; k < N-1 ; ++k)
        // k is indexing of time
        // Integration by trapezoid rule
        {
            gamma_1[i] = gamma_1[i] + (dt / 2) * (ker[k] * cos((i * dw) * (k * dt)) + ker[k+1] * cos((i * dw) * ((k+1) * dt)));
            gamma_2[i] = gamma_2[i] + (dt / 2) * (ker[k] * sin((i * dw) * (k * dt)) + ker[k+1] * sin((i * dw) * ((k+1) * dt)));
        }
        cout << "Fourier transform " << static_cast<double>(i+1)*100/static_cast<double>(N_abs) << "% done..." << endl;
    }
    
    double *absorption = new double[N_abs];
    // Not multiplied by C(0), so unit in 1/fs
    for (int k = 0 ; k < N_abs ; ++k)
    {
        absorption[k] = (k * k * dw * dw) * (w_bar_squared * gamma_1[k] / ((w_bar_squared - k * k * dw * dw - k * dw * gamma_2[k]) * (w_bar_squared - k * k * dw * dw - k * dw * gamma_2[k]) + (k * dw * gamma_1[k]) * (k * dw * gamma_1[k])));
        cout << "Calculating absorption spectrum " << static_cast<double>(k+1)*100/static_cast<double>(N_abs) << "% done..." << endl;
    }
    cout << endl;
    
    cout << "Exporting absorption data..." << endl;
    
    ofstream myfile;
    myfile.open ("absorbance.txt");
    
    for (int k =0 ; k < N_abs ; ++k)
    {
        myfile << k * dw << " " << absorption[k] << endl;
    }
    
    myfile.close();
     
    delete[] correlation;
    delete[] time;
    delete[] ker;
    delete[] gamma_1;
    delete[] gamma_2;
    delete[] absorption;
    
    return 0;
    
}
