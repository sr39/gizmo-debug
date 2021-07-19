#pragma once


// Collisional Excitation

double g0(double x);

double g(double u, void* params);

double g_integral(double y2);

// Collisional Ionisation

double f(double u, void* params);

double f_integral(double y2);

// Recombination

double recomb_rate_f(double u, void* params);

double recomb_rate_integral(double y2);

double recomb_cool_f(double u, void* params);

double recomb_cool_integral(double y2);

