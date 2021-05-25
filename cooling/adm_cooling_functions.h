#pragma once


// Collisional Excitation

double g0(double x);

double g(double u, void* params);

double g_integral(double y2);

// Collisional Ionisation

double f(double u, void* params);

double f_integral(double y2);
