#pragma once
#include <vector>
#include <iostream>
#include "Struct.h"
int Nx = 25;
int Ny = 16;
int Nz = 16;
int Nr = 4;
int CtrX = 8;
int CtrY = 8;
int CtrZ = 8;
double U = 0.1;
double Re = 100;
double niu = U * Nx / Re;
double tau = 3.0 * niu + 0.5;

std::vector<int> gridSize{ Nx, Ny, Nz,Nr,CtrX,CtrY,CtrZ };
std::vector<double> inVel{ 1.0,0.0,0.0 };
double pressure = 1.0;