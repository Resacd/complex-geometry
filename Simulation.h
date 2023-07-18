#pragma once

#include "VelocitySet.h"
#include "Subdomain.h"
#include "OutputFile.h"
#include "Boundary.h"

class Simulation
{

public:
	int myrank_, num_procs_;
	SubDomain SubDomain_;
	Boundary Boundary_;
	double pressure_;
	double m_relaxation;
	std::vector<double> InletVelocity_;
	Simulation(SubDomain& SubDomain, double tau, std::vector<double>& inVel, double pressure);
	~Simulation();
	void Collision();
	void Streaming();
	void ApplyInletBc();
	void ApplyOutletBc();
	void ApplyFreeSlipBC();
	void ApplyBounceBackBc();
	void ApplyImmersedBc();
	void Communicate();
	void printData(const std::vector<double>& data, const std::string& label);
	void Iterate();
	void SetRelaxation(double tau);
	void WriteOutput(int iter);
};

