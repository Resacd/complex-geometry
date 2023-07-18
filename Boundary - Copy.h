#pragma once
#include "Node.h"

class Boundary
{

private:
	std::vector<double> BoundaryVel_;

public:
	Boundary();
	~Boundary();
	double Density();
	std::vector<double> Velocity();
	void SetVel(double ux, double uy, double uz);
	double Apply_Equilibrium_Bc(std::shared_ptr<VelocitySet> velSet, int dir);
	void ApplyZouHeBC(std::shared_ptr<VelocitySet> velSet, Node &lat, std::vector<double>& velocity, int dir);
	void Apply_BounceBack(Node& lat, int dir, int opp_dir);
	void ApplyFreeSlipBC(std::shared_ptr<VelocitySet> velSet, Node& lat, int dir, int opp_dir);
	void ApplyPressureDrivenBC(std::shared_ptr<VelocitySet> velSet, Node& lat, int dir, double pressure);
	void calculateForceOnImmersedBoundary(std::shared_ptr<VelocitySet> velSet,
		std::vector<std::shared_ptr<Node>> &immersedBoundaryNodes);
	void ApplyIBM(std::shared_ptr<VelocitySet> velSet,
		const std::vector<std::shared_ptr<Node>>& immersedBoundaryNodes);
};

