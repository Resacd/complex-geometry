#include "Boundary.h"
#include <cmath>
#include <iostream>

Boundary::Boundary() {



}

Boundary::~Boundary() {


}

double Boundary::Density() {

	return 1;
}

void Boundary::SetVel(double ux, double uy, double uz)
{
	BoundaryVel_.push_back(ux);
	BoundaryVel_.push_back(uy);
	BoundaryVel_.push_back(uz);
}

std::vector<double> Boundary::Velocity() {

	return BoundaryVel_;
}

double Boundary::Apply_Equilibrium_Bc(std::shared_ptr<VelocitySet> velSet, int dir) {

	Velocity();
	double du = BoundaryVel_[0] * velSet->GetDirection(dir)[0] +
                BoundaryVel_[1] * velSet->GetDirection(dir)[1] +
                BoundaryVel_[2] * velSet->GetDirection(dir)[2];
	double u_sqr = BoundaryVel_[0] * BoundaryVel_[0] + 
                   BoundaryVel_[1] * BoundaryVel_[1] + 
                   BoundaryVel_[2] * BoundaryVel_[2];

	return velSet->GetWeight(dir) * Density() * (1 + 3 * du + 9 / 2 * du * du - 3 / 2 * u_sqr);

}

void Boundary::ApplyZouHeBC(std::shared_ptr<VelocitySet> velSet, Node& lat, std::vector<double>& velocity, int dir)
{
    double du = velocity[0] * velSet->GetDirection(dir)[0] +
        velocity[1] * velSet->GetDirection(dir)[1] +
        velocity[2] * velSet->GetDirection(dir)[2];

    double u_sqr = velocity[0] * velocity[0] +
        velocity[1] * velocity[1] +
        velocity[2] * velocity[2];

    double equilibriumDist = velSet->GetWeight(dir) * Density() * (1 + 3 * du + 9.0 / 2.0 * du * du - 3.0 / 2.0 * u_sqr);
    double distributionDiff = equilibriumDist - lat.m_distributions[dir];
    double bounceBackDist = lat.m_distributions[dir] - distributionDiff;

    lat.m_distributions[dir] = bounceBackDist;

    int opp_dir = velSet->OppositeDirection(dir);
    if (lat.m_neighbours[dir] != nullptr) {
        int opp_dir = velSet->OppositeDirection(dir);
        lat.m_neighbours[dir]->m_distributions[opp_dir] = bounceBackDist;
    }
    else {

    }
}


void Boundary::Apply_BounceBack(Node& lat, int dir,int opp_dir)
{

    lat.m_distributions[opp_dir] = lat.m_newDistributions[dir];
}

void Boundary::ApplyFreeSlipBC(std::shared_ptr<VelocitySet> velSet, Node& lat, int dir, int opp_dir)
{
    if (lat.m_neighbours[dir] != nullptr) {
        std::swap(lat.m_distributions[dir], lat.m_neighbours[dir]->m_distributions[opp_dir]);
    }
    else {

    }
}


void Boundary::ApplyPressureDrivenBC(std::shared_ptr<VelocitySet> velSet, Node& lat, int dir, double pressure)
{
    int opp_dir = velSet->OppositeDirection(dir);
    double equilibriumDist = pressure / (velSet->GetWeight(dir) * lat.Density());

    double distributionDiff = equilibriumDist - lat.m_distributions[dir];
    double bounceBackDist = equilibriumDist - distributionDiff;

    lat.m_distributions[dir] += distributionDiff;

    if (lat.m_neighbours[dir] != nullptr) {
        lat.m_neighbours[dir]->m_distributions[opp_dir] = bounceBackDist;
    }
}


void Boundary::calculateForceOnImmersedBoundary(std::shared_ptr<VelocitySet> velSet, 
    std::vector<std::shared_ptr<Node>>& immersedBoundaryNodes) {

    for (const auto& node : immersedBoundaryNodes) {
        double forceX = 0.0;
        double forceY = 0.0;
        double forceZ = 0.0;

        for (int dir = 0; dir < _nLatNodes; dir++) {
            if (node->m_neighbours[dir]->nodeType != NodeType::Solid) {
                forceX += node->m_neighbours[dir]->Density() * 
                    node->m_neighbours[dir]->Velocity(velSet)[0] * velSet->GetDirection(dir)[0];

                forceY += node->m_neighbours[dir]->Density() *
                    node->m_neighbours[dir]->Velocity(velSet)[1] * velSet->GetDirection(dir)[1];

                forceZ += node->m_neighbours[dir]->Density() * 
                    node->m_neighbours[dir]->Velocity(velSet)[2] * velSet->GetDirection(dir)[2];
            }
        }

        node->forceX = forceX;
        node->forceY = forceY;
        node->forceZ = forceZ;
    }
}

void Boundary::ApplyIBM(std::shared_ptr<VelocitySet> velSet,
    const std::vector<std::shared_ptr<Node>>& immersedBoundaryNodes) {

    for (const auto& node : immersedBoundaryNodes) {
        for (int dir = 0; dir < _nLatNodes; dir++) {
            if (node->m_neighbours[dir]->nodeType == NodeType::Fluid) {
                node->m_neighbours[dir]->m_distributions[dir] += 
                    velSet->GetWeight(dir) * (node->forceX * velSet->GetDirection(dir)[0] +
                        node->forceY * velSet->GetDirection(dir)[1] + node->forceZ * velSet->GetDirection(dir)[2]);
                node->m_neighbours[dir]->m_newDistributions[dir] = node->m_neighbours[dir]->m_distributions[dir];
            }
        }
    }
}

