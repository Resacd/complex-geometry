#pragma once
#include <iostream>
#include "VelocitySet.h"
#include "include.h"
#include "Struct.h"

class Node
{



public:
	std::vector<double> m_distributions;
	std::vector<double> m_newDistributions;
	std::shared_ptr<Node> m_neighbours[_nLatNodes];
	int x_position;
	int y_position;
	int z_position;
	double forceX = 0.0;
	double forceY = 0.0;
	double forceZ = 0.0;
	NodeType nodeType;

	static const double PI;
	Node() = default;
	Node(int coord_x, int coord_y,int coord_z,NodeType type);
	~Node();

	double Equilibrium(std::shared_ptr<VelocitySet> velSet, int dir);
	virtual double Density();
	virtual std::vector<double> Velocity(std::shared_ptr<VelocitySet> velSet);
	double Pressure(std::shared_ptr<VelocitySet> velSet);
	void Stream(int dir);
	void BounceBack(int dir, int opp_dir);
	virtual void ApplyBC();
	void Initialize(std::shared_ptr<VelocitySet> velSet, int dir);

};

