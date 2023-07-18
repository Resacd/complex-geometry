#pragma once
#ifndef _MESH_GENERATOR_H_
#define _MESH_GENERATOR_H_

#include "Struct.h"
#include <tuple>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <fstream>
#include "Struct.h"
#include <unordered_set>
#include <omp.h>
#include <map>
#include "GeometryWriter.h"

class Geometry {

public:

	Writer Writer;
	std::vector<std::vector<std::vector<Point>>> pointsBox;
	std::vector<std::tuple<int, int, int>> closestIndices;
	std::vector<Triangle> triangles;
	double spacing;
	Geometry(std::vector<Triangle>& tri, double dx);
	~Geometry();
	void generatePointsInBoundingBox();
	Line generateLine(const Point& startPoint, double lineLength);
	bool doLinesIntersect(const Line& line1, const Line& line2);
	bool isPointInsideTriangle(double x, double y, const Triangle& triangle);
	int countIntersectingTriangles(Line& line);
	void labelPoints();

	/*std::vector<std::vector<std::vector<Point>>> generatePointsInBoundingBox(const std::vector<Point>&
	points, double spacing);*/
	std::vector<std::tuple<int, int, int>> findClosestPoints(std::vector<Point>& PointsOff);
	void fillMissingPoints();
	std::map<std::pair<int, int>, std::vector<int>> extractThirdElements();
	std::vector<int> getKVector(const std::map<std::pair<int, int>, std::vector<int>>& thirdElements,
		int targetI, int targetJ);
	std::vector<std::vector<std::vector<bool>>> labelPointsInComplexGeometry(const std::map<std::pair<int, int>, std::vector<int>>& thirdElements);
	std::vector<std::vector<std::vector<bool>>> getPointLabels(std::vector<Point>& PointsOff,double spacing);
	void writePointsToFile(const std::vector<std::vector<std::vector<bool>>>& pointLabels, const std::string& filename);
};

#endif