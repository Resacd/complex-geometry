#include "ComplexGeometry.h"
#include <queue>
#include <unordered_set>
#include <map>

Geometry::Geometry(std::vector<Triangle>& tri, double dx) {

    triangles = tri;
    spacing = dx;
    generatePointsInBoundingBox();
    Writer.writeBoxToFile(pointsBox, "box.vtk");
    labelPoints();
    Writer.writeLabeledPointsToFile(pointsBox, "labeled.vtk");

}

Geometry::~Geometry() {


}


void Geometry::generatePointsInBoundingBox() {

    Point minPoint = { std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
    Point maxPoint = { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };

    for (const auto& triangle : triangles) {
        minPoint.x = std::min({ minPoint.x, triangle.x1, triangle.x2, triangle.x3 });
        minPoint.y = std::min({ minPoint.y, triangle.y1, triangle.y2, triangle.y3 });
        minPoint.z = std::min({ minPoint.z, triangle.z1, triangle.z2, triangle.z3 });
        maxPoint.x = std::max({ maxPoint.x, triangle.x1, triangle.x2, triangle.x3 });
        maxPoint.y = std::max({ maxPoint.y, triangle.y1, triangle.y2, triangle.y3 });
        maxPoint.z = std::max({ maxPoint.z, triangle.z1, triangle.z2, triangle.z3 });
    }

    int nx = std::ceil((maxPoint.x - minPoint.x) / spacing);
    int ny = std::ceil((maxPoint.y - minPoint.y) / spacing);
    int nz = std::ceil((maxPoint.z - minPoint.z) / spacing);

    pointsBox.resize(nx + 1);
    for (int i = 0; i <= nx; ++i) {
        pointsBox[i].resize(ny + 1);
        for (int j = 0; j <= ny; ++j) {
            pointsBox[i][j].resize(nz + 1);
        }
    }

    double x, y, z;
    for (int k = 0; k <= nz; ++k) {
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                x = minPoint.x + i * spacing;
                y = minPoint.y + j * spacing;
                z = minPoint.z + k * spacing;
                Point p{ x, y, z };
                pointsBox[i][j][k] = p;
            }
        }
    }
}

Line Geometry::generateLine(const Point& startPoint, double lineLength) {
    Line line;
    line.start = startPoint;
    line.start.y = startPoint.y;
    line.end.x = startPoint.x;
    line.end.y = startPoint.y;
    line.end.z = startPoint.z + lineLength;
    return line;
}

bool Geometry::doLinesIntersect(const Line& line1, const Line& line2) {

    double x1 = line1.start.x;
    double y1 = line1.start.y;
    double z1 = line1.start.z;
    double x2 = line1.end.x;
    double y2 = line1.end.y;
    double z2 = line1.end.z;

    double x3 = line2.start.x;
    double y3 = line2.start.y;
    double z3 = line2.start.z;
    double x4 = line2.end.x;
    double y4 = line2.end.y;
    double z4 = line2.end.z;

    double ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1));
    double ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1));

    if (ua >= 0 && ua <= 1 && ub >= 0 && ub <= 1) {
        return true;
    }

    return false;
}

bool Geometry::isPointInsideTriangle(double x, double y, const Triangle& triangle) {

    double x1 = triangle.x1;
    double y1 = triangle.y1;
    double x2 = triangle.x2;
    double y2 = triangle.y2;
    double x3 = triangle.x3;
    double y3 = triangle.y3;

    double d1 = ((x - x2) * (y1 - y2) - (y - y2) * (x1 - x2));
    double d2 = ((x - x3) * (y2 - y3) - (y - y3) * (x2 - x3));
    double d3 = ((x - x1) * (y3 - y1) - (y - y1) * (x3 - x1));

    bool hasNegative = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool hasPositive = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(hasNegative && hasPositive);
}

int Geometry::countIntersectingTriangles(Line& line) {
    const double threshold = 0.0000001;
    int count = 0;
    for (const auto& triangle : triangles) {
        Line edge1, edge2, edge3;
        edge1.start = { triangle.x1, triangle.y1, triangle.z1 };
        edge1.end = { triangle.x2, triangle.y2, triangle.z2 };
        edge2.start = { triangle.x2, triangle.y2, triangle.z2 };
        edge2.end = { triangle.x3, triangle.y3, triangle.z3 };
        edge3.start = { triangle.x3, triangle.y3, triangle.z3 };
        edge3.end = { triangle.x1, triangle.y1, triangle.z1 };

        if (doLinesIntersect(line, edge1) || doLinesIntersect(line, edge2) || doLinesIntersect(line, edge3)) {

            count++;
            continue;
        }

        double px = line.start.x;
        double py = line.start.y;
        double pz = line.start.z;

        double nx = triangle.nx;
        double ny = triangle.ny;
        double nz = triangle.nz;

        double qx = triangle.x1;
        double qy = triangle.y1;
        double qz = triangle.z1;

        double t = ((qx - px) * nx + (qy - py) * ny + (qz - pz) * nz) /
            ((line.end.x - line.start.x) * nx + (line.end.y - line.start.y) * ny + (line.end.z - line.start.z) * nz);

        if (t >= 0 && t <= 1) {

            double ix = px + t * (line.end.x - line.start.x);
            double iy = py + t * (line.end.y - line.start.y);
            double iz = pz + t * (line.end.z - line.start.z);

            if (isPointInsideTriangle(ix, iy, triangle) && t > threshold && t < 1 - threshold) {
                count++;
            }
        }
    }
    return count;
}

void Geometry::labelPoints() {
#pragma omp parallel for collapse(3)
    for (int i = 0; i < pointsBox.size(); i++) {
        for (int j = 0; j < pointsBox[i].size(); j++) {
            for (int k = 0; k < pointsBox[i][j].size(); k++) {
                Point& point = pointsBox[i][j][k];

                Line line;
                line.start = point;
                line.end.x = point.x;
                line.end.y = point.y;
                line.end.z = point.z + 2.0;

                int intersectingCount = countIntersectingTriangles(line);
                if (intersectingCount % 2 == 1) {
                    point.Type = NodeType::Fluid;
                }
                else {
                    point.Type = NodeType::Solid;
                }
            }
        }
    }
}


//void Geometry::labelPoints() {
//    for (int i = 0; i < pointsBox.size(); i++) {
//        for (int j = 0; j < pointsBox[i].size(); j++) {
//            for (int k = 0; k < pointsBox[i][j].size(); k++) {
//                Point& point = pointsBox[i][j][k];
//
//                Line line;
//                line.start = point;
//                line.end.x = point.x; 
//                line.end.y = point.y; 
//                line.end.z = point.z + 2.0; 
//
//                int intersectingCount = countIntersectingTriangles(line);
//                if (intersectingCount % 2 == 1) {
//                    point.Type = NodeType::Fluid;
//                }
//                else {
//                    point.Type = NodeType::Solid;
//                }
//            }
//        }
//    }
//}


//std::vector<std::vector<std::vector<Point>>> Geometry::generatePointsInBoundingBox(const std::vector<Point>& points, double spacing) {
//    Point minPoint = { std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
//    Point maxPoint = { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };
//
//    for (const Point& p : points) {
//        minPoint.x = std::min(p.x, minPoint.x);
//        minPoint.y = std::min(p.y, minPoint.y);
//        minPoint.z = std::min(p.z, minPoint.z);
//        maxPoint.x = std::max(p.x, maxPoint.x);
//        maxPoint.y = std::max(p.y, maxPoint.y);
//        maxPoint.z = std::max(p.z, maxPoint.z);
//    }
//
//    int nx = std::ceil((maxPoint.x - minPoint.x) / spacing);
//    int ny = std::ceil((maxPoint.y - minPoint.y) / spacing);
//    int nz = std::ceil((maxPoint.z - minPoint.z) / spacing);
//
//    pointsBox.resize(nx + 1);
//    for (int i = 0; i < nx + 1; i++) {
//        pointsBox[i].resize(ny + 1);
//        for (int j = 0; j < ny + 1; j++) {
//            pointsBox[i][j].resize(nz + 1);
//        }
//    }
//    double x, y, z;
//    for (int k = 0; k < nz + 1; k++) {
//        for (int j = 0; j < ny + 1; j++) {
//            for (int i = 0; i < nx + 1; i++) {
//                x = minPoint.x + i * spacing;
//                y = minPoint.y + j * spacing;
//                z = minPoint.z + k * spacing;
//                Point p{ x, y, z };
//                pointsBox[i][j][k] = p;
//            }
//        }
//    }
//
//    return pointsBox;
//}

std::vector<std::tuple<int, int, int>> Geometry::findClosestPoints(std::vector<Point>& PointsOff) {
    std::vector<std::vector<std::vector<bool>>> pointMask(pointsBox.size(), std::vector<std::vector<bool>>(pointsBox[0].size(), std::vector<bool>(pointsBox[0][0].size(), false)));

    std::vector<std::tuple<int, int, int>> closestIndices;

    for (const auto& pOff : PointsOff) {
        float minDist = std::numeric_limits<float>::max();
        int iMin = -1, jMin = -1, kMin = -1;
        for (int i = 0; i < pointsBox.size(); ++i) {
            for (int j = 0; j < pointsBox[i].size(); ++j) {
                for (int k = 0; k < pointsBox[i][j].size(); ++k) {
                    if (!pointMask[i][j][k]) {
                        const auto& pBox = pointsBox[i][j][k];
                        float dist = std::sqrt(std::pow(pOff.x - pBox.x, 2) + std::pow(pOff.y - pBox.y, 2) + std::pow(pOff.z - pBox.z, 2));
                        if (dist < minDist) {
                            minDist = dist;
                            iMin = i;
                            jMin = j;
                            kMin = k;
                        }
                    }
                }
            }
        }

        if (iMin != -1 && jMin != -1 && kMin != -1) {
            closestIndices.emplace_back(iMin, jMin, kMin);
            pointMask[iMin][jMin][kMin] = true;
        }
        else {

        }
    }

    return closestIndices;
}

void Geometry::fillMissingPoints() {
    std::vector<std::tuple<int, int, int>> filledIndices;
    for (int i = 0; i < closestIndices.size(); ++i) {
        filledIndices.push_back(closestIndices[i]);
        if (i < closestIndices.size() - 1) {
            int currentIndex = std::get<2>(closestIndices[i]);
            int nextIndex = std::get<2>(closestIndices[i + 1]);
            if (nextIndex - currentIndex > 1) {
                // Fill in missing points with regular spacing
                int missingCount = nextIndex - currentIndex - 1;
                int currentI = std::get<0>(closestIndices[i]);
                int currentJ = std::get<1>(closestIndices[i]);
                for (int k = 1; k <= missingCount; ++k) {
                    filledIndices.push_back(std::make_tuple(currentI, currentJ, currentIndex + k));
                }
            }
        }
    }
    closestIndices = filledIndices;
}



std::map<std::pair<int, int>, std::vector<int>> Geometry::extractThirdElements() {
    std::map<std::pair<int, int>, std::vector<int>> thirdElements;

    for (const auto& tuple : closestIndices) {
        int i, j, k;
        std::tie(i, j, k) = tuple;
        thirdElements[std::make_pair(i, j)].push_back(k);
    }

    return thirdElements;
}

std::vector<int> Geometry::getKVector(const std::map<std::pair<int, int>, std::vector<int>>& thirdElements,
    int targetI, int targetJ) {
    std::pair<int, int> key = std::make_pair(targetI, targetJ);
    auto it = thirdElements.find(key);
    if (it != thirdElements.end()) {
        return it->second;
    }
    else {
        return {}; 
    }
}

std::vector<std::vector<std::vector<bool>>> Geometry::labelPointsInComplexGeometry(const std::map<std::pair<int, int>, std::vector<int>>& thirdElements) {
    std::vector<std::vector<std::vector<bool>>> pointLabels(pointsBox.size(), std::vector<std::vector<bool>>(pointsBox[0].size(), std::vector<bool>(pointsBox[0][0].size(), false)));

    for (const auto& entry : thirdElements) {
        const std::pair<int, int>& key = entry.first;
        const std::vector<int>& kValues = entry.second;

        if (kValues.size() == 1) {
            int i = key.first;
            int j = key.second;
            int k = kValues[0];
            if (k != -1) {  // Skip labeling if the index is a dummy index
                pointLabels[i][j][k] = true;
            }
        }
        else if (kValues.size() > 1) {
            std::vector<int> sortedKValues = kValues;  // Create a copy of kValues
            std::sort(sortedKValues.begin(), sortedKValues.end());
            for (size_t idx = 0; idx < sortedKValues.size() - 1; idx += 2) {
                int i = key.first;
                int j = key.second;
                int k1 = sortedKValues[idx];
                int k2 = sortedKValues[idx + 1];
                if (k1 != -1 && k2 != -1) {  // Skip labeling if any of the indices are dummy indices
                    for (int k = k1; k <= k2; ++k) {
                        pointLabels[i][j][k] = true;
                    }
                }
            }
        }
    }

    return pointLabels;
}


//std::vector<std::vector<std::vector<bool>>> Geometry::getPointLabels(std::vector<Point>& PointsOff, double spacing) {
//    pointsBox = generatePointsInBoundingBox(PointsOff, spacing);
//    closestIndices = findClosestPoints(PointsOff);
//    std::map<std::pair<int, int>, std::vector<int>> thirdElements = extractThirdElements();
//    std::vector<std::vector<std::vector<bool>>> pointLabels = labelPointsInComplexGeometry(thirdElements);
//
//    return pointLabels;
//}


void Geometry::writePointsToFile( const std::vector<std::vector<std::vector<bool>>>& pointLabels,const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cout << "Failed to open file: " << filename << std::endl;
        return;
    }

    file << "x,y,z,label" << std::endl;

    for (int k = 0; k < pointsBox[0][0].size(); ++k) {
        for (int j = 0; j < pointsBox[0].size(); ++j) {
            for (int i = 0; i < pointsBox.size(); ++i) {
                const Point& point = pointsBox[i][j][k];
                bool label = pointLabels[i][j][k];

                if (label) {

                    file << point.x << "," << point.y << "," << point.z << ",true" << std::endl;
                }
            }
        }
    }
    file.close();
}
