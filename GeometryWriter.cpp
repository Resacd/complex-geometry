#include "GeometryWriter.h"


Writer::Writer() {


}

Writer::~Writer() {


}

void Writer::writePointoffToFile(const std::vector<Point>& points, const std::string& filename) {
    // Open the output file for writing
    std::ofstream outFile(filename);

    // Write the header information for the VTK file format
    outFile << "# vtk DataFile Version 3.0" << std::endl;
    outFile << "Points data" << std::endl;
    outFile << "ASCII" << std::endl;
    outFile << "DATASET POLYDATA" << std::endl;

    // Write the points data
    outFile << "POINTS " << points.size() << " double" << std::endl;
    for (const auto& point : points) {
        outFile << point.x << " " << point.y << " " << point.z << std::endl;
    }

    // Close the output file
    outFile.close();
}

void Writer::writeBoxToFile(const std::vector<std::vector<std::vector<Point>>>& points, const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Unable to open file " << filename << " for writing" << std::endl;
        return;
    }

    int nx = points.size();
    int ny = points[0].size();
    int nz = points[0][0].size();

    int numPoints = nx * ny * nz;

    file << "# vtk DataFile Version 3.0\n";
    file << "Points\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << numPoints << " float\n";

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                file << (points[i][j][k].x) << " " << (points[i][j][k].y) << " " << (points[i][j][k].z) << "\n";
            }
        }
    }


    file.close();
    std::cout << "Points written to file " << filename << std::endl;
}

void Writer::writeClosestPointsToVTK(const std::string& filename, 
    std::vector<std::vector<std::vector<Point>>>& pointsBox,
    const std::vector<std::tuple<int, int, int>>& closestIndices) {
    std::ofstream outfile(filename);

    outfile << "# vtk DataFile Version 4.0\n";
    outfile << "vtk output\n";
    outfile << "ASCII\n";
    outfile << "DATASET UNSTRUCTURED_GRID\n";

    // Write the POINTS section
    outfile << "POINTS " << closestIndices.size() << " float\n";
    for (const auto& index : closestIndices) {
        const auto& point = pointsBox[std::get<0>(index)][std::get<1>(index)][std::get<2>(index)];
        outfile << point.x << " " << point.y << " " << point.z << "\n";
    }

    // Write the CELLS section
    outfile << "CELLS " << closestIndices.size() << " " << closestIndices.size() * 2 << "\n";
    for (int i = 0; i < closestIndices.size(); ++i) {
        outfile << "1 " << i << "\n";
    }

    // Write the CELL_TYPES section
    outfile << "CELL_TYPES " << closestIndices.size() << "\n";
    for (int i = 0; i < closestIndices.size(); ++i) {
        outfile << "1\n";
    }

    outfile.close();
}


void Writer::writeLabeledPointsToFile(const std::vector<std::vector<std::vector<Point>>>& pointsBox, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Write the header information
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "Points Inside Complex Geometry" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Collect the fluid points
    std::vector<Point> fluidPoints;
    for (const auto& plane : pointsBox) {
        for (const auto& row : plane) {
            for (const auto& point : row) {
                if (point.Type == NodeType::Fluid) {
                    fluidPoints.push_back(point);
                }
            }
        }
    }

    // Write the POINTS section
    file << "POINTS " << fluidPoints.size() << " double" << std::endl;
    for (const auto& point : fluidPoints) {
        file << point.x << " " << point.y << " " << point.z << std::endl;
    }

    // Write the CELLS section
    int numCells = fluidPoints.size();
    file << "CELLS " << numCells << " " << numCells * 2 << std::endl;
    for (int i = 0; i < numCells; ++i) {
        file << "1 " << i << std::endl;
    }

    // Write the CELL_TYPES section
    file << "CELL_TYPES " << numCells << std::endl;
    for (int i = 0; i < numCells; ++i) {
        file << "1" << std::endl;
    }

    file.close();
}
