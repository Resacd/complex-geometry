#include <iostream>
#include "LatticeModel.h"
#include "Simulation.h"
#include "OutputFile.h"
#include <mpi.h>
#include "Subdomain.h"
#include "Subdomain.h"
#include "Input.h"
#include "GeometryReader.h"
#include "GeometryWriter.h"
#include "ComplexGeometry.h"



int main(int argc, char* argv[])
{
	auto start = std::chrono::high_resolution_clock::now();
	int myrank, num_procs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	omp_set_num_threads(omp_get_max_threads());
	int numThreads = omp_get_max_threads();
	omp_set_num_threads(numThreads);
	std::cout << "Number of OpenMP threads: " << numThreads << std::endl;
	std::shared_ptr<LatticeModel> velSet = std::make_shared<LatticeModel>();
	Reader Reader;
	Writer Writer;
	double spacing = 0.1;

	std::vector<Triangle> tris = Reader.readSTLASCII("Vessel.stl");
	//std::vector<Triangle> tris = Reader.readSTLASCII("CylinderTest.stl");
	Geometry Geometry(tris,spacing);

	MPI_Finalize();
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "Elapsed time: " << elapsed.count() << " seconds." << std::endl;
	return 0;
}




