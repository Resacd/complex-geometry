#include "Subdomain.h"

SubDomain::SubDomain(std::shared_ptr<VelocitySet> velSet, std::vector<int> &gridSize) {

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs_);
    velSet_ = velSet;
    gridSize_ = gridSize;
    Create_cart_comm();
    CreateNodes();
    Connectivity();
    Initialize();
    GetInletBoundary();
    GetOutletBoundary();
    //GetBounceBackBoundary();
    GetImmersedBoundary();
    PrintBoundaryNodes();
    MPI_Comm_free(&cart_comm_);
   
}

SubDomain::~SubDomain() {

}

void SubDomain::Create_cart_comm() {
    
    int const ndims = 3;
    dims_ = { 0, 0, 0 };
    int periods[ndims] = { 0, 0, 0 };
    int reorder = 1;
    coords_ = { 0, 0, 0 };
    int remain_dims[3] = { 1, 1, 1 };
    int Nx = gridSize_[0];
    int Ny = gridSize_[1];
    int Nz = gridSize_[2];
    MPI_Dims_create(num_procs_, ndims, dims_.data());
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims_.data(), periods, reorder, &cart_comm_);
    MPI_Cart_shift(cart_comm_, 0, 1, &my_left_, &my_right_);
    MPI_Cart_shift(cart_comm_, 1, 1, &my_down_, &my_up_);
    MPI_Cart_shift(cart_comm_, 2, 1, &my_front_, &my_back_);
    MPI_Comm_rank(cart_comm_, &myrank_);
    MPI_Cart_coords(cart_comm_, myrank_, ndims, coords_.data());
    int cart_coords[ndims];
    MPI_Cart_get(cart_comm_, ndims, dims_.data(), periods, cart_coords);

    int halo_size = 1;
    my_Nx_ = Nx / dims_[0];
    my_Ny_ = Ny / dims_[1];
    my_Nz_ = Nz / dims_[2];
    start_x_ = coords_[0] * my_Nx_;
    start_y_ = coords_[1] * my_Ny_;
    start_z_ = coords_[2] * my_Nz_;
    end_x_ = start_x_ + my_Nx_ - 1;
    end_y_ = start_y_ + my_Ny_ - 1;
    end_z_ = start_z_ + my_Nz_ - 1;

    if (Nx % dims_[0] > 0 && coords_[0] == dims_[0] - 1) {
        my_Nx_ += Nx % dims_[0];
    }
    if (Ny % dims_[1] > 0 && coords_[1] == dims_[1] - 1) {
        my_Ny_ += Ny % dims_[1];
    }
    if (Nz % dims_[2] > 0 && coords_[2] == dims_[2] - 1) {
        my_Nz_ += Nz % dims_[2];
    }

    if (dims_[0] > 1) {
        if (coords_[0] == 0) {
            my_Nx_ += halo_size;
            end_x_ += halo_size;
        }
        else if (coords_[0] == dims_[0] - 1) {
            my_Nx_ += halo_size;
            start_x_ += 2 * halo_size * coords_[0] - halo_size;
            end_x_ = start_x_ + my_Nx_ - 1;
        }
        else {
            my_Nx_ += 2 * halo_size;
            start_x_ += coords_[0] * 2 * halo_size - halo_size;
            end_x_ = start_x_ + my_Nx_ - 1;
        }
    }

    if (dims_[1] > 1) {
        if (coords_[1] == 0) {
            my_Ny_ += halo_size;
            end_y_ += halo_size;
        }
        else if (coords_[1] == dims_[1] - 1) {
            my_Ny_ += halo_size;
            start_y_ += 2 * halo_size * coords_[1] - halo_size;
            end_y_ = start_y_ + my_Ny_ - 1;
        }
        else {
            my_Ny_ += 2 * halo_size;
            start_y_ += coords_[1] * 2 * halo_size - halo_size;
            end_y_ = start_y_ + my_Ny_ - 1;
        }
    }

    if (dims_[2] > 1) {
        if (coords_[2] == 0) {
            my_Nz_ += halo_size;
            end_z_ += halo_size;
        }
        else if (coords_[2] == dims_[2] - 1) {
            my_Nz_ += halo_size;
            start_z_ += 2 * halo_size * coords_[2] - halo_size;
            end_z_ = start_z_ + my_Nz_ - 1;
        }
        else {
            my_Nz_ += 2 * halo_size;
            start_z_ += coords_[2] * 2 * halo_size - halo_size;
            end_z_ = start_z_ + my_Nz_ - 1;
        }
    }

    NewgridSize_.resize(3);

    NewgridSize_[0] = (dims_[0] - 1) * 2 * halo_size + gridSize_[0];
    NewgridSize_[1] = (dims_[1] - 1) * 2 * halo_size + gridSize_[1];
    NewgridSize_[2] = (dims_[2] - 1) * 2 * halo_size + gridSize_[2];


    for (int p = 0; p < num_procs_; p++) {
        if (myrank_ == p) {

            printf("Processor %d:\n", myrank_);
            printf("my_Nx_ = %d, my_Ny_ = %d, my_Nz_ = %d\n", my_Nx_, my_Ny_, my_Nz_);
            printf("NewgridSize_[0] = %d, NewgridSize_[1] = %d, NewgridSize_[2] = %d\n",
                NewgridSize_[0], NewgridSize_[1], NewgridSize_[2]);
            printf("start_x_ = %d, end_x_ = %d\n", start_x_, end_x_);
            printf("start_y_ = %d, end_y_ = %d\n", start_y_, end_y_);
            printf("start_z_ = %d, end_z_ = %d\n", start_z_, end_z_);

            for (int k = 0; k < my_Nz_; k++) {
                for (int j = 0; j < my_Ny_; j++) {
                    for (int i = 0; i < my_Nx_; i++) {
                        int global_i = i + start_x_;
                        int global_j = j + start_y_;
                        int global_k = k + start_z_;
                        std::tuple<int, int, int> t1 = std::make_tuple(i, j, k);
                        std::tuple<int, int, int> t2 = std::make_tuple(global_i, global_j, global_k);
                        std::pair<std::tuple<int, int, int>, std::tuple<int, int, int>> p = std::make_pair(t1, t2);
                        local2global_.push_back(p);
                        printf("Process %d [%d,%d,%d]: local index = (%d, %d, %d), global index = (%d, %d, %d)\n",
                            myrank_, coords_[0], coords_[1], coords_[2], i, j, k,
                            global_i, global_j, global_k);
                    }
                }
            }
            fflush(stdout);
        }
        MPI_Barrier(cart_comm_); 
    }
}

bool SubDomain::IsInSubDomain(int coord_x, int coord_y, int coord_z) {

    if (coord_x >= 0 && coord_x < my_Nx_
        && coord_y >= 0 && coord_y < my_Ny_
        && coord_z >= 0 && coord_z < my_Nz_)
        return true;
    else
        return false;

}

bool SubDomain::IsSolidInSubDomain(int coord_x, int coord_y, int coord_z)
{
    int global_coordX = coord_x + start_x_;
    int global_coordY = coord_y + start_y_;
    int a = pow(global_coordX - gridSize_[4], 2);
    int b = pow(global_coordY - gridSize_[5], 2);
    int tmp = pow(a, 2) + pow(b, 2);

    if (IsInSubDomain(coord_x, coord_y, coord_z) && tmp < pow(gridSize_[3], 2))
        return true;
    else
        return false;
}


bool SubDomain::IsFluidInSubDomain(int coord_x, int coord_y, int coord_z)
{
    int global_coordX = coord_x + start_x_;
    int global_coordY = coord_y + start_y_;
    int a = pow(global_coordX - gridSize_[4], 2);
    int b = pow(global_coordY - gridSize_[5], 2);
    int tmp = pow(a, 2) + pow(b, 2);

    if (IsInSubDomain(coord_x, coord_y, coord_z) && tmp >= pow(gridSize_[3], 2))
        return true;

    else
        return false;
}

void SubDomain::CreateNodes()
{
    std::vector<std::shared_ptr<Node>> tmp1;
    std::vector<std::vector<std::shared_ptr<Node>>> tmp2;

    for (int k = 0; k < my_Nz_; k++)
    {
        tmp2.resize(0);
        for (int j = 0; j < my_Ny_; j++)
        {
            tmp1.resize(0);
            for (int i = 0; i < my_Nx_; i++)
            {
                if (IsSolidInSubDomain(i, j, k)) {
                    tmp1.emplace_back(std::make_shared<Node>(i, j, k, NodeType::Solid));
                }
                else if (IsFluidInSubDomain(i, j, k)) {
                    tmp1.emplace_back(std::make_shared<Node>(i, j, k, NodeType::Fluid));
                }
                else {
                    tmp1.emplace_back(std::shared_ptr<Node>());
                }
            }
            tmp2.push_back(tmp1);
        }
        lattice_.push_back(tmp2);
    }
}

void SubDomain::Connectivity()
{
    for (int k = 0; k < my_Nz_; k++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int i = 0; i < my_Nx_; i++) {
                lattice_[k][j][i]->m_neighbours[0] = nullptr;
                for (int dir = 1; dir < _nLatNodes; dir++) {
                    int x_neighbor = lattice_[k][j][i]->x_position + velSet_->GetDirection(dir)[0];
                    int y_neighbor = lattice_[k][j][i]->y_position + velSet_->GetDirection(dir)[1];
                    int z_neighbor = lattice_[k][j][i]->z_position + velSet_->GetDirection(dir)[2];
                    if (IsInSubDomain(x_neighbor, y_neighbor, z_neighbor)) {
                        lattice_[k][j][i]->m_neighbours[dir] = lattice_[z_neighbor][y_neighbor][x_neighbor];
                    }
                    else
                        lattice_[k][j][i]->m_neighbours[dir] = nullptr;
                }
            }
        }
    }
}

void SubDomain::WriteGeometryToFile(const std::string& filename)
{
    std::ofstream file(filename);
    if (!file)
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    file << "# vtk DataFile Version 2.0\n";
    file << "Geometry\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << my_Nx_ << " " << my_Ny_ << " " << my_Nz_ << "\n";
    file << "ORIGIN 0 0 0\n";
    file << "SPACING 1 1 1\n";
    file << "POINT_DATA " << my_Nx_ * my_Ny_ * my_Nz_ << "\n";

    file << "SCALARS SolidNode unsigned_char 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int k = 0; k < my_Nz_; k++)
    {
        for (int j = 0; j < my_Ny_; j++)
        {
            for (int i = 0; i < my_Nx_; i++)
            {
                if (IsSolidInSubDomain(i, j, k))
                    file << "255\n";
                else
                    file << "0\n";
            }
        }
    }

    file << "SCALARS FluidNode unsigned_char 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int k = 0; k < my_Nz_; k++)
    {
        for (int j = 0; j < my_Ny_; j++)
        {
            for (int i = 0; i < my_Nx_; i++)
            {
                if (IsFluidInSubDomain(i, j, k))
                    file << "255\n";
                else
                    file << "0\n";
            }
        }
    }

    file.close();
}


void SubDomain::Initialize()
{
    for (int k = 0; k < my_Nz_; k++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int i = 0; i < my_Nx_; i++) {

                for (int dir = 0; dir < _nLatNodes; dir++) {
                  lattice_[k][j][i]->Initialize(velSet_, dir);
                }
            }
        }
    }
}


bool SubDomain::IfProcContainsInletBoundary() {

    if (coords_[0] == 0)
        return 1;
    else
        return 0;
}

bool SubDomain::IfProcContainsOutletBoundary() {

    if (coords_[0] == dims_[0]-1)
        return 1;
    else
        return 0;
}

bool SubDomain::IfProcContainsBounceBackBoundary() {

    if (coords_[1]==0 || coords_[1]==dims_[1]-1 || coords_[2] == 0 || coords_[2] == dims_[2] - 1)
        return 1;
    else
        return 0;
}


void SubDomain::GetInletBoundary() {
    if (IfProcContainsInletBoundary()) {

        for (int j = start_y_; j <= end_y_; j++) {
            for (int k = start_z_; k <= end_z_; k++) {
                if (k != NewgridSize_[2] - 1) {

                    InletBoundaryNode_.push_back(lattice_[k - start_z_][j - start_y_][start_x_ - start_x_]);
                }
            }
        }
    }
}

void SubDomain::GetOutletBoundary() {
    if (IfProcContainsOutletBoundary()) {

        for (int j = start_y_; j <= end_y_; j++) {
            for (int k = start_z_; k <= end_z_; k++) {
                if (k != NewgridSize_[2] - 1) {

                    OutletBoundaryNode_.push_back(lattice_[k - start_z_][j - start_y_][end_x_ - start_x_]);
                }
            }
        }
    }
}


void SubDomain::GetFreeSlipBoundary() {
    if (IfProcContainsBounceBackBoundary()) {

        //front side
        if (coords_[1] == 0) {
            for (int i = start_x_; i <= end_x_; i++) {
                for (int k = start_z_; k <= end_z_; k++) {
                    if (k != NewgridSize_[2] - 1) {

                        FreeSlipBoundaryNode_.push_back(lattice_[k - start_z_][start_y_ - start_y_][i - start_x_]);
                    }
                }
            }
        }
        //back side 
        if (coords_[1] == dims_[1] - 1) {
            for (int i = start_x_; i <= end_x_; i++) {
                for (int k = start_z_; k <= end_z_; k++) {
                    if (k != NewgridSize_[2] - 1) {

                        FreeSlipBoundaryNode_.push_back(lattice_[k - start_z_][end_y_ - start_y_][i - start_x_]);
                    }
                }
            }
        }
        //down side
        if (coords_[2] == 0) {
            for (int i = start_x_; i <= end_x_; i++) {
                for (int j = start_y_; j <= end_y_; j++) {

                    FreeSlipBoundaryNode_.push_back(lattice_[start_z_ - start_z_][j - start_y_][i - start_x_]);
                }
            }
        }

        //up side
        if (coords_[2] == dims_[2] - 1) {
            for (int i = start_x_; i <= end_x_; i++) {
                for (int j = start_y_; j <= end_y_; j++) {
                    if (end_z_ != NewgridSize_[2] - 1) {

                        FreeSlipBoundaryNode_.push_back(lattice_[end_z_-start_z_][j-start_y_][i-start_x_]);
                    }
                }
            }
        }
    }
}

void SubDomain::GetImmersedBoundary() {

    for (int k = 0; k < my_Nz_; k++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int i = 0; i < my_Nx_; i++) {

                if (lattice_[k][j][i]->nodeType == NodeType::Solid) {

                    ImmersedBoundaryNode_.push_back(lattice_[k][j][i]);

                }
            }
        }
    }
}

void SubDomain::PrintBoundaryNodes() {


    std::cout << "Inlet Boundary Nodes (Processor " << myrank_ << "):\n";
    for (const auto& node : InletBoundaryNode_) {
        std::cout << "(" << node->x_position << ", " << node->y_position << ", " << node->z_position << ")\n";
    }
    std::cout << "\n";

    std::cout << "Outlet Boundary Nodes (Processor " << myrank_ << "):\n";
    for (const auto& node : OutletBoundaryNode_) {
        std::cout << "(" << node->x_position << ", " << node->y_position << ", " << node->z_position << ")\n";
    }
    std::cout << "\n";

    std::cout << "BounceBack Boundary Nodes (Processor " << myrank_ << "):\n";
    for (const auto& node : BounceBackBoundaryNode_) {
        std::cout << "(" << node->x_position << ", " << node->y_position << ", " << node->z_position << ")\n";
    }
    std::cout << "\n";

    std::cout << "Immersed Boundary Nodes (Processor " << myrank_ << "):\n";
    for (const auto& node : ImmersedBoundaryNode_) {
        std::cout << "(" << node->x_position << ", " << node->y_position << ", " << node->z_position << ")\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);
}




const std::shared_ptr<VelocitySet> SubDomain::GetVelSet()
{
    return velSet_;
}



