#include "Simulation.h"


Simulation::Simulation(SubDomain& SubDomain,double tau,std::vector<double> &inVel,double pressure)
	: SubDomain_(SubDomain)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs_);
	SetRelaxation(tau);
    InletVelocity_ = inVel;
    pressure_ = pressure;
}

Simulation::~Simulation()
{
}


void Simulation::Collision()
{
    for (int k = 0; k < SubDomain_.my_Nz_; k++) {
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {

                    SubDomain_.lattice_[k][j][i]->m_newDistributions[dir] =
                        SubDomain_.lattice_[k][j][i]->m_distributions[dir] -
                        1 / m_relaxation * (SubDomain_.lattice_[k][j][i]->m_distributions[dir] -
                            SubDomain_.lattice_[k][j][i]->Equilibrium(SubDomain_.GetVelSet(), dir));

                }
            }
        }
    }
}


void Simulation::Streaming()
{
    for (int k = 0; k < SubDomain_.my_Nz_; k++) {
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 1; dir < _nLatNodes; dir++) {
                    if (SubDomain_.lattice_[k][j][i]->m_neighbours[dir] != nullptr) {
                        SubDomain_.lattice_[k][j][i]->Stream(dir);
                    }
                }
            }
        }
    }
}


void Simulation::printData(const std::vector<double>& data, const std::string& label) {
    std::cout << label << " data:"<<std::endl;
    for (const auto& value : data) {
        std::cout << value<<std::endl;
    }
    std::cout << std::endl;
}

void Simulation::Communicate() {
    int tag_xp = 0;
    int tag_xm = 1;
    int tag_yp = 2;
    int tag_ym = 3;
    int tag_zp = 4;
    int tag_zm = 5;

    MPI_Datatype Distribution;
    MPI_Type_contiguous(_nLatNodes, MPI_DOUBLE, &Distribution);
    MPI_Type_commit(&Distribution);

    MPI_Status status;

    if (SubDomain_.my_right_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    send_data.push_back(SubDomain_.lattice_[k][j][SubDomain_.my_Nx_ - 1]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size());

        //printData(send_data, "Send (right)");

        MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_right_, tag_xp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_right_, tag_xm,
            MPI_COMM_WORLD, &status);

        //printData(recv_data, "Receive (right)");

        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.lattice_[k][j][0]->m_distributions[dir] = recv_data[index++];
                }
            }
        }
    }


    if (SubDomain_.my_left_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    send_data.push_back(SubDomain_.lattice_[k][j][0]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size());
        //printData(send_data, "Send (left)");
        MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_left_, tag_xm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_left_, tag_xp,
            MPI_COMM_WORLD, &status);
        //printData(recv_data, "Receive (left)");
        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.lattice_[k][j][SubDomain_.my_Nx_ - 1]->m_distributions[dir] = recv_data[index++];
                }
            }
        }
    }


    if (SubDomain_.my_front_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    send_data.push_back(SubDomain_.lattice_[k][SubDomain_.my_Ny_ - 1][i]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        //printData(send_data, "Send (front)");
        MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_front_, tag_xm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_front_, tag_xp,
            MPI_COMM_WORLD, &status);
        //printData(recv_data, "Receive (front)");
        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.lattice_[k][0][i]->m_distributions[dir] = recv_data[index++];
                }
            }
        }
    }


    if (SubDomain_.my_back_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    send_data.push_back(SubDomain_.lattice_[k][0][i]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        //printData(send_data, "Send (back)");
        MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_back_, tag_xm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_back_, tag_xp,
            MPI_COMM_WORLD, &status);
        //printData(recv_data, "Receive (back)");
        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.lattice_[k][SubDomain_.my_Ny_ - 1][i]->m_distributions[dir] = recv_data[index++];
                }
            }
        }
    }


    if (SubDomain_.my_up_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    send_data.push_back(SubDomain_.lattice_[SubDomain_.my_Nz_ - 1][j][i]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        //printData(send_data, "Send (up)");
        MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_up_, tag_xm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_up_, tag_xp,
            MPI_COMM_WORLD, &status);
        //printData(recv_data, "Receive (up)");
        int index = 0;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.lattice_[0][j][i]->m_distributions[dir] = recv_data[index++];
                }
            }
        }
    }


    if (SubDomain_.my_down_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    send_data.push_back(SubDomain_.lattice_[0][j][i]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        //printData(send_data, "Send (down)");
        MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_down_, tag_xm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_down_, tag_xp,
            MPI_COMM_WORLD, &status);
        //printData(recv_data, "Receive (down)");
        int index = 0;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.lattice_[SubDomain_.my_Nz_ - 1][j][i]->m_distributions[dir] = recv_data[index++];
                }
            }
        }
    }

}


void Simulation::Iterate()
{
    for (int iter = 0; iter < 50; iter++) {
        Collision();
        Streaming();
        ApplyInletBc();
        ApplyOutletBc();
        ApplyFreeSlipBC();
        ApplyBounceBackBc();
        Communicate();
        WriteOutput(iter);
        if (myrank_ == 0) {
            std::cout << "Iteration" << " " << iter << " " << "has been performed successfully " << std::endl;
        }
    }
}

void Simulation::ApplyInletBc() {

    for (int lat = 0; lat < SubDomain_.InletBoundaryNode_.size(); lat++) {
        Node& node = *(SubDomain_.InletBoundaryNode_[lat]);
        for (int dir = 1; dir < _nLatNodes; dir++)
        {
            Boundary_.ApplyZouHeBC(SubDomain_.velSet_, node, InletVelocity_, dir);
        }
    }
}

void Simulation::ApplyOutletBc() {

    for (int lat = 0; lat < SubDomain_.OutletBoundaryNode_.size(); lat++) {
        Node& node = *(SubDomain_.OutletBoundaryNode_[lat]);
        for (int dir = 1; dir < _nLatNodes; dir++)
        {
            int opp_dir = SubDomain_.GetVelSet()->OppositeDirection(dir);
            Boundary_.ApplyPressureDrivenBC(SubDomain_.velSet_, node, dir, pressure_);
         
        }
    }

}

void Simulation::ApplyFreeSlipBC() {

    for (int lat = 0; lat < SubDomain_.FreeSlipBoundaryNode_.size(); lat++) {
        Node& node = *(SubDomain_.OutletBoundaryNode_[lat]);
        for (int dir = 1; dir < _nLatNodes; dir++)
        {
            int opp_dir = SubDomain_.GetVelSet()->OppositeDirection(dir);
            Boundary_.ApplyFreeSlipBC(SubDomain_.velSet_, node, dir, opp_dir);
        }
    }
}

void Simulation::ApplyBounceBackBc() {

    for (int lat = 0; lat < SubDomain_.ImmersedBoundaryNode_.size(); lat++) {
        for (int dir = 1; dir < _nLatNodes; dir++)
        {
            int opp_dir = SubDomain_.GetVelSet()->OppositeDirection(dir);
            SubDomain_.ImmersedBoundaryNode_[lat]->BounceBack(dir, opp_dir);
        }
    }

}

void Simulation::ApplyImmersedBc() {

    //Boundary_.ApplyIBM(SubDomain_.velSet_, SubDomain_.ImmersedBoundaryNode_);
    for (int lat = 0; lat < SubDomain_.ImmersedBoundaryNode_.size(); lat++) {
        for (int dir = 1; dir < _nLatNodes; dir++)
        {
            int opp_dir = SubDomain_.GetVelSet()->OppositeDirection(dir);
            SubDomain_.ImmersedBoundaryNode_[lat]->BounceBack(dir, opp_dir);
        }
    }
}

void Simulation::SetRelaxation(double tau)
{
	m_relaxation = tau;
}

void Simulation::WriteOutput(int iter)
{
	OutputFile file((iter));
	file.WriteCSV(SubDomain_.lattice_, SubDomain_.velSet_,
        SubDomain_.start_x_,SubDomain_.start_y_,SubDomain_.start_z_,
        SubDomain_.end_x_,SubDomain_.end_y_,SubDomain_.end_z_);
}