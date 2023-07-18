#include "OutputFile.h"
#include "Struct.h"

OutputFile::OutputFile(int outputIdx)
    : m_outputIdx(outputIdx)
{
}

OutputFile::~OutputFile()
{
}

void OutputFile::WriteCSV(std::vector<std::vector<std::vector<std::shared_ptr<Node>>>>& lattice,
    std::shared_ptr<VelocitySet> velSet, int start_x, int start_y, int start_z, int end_x, int end_y, int end_z)
{
    std::unique_lock<std::mutex> lock(m_fileMutex);

    std::string file_name = std::to_string(m_outputIdx) + ".csv";
    std::ofstream velocityCSV(file_name, std::ios_base::trunc);

    if (!velocityCSV.is_open()) {
        return;
    }

    velocityCSV << "coordX,coordY,coordZ,velX,velY,velZ,pressure,velMag" << std::endl;

    std::vector<std::vector<double>> all_data; 

    for (int z = start_z; z <= end_z; ++z) {
        for (int y = start_y; y <= end_y; ++y) {
            for (int x = start_x; x <= end_x; ++x) {
                int global_x = x;
                int global_y = y;
                int global_z = z;
                auto& node = lattice[z - start_z][y - start_y][x - start_x];
                if (node != nullptr) {
                    double velMag = sqrt(node->Velocity(velSet)[0] * node->Velocity(velSet)[0] +
                        node->Velocity(velSet)[1] * node->Velocity(velSet)[1] + node->Velocity(velSet)[2] *
                        node->Velocity(velSet)[2]);

                    // Store data in the vector
                    all_data.push_back({ static_cast<double>(global_x),
                        static_cast<double>(global_y),
                        static_cast<double>(global_z),
                        node->Velocity(velSet)[0],
                        node->Velocity(velSet)[1],
                        node->Velocity(velSet)[2],
                        node->Pressure(velSet),
                        velMag });
                }
            }
        }
    }

    std::sort(all_data.begin(), all_data.end(),
        [](const std::vector<double>& a, const std::vector<double>& b) {
            return std::tie(a[0], a[1], a[2]) < std::tie(b[0], b[1], b[2]);
        });

    for (const auto& data : all_data) {
        velocityCSV << std::fixed << std::setprecision(6);

        velocityCSV << static_cast<int>(data[0]) << ","
            << static_cast<int>(data[1]) << ","
            << static_cast<int>(data[2]) << ","
            << data[3] << ","
            << data[4] << ","
            << data[5] << ","
            << data[6] << ","
            << data[7] << std::endl;
    }

    lock.unlock();
}