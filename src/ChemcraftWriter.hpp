#pragma once

#include "fstream"
#include <string>
#include "Structure.hpp"

class ChemcraftWriter {
private:
    std::ofstream _fp;
    size_t _clusterIdx = 0;
    size_t _singleIdx = 0;
    size_t _tempStructIdx = 0;

    void startJob(const std::string& type, const std::string& comment);

public:
    explicit ChemcraftWriter(const std::string& filepath);
    ~ChemcraftWriter();

    ChemcraftWriter(const ChemcraftWriter&) = delete;
    ChemcraftWriter& operator=(const ChemcraftWriter&) = delete;

    void startCluster(const std::string& comment = "");
    void startSingle(const std::string& comment = "");

    void writeStructure(
        const Structure& structure,
        std::optional<double> energy_eV = std::nullopt
    );

    void close();
};
