#ifndef __RESULT__
#define __RESULT__

#include <string>
#include <vector>

#include "Structure.hpp"

// indices (clusters, singles) represent the indices from the initially provided structure list
class Result {
    std::vector<std::vector<int>> _clusters;
    std::vector<int> _singles;
    std::vector<Structure> _structures;

public:
    Result(std::vector<Structure> &structures);

    void addSingle(const int structureIdx);
    void addCluster(const std::vector<int> structureIdxs);

    void exportDat(const std::string filepath) const;
    void exportChemcraft(const std::string filepath) const;

    std::vector<std::vector<int>> getClusters() const;
    std::vector<int> getSingles() const;
    std::vector<Structure> getStructures() const;
};

#endif
