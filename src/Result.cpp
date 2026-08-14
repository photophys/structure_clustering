#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <iostream>

#include "Result.hpp"
#include "constants.hpp"
#include "ChemcraftWriter.hpp"

Result::Result(std::vector<Structure> &structures)
    : _clusters({}), _singles({}), _structures(structures) {};

void Result::addSingle(const int structureIdx) { _singles.push_back(structureIdx); }

void Result::addCluster(const std::vector<int> structureIdxs) {
    _clusters.push_back(structureIdxs);
}

std::string getElementSymbol(int atomicNumber) {
    const auto &elements = CHEMICAL_ELEMENTS;

    if (atomicNumber < 1 || atomicNumber > elements.size()) {
        throw std::out_of_range("Invalid atomic number");
    }

    return elements[atomicNumber - 1];
}

void Result::exportDat(const std::string filepath) const {
    std::ofstream file;
    file.open(filepath);

    // write groups
    file << "@GROUPS" << std::endl;
    for (auto &cluster : _clusters) {
        for (auto &id : cluster) {
            file << id + 1 << " ";
        }
        file << std::endl;
    }

    // write singles
    file << "@UNIQUES" << std::endl;
    for (auto &id : _singles) {
        file << id + 1 << " ";
    }
    file << std::endl;

    // write graphs
    typedef boost::graph_traits<Graph>::edge_iterator edge_iterator;
    for (size_t i = 0; i < _structures.size(); i++) {
        auto &graph = _structures[i].getGraph();
        file << "@GRAPH-" << i + 1 << std::endl;
        std::pair<edge_iterator, edge_iterator> ei = edges(graph);
        for (edge_iterator it = ei.first; it != ei.second; ++it) {
            file << source(*it, graph) + 1 << " " << target(*it, graph) + 1 << std::endl;
        }
    }

    // write geometries
    for (size_t i = 0; i < _structures.size(); i++) {
        file << "@STRUCTURE-" << i + 1 << std::endl;
        auto &structure = _structures[i];
        for (int j = 0; j < structure.numAtoms(); j++) {
            auto &atom = structure.getAtom(j);
            file << getElementSymbol(atom.atomicNumber()) << " " << atom.position().x() << " "
                 << atom.position().y() << " " << atom.position().z() << std::endl;
        }
    }

    file.close();
}

void Result::exportChemcraft(const std::string filepath) const {
    auto chemcraftFile = ChemcraftWriter(filepath);
    // write groups
    for (auto &cluster : _clusters) {
        chemcraftFile.startCluster(_structures[cluster[0]].getHash());
        for (auto &id : cluster) {
            auto &structure = _structures[id];
            chemcraftFile.writeStructure(structure);
        }
    }

    // write singles
    for (auto &id : _singles) {
        auto &structure = _structures[id];
        chemcraftFile.startSingle(structure.getHash());
        chemcraftFile.writeStructure(structure);
    }

    chemcraftFile.close();
};

std::vector<std::vector<int>> Result::getClusters() const { return _clusters; };
std::vector<int> Result::getSingles() const { return _singles; };
std::vector<Structure> Result::getStructures() const { return _structures; };
