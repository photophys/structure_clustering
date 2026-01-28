#include <iostream>
#include <stdexcept>
#include <string>

#include "Machine.hpp"
#include "Result.hpp"
#include "Structure.hpp"
#include "constants.hpp"

Machine::Machine() : _onlyConnectedGraphs(true), _covalentRadii(DEFAULT_COVALENT_RADII) {}

void Machine::setOnlyConnectedGraphs(bool onlyConnectedGraphs) {
    _onlyConnectedGraphs = onlyConnectedGraphs;
}

void Machine::setCovalentRadius(int atomicNumber, double radius) {
    if (atomicNumber <= 0 || atomicNumber > _covalentRadii.size()) {
        throw std::out_of_range("Invalid atomic number");
    }
    _covalentRadii[atomicNumber - 1] = radius;
}

void Machine::addPairDistance(int atomicNumberA, int atomicNumberB, double maxDistance) {
    _pairDistances[{atomicNumberA, atomicNumberB}] = maxDistance;
    _pairDistances[{atomicNumberB, atomicNumberA}] = maxDistance;
}

bool Machine::isOnlyConnectedGraphs() const { return _onlyConnectedGraphs; }

double Machine::getCovalentRadius(int atomicNumber) const {
    if (atomicNumber <= 0 || atomicNumber > _covalentRadii.size()) {
        throw std::out_of_range("Invalid atomic number");
    }

    double radius = _covalentRadii[atomicNumber - 1];

    if (radius == -1) {
        throw std::out_of_range("No covalent radius specified for atomic number " +
                                std::to_string(atomicNumber));
    } else {
        return radius;
    }
}

double Machine::getMaxPairDistance(int atomicNumberA, int atomicNumberB) const {
    std::pair<int, int> pair1 = std::make_pair(atomicNumberA, atomicNumberB);
    if (_pairDistances.find(pair1) != _pairDistances.end()) {
        return _pairDistances.at(pair1);
    }

    std::pair<int, int> pair2 = std::make_pair(atomicNumberB, atomicNumberA);
    if (_pairDistances.find(pair2) != _pairDistances.end()) {
        return _pairDistances.at(pair2);
    }

    return -1;
}

namespace {
    struct VertexInvariant {
        using Map = std::map<std::string, size_t>;
        Graph const &_graph;
        Map &_mappings;

        using result_type = size_t;
        using argument_type = Graph::vertex_descriptor;

        size_t operator()(argument_type u) const {
            return _mappings.at(boost::get(boost::vertex_name, _graph, u));
        }
        size_t max() const { return _mappings.size(); }

        void collect_names() {
            for (auto vd : boost::make_iterator_range(boost::vertices(_graph))) {
                size_t next_id = _mappings.size();
                _mappings.insert({boost::get(boost::vertex_name, _graph, vd), next_id});
            }
        }
    };
} // namespace

template <typename graph1, typename graph2>
bool is_named_vertices_isomorphic(const graph1 &g, const graph2 &h) noexcept {
    auto ref_index_map = get(boost::vertex_index, g);
    using vd = typename boost::graph_traits<graph1>::vertex_descriptor;
    std::vector<vd> iso(boost::num_vertices(g));

    VertexInvariant::Map shared_names;
    VertexInvariant inv1{g, shared_names};
    VertexInvariant inv2{h, shared_names};

    inv1.collect_names();
    inv2.collect_names();

    return boost::isomorphism(
        g, h,
        boost::isomorphism_map(make_iterator_property_map(iso.begin(), ref_index_map))
            .vertex_invariant1(inv1)
            .vertex_invariant2(inv2));
}

Result Machine::cluster(std::vector<Structure> &structures) {
    // construct graphs for each structure
    for (auto &structure : structures) {
        structure.constructGraph(*this);
    }

    // temporarily write clustered indices here
    std::vector<std::vector<int>> clusterIndices;

    bool found;

    // For each structure, check if its graph is isomorphic to any of the other graphs.
    // The graphs can only be isomorphic if they have the same number of bonds and the same vertex
    // labels. If specified in the options, the graphs are also required to be fully connected.
    for (int i = 0; i < (int)structures.size(); i++) {
        found = false;
        for (int j = 0; j < (int)clusterIndices.size(); j++) {
            if (structures[i].getNumConnections() ==
                structures[clusterIndices[j][0]].getNumConnections()) {
                found = is_named_vertices_isomorphic(structures[clusterIndices[j][0]].getGraph(),
                                                     structures[i].getGraph());
                if (found && (!_onlyConnectedGraphs || structures[i].isGraphFullyConnected())) {
                    clusterIndices[j].push_back(i);
                }
                if (found)
                    break;
            }
        }
        if (!found && (!_onlyConnectedGraphs || structures[i].isGraphFullyConnected()))
            clusterIndices.push_back({i});
    }

    // cleanup result
    Result result = Result(structures);
    for (const auto &cluster : clusterIndices) {
        if (cluster.size() > 1) {
            result.addCluster(cluster);
        } else {
            result.addSingle(cluster[0]);
        }
    }

    return result;
}
