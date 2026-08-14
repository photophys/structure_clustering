#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/graph/connected_components.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/range/iterator_range.hpp>

#include "Machine.hpp"
#include "Result.hpp"
#include "Structure.hpp"
#include "constants.hpp"
#include "hash.hpp"


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
    const auto pair = std::minmax(atomicNumberA, atomicNumberB);
    _pairDistances[{pair.first, pair.second}] = maxDistance;
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
    const auto pair = std::minmax(atomicNumberA, atomicNumberB);
    const auto it = _pairDistances.find({pair.first, pair.second});
    return it != _pairDistances.end() ? it->second : -1;
}

namespace {
    std::vector<int> component_sizes(const Graph &g) {
        const std::size_t n = boost::num_vertices(g);
        std::vector<int> component(n);
        const int num = boost::connected_components(g, component.data());

        std::vector<int> sizes(num, 0);
        for (std::size_t i = 0; i < component.size(); ++i) {
            ++sizes[component[i]];
        }
        std::sort(sizes.begin(), sizes.end());
        return sizes;
    }

    std::string cheap_graph_signature(const Graph &g) {
        std::map<std::string, int> atom_counts;
        std::vector<std::pair<std::string, std::size_t>> atom_degree_pairs;
        atom_degree_pairs.reserve(boost::num_vertices(g));

        for (auto v : boost::make_iterator_range(boost::vertices(g))) {
            const std::string label = graph_hash::vertex_label(g, v);
            ++atom_counts[label];
            atom_degree_pairs.push_back(std::make_pair(label, boost::out_degree(v, g)));
        }
        std::sort(atom_degree_pairs.begin(), atom_degree_pairs.end());

        std::ostringstream key;
        key << "n=" << boost::num_vertices(g)
            << "|m=" << boost::num_edges(g);

        key << "|atoms=";
        for (const auto &entry : atom_counts) {
            key << entry.first << ':' << entry.second << ',';
        }

        key << "|deg=";
        for (const auto &entry : atom_degree_pairs) {
            key << entry.first << ':' << entry.second << ',';
        }

        key << "|components=";
        const auto sizes = component_sizes(g);
        for (int size : sizes) {
            key << size << ',';
        }

        return key.str();
    }

    struct ColorInvariant {
        const Graph &_graph;
        const std::vector<std::size_t> &_color_ids;
        std::size_t _max;

        using result_type = std::size_t;
        using argument_type = Graph::vertex_descriptor;

        std::size_t operator()(argument_type u) const {
            auto index = boost::get(boost::vertex_index, _graph);
            return _color_ids[index[u]];
        }

        std::size_t max() const { return _max; }
    };

    bool same_color_multiset(std::vector<std::uint64_t> a, std::vector<std::uint64_t> b) {
        if (a.size() != b.size()) {
            return false;
        }
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        return a == b;
    }

    std::vector<std::size_t> compress_colors(
        const std::vector<std::uint64_t> &colors,
        const std::unordered_map<std::uint64_t, std::size_t> &mapping) {
        std::vector<std::size_t> ids(colors.size());
        for (std::size_t i = 0; i < colors.size(); ++i) {
            ids[i] = mapping.at(colors[i]);
        }
        return ids;
    }

    bool is_refined_vertices_isomorphic(const Graph &g,
                                        const std::vector<std::uint64_t> &colors_g,
                                        const Graph &h,
                                        const std::vector<std::uint64_t> &colors_h) noexcept {
        if (boost::num_vertices(g) != boost::num_vertices(h) ||
            boost::num_edges(g) != boost::num_edges(h)) {
            return false;
        }

        if (!same_color_multiset(colors_g, colors_h)) {
            return false;
        }

        std::vector<std::uint64_t> unique_colors;
        unique_colors.reserve(colors_g.size() + colors_h.size());
        unique_colors.insert(unique_colors.end(), colors_g.begin(), colors_g.end());
        unique_colors.insert(unique_colors.end(), colors_h.begin(), colors_h.end());
        std::sort(unique_colors.begin(), unique_colors.end());
        unique_colors.erase(std::unique(unique_colors.begin(), unique_colors.end()), unique_colors.end());

        std::unordered_map<std::uint64_t, std::size_t> mapping;
        mapping.reserve(unique_colors.size());
        for (std::size_t i = 0; i < unique_colors.size(); ++i) {
            mapping[unique_colors[i]] = i;
        }

        const auto color_ids_g = compress_colors(colors_g, mapping);
        const auto color_ids_h = compress_colors(colors_h, mapping);

        auto ref_index_map = boost::get(boost::vertex_index, g);
        using vd = boost::graph_traits<Graph>::vertex_descriptor;
        std::vector<vd> iso(boost::num_vertices(g));

        ColorInvariant inv_g{g, color_ids_g, unique_colors.size()};
        ColorInvariant inv_h{h, color_ids_h, unique_colors.size()};

        return boost::isomorphism(
            g, h,
            boost::isomorphism_map(boost::make_iterator_property_map(iso.begin(), ref_index_map))
                .vertex_invariant1(inv_g)
                .vertex_invariant2(inv_h));
    }
} // namespace

Result Machine::cluster(std::vector<Structure> &structures) {
    for (auto &structure : structures) {
        structure.constructGraph(*this);
    }

    const int n = static_cast<int>(structures.size());
    std::vector<std::vector<int>> clusterIndices;
    std::unordered_map<std::string, std::vector<int>> buckets;
    buckets.reserve(structures.size() * 2 + 1);

    std::vector<std::vector<std::uint64_t>> wlColors(n);
    // Build one key per graph: (This removes the global pairwise comparison loop.)
    for (int i = 0; i < n; ++i) {
        const Graph &g = structures[i].getGraph();

        if (_onlyConnectedGraphs && structures[i].getNumFragments() != 1) {
            // Disconnected graphs are sorted out when only connected graphs are requested.
            continue;
        }

        wlColors[i] = graph_hash::wl_vertex_colors(g);
        const std::uint64_t wlHash = graph_hash::wl_multiset_hash(g, wlColors[i]);

        std::string key = cheap_graph_signature(g);
        key += "|wl=";
        key += to_hex_string(wlHash);

        buckets[key].push_back(i);
    }

    // Most buckets should already be very close to the final clusters.
    // Exact isomorphism checks are performed only inside buckets that
    // survived the cheap signature + integer WL filters.
    for (const auto &bucket : buckets) {
        const auto &indices = bucket.second;
        if (indices.size() == 1) {
            clusterIndices.push_back({indices[0]});
            continue;
        }

        std::vector<std::vector<int>> localClusters;
        for (int idx : indices) {
            bool found = false;

            for (auto &cluster : localClusters) {
                const int rep = cluster.front();
                if (is_refined_vertices_isomorphic(structures[rep].getGraph(), wlColors[rep],
                                                   structures[idx].getGraph(), wlColors[idx])) {
                    cluster.push_back(idx);
                    found = true;
                    break;
                }
            }

            if (!found) {
                localClusters.push_back({idx});
            }
        }

        clusterIndices.insert(clusterIndices.end(), localClusters.begin(), localClusters.end());
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
