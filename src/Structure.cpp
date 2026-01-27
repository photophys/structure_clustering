#include <iostream>
#include <string>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

#include "Atom.hpp"
#include "Machine.hpp"
#include "Structure.hpp"
#include "hash.hpp"

double distance(const Atom &atom_a, const Atom &atom_b) {
    return sqrt(pow(atom_b.position().x() - atom_a.position().x(), 2) +
                pow(atom_b.position().y() - atom_a.position().y(), 2) +
                pow(atom_b.position().z() - atom_a.position().z(), 2));
}

bool isConnected(const double &r_a, const double &r_b, const double &distance_ab) {
    return .8 < distance_ab && distance_ab < r_a + r_b + .4;
}

Structure::Structure(const int id) : _id(id), _graph({}) {}

void Structure::addAtom(const Atom &atom) { _atoms.push_back(atom); }

int Structure::numAtoms() const { return _atoms.size(); }
const Atom &Structure::getAtom(int index) const { return _atoms[index]; }
const int Structure::getNumConnections() const { return boost::num_edges(_graph); }
const Graph &Structure::getGraph() const { return _graph; }
const bool Structure::isGraphFullyConnected() const {
    std::vector<int> component(boost::num_vertices(_graph));
    int num = boost::connected_components(_graph, &component[0]);
    return num == 1;
}

void Structure::constructGraph(const Machine &machine) {
    // add graph vertices
    auto vertex_name_map = get(boost::vertex_name, _graph);
    std::vector<Graph::vertex_descriptor> vertex_descriptors(this->numAtoms());

    for (std::size_t i = 0; i < this->numAtoms(); ++i) {
        auto vd = boost::add_vertex(_graph);
        vertex_descriptors[i] = vd;
        vertex_name_map[vd] = this->getAtom(i).atomicNumber();
    }

    // check connectivity and add graph edges
    for (int i = 0; i < this->numAtoms(); i++) {
        for (int j = i + 1; j < this->numAtoms(); j++) {
            const Atom &atom_a = this->getAtom(i);
            const Atom &atom_b = this->getAtom(j);

            const int &atomicNumberA = atom_a.atomicNumber();
            const int &atomicNumberB = atom_b.atomicNumber();

            double distance_ab = distance(atom_a, atom_b);
            double r_a = machine.getCovalentRadius(atomicNumberA);
            double r_b = machine.getCovalentRadius(atomicNumberB);
            double maxPairDistance = machine.getMaxPairDistance(atomicNumberA, atomicNumberB);

            if (isConnected(r_a, r_b, distance_ab) || (distance_ab <= maxPairDistance)) {

                boost::add_edge(vertex_descriptors[i], vertex_descriptors[j], _graph);
            }
        }
    }
};

/**
 * Compute a hash value that represents the structure of a graph.
 *
 * The hash is intended to be invariant to vertex ordering and
 * sensitive to local neighborhood structure.
 *
 * Weisfeiler-Lehman (WL) refinement:
 * Graphs do not have a natural ordering of vertices. WL refinement creates a canonical,
 * order-independent description of a graph’s structure.
 * 1) Start with simple labels (element names, not unique).
 * 2) Repeatedly update each label using:
 *    - the current label of the vertex
 *    - the multiset of neighbor labels (https://en.wikipedia.org/wiki/Multiset)
 * 3) After several iterations, vertices with different local structures almost always
 * have different labels.
 * https://en.wikipedia.org/wiki/Weisfeiler_Leman_graph_isomorphism_test
 *
 * @param g           The input graph
 * @param iterations  Number of WL refinement iterations (controls locality depth)
 * @return            A hash representing the graph structure
 */
std::size_t hash_graph(const Graph& g, std::size_t iterations = 10)
{
    // Access the vertex "name" property
    auto name = get(boost::vertex_name, g);
    
    // One label per vertex (labels are refined over iterations)
    std::vector<std::string> labels(num_vertices(g));

    // --- Initial labeling ---
    // Each vertex starts with the element name
    for (auto v : boost::make_iterator_range(vertices(g)))
        labels[v] = name[v];

    // --- Weisfeiler-Lehman refinement ---
    // Each iteration updates vertex labels based on:
    //   - the current label of the vertex
    //   - the multiset of labels of its neighbors
    // (That means: Two vertices become distinguishable if
    // their neighborhoods differ.)
    for (std::size_t it = 0; it < iterations; ++it) {
        std::vector<std::string> new_labels(labels.size());

        for (auto v : boost::make_iterator_range(vertices(g))) {
            std::vector<std::string> neigh;

            // Collect labels of neighboring vertices
            for (auto u : boost::make_iterator_range(adjacent_vertices(v, g)))
                neigh.push_back(labels[u]);

            // Sort to make the neighborhood representation independent of order
            std::sort(neigh.begin(), neigh.end());

            // Combine the vertex label and its neighborhood into a single string
            // Example: C(H,H,O)
            std::ostringstream oss;
            oss << labels[v] << "(";
            for (auto& n : neigh)
                oss << n << ",";
            oss << ")";

            new_labels[v] = oss.str();
        }

        // Replace old labels with refined labels
        labels.swap(new_labels);
    }

    // --- Canonization ---
    // The final graph representation is a sorted multiset of vertex labels.
    // Sorting removes dependence on vertex indices.
    // https://en.wikipedia.org/wiki/Graph_canonization
    std::sort(labels.begin(), labels.end());

    std::ostringstream canonical;
    for (auto& l : labels)
        canonical << l << ";";

    // Get the hash (this hash is of type size_t - useful for fast lookup)
    return std::hash<std::string>{}(canonical.str());
}

/**
 * Return a hexadecimal string representing the structure hash.
 *
 * Internally, it computes a graph hash using WL refinement and converts
 * it to a fixed-width hexadecimal string.
 * 
 * @return Structure hash
 */
const std::string Structure::getHash() const {
    std::uint64_t h = static_cast<std::uint64_t>(hash_graph(_graph));
    return to_hex_string(h);
}

/**
 * Compute the number of fragments (Boost: connected components) in the graph.
 *
 * @return Number of connected components
 */
const int Structure::getNumFragments() const {
    std::vector<int> component(boost::num_vertices(_graph));
    return boost::connected_components(_graph, &component[0]);
}

/**
 * Return the vertex indices belonging to a specific fragment
 * (Boost: connected component).
 *
 * @param fragmentIndex Index of the fragment
 * @return              List of vertex indices in that fragment
 */
const std::vector<int> Structure::getFragmentAtomIndices(int fragmentIndex) const
{
    const std::size_t n = boost::num_vertices(_graph);

    // Component index for each vertex
    std::vector<int> component(n);
    boost::connected_components(_graph, component.data());

    std::vector<int> atomIndices;
    atomIndices.reserve(n);

    // Collect vertices that belong to the requested component
    for (std::size_t v = 0; v < n; ++v) {
        if (component[v] == fragmentIndex) {
            atomIndices.push_back(static_cast<int>(v));
        }
    }

    return atomIndices;
}
