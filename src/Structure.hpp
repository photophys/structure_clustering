#ifndef __STRUCTURE__
#define __STRUCTURE__

#include <iostream>
#include <string>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

#include "Atom.hpp"

class Machine; // forward declaration

using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                    boost::property<boost::vertex_name_t, std::string>>;

class Structure {
    int _id;
    std::vector<Atom> _atoms;
    Graph _graph;

public:
    Structure(const int id);

    void addAtom(const Atom &atom);

    int numAtoms() const;
    const Atom &getAtom(int index) const;
    const int getNumConnections() const;
    const Graph &getGraph() const;
    const bool isGraphFullyConnected() const;

    void constructGraph(const Machine &machine);

    const std::string getHash() const;
    const std::string getCheapHash() const;

    const int getNumFragments() const;
    const std::vector<int> getFragmentAtomIndices(int fragmentIndex) const;

    int getId() const;
};

#endif
