# structure-clustering – Fast Exact Clustering of Molecular Structures through Parameterized Connectivity Graphs

**structure-clustering** groups molecular geometries by their distance-derived connectivity.

For each structure, the algorithm builds an undirected graph in which atoms are element-labelled vertices and connected atoms are edges. Structures are assigned the same cluster only when their graphs are exactly isomorphic. Fast graph signatures and Weisfeiler–Lehman refinement are used to reduce the number of expensive isomorphism checks. The edges describe the molecule's connectivity, not it's bonds.

<img src="https://github.com/user-attachments/assets/fef206d6-e039-49ce-911d-627068841853" width="50%" alt="Example clusters of Ag+(H2O)4 structures" />


## Installation

Install the package from PyPI:
```bash
pip install structure_clustering
```

Upgrade an existing installation with:
```bash
pip install --upgrade structure_clustering
```

Prebuilt wheels are available for common platforms. Building from source requires a C++17 compiler and the [Boost Graph Library](https://github.com/boostorg/graph) headers.

## Quick start

Cluster a multi-XYZ file:
```bash
structure_clustering structures.xyz
```

By default, only structures with fully connected graphs are clustered; structures with disconnected graphs are sorted out. The command writes the clustered structures to the [Chemcraft](https://www.chemcraftprog.com/)-compatible file `sc.chemcraft.chd`.

<img width="427" height="325" alt="Chemcraft visualization demo" src="https://github.com/user-attachments/assets/0a75b40f-b593-45ce-be1a-2e2e8ce6b77c" /><br/>


To choose another Chemcraft output path:
```bash
structure_clustering structures.xyz -ec results.chemcraft.chd
```

To write only a deduplicated set of structures, use `--representatives-only`:
```bash
structure_clustering structures.xyz --representatives-only
```

This keeps the first structure from each cluster and all unique singles in the exported file. Clustering itself is unchanged; the option only filters the written output. If native output is requested with `-e`, the same filtering is applied there as well.

The native clustering format is optional and can be written in addition:
```bash
structure_clustering structures.xyz -e
```

This writes `sc.dat` alongside the default `sc.chemcraft.chd`. To choose paths for both outputs:
```bash
structure_clustering structures.xyz -e results.dat -ec results.chemcraft.chd
```

The native `.dat` output can be visualized with [cluster-vis](https://photophys.github.io/cluster-vis/).

## Command-line interface

```text
usage: structure_clustering [-h] [--disconnected] [--config CONFIG]
                            [-e [FILE]] [-ec [FILE]]
                            [--representatives-only]
                            xyz_file
```

Options:

- `xyz_file`: multi-XYZ file containing the input structures.
- `--config CONFIG`: TOML file with connectivity settings.
- `--disconnected`: include disconnected graphs in clustering. By default, disconnected structures are sorted out. An explicit CLI flag overrides `options.only_connected_graphs` from the TOML file.
- `-e [FILE]`, `--export [FILE]`: additionally write the native clustering output. Without a filename, uses `sc.dat`.
- `-ec [FILE]`, `--export-chemcraft [FILE]`: set the Chemcraft-compatible output path. Chemcraft output is written by default to `sc.chemcraft.chd`. Note that the file extension needs to be **`.chd`**, otherwise Chemcraft will not read it properly.
- `--representatives-only`: export only the first structure from each cluster, plus all unique single structures. This affects both Chemcraft and native output and does not change clustering itself.
- `-h`, `--help`: show the command-line help.

If the installed `structure_clustering` command is not on your `PATH`, use:

```bash
python -m structure_clustering structures.xyz
```

### Example

```bash
structure_clustering structures.xyz -ec clusters.chd --representatives-only
```

A typical run reports the number of clusters and singles, the number of structures sorted out, cluster-size statistics, and connectivity statistics. In the default connected-only mode, the sorted-out count includes disconnected structures as well as redundant members removed when one representative per cluster is retained.

## Connectivity model and configuration

Connectivity is derived from interatomic distances.

For atom pairs without a dedicated pair-distance rule, atoms are connected when
```math
0.8\,\text{Å} < d < r(A) + r(B) + 0.4\,\text{Å,}
```
where $d$ is the distance between the atoms $A$ and $B$, with their covalent radii $r(A)$ and $r(B)$.

A pair-specific distance replaces that rule for the selected element pair. The CLI defines an O–H maximum distance of 2.3 Å by default; it can be replaced in the configuration file.

Example `sc_config.toml`:

```toml
[covalent]
He = 0.90
Ag = 1.59

[pair]
O-H = 2.30

[options]
only_connected_graphs = true
```

All sections and entries are optional:

- `[covalent]` overrides the default covalent radius of an element.
- `[pair]` defines a maximum distance for a specific unordered atom pair. `O-H` and `H-O` therefore refer to the same pair.
- `[options].only_connected_graphs` controls whether disconnected graphs are eligible for clustering.

Distances are in Ångström and element symbols are case-sensitive.

### Connected vs. disconnected graphs

Connected-only clustering is the default. In this mode, a disconnected input structure is **sorted out entirely**: it is not compared with other structures and does not appear in either `result.clusters` or `result.singles`.

Use `--disconnected` to include disconnected structures in clustering by graph topology. With that option enabled, disconnected structures can appear in clusters or as singles just like connected structures.

## How clustering works

For each input structure, `structure_clustering`:

1. Builds an undirected, element-labelled connectivity graph from the geometry.
2. Computes inexpensive graph descriptors such as vertex/edge counts, element counts, element-degree combinations, and connected-component sizes.
3. Applies WL refinement to obtain order-independent vertex colours and a compact graph signature.
4. Buckets structures with matching signatures.
5. Runs exact Boost graph-isomorphism checks only inside the surviving buckets.

The signatures and WL hashes are **filters, not proofs of graph identity**. Final cluster membership still requires exact labelled graph isomorphism.

This makes the current implementation conceptually different from the previous version, which compared each structure against existing cluster representatives much more broadly. The new approach moves most rejection work into reusable signatures and reserves exact isomorphism for plausible matches.

Graphs do not have a natural ordering of vertices. [Weisfeiler-Lehman](https://en.wikipedia.org/wiki/Weisfeiler_Leman_graph_isomorphism_test) (WL) refinement creates a canonical, order-independent description of a graph’s structure.

## Python API

### Cluster structures from a multi-XYZ file

```python
import structure_clustering

machine = structure_clustering.Machine()

# Optional connectivity overrides.
machine.setCovalentRadius(1, 0.42)
machine.addPairDistance(8, 1, 2.30)
machine.setOnlyConnectedGraphs(True)

structures = structure_clustering.import_multi_xyz("structures.xyz")
result = machine.cluster(structures)

print("clusters:", result.clusters)
print("singles:", result.singles)
```

`result.clusters` contains the clusters (groups), `result.singles` contains eligible structures that have no clustered equivalent. When connected-only mode is active, disconnected structures are excluded from both collections.

The indices in `clusters` and `singles` are **zero-based indices into the original input structure list**, so skipped disconnected structures can create gaps in the reported indices. Do not rely on the ordering of the cluster list itself.

### Create a structure programmatically

```python
from structure_clustering import Atom, Structure

structure = Structure(0)
structure.addAtom(Atom(8, -1.674872668,  0.000000, -0.984966492))
structure.addAtom(Atom(1, -1.674872668,  0.759337, -0.388923492))
structure.addAtom(Atom(1, -1.674872668, -0.759337, -0.388923492))
```

The integer passed to `Structure(...)` is the structure's stored ID; cluster membership is still reported using the structure's position in the input list.

### Inspect clustered structures

After clustering, graph-derived information is available through `result.structures`:

```python
structure = result.structures[5]

print("atoms:", structure.numAtoms)
print("connections:", structure.numConnections)
print("fragments:", structure.numFragments)
print("hash:", structure.getHash())

first_atom = structure.getAtom(0)
print("atomic number:", first_atom.atomic_number)
print("x coordinate:", first_atom.position.x)

print("fragment 0 atom indices:", structure.getFragmentAtomIndices(0))
```

### Export from Python

```python
result.export("clusters.dat")
result.exportChemcraft("clusters.chemcraft.xyz")
```

The native `.dat` format contains cluster groups, singles, graph edges, and geometries. Group and structure numbers written to that format are one-based, while the Python result indices are zero-based.

## Structure hashes

`Structure.getHash()` returns a 16-character hexadecimal WL-style graph hash. The implementation uses deterministic 64-bit [FNV-1a-based](https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function) mixing, so the hash does not depend on Python's process-randomized hashing or C++ `std::hash` implementation details.

The hash is independent of vertex ordering and useful for keeping track of clusters across multiple program runs (e.g., molecular dynamics simulations). It is not a substitute for exact graph isomorphism: hash collisions and WL-indistinguishable non-isomorphic graphs are possible but very unlikely. 

Example:

```text
0504d8ff3dc965c0
```

Example structure with two connected components:

![Structure Clustering example with two fragments](https://github.com/user-attachments/assets/8b3560e8-c334-4ac5-beac-e7ee47e2633d)

## Changes to v1.1.6

- **Significantly improved performance:** graph signatures and integer WL refinement now bucket plausible matches before exact isomorphism, instead of broadly comparing each structure with cluster representatives. Graph isomorphism remains the final equivalence test.
- **Deterministic hashing:** the previous string-WL approach is replaced by shared integer-WL/FNV-based hashing.
- **Chemcraft CLI output:** the CLI now writes Chemcraft-compatible output by default.

*No breaking changes compared to v1.1.6. This version is confirmed to produce the exact same results as v1.1.6.*

## Development

Local development requires Python, CMake, a C++17 compiler, and Boost headers.

Create and activate a virtual environment, then install the project from the repository root:

```bash
python -m venv ~/venvs/structure_clustering_dev
source ~/venvs/structure_clustering_dev/bin/activate
pip install .
```

Python bindings are defined in `src/main.cpp`. When adding a new C++ API method or property that should be available from Python, expose it there as well.

The GitHub Actions workflow builds wheels for the configured Python/platform matrix when the corresponding workflow is triggered.

## License

**structure-clustering** is distributed under the MIT License. See [LICENSE](LICENSE).
