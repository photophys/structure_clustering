# structure_clustering &ndash; Cluster Molecular Structures Into Groups of Similar Ones

**structure_clustering** is a Python package to cluster molecular structures into groups of similar ones. Our approach involves analysing the intermolecular distances to represent each structure's connectivity as an undirected, vertex-labelled graph. It then uses graph isomorphism to identify structures that belong to the same group. The package offers a command-line interface for clustering a multi-XYZ file or can be used within your Python code.

<img src="https://github.com/user-attachments/assets/fef206d6-e039-49ce-911d-627068841853" width="50%" />[^1]

[^1]: The figure shows exemplary clusters from Ag⁺(H₂O)₄ structures.


## Installation

You can install structure_clustering via pip:
```bash
pip install structure_clustering
```

Prebuilt wheels are available for most platforms (Windows, Linux, MacOS). If you prefer to compile and build the wheel yourself, ensure that the [Boost Graph Library](https://www.boost.org/doc/libs/release/libs/graph/doc/index.html) is installed system-wide.


## Using the Command-Line Interface
You can invoke the structure_clustering script using the `structure_clustering` command.

<details>
  <summary>Use this method if the command does not work</summary>

  On some systems, scripts installed via pip are not added to the system's `PATH`. You can either [add](https://stackoverflow.com/a/70680333/17726525) them to your `PATH`, or run the script directly by invoking `python3 -m structure_clustering`.
</details>

```bash
usage: structure_clustering <xyz_file> [--config CONFIG] [--output OUTPUT] [--disconnected]

Cluster molecular structures into groups.

positional arguments:
  xyz_file         path of the multi-xyz-file containing the structures

options:
  --config CONFIG  path of the config TOML file
  --output OUTPUT  path of the resulting output file, defaults to <xyz_file>.sc.dat
  --disconnected   if you want to include disconnected graphs
  -h, --help       show this help message and exit
```

For example, to cluster an xyz file:
```bash
structure_clustering my_structures.xyz
```
To specify a custom distance for recognising O-H connectivity (see the next section), use a TOML config file:
```bash
structure_clustering my_structures.xyz --config sc_config.toml
```

In both cases, a file named `my_structures.xyz.sc.dat` will be created, which you can import at <a href="https://photophys.github.io/cluster-vis/"><img src="https://raw.githubusercontent.com/photophys/MOLGA.jl/refs/heads/main/docs/src/assets/logo.svg" height="15px" /> https://photophys.github.io/cluster-vis/</a> to visualise the results of your clustering process.

The terminal output will look like this:
```
Loading configuration from demo_config.toml
Using covalent radius of 1.59 for Ag
Using pair distance of 2.3 for O-H
Clustering does not include disconnected graphs

Using 437 structures from structures.xyz
Clustering finished <structure_clustering._core.Result object at 0x7f7c949c37b0>
  14 clusters (total 318 structures)
  13 unique single structures
  132 (30.21%) structures sorted out (305 remaining)
  cluster size: Avg=22.7 Med=4.5 Q1=2.2 Q3=23.5
  connections/structure: Avg=12.2 Med=12.0 Q1=12.0 Q3=12.0 (all 437)
  connections/structure: Avg=12.4 Med=12.0 Q1=12.0 Q3=12.0 (remaining 305)
Writing output file to structures.xyz.sc.dat ...

🚀 Open https://photophys.github.io/cluster-vis/ to visualize your results
```


## Configuration File
You can use a TOML file to control the parameters of the command-line interface. The `[covalent]` section allows you to override the algorithm's default covalent radii. In the `[pair]` section, you can specify a maximum distance for pairs of atoms.
```toml
[covalent]
He = 0.9
Ag = 1.59

[pair]
O-H = 2.3

[options]
only_connected_graphs = true
```
All settings are optional. Distances are given in Angstrom. Elements are case-sensitive. If you specify `only_connected_graphs` in the config file, this will overwrite your setting from the command-line switch.


## Example Code
```py
import structure_clustering
from structure_clustering import Structure, Atom

sc_machine = structure_clustering.Machine()

sc_machine.setCovalentRadius(1, 0.42)  # change hydrogen covalent radius to 0.42
sc_machine.addPairDistance(8, 1, 2.3)  # extend max distance for O-H pairs to 2.3 Ang

sc_machine.setOnlyConnectedGraphs(True)  # only include fully connected graphs (default)

# you will need some structures
population = structure_clustering.import_multi_xyz("structs.xyz")

# you can also create your structures programmatically
structure = Structure()
structure.addAtom(Atom(8, -1.674872668, 0.0, -0.984966492))
structure.addAtom(Atom(1, -1.674872668, 0.759337, -0.388923492))
structure.addAtom(Atom(1, -1.674872668, -0.759337, -0.388923492))
population += [structure]  # add this structure to our population

sc_result = sc_machine.cluster(population)

print("clusters", sc_result.clusters)
print("singles", sc_result.singles)

# Output (indices from the original structure list):
# clusters [[0, 11], [1, 2, 4, 6, 12, 13, 14, 15, 19], [3, 17, 18, 23]]
# singles [9, 16, 22]
```


## License
The structure_clustering package is licensed under the MIT License. See the [LICENSE file](LICENSE) for more details.
