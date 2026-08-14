import argparse
import re
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
import toml

import structure_clustering
from structure_clustering import element_to_atomic_number


DEFAULT_OH_MAX_DISTANCE = 2.3
DEFAULT_NATIVE_OUTPUT = Path("sc.dat")
DEFAULT_CHEMCRAFT_OUTPUT = Path("sc.chemcraft.chd")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="structure_clustering",
        description=(
            "Cluster molecular structures by user-defined, distance-derived "
            "connectivity graphs."
        ),
    )

    parser.add_argument(
        "xyz_file",
        type=Path,
        help="Path to the multi-XYZ file containing the structures.",
    )

    parser.add_argument(
        "--disconnected",
        action="store_true",
        default=None,
        help=(
            "Include disconnected graphs in clustering. "
            "This overrides options.only_connected_graphs in the TOML configuration."
        ),
    )

    parser.add_argument(
        "--config",
        type=Path,
        help="Path to a TOML configuration file.",
    )

    parser.add_argument(
        "-e",
        "--export",
        nargs="?",
        const=DEFAULT_NATIVE_OUTPUT,
        default=None,
        type=Path,
        metavar="FILE",
        help=(
            "Export native clustering output. "
            f"Without FILE, writes to {DEFAULT_NATIVE_OUTPUT}."
        ),
    )

    parser.add_argument(
        "-ec",
        "--export-chemcraft",
        nargs="?",
        const=DEFAULT_CHEMCRAFT_OUTPUT,
        default=DEFAULT_CHEMCRAFT_OUTPUT,
        type=Path,
        metavar="FILE",
        help=(
            "Export Chemcraft-compatible XYZ output. "
            f"Defaults to {DEFAULT_CHEMCRAFT_OUTPUT}."
        ),
    )

    parser.add_argument(
        "--representatives-only",
        action="store_true",
        help=(
            "Export only one representative from each cluster, "
            "plus all unique single structures."
        ),
    )

    return parser


def summarize(values: Sequence[int]) -> str:
    if not values:
        return "n/a"

    array = np.asarray(values, dtype=float)

    return (
        f"Avg={np.mean(array):.1f} "
        f"Med={np.median(array):.1f} "
        f"Q1={np.percentile(array, 25):.1f} "
        f"Q3={np.percentile(array, 75):.1f}"
    )


def parse_pair_name(pair_name: str) -> Tuple[int, int]:
    """Convert a pair key such as 'C-O' to two atomic numbers."""
    elements = [element.strip() for element in pair_name.split("-")]

    if len(elements) != 2 or not all(elements):
        raise ValueError(
            f"Invalid pair key {pair_name!r}; expected exactly two elements, "
            "for example 'C-O'."
        )

    try:
        atomic_number_a = element_to_atomic_number(elements[0])
        atomic_number_b = element_to_atomic_number(elements[1])
    except Exception as exc:
        raise ValueError(
            f"Could not interpret pair key {pair_name!r} as element symbols."
        ) from exc

    return atomic_number_a, atomic_number_b


def load_configuration(
    config_path: Path,
    machine: structure_clustering.Machine,
    only_connected_graphs: bool,
) -> bool:
    try:
        config: Dict[str, Any] = toml.load(str(config_path))
    except (OSError, toml.TomlDecodeError) as exc:
        raise ValueError(
            f"Could not read TOML configuration {config_path}: {exc}"
        ) from exc

    options = config.get("options", {})
    if not isinstance(options, dict):
        raise ValueError("'options' in the TOML configuration must be a table.")

    if "only_connected_graphs" in options:
        configured_value = options["only_connected_graphs"]

        if not isinstance(configured_value, bool):
            raise ValueError(
                "'options.only_connected_graphs' must be true or false."
            )

        only_connected_graphs = configured_value

    covalent_settings = config.get("covalent", {})
    if not isinstance(covalent_settings, dict):
        raise ValueError("'covalent' in the TOML configuration must be a table.")

    for element, radius in covalent_settings.items():
        try:
            atomic_number = element_to_atomic_number(element.strip())
            radius = float(radius)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid covalent-radius entry for {element!r}: {radius!r}"
            ) from exc

        if radius <= 0.0:
            raise ValueError(
                f"Covalent radius for {element!r} must be positive; got {radius}."
            )

        print(f"Using covalent radius {radius:.4g} Å for {element}")
        machine.setCovalentRadius(atomic_number, radius)

    pair_settings = config.get("pair", {})
    if not isinstance(pair_settings, dict):
        raise ValueError("'pair' in the TOML configuration must be a table.")

    for pair_name, distance in pair_settings.items():
        try:
            atomic_number_a, atomic_number_b = parse_pair_name(pair_name)
            distance = float(distance)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid pair-distance entry for {pair_name!r}: {distance!r}"
            ) from exc

        if distance <= 0.0:
            raise ValueError(
                f"Pair distance for {pair_name!r} must be positive; got {distance}."
            )

        print(f"Using pair distance {distance:.4g} Å for {pair_name}")
        machine.addPairDistance(atomic_number_a, atomic_number_b, distance)

    return only_connected_graphs


def validate_output_paths(
    parser: argparse.ArgumentParser,
    native_output: Optional[Path],
    chemcraft_output: Optional[Path],
) -> None:
    output_paths = [
        path
        for path in (native_output, chemcraft_output)
        if path is not None
    ]

    if (
        native_output is not None
        and chemcraft_output is not None
        and native_output.resolve() == chemcraft_output.resolve()
    ):
        parser.error(
            "Native and Chemcraft outputs must use different filenames."
        )

    for output_path in output_paths:
        if output_path.exists() and output_path.is_dir():
            parser.error(f"Output path is a directory: {output_path}")

        if not output_path.parent.exists():
            parser.error(
                f"Output directory does not exist: {output_path.parent}"
            )


def export_representative_chemcraft(
    result,
    output_path: Path,
) -> None:
    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary_path = Path(temporary_directory) / "full.chemcraft.xyz"
        result.exportChemcraft(str(temporary_path))
        lines = temporary_path.read_text().splitlines(keepends=True)

    filtered_lines = []
    job_pattern = re.compile(r"^\[Job \d+\]\s*$")
    index = 0

    while index < len(lines) and job_pattern.match(lines[index].strip()) is None:
        filtered_lines.append(lines[index])
        index += 1

    while index < len(lines):
        block_start = index
        index += 1
        while index < len(lines) and job_pattern.match(lines[index].strip()) is None:
            index += 1

        block = lines[block_start:index]
        is_cluster = len(block) > 1 and block[1].startswith("Cluster ")

        if not is_cluster:
            filtered_lines.extend(block)
            continue

        first_geometry_end = next(
            (
                line_index
                for line_index, line in enumerate(block)
                if line.strip() == "[/geometry]"
            ),
            None,
        )
        if first_geometry_end is None:
            raise RuntimeError(
                "Could not identify the first geometry in Chemcraft cluster output."
            )

        filtered_lines.extend(block[: first_geometry_end + 1])

        trailing_geometry = next(
            (
                line
                for line in block[first_geometry_end + 1 :]
                if line.strip() == "[Geometry]"
            ),
            None,
        )
        if trailing_geometry is not None:
            filtered_lines.append(trailing_geometry)

        filtered_lines.extend(
            line for line in block[first_geometry_end + 1 :]
            if line.strip() == "[/Job]"
        )

    output_path.write_text("".join(filtered_lines))


def export_representative_native(
    result,
    output_path: Path,
    clusters: Sequence[Sequence[int]],
    singles: Sequence[int],
    representative_indices: Sequence[int],
) -> None:
    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary_path = Path(temporary_directory) / "full.dat"
        result.export(str(temporary_path))
        lines = temporary_path.read_text().splitlines(keepends=True)

    section_pattern = re.compile(r"^@(GRAPH|STRUCTURE)-(\d+)\s*$")
    sections = {"GRAPH": {}, "STRUCTURE": {}}

    index = 0
    while index < len(lines):
        match = section_pattern.match(lines[index].strip())
        if match is None:
            index += 1
            continue

        section_type, old_number_text = match.groups()
        old_index = int(old_number_text) - 1
        section_start = index + 1
        index = section_start
        while index < len(lines) and section_pattern.match(lines[index].strip()) is None:
            index += 1
        sections[section_type][old_index] = lines[section_start:index]

    old_to_new = {
        old_index: new_index
        for new_index, old_index in enumerate(representative_indices, start=1)
    }

    output_lines = ["@GROUPS\n"]
    for cluster in clusters:
        if cluster:
            output_lines.append(f"{old_to_new[cluster[0]]}\n")

    output_lines.append("@UNIQUES\n")
    output_lines.append(
        " ".join(str(old_to_new[single]) for single in singles) + "\n"
    )

    for section_type in ("GRAPH", "STRUCTURE"):
        for new_index, old_index in enumerate(representative_indices, start=1):
            if old_index not in sections[section_type]:
                raise RuntimeError(
                    f"Could not find @{section_type}-{old_index + 1} in native output."
                )
            output_lines.append(f"@{section_type}-{new_index}\n")
            output_lines.extend(sections[section_type][old_index])

    output_path.write_text("".join(output_lines))


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.xyz_file.is_file():
        parser.error(f"Input XYZ file does not exist: {args.xyz_file}")

    if args.config is not None and not args.config.is_file():
        parser.error(f"Configuration file does not exist: {args.config}")

    validate_output_paths(
        parser,
        native_output=args.export,
        chemcraft_output=args.export_chemcraft,
    )

    # Default behavior: cluster only fully connected graphs.
    only_connected_graphs = True

    machine = structure_clustering.Machine()

    # Default O-H connectivity threshold. A TOML [pair] entry for O-H or H-O
    # replaces this setting.
    machine.addPairDistance(1, 8, DEFAULT_OH_MAX_DISTANCE)

    if args.config is not None:
        print(f"Loading configuration from {args.config}")

        try:
            only_connected_graphs = load_configuration(
                args.config,
                machine,
                only_connected_graphs,
            )
        except ValueError as exc:
            parser.error(str(exc))

    # An explicit command-line option takes precedence over TOML.
    if args.disconnected is not None:
        only_connected_graphs = not args.disconnected

    machine.setOnlyConnectedGraphs(only_connected_graphs)

    if only_connected_graphs:
        print("Clustering only connected graphs")
    else:
        print("Clustering includes disconnected graphs")

    structures = structure_clustering.import_multi_xyz(str(args.xyz_file))

    if not structures:
        parser.error(f"No structures were read from {args.xyz_file}")

    num_input_structures = len(structures)
    print(f"\nUsing {num_input_structures} structures from {args.xyz_file}")

    result = machine.cluster(structures)

    # Result owns the graph-constructed structures returned from C++.
    clustered_structures = result.structures
    clusters = [list(cluster) for cluster in result.clusters]
    singles = list(result.singles)

    classified_indices = singles + [
        structure_index
        for cluster in clusters
        for structure_index in cluster
    ]

    expected_indices = {
        index
        for index, structure in enumerate(clustered_structures)
        if not only_connected_graphs or structure.numFragments == 1
    }
    actual_indices = set(classified_indices)

    if len(classified_indices) != len(actual_indices) or actual_indices != expected_indices:
        raise RuntimeError(
            "Internal clustering error: result indices do not form a complete, "
            "non-overlapping partition of the eligible input structures."
        )

    cluster_sizes = [len(cluster) for cluster in clusters]
    num_clustered_structures = sum(cluster_sizes)

    # Each cluster with more than one member retains one representative.
    remaining_indices = sorted(
        singles + [cluster[0] for cluster in clusters if cluster]
    )

    num_remaining_structures = len(remaining_indices)
    num_sorted_out = num_input_structures - num_remaining_structures
    percentage_sorted_out = 100.0 * num_sorted_out / num_input_structures

    all_connection_counts = [
        structure.numConnections
        for structure in clustered_structures
    ]

    remaining_connection_counts = [
        all_connection_counts[index]
        for index in remaining_indices
    ]

    print("\nClustering finished")
    print(
        f"  {len(clusters)} clusters "
        f"({num_clustered_structures} structures in non-singleton clusters)"
    )
    print(f"  {len(singles)} unique single structures")
    print(
        f"  {num_sorted_out} ({percentage_sorted_out:.2f}%) structures removed "
        f"({num_remaining_structures} representatives remain)"
    )
    print(f"  Cluster size: {summarize(cluster_sizes)}")
    print(
        f"  Connections/structure: {summarize(all_connection_counts)} "
        f"(all {num_input_structures})"
    )
    print(
        f"  Connections/structure: {summarize(remaining_connection_counts)} "
        f"(remaining {num_remaining_structures})"
    )

    if args.export is not None:
        print(f"\nWriting native output file to {args.export} ...")
        if args.representatives_only:
            export_representative_native(
                result,
                args.export,
                clusters,
                singles,
                remaining_indices,
            )
        else:
            result.export(str(args.export))

    if args.export_chemcraft is not None:
        print(
            f"\nWriting Chemcraft output file to "
            f"{args.export_chemcraft} ..."
        )
        if args.representatives_only:
            export_representative_chemcraft(result, args.export_chemcraft)
        else:
            result.exportChemcraft(str(args.export_chemcraft))

    if args.export is not None:
        print("\n🚀 \033[1mOpen https://photophys.github.io/cluster-vis/ to visualize the native output\033[0m")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
