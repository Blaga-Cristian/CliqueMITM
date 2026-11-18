import os
import random
import re
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
EXE_SUFFIX = ".exe" if os.name == "nt" else ""
SOLVER = ROOT / "build" / f"balanced_mitm{EXE_SUFFIX}"

SEED = 20261118
NUM_GRAPHS = 180
DENSITIES = (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0)
ORDERINGS = ("input", "high", "degen", "clique3", "clique4")


def generate_graph(n, probability, rng):
    return [
        (u, v)
        for u in range(n)
        for v in range(u)
        if rng.random() < probability
    ]


def write_graph(path, n, edges):
    with path.open("w") as file:
        file.write(f"{n} {len(edges)}\n")
        for u, v in edges:
            file.write(f"{u} {v}\n")


def run_solver(graph_path, algorithm, split, ordering):
    command = [
        str(SOLVER),
        str(graph_path),
        algorithm,
        "--order", ordering,
        "--split", str(split),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=True)

    match = re.search(r"counts=(\S+)", result.stdout)
    if match is None:
        raise RuntimeError(f"Missing histogram in output: {result.stdout}")

    return tuple(int(value) for value in match.group(1).split(",") if value)


def main():
    rng = random.Random(SEED)

    with tempfile.TemporaryDirectory() as directory:
        graph_path = Path(directory) / "test.edges"

        for test_id in range(NUM_GRAPHS):
            n = rng.randint(1, 14)
            probability = rng.choice(DENSITIES)
            edges = generate_graph(n, probability, rng)
            write_graph(graph_path, n, edges)

            expected = run_solver(graph_path, "enumerate", 0, "input")
            splits = sorted({0, min(4, n), n // 2})

            for split in splits:
                for ordering in ORDERINGS:
                    result = run_solver(graph_path, "mitm-dp", split, ordering)
                    if result != expected:
                        raise AssertionError(
                            f"Mismatch for graph {test_id}: "
                            f"n={n}, p={probability}, split={split}, order={ordering}"
                        )

    print(f"PASS: {NUM_GRAPHS} graphs; all MITM configurations match enumeration")


if __name__ == "__main__":
    main()
