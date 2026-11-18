import os
import random
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
EXE_SUFFIX = ".exe" if os.name == "nt" else ""
SOLVER = ROOT / "build" / f"block_mitm{EXE_SUFFIX}"
TEST_DIR = ROOT / "data" / "small_tests"

SEED = 20261118
DENSITIES = (0.10, 0.25, 0.50, 0.75, 0.90, 0.97)
ORDERINGS = ("high", "low", "peel", "clique3")


def write_graph(path, n, edges):
    with path.open("w") as file:
        file.write(f"{n} {len(edges)}\n")
        for u, v in edges:
            file.write(f"{u} {v}\n")


def brute_force_histogram(n, edges):
    neighbors = [0] * n
    for u, v in edges:
        neighbors[u] |= 1 << v
        neighbors[v] |= 1 << u

    is_clique = bytearray(1 << n)
    is_clique[0] = 1

    counts = [0] * (n + 1)
    counts[0] = 1

    for mask in range(1, 1 << n):
        lowest_bit = mask & -mask
        vertex = lowest_bit.bit_length() - 1
        rest = mask ^ lowest_bit

        if is_clique[rest] and (rest & ~neighbors[vertex]) == 0:
            is_clique[mask] = 1
            counts[mask.bit_count()] += 1

    return counts


def run_solver(graph_path, ordering, budget):
    command = [
        str(SOLVER),
        str(graph_path),
        "--order", ordering,
        "--budget", str(budget),
        "--seconds", "2",
        "--node-limit", "10000000",
    ]
    output = subprocess.check_output(command, text=True)
    fields = dict(item.split("=", 1) for item in output.split())

    if fields.get("status") != "COMPLETE":
        raise AssertionError(f"Solver did not finish: {fields}")

    histogram = [int(value) for value in fields["counts"].split(",")]
    return histogram


def main():
    TEST_DIR.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    checks = 0

    for n in range(1, 15):
        for test_id in range(8):
            probability = rng.choice(DENSITIES)
            edges = [
                (u, v)
                for u in range(n)
                for v in range(u + 1, n)
                if rng.random() < probability
            ]
            graph_path = TEST_DIR / f"g{n}_{test_id}.edges"
            write_graph(graph_path, n, edges)

            expected = brute_force_histogram(n, edges)
            for ordering in ORDERINGS:
                result = run_solver(graph_path, ordering, min(4, n))
                if result != expected:
                    raise AssertionError(f"Mismatch: {graph_path.name}, {ordering}")
                checks += 1

    print(f"PASS: {checks} configurations match exhaustive clique counts")


if __name__ == "__main__":
    main()
