import os
import random
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
EXE_SUFFIX = ".exe" if os.name == "nt" else ""
SOLVER = ROOT / "build" / f"block_mitm_compact{EXE_SUFFIX}"

SEED = 5181
DENSITIES = (0.0, 0.2, 0.5, 0.85, 1.0)


def generate_graph(n, probability, rng):
    neighbors = [0] * n
    edges = []

    for u in range(n):
        for v in range(u + 1, n):
            if rng.random() < probability:
                neighbors[u] |= 1 << v
                neighbors[v] |= 1 << u
                edges.append((u, v))

    return neighbors, edges


def brute_force_histogram(n, neighbors):
    histogram = [0] * (n + 1)
    histogram[0] = 1

    is_clique = bytearray(1 << n)
    is_clique[0] = 1

    for mask in range(1, 1 << n):
        bit = mask & -mask
        vertex = bit.bit_length() - 1
        rest = mask ^ bit

        if is_clique[rest] and (rest & neighbors[vertex]) == rest:
            is_clique[mask] = 1
            histogram[mask.bit_count()] += 1

    return histogram


def write_graph(path, n, edges):
    with path.open("w") as file:
        file.write(f"{n} {len(edges)}\n")
        for u, v in edges:
            file.write(f"{u} {v}\n")


def run_solver(path, budget):
    command = [
        str(SOLVER),
        str(path),
        "--budget", str(budget),
        "--order", "clique3",
        "--seconds", "30",
        "--node-limit", "0",
    ]
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=30,
    )

    if result.returncode != 0 or "status=COMPLETE" not in result.stdout:
        raise AssertionError((result.returncode, result.stdout, result.stderr))

    return [int(value) for value in result.stdout.split("counts=", 1)[1].split(",")]


def main():
    rng = random.Random(SEED)
    checks = 0

    with tempfile.TemporaryDirectory() as directory:
        graph_path = Path(directory) / "test.edges"

        for n in range(1, 17):
            for test_id in range(8):
                probability = DENSITIES[test_id % len(DENSITIES)]
                neighbors, edges = generate_graph(n, probability, rng)
                write_graph(graph_path, n, edges)

                expected = brute_force_histogram(n, neighbors)
                budgets = sorted({1, min(n, 5), min(n, 11), n})

                for budget in budgets:
                    actual = run_solver(graph_path, budget)
                    if actual != expected:
                        raise AssertionError(
                            f"Histogram mismatch: n={n}, p={probability}, budget={budget}"
                        )
                    checks += 1

    print(f"PASS: {checks} exact clique-histogram comparisons")


if __name__ == "__main__":
    main()
