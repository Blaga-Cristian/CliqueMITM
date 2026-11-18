import csv
import os
import re
import statistics
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GRAPHS_DIR = ROOT / "data" / "graphs"
RESULTS_DIR = ROOT / "results"
SOURCE_DIR = ROOT / "src"
BUILD_DIR = ROOT / "build"
PIVOTER_DIR = ROOT.parent / "Pivoter"

EXE = ".exe" if os.name == "nt" else ""
BALANCED = BUILD_DIR / f"balanced_mitm{EXE}"
BLOCK = BUILD_DIR / f"block_mitm{EXE}"
PIVOTER = PIVOTER_DIR / "bin" / f"degeneracy_cliques{EXE}"

REPETITIONS = 3
PROCESS_TIMEOUT_SECONDS = 8
BLOCK_TIMEOUT_SECONDS = 5
BLOCK_NODE_LIMIT = 100_000_000


FIELDS = [
    "graph", "solver", "median_wall_ms", "status",
    "nodes", "histogram", "input_file",
]


def generate_hamming_graph():
    """Construct the DIMACS-style Hamming(6,4) graph in local data/graphs."""
    n = 64
    edges = [
        (u, v)
        for u in range(n)
        for v in range(u + 1, n)
        if (u ^ v).bit_count() >= 4
    ]

    path = GRAPHS_DIR / "hamming6-4_constructed.edges"
    with path.open("w", encoding="utf-8") as file:
        file.write(f"{n} {len(edges)}\n")
        for u, v in edges:
            file.write(f"{u} {v}\n")
    return path


def compile_cpp(source, executable):
    if not source.is_file():
        raise FileNotFoundError(f"Missing C++ source: {source}")

    command = [
        "g++", "-std=c++20", "-O3", "-march=native",
        str(source), "-o", str(executable),
    ]
    print(f"Building {source.name} ...", flush=True)
    subprocess.run(command, check=True)


def compile_solvers():
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    compile_cpp(SOURCE_DIR / "balanced_mitm.cpp", BALANCED)
    compile_cpp(SOURCE_DIR / "block_mitm.cpp", BLOCK)


def solver_command(solver, graph):
    if solver == "balanced-DP":
        return [
            str(BALANCED), str(graph), "mitm-dp",
            "--split", "14", "--order", "clique3",
        ], ROOT

    if solver == "block-clique3":
        return [
            str(BLOCK), str(graph),
            "--budget", "16", "--order", "clique3",
            "--seconds", str(BLOCK_TIMEOUT_SECONDS),
            "--node-limit", str(BLOCK_NODE_LIMIT),
        ], ROOT

    if solver == "Pivoter":
        return [
            str(PIVOTER), "-i", str(graph),
            "-t", "A", "-k", "0", "-d", "0",
        ], PIVOTER_DIR

    raise ValueError(f"Unknown solver: {solver}")


def parse_histogram(output, solver, n):
    if solver == "Pivoter":
        counts = [0] * (n + 1)
        # Pivoter prints rows like: 3, 960.000000
        matches = re.findall(
            r"^\s*(\d+)\s*,\s*(\d+)(?:\.(\d+))?\s*$",
            output,
            flags=re.MULTILINE,
        )
        if not matches:
            raise ValueError("Could not parse Pivoter's clique counts")

        for size_text, value_text, decimals in matches:
            size = int(size_text)
            if size > n or (decimals and set(decimals) != {"0"}):
                raise ValueError("Unexpected Pivoter histogram entry")
            counts[size] = int(value_text)
        if counts[0] != 1:
            raise ValueError("Pivoter did not report the empty clique")
        return counts

    match = re.search(r"counts=([\d,]+)", output)
    if not match:
        raise ValueError(f"Could not parse {solver}'s clique counts")

    counts = [int(x) for x in match.group(1).strip(",").split(",")]
    if len(counts) > n + 1 or not counts or counts[0] != 1:
        raise ValueError(f"Invalid histogram from {solver}")
    return counts + [0] * (n + 1 - len(counts))


def graph_vertex_count(path):
    with path.open(encoding="utf-8") as file:
        for line in file:
            fields = line.split()
            if not fields or fields[0] in ("#", "c"):
                continue
            return int(fields[2] if fields[0] == "p" else fields[0])
    raise ValueError(f"No graph header in {path}")


def benchmark(graph, solver):
    command, cwd = solver_command(solver, graph)
    n = graph_vertex_count(graph)
    times_ms = []
    histogram = None
    nodes = ""
    status = "complete"

    for run_index in range(1, REPETITIONS + 1):
        started = time.perf_counter()
        try:
            result = subprocess.run(
                command,
                cwd=cwd,
                capture_output=True,
                text=True,
                timeout=PROCESS_TIMEOUT_SECONDS,
            )
        except subprocess.TimeoutExpired:
            status = "TIMEOUT"
            break
        except OSError as error:
            status = f"EXEC_ERROR: {error}"
            break

        elapsed_ms = (time.perf_counter() - started) * 1000

        if result.returncode != 0:
            message = result.stderr.strip().splitlines()
            status = "ERROR: " + (message[-1] if message else str(result.returncode))
            break

        if solver == "block-clique3" and "status=COMPLETE" not in result.stdout:
            status = "INTERNAL_LIMIT"
            break

        try:
            current = parse_histogram(result.stdout, solver, n)
        except ValueError as error:
            status = f"PARSE_ERROR: {error}"
            break

        if histogram is not None and histogram != current:
            raise RuntimeError(f"{solver} returned different counts on repeated runs")

        histogram = current
        times_ms.append(elapsed_ms)

        match = re.search(r"(?:nodes|states)=(\d+)", result.stdout)
        if match:
            nodes = match.group(1)

        print(f"    Run {run_index}/{REPETITIONS}: {elapsed_ms:.2f} ms", flush=True)

    return {
        "solver": solver,
        "median_wall_ms": round(statistics.median(times_ms), 3) if status == "complete" else "",
        "status": status,
        "nodes": nodes,
        "histogram": ";".join(map(str, histogram)) if status == "complete" else "",
    }


def validate_graph_results(results):
    completed = [r for r in results if r["status"] == "complete"]
    if len(completed) < 2:
        return

    expected = completed[0]
    for other in completed[1:]:
        if other["histogram"] != expected["histogram"]:
            raise RuntimeError(
                f"Mismatch: {expected['solver']} and {other['solver']} "
                f"on {expected['graph']}"
            )
    print("  Exact histogram comparison: OK", flush=True)


def save_results(rows):
    path = RESULTS_DIR / "short_runs.csv"
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nSaved results to: {path}")


def main():
    GRAPHS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    hamming = generate_hamming_graph()
    graphs = [
        ("Moon-Moser 39 shuffled", GRAPHS_DIR / "MM39_shuffle29.edges"),
        ("ER40 p≈0.94", GRAPHS_DIR / "ER40_p94_s11.edges"),
        ("Hamming6-4 (constructed)", hamming),
        ("Random 80 p≈0.7", GRAPHS_DIR / "rand_80_70.edges"),
    ]

    for name, path in graphs:
        if not path.is_file():
            raise FileNotFoundError(f"Graph missing ({name}): {path}")

    compile_solvers()
    pivoter_available = PIVOTER.is_file()
    if not pivoter_available:
        print(
            f"\nPivoter is not built for this environment: {PIVOTER}\n"
            "Skipping Pivoter runs; compile it separately to include them.\n",
            flush=True,
        )

    all_rows = []
    for name, graph in graphs:
        print(f"\n{name}", flush=True)
        solvers = ["block-clique3"]
        if name in ("Moon-Moser 39 shuffled", "ER40 p≈0.94"):
            solvers.insert(0, "balanced-DP")
        if pivoter_available:
            solvers.append("Pivoter")

        current_rows = []
        for solver in solvers:
            print(f"  {solver}", flush=True)
            result = benchmark(graph, solver)
            row = {
                "graph": name,
                **result,
                "input_file": str(graph.relative_to(ROOT)),
            }
            current_rows.append(row)
            all_rows.append(row)
            print(f"    Result: {result['status']}", flush=True)

        validate_graph_results(current_rows)
        save_results(all_rows)

    print(f"\nFinished: {len(all_rows)} benchmark records.")


if __name__ == "__main__":
    main()

