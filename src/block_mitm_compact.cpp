#include <bits/stdc++.h>

using namespace std;

using u128 = unsigned __int128;
using u64 = uint64_t;
using u32 = uint32_t;

string toDecimal(u128 value) {
    if (value == 0) return "0";

    string result;
    while (value > 0) {
        result.push_back('0' + value % 10);
        value /= 10;
    }
    reverse(result.begin(), result.end());
    return result;
}

int firstBit(u128 mask) {
    u64 low = static_cast<u64>(mask);
    return low ? __builtin_ctzll(low)
               : 64 + __builtin_ctzll(static_cast<u64>(mask >> 64));
}

int popcount(u128 mask) {
    return __builtin_popcountll(static_cast<u64>(mask)) +
           __builtin_popcountll(static_cast<u64>(mask >> 64));
}

struct Graph {
    int n;
    vector<u128> adj;
    long long m = 0;

    explicit Graph(int size = 0) : n(size), adj(size, 0) {
        if (n > 120) {
            throw runtime_error("Maximum 120 vertices (128-bit clique counter)");
        }
    }

    void addEdge(int u, int v) {
        if (u < 0 || v < 0 || u >= n || v >= n) {
            throw runtime_error("edge index out of range");
        }
        if (u == v) return;

        u128 bit = static_cast<u128>(1) << v;
        if ((adj[u] & bit) == 0) {
            adj[u] |= bit;
            adj[v] |= static_cast<u128>(1) << u;
            ++m;
        }
    }

    u128 allVertices() const {
        return n ? (static_cast<u128>(1) << n) - 1 : 0;
    }
};

Graph readGraph(const string& path) {
    ifstream input(path);
    if (!input) throw runtime_error("Can't open " + path);

    string line;
    int n = -1, m = -1;

    while (getline(input, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'c') continue;

        istringstream stream(line);
        if (line[0] == 'p') {
            string marker, format;
            stream >> marker >> format >> n >> m;
        } else {
            stream >> n >> m;
        }
        break;
    }

    if (n < 0 || m < 0) throw runtime_error("Bad graph header");

    Graph graph(n);
    for (int i = 0; i < m; ++i) {
        int u, v;
        if (!(input >> u >> v)) throw runtime_error("Bad edge list");
        graph.addEdge(u, v);
    }
    return graph;
}

struct Timeout {};

struct Solver {
    const Graph& graph;
    int budget;
    string order;
    double cutoffSeconds;
    u64 nodeLimit;

    vector<u128> counts;
    u128 remaining;
    u64 nodes = 0;
    u64 blocks = 0;
    u64 preprocessEntries = 0;
    bool timedOut = false;
    double selectionMs = 0;
    double dpMs = 0;
    double dfsMs = 0;

    chrono::steady_clock::time_point startedAt;

    vector<u32> crossNeighbors;
    vector<u32> internalNeighbors;
    vector<u32> dp;
    vector<u32> offsets;
    int supportSize = 0;
    u32 fullMask = 0;

    Solver(const Graph& graph, int budget, string order, double seconds, u64 limit)
        : graph(graph),
          budget(budget),
          order(move(order)),
          cutoffSeconds(seconds),
          nodeLimit(limit),
          counts(graph.n + 1, 0),
          remaining(graph.allVertices()) {
        if (budget < 1 || budget > 25) {
            throw runtime_error("--budget must be 1..25 for this implementation");
        }
    }

    void checkLimits() const {
        if (nodeLimit && nodes > nodeLimit) throw Timeout{};
        if (cutoffSeconds > 0 &&
            chrono::duration<double>(chrono::steady_clock::now() - startedAt).count() > cutoffSeconds) {
            throw Timeout{};
        }
    }

    // Count cliques of sizes 0, 1, 2, and 3 in an induced vertex subset.
    u64 countUpToThree(u128 vertices) const {
        u64 count = 1;
        u128 a = vertices;

        while (a) {
            int v = firstBit(a);
            a &= a - 1;
            ++count;

            u128 b = a & graph.adj[v];
            while (b) {
                int w = firstBit(b);
                b &= b - 1;
                ++count;
                count += popcount(b & graph.adj[w]);
            }
        }
        return count;
    }

    vector<int> chooseSupport() const {
        vector<int> vertices;
        u128 candidates = remaining;
        while (candidates) {
            int v = firstBit(candidates);
            candidates &= candidates - 1;
            vertices.push_back(v);
        }

        if (order == "peel" || order == "clique3") {
            int leftTarget = max(0, static_cast<int>(vertices.size()) - budget);
            u128 left = 0;
            u128 pool = remaining;

            while (popcount(left) < leftTarget) {
                u64 bestScore = ULLONG_MAX;
                int best = -1;
                int bestDegree = INT_MAX;
                u128 choices = pool;

                while (choices) {
                    int v = firstBit(choices);
                    choices &= choices - 1;

                    int degree = popcount(graph.adj[v] & remaining);
                    u64 score = order == "peel"
                        ? popcount(graph.adj[v] & pool)
                        : countUpToThree(graph.adj[v] & left);

                    if (score < bestScore ||
                        (score == bestScore && degree < bestDegree) ||
                        (score == bestScore && degree == bestDegree && (best == -1 || v < best))) {
                        bestScore = score;
                        bestDegree = degree;
                        best = v;
                    }
                }

                left |= static_cast<u128>(1) << best;
                pool &= ~(static_cast<u128>(1) << best);
            }

            vector<int> support;
            while (pool) {
                int v = firstBit(pool);
                pool &= pool - 1;
                support.push_back(v);
            }
            return support;
        }

        vector<pair<int, int>> ranked;
        for (int v : vertices) {
            ranked.push_back({popcount(remaining & graph.adj[v]), v});
        }
        stable_sort(ranked.begin(), ranked.end(), [&](auto a, auto b) {
            if (a.first != b.first) {
                return order == "low" ? a.first < b.first : a.first > b.first;
            }
            return a.second < b.second;
        });

        vector<int> support;
        for (int i = 0; i < min(budget, static_cast<int>(ranked.size())); ++i) {
            support.push_back(ranked[i].second);
        }
        return support;
    }

    void buildDP(const vector<int>& support) {
        supportSize = static_cast<int>(support.size());
        fullMask = (1U << supportSize) - 1;

        internalNeighbors.assign(supportSize, 0);
        crossNeighbors.assign(graph.n, 0);

        for (int j = 0; j < supportSize; ++j) {
            for (int i = 0; i < supportSize; ++i) {
                if (graph.adj[support[j]] & (static_cast<u128>(1) << support[i])) {
                    internalNeighbors[j] |= 1U << i;
                }
            }
            for (int v = 0; v < graph.n; ++v) {
                if (graph.adj[v] & (static_cast<u128>(1) << support[j])) {
                    crossNeighbors[v] |= 1U << j;
                }
            }
        }

        size_t subsetCount = 1ULL << supportSize;

        // Each mask only needs coefficients of degree 0..popcount(mask).
        offsets.resize(subsetCount + 1);
        offsets[0] = 0;
        for (size_t mask = 0; mask < subsetCount; ++mask) {
            offsets[mask + 1] = offsets[mask] + 1 + __builtin_popcountll(mask);
        }

        dp.resize(offsets[subsetCount]);
        dp[0] = 1;

        for (size_t mask = 1; mask < subsetCount; ++mask) {
            int v = __builtin_ctzll(mask);
            size_t smaller = mask & (mask - 1);
            size_t compatible = smaller & internalNeighbors[v];

            const u32* without = &dp[offsets[smaller]];
            const u32* with = &dp[offsets[compatible]];
            u32* dest = &dp[offsets[mask]];

            int size = __builtin_popcountll(mask);
            int compatibleSize = __builtin_popcountll(compatible);
            dest[0] = 1;

            for (int k = 1; k <= size; ++k) {
                u32 skip = k < size ? without[k] : 0;
                u32 take = k - 1 <= compatibleSize ? with[k - 1] : 0;
                dest[k] = skip + take;
            }
        }
        preprocessEntries += subsetCount;
    }

    void dfs(u128 candidates, u32 compatible, int selected) {
        if (compatible == 0) return;

        ++nodes;
        if ((nodes & 4095) == 0) checkLimits();

        const u32* histogram = &dp[offsets[compatible]];
        int maxSize = __builtin_popcount(compatible);
        for (int k = 1; k <= maxSize; ++k) {
            counts[selected + k] += histogram[k];
        }

        while (candidates) {
            int v = firstBit(candidates);
            candidates &= candidates - 1;

            u32 nextCompatible = compatible & crossNeighbors[v];
            if (nextCompatible == 0) continue;

            dfs(candidates & graph.adj[v], nextCompatible, selected + 1);
        }
    }

    void run() {
        startedAt = chrono::steady_clock::now();
        counts[0] = 1;

        try {
            while (remaining) {
                checkLimits();

                auto selectStart = chrono::steady_clock::now();
                vector<int> support = chooseSupport();
                auto selectEnd = chrono::steady_clock::now();
                selectionMs += chrono::duration<double, milli>(selectEnd - selectStart).count();

                u128 supportMask = 0;
                for (int v : support) {
                    supportMask |= static_cast<u128>(1) << v;
                }

                auto dpStart = chrono::steady_clock::now();
                buildDP(support);
                auto dpEnd = chrono::steady_clock::now();
                dpMs += chrono::duration<double, milli>(dpEnd - dpStart).count();

                u128 left = remaining & ~supportMask;
                ++blocks;
                try {
                    dfs(left, fullMask, 0);
                } catch (const Timeout&) {
                    dfsMs += chrono::duration<double, milli>(
                        chrono::steady_clock::now() - dpEnd
                    ).count();
                    throw;
                }
                dfsMs += chrono::duration<double, milli>(
                    chrono::steady_clock::now() - dpEnd
                ).count();

                remaining = left;
            }
        } catch (const Timeout&) {
            timedOut = true;
        }
    }

    double elapsedMs() const {
        return chrono::duration<double, milli>(chrono::steady_clock::now() - startedAt).count();
    }
};

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            cerr << "Usage: block_mitm graph.edges [--budget 16] "
                 << "[--order high|low|peel|clique3] [--seconds 3] "
                 << "[--node-limit 10000000]\n";
            return 1;
        }

        string file = argv[1];
        string order = "high";
        int budget = 16;
        double seconds = 5;
        u64 nodeLimit = 10000000;

        for (int i = 2; i < argc; ++i) {
            string option = argv[i];
            if (i + 1 >= argc) throw runtime_error("missing value for " + option);
            string value = argv[++i];

            if (option == "--budget") budget = stoi(value);
            else if (option == "--order") order = value;
            else if (option == "--seconds") seconds = stod(value);
            else if (option == "--node-limit") nodeLimit = stoull(value);
            else throw runtime_error("bad argument " + option);
        }

        if (order != "high" && order != "low" && order != "peel" && order != "clique3") {
            throw runtime_error("order must be high, low, peel or clique3");
        }

        Graph graph = readGraph(file);
        Solver solver(graph, budget, order, seconds, nodeLimit);
        solver.run();

        cout << "status=" << (solver.timedOut ? "LIMIT" : "COMPLETE")
             << " n=" << graph.n << " m=" << graph.m
             << " ms=" << fixed << setprecision(4) << solver.elapsedMs()
             << " nodes=" << solver.nodes << " blocks=" << solver.blocks
             << " dp_states=" << solver.preprocessEntries
             << " budget=" << budget << " order=" << order
             << " selection_ms=" << solver.selectionMs
             << " dp_ms=" << solver.dpMs
             << " dfs_ms=" << solver.dfsMs << " counts=";

        for (int i = 0; i <= graph.n; ++i) {
            if (i) cout << ',';
            cout << toDecimal(solver.counts[i]);
        }
        cout << '\n';
        return 0;
    } catch (const exception& error) {
        cerr << error.what() << '\n';
        return 2;
    }
}
