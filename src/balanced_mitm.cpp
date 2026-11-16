#include <bits/stdc++.h>

using namespace std;

using u128 = unsigned __int128;
using u64 = uint64_t;

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

struct Graph {
    int n;
    vector<u64> adj;
    long long m = 0;

    explicit Graph(int size = 0) : n(size), adj(size, 0) {
        if (n > 63) {
            throw runtime_error("This implementation supports n <= 63");
        }
    }

    void addEdge(int u, int v) {
        if (u == v) return;
        if (u < 0 || v < 0 || u >= n || v >= n) {
            throw runtime_error("Bad edge");
        }
        if ((adj[u] & (1ULL << v)) == 0) {
            adj[u] |= 1ULL << v;
            adj[v] |= 1ULL << u;
            ++m;
        }
    }

    u64 allVertices() const {
        return n == 0 ? 0 : (1ULL << n) - 1;
    }
};

Graph readGraph(const string& path) {
    ifstream input(path);
    if (!input) throw runtime_error("Cannot read graph file " + path);

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

    if (n < 0 || m < 0) throw runtime_error("Invalid header");

    Graph graph(n);
    for (int i = 0; i < m; ++i) {
        int u, v;
        if (!(input >> u >> v)) throw runtime_error("Missing edges");
        graph.addEdge(u, v);
    }
    return graph;
}

struct Stats {
    vector<u128> counts;
    u64 states = 0;

    explicit Stats(int n) : counts(n + 1, 0) {}
};

void enumerateCliques(const Graph& graph, Stats& stats, u64 candidates, int size) {
    ++stats.counts[size];
    ++stats.states;

    while (candidates) {
        int v = __builtin_ctzll(candidates);
        candidates &= candidates - 1;
        enumerateCliques(graph, stats, candidates & graph.adj[v], size + 1);
    }
}

Stats exhaustiveEnumeration(const Graph& graph) {
    Stats stats(graph.n);
    enumerateCliques(graph, stats, graph.allVertices(), 0);
    return stats;
}

// A compact Pivoter-style counting baseline using a pivoted search tree.
struct Pivoter {
    const Graph& graph;
    Stats stats;
    vector<vector<u128>> choose;

    explicit Pivoter(const Graph& graph)
        : graph(graph),
          stats(graph.n),
          choose(graph.n + 1, vector<u128>(graph.n + 1, 0)) {
        for (int i = 0; i <= graph.n; ++i) {
            choose[i][0] = choose[i][i] = 1;
            for (int j = 1; j < i; ++j) {
                choose[i][j] = choose[i - 1][j - 1] + choose[i - 1][j];
            }
        }
    }

    void dfs(u64 candidates, int held, int optional) {
        ++stats.states;
        if (candidates == 0) {
            for (int k = 0; k <= optional; ++k) {
                stats.counts[held + k] += choose[optional][k];
            }
            return;
        }

        int pivot = -1;
        int bestDegree = -1;
        u64 remaining = candidates;

        while (remaining) {
            int v = __builtin_ctzll(remaining);
            remaining &= remaining - 1;
            int degree = __builtin_popcountll(candidates & graph.adj[v]);

            if (degree > bestDegree) {
                bestDegree = degree;
                pivot = v;
            }
        }

        dfs(candidates & graph.adj[pivot], held, optional + 1);

        u64 excluded = candidates & ~(graph.adj[pivot] | (1ULL << pivot));
        remaining = candidates;

        while (excluded) {
            int v = __builtin_ctzll(excluded);
            excluded &= excluded - 1;
            remaining &= ~(1ULL << v);
            dfs(remaining & graph.adj[v], held + 1, optional);
        }
    }

    Stats run() {
        stats.counts[0] = 1;

        // Peel vertices by minimum degree to obtain a degeneracy ordering.
        u64 remaining = graph.allVertices();
        vector<int> ordering;
        ordering.reserve(graph.n);

        while (remaining) {
            int best = -1;
            int minDegree = INT_MAX;
            u64 candidates = remaining;

            while (candidates) {
                int v = __builtin_ctzll(candidates);
                candidates &= candidates - 1;
                int degree = __builtin_popcountll(graph.adj[v] & remaining);

                if (degree < minDegree) {
                    minDegree = degree;
                    best = v;
                }
            }

            ordering.push_back(best);
            remaining &= ~(1ULL << best);
        }

        remaining = graph.allVertices();
        for (int v : ordering) {
            remaining &= ~(1ULL << v);
            dfs(graph.adj[v] & remaining, 1, 0);
        }
        return move(stats);
    }
};

struct MITM {
    const Graph& graph;
    Stats stats;
    int leftSize, rightSize;
    bool useDirectDP;
    size_t rowWidth;
    vector<u128> table;
    vector<uint32_t> rightNeighbors;
    vector<uint32_t> crossNeighbors;

    MITM(const Graph& graph, int split = -1, bool directDP = false)
        : graph(graph),
          stats(graph.n),
          rightSize(split < 0 ? graph.n / 2 : split),
          useDirectDP(directDP) {
        leftSize = graph.n - rightSize;
        if (rightSize < 0 || rightSize > 21 || leftSize > 63) {
            throw runtime_error("MITM: support half must have <=21 vertices");
        }
        rowWidth = rightSize + 1;
    }

    void buildTable() {
        size_t subsetCount = 1ULL << rightSize;
        table.assign(subsetCount * rowWidth, 0);
        rightNeighbors.resize(rightSize);
        crossNeighbors.resize(leftSize);

        uint32_t full = rightSize ? (1U << rightSize) - 1 : 0;
        for (int i = 0; i < rightSize; ++i) {
            rightNeighbors[i] = (graph.adj[leftSize + i] >> leftSize) & full;
        }
        for (int i = 0; i < leftSize; ++i) {
            crossNeighbors[i] = (graph.adj[i] >> leftSize) & full;
        }

        table[0] = 1;
        if (useDirectDP) {
            for (size_t mask = 1; mask < subsetCount; ++mask) {
                int v = __builtin_ctzll(mask);
                size_t smaller = mask & (mask - 1);

                u128* dest = &table[mask * rowWidth];
                const u128* without = &table[smaller * rowWidth];
                const u128* with = &table[(smaller & rightNeighbors[v]) * rowWidth];
                int size = __builtin_popcountll(mask);

                dest[0] = without[0];
                for (int k = 1; k <= size; ++k) {
                    dest[k] = without[k] + with[k - 1];
                }
            }
            return;
        }

        // Alternatively, mark clique subsets and apply the subset zeta transform.
        vector<uint8_t> isClique(subsetCount);
        isClique[0] = 1;

        for (size_t mask = 1; mask < subsetCount; ++mask) {
            int v = __builtin_ctzll(mask);
            size_t smaller = mask & (mask - 1);
            isClique[mask] = isClique[smaller] && ((smaller & ~rightNeighbors[v]) == 0);
            if (isClique[mask]) {
                table[mask * rowWidth + __builtin_popcountll(mask)] = 1;
            }
        }

        for (int bit = 0; bit < rightSize; ++bit) {
            for (size_t mask = 0; mask < subsetCount; ++mask) {
                if ((mask & (1ULL << bit)) == 0) continue;

                u128* dest = &table[mask * rowWidth];
                const u128* source = &table[(mask ^ (1ULL << bit)) * rowWidth];
                int maxSize = __builtin_popcountll(mask) - 1;
                for (int k = 0; k <= maxSize; ++k) {
                    dest[k] += source[k];
                }
            }
        }
    }

    void dfs(u64 candidates, uint32_t allowed, int selected) {
        ++stats.states;
        const u128* counts = &table[static_cast<size_t>(allowed) * rowWidth];

        for (int k = 0; k <= rightSize; ++k) {
            stats.counts[selected + k] += counts[k];
        }

        while (candidates) {
            int v = __builtin_ctzll(candidates);
            candidates &= candidates - 1;
            dfs(candidates & graph.adj[v], allowed & crossNeighbors[v], selected + 1);
        }
    }

    Stats run() {
        buildTable();
        u64 left = leftSize ? (1ULL << leftSize) - 1 : 0;
        uint32_t right = rightSize ? (1U << rightSize) - 1 : 0;
        dfs(left, right, 0);
        return move(stats);
    }
};

// Counts all cliques of sizes 0..maxSize in the induced graph on `vertices`.
array<uint64_t, 5> smallCliqueCounts(const Graph& graph, u64 vertices, int maxSize) {
    array<uint64_t, 5> counts{};
    counts[0] = 1;

    u64 a = vertices;
    while (a) {
        int v = __builtin_ctzll(a);
        a &= a - 1;
        ++counts[1];
        if (maxSize == 1) continue;

        u64 b = a & graph.adj[v];
        while (b) {
            int w = __builtin_ctzll(b);
            b &= b - 1;
            ++counts[2];
            if (maxSize == 2) continue;

            u64 c = b & graph.adj[w];
            if (maxSize == 3) {
                counts[3] += __builtin_popcountll(c);
                continue;
            }

            while (c) {
                int x = __builtin_ctzll(c);
                c &= c - 1;
                ++counts[3];
                counts[4] += __builtin_popcountll(c & graph.adj[x]);
            }
        }
    }
    return counts;
}

// Build the explicit side greedily; the remaining vertices form the DP side.
vector<int> cliqueAwareLeft(const Graph& graph, int targetSize, int maxSize) {
    vector<int> ordering;
    ordering.reserve(graph.n);

    u64 chosen = 0;
    u64 remaining = graph.allVertices();
    vector<int> degrees(graph.n);
    for (int v = 0; v < graph.n; ++v) {
        degrees[v] = __builtin_popcountll(graph.adj[v]);
    }

    for (int i = 0; i < targetSize; ++i) {
        int best = -1;
        int bestDegree = INT_MAX;
        uint64_t bestCost = ULLONG_MAX;
        u64 candidates = remaining;

        while (candidates) {
            int v = __builtin_ctzll(candidates);
            candidates &= candidates - 1;

            auto counts = smallCliqueCounts(graph, chosen & graph.adj[v], maxSize);
            uint64_t cost = 0;
            for (int size = 0; size <= maxSize; ++size) {
                cost += counts[size];
            }

            if (cost < bestCost ||
                (cost == bestCost && degrees[v] < bestDegree) ||
                (cost == bestCost && degrees[v] == bestDegree && v < best)) {
                best = v;
                bestCost = cost;
                bestDegree = degrees[v];
            }
        }

        ordering.push_back(best);
        chosen |= 1ULL << best;
        remaining &= ~(1ULL << best);
    }

    while (remaining) {
        int v = __builtin_ctzll(remaining);
        remaining &= remaining - 1;
        ordering.push_back(v);
    }
    return ordering;
}

vector<int> degeneracyOrdering(const Graph& graph) {
    vector<int> ordering;
    ordering.reserve(graph.n);
    u64 remaining = graph.allVertices();

    while (remaining) {
        int best = -1;
        int minDegree = INT_MAX;
        u64 candidates = remaining;

        while (candidates) {
            int v = __builtin_ctzll(candidates);
            candidates &= candidates - 1;
            int degree = __builtin_popcountll(graph.adj[v] & remaining);

            if (degree < minDegree) {
                minDegree = degree;
                best = v;
            }
        }

        ordering.push_back(best);
        remaining &= ~(1ULL << best);
    }
    return ordering;
}

Graph reorderGraph(const Graph& graph, const string& mode, int split) {
    if (mode == "input") return graph;

    vector<int> ordering(graph.n);
    iota(ordering.begin(), ordering.end(), 0);

    vector<int> degrees(graph.n);
    for (int v = 0; v < graph.n; ++v) {
        degrees[v] = __builtin_popcountll(graph.adj[v]);
    }

    if (mode == "high" || mode == "low") {
        stable_sort(ordering.begin(), ordering.end(), [&](int u, int v) {
            if (degrees[u] != degrees[v]) {
                return mode == "high" ? degrees[u] < degrees[v]
                                      : degrees[u] > degrees[v];
            }
            return u < v;
        });
    } else if (mode == "degen") {
        ordering = degeneracyOrdering(graph);
    } else if (mode == "clique3" || mode == "clique4") {
        if (split < 0) {
            throw runtime_error("Clique-aware ordering requires --split");
        }
        ordering = cliqueAwareLeft(graph, graph.n - split, mode == "clique3" ? 3 : 4);
    } else {
        throw runtime_error("order must be input, high, low, degen, clique3 or clique4");
    }

    Graph reordered(graph.n);
    for (int i = 0; i < graph.n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (graph.adj[ordering[i]] & (1ULL << ordering[j])) {
                reordered.addEdge(i, j);
            }
        }
    }
    return reordered;
}

int main(int argc, char** argv) {
    try {
        if (argc < 3) {
            cerr << "Usage: clique_stats graph.edges pivoter|mitm|mitm-dp|enumerate "
                 << "[--split num] [--order input|high|low|degen|clique3|clique4]\n";
            return 1;
        }

        string method = argv[2];
        string order = "input";
        int split = -1;

        for (int i = 3; i < argc; ++i) {
            string argument = argv[i];
            if (argument == "--split" && i + 1 < argc) {
                split = stoi(argv[++i]);
            } else if (argument == "--order" && i + 1 < argc) {
                order = argv[++i];
            } else {
                throw runtime_error("Bad argument: " + argument);
            }
        }

        Graph loaded = readGraph(argv[1]);
        auto start = chrono::steady_clock::now();
        Graph graph = reorderGraph(loaded, order, split);
        auto orderedAt = chrono::steady_clock::now();

        Stats stats(graph.n);
        if (method == "pivoter") {
            stats = Pivoter(graph).run();
        } else if (method == "mitm") {
            stats = MITM(graph, split, false).run();
        } else if (method == "mitm-dp") {
            stats = MITM(graph, split, true).run();
        } else if (method == "enumerate") {
            stats = exhaustiveEnumeration(graph);
        } else {
            throw runtime_error("Unknown method");
        }

        auto finish = chrono::steady_clock::now();
        double elapsedMs = chrono::duration<double, milli>(finish - start).count();
        double orderingMs = chrono::duration<double, milli>(orderedAt - start).count();

        cout << "method=" << method << " order=" << order << " split=" << split
             << " n=" << graph.n << " m=" << graph.m
             << " ms=" << setprecision(9) << elapsedMs
             << " ordering_ms=" << orderingMs
             << " states=" << stats.states << " counts=";
        for (u128 count : stats.counts) cout << toDecimal(count) << ',';
        cout << '\n';
        return 0;
    } catch (const exception& error) {
        cerr << error.what() << '\n';
        return 2;
    }
}
