// Algorithms computating the temporal order for various duplication-divergence models.
// Compile: g++ dd_temporal_algorithms.cpp -O3 -o ./dd_temporal_algorithms
// Run: ./dd_temporal_algorithms synthetic all MODE n n0 PARAMETERS
//   or ./dd_temporal_algorithms synthetic ALGORITHM MODE n n0 PARAMETERS

#include "./dd_temporal.h"

#include <deque>
#include <queue>

using namespace std;

typedef set<int> Bin;
typedef deque<Bin> BinningScheme;
typedef vector<pair<int, int>> PairingScheme;

const int G_TRIES = 100, SIGMA_TRIES = 10000;

enum TemporalAlgorithm {
  DEGREE_SORT, DEGREE_PEEL, NEIGHBORHOOD_SORT, NEIGHBORHOOD_PEEL, NEIGHBORHOOD_RANK,
  PROBABILITY_SORT, PROBABILITY_SUM_SORT
};

const std::map<TemporalAlgorithm, std::string> SHORT_ALGORITHM_NAME = {
  { TemporalAlgorithm::DEGREE_SORT, "sort_by_degree" },
  { TemporalAlgorithm::DEGREE_PEEL, "peel_by_degree" },
  { TemporalAlgorithm::NEIGHBORHOOD_SORT, "sort_by_neighborhood" },
  { TemporalAlgorithm::NEIGHBORHOOD_PEEL, "peel_by_neighborhood" },
  { TemporalAlgorithm::NEIGHBORHOOD_RANK, "rank_by_neighborhood" },
  { TemporalAlgorithm::PROBABILITY_SORT, "sort_by_probability" },
  { TemporalAlgorithm::PROBABILITY_SUM_SORT, "sort_by_probability_sum" },
};

const std::map<TemporalAlgorithm, std::string> LONG_ALGORITHM_NAME = {
  { TemporalAlgorithm::DEGREE_SORT, "Sort vertices by degree descending" },
  { TemporalAlgorithm::DEGREE_PEEL, "Peel vertices by degree (smallest last)" },
  { TemporalAlgorithm::NEIGHBORHOOD_SORT, "Sort vertices by neighborhood subset descending" },
  { TemporalAlgorithm::NEIGHBORHOOD_PEEL,
      "Peel vertices by neighborhood subset relation (smallest last)" },
  { TemporalAlgorithm::NEIGHBORHOOD_RANK, "Rank vertices in DAG of neighborhood subset relation" },
  { TemporalAlgorithm::PROBABILITY_SORT, "Sort vertices if p_uv > threshold" },
  { TemporalAlgorithm::PROBABILITY_SUM_SORT, "Sort vertices by sum of p_uv" },
};

const std::map<std::string, TemporalAlgorithm> REVERSE_ALGORITHM_NAME = {
  { "sort_by_degree", TemporalAlgorithm::DEGREE_SORT },
  { "peel_by_degree", TemporalAlgorithm::DEGREE_PEEL },
  { "sort_by_neighborhood", TemporalAlgorithm::NEIGHBORHOOD_SORT },
  { "peel_by_neighborhood", TemporalAlgorithm::NEIGHBORHOOD_PEEL },
  { "rank_by_neighborhood", TemporalAlgorithm::NEIGHBORHOOD_RANK },
  { "sort_by_probability", TemporalAlgorithm::PROBABILITY_SORT },
  { "sort_by_probability_sum", TemporalAlgorithm::PROBABILITY_SUM_SORT },
};

class DAG {
 private:
  vector<set<int>> H;
  vector<int> sources;
  int n0;

 public:
  DAG(Graph &G, const int &n0_value)
      : H(get_graph_size(G)), sources(get_graph_size(G)), n0(n0_value) { }

  void add_edge(const int &u, const int &v) {
    H[u].insert(v), ++sources[v];
  }

  set<int> get_neighbors(const int &v) {
    return H[v];
  }

  set<int> get_sources() {
    set<int> out;
    for (size_t i = n0; i < sources.size(); i++) {
      if (is_source(i)) {
        out.insert(i);
      }
    }
    return out;
  }

  bool decrement_source(const int &v) {
    return --sources[v];
  }

  bool is_source(const int &v) {
    return sources[v] == 0;
  }
};

BinningScheme sort_by_degree(const Graph &G, const int &n0) {
  int n = get_graph_size(G);
  BinningScheme out(n);
  for (const auto &v : get_vertices(G)) {
    int index = get_index(G, v);
    if (index < n0) {
      continue;
    }
    out[n - 1 - get_degree(G, v)].insert(index);
  }
  return out;
}

BinningScheme peel_by_degree(Graph &G, const int &n0) {
  BinningScheme out;
  while (get_graph_size(G) > n0) {
    set<Vertex> best_vertices;
    int min_degree = numeric_limits<int>::max();
    for (const auto &v : get_vertices(G)) {
      if (get_index(G, v) < n0) {
        continue;
      }
      if (get_degree(G, v) < min_degree) {
        best_vertices.clear(), min_degree = get_degree(G, v);
      }
      if (get_degree(G, v) == min_degree) {
        best_vertices.insert(v);
      }
    }
    out.push_front(Bin());
    Bin &current_bin = out.front();
    for (auto v : best_vertices) {
      current_bin.insert(get_index(G, v)), delete_vertex(G, v);
    }
  }
  return out;
}

PairingScheme sort_by_neighborhood(Graph &G, const int &n0) {
  PairingScheme out;
  auto V(get_vertices(G));
  for (size_t i = 0; i < V.size(); i++) {
    int u = get_index(G, V[i]);
    if (u < n0) {
      continue;
    }
    auto Nu = get_neighbors(G, V[i]);
    for (size_t j = i + 1; j < V.size(); j++) {
      int v = get_index(G, V[j]);
      if (v < n0) {
        continue;
      }
      auto Nv = get_neighbors(G, V[j]);
      if (Nu.size() < Nv.size() && includes(Nv.begin(), Nv.end(), Nu.begin(), Nu.end())) {
        out.push_back(make_pair(v, u));
      } else if (Nv.size() < Nu.size() && includes(Nu.begin(), Nu.end(), Nv.begin(), Nv.end())) {
        out.push_back(make_pair(u, v));
      }
    }
  }
  return out;
}

BinningScheme peel_by_neighborhood(Graph &G, const int &n0) {
  BinningScheme out;
  while (get_graph_size(G) > n0) {
    set<Vertex> best_vertices;
    for (const auto &u : get_vertices(G)) {
      if (get_index(G, u) < n0) {
        continue;
      }
      auto Nu = get_neighbors(G, u);
      bool has_descendants = false;
      for (const auto &v : get_vertices(G)) {
        if (u == v || get_index(G, v) < n0) {
          continue;
        }
        auto Nv = get_neighbors(G, v);
        if (Nv.size() < Nu.size() && includes(Nu.begin(), Nu.end(), Nv.begin(), Nv.end())) {
          has_descendants = true;
          break;
        }
      }
      if (!has_descendants) {
        best_vertices.insert(u);
      }
    }
    out.push_front(Bin());
    Bin &current_bin = out.front();
    for (auto v : best_vertices) {
      current_bin.insert(get_index(G, v)), delete_vertex(G, v);
    }
  }
  return out;
}

BinningScheme get_rank_from_DAG(DAG &H) {
  BinningScheme out;
  Bin sources = H.get_sources();
  out.push_front(sources);
  while (!out.front().empty()) {
    Bin current_bin;
    for (const auto &v : out.front()) {
      for (const auto &u : H.get_neighbors(v)) {
        H.decrement_source(u);
        if (H.is_source(u)) {
          current_bin.insert(u);
        }
      }
    }
    out.push_front(current_bin);
  }
  return out;
}

BinningScheme rank_by_neighborhood(Graph &G, const int &n0) {
  DAG H(G, n0);
  auto V(get_vertices(G));
  for (size_t i = 0; i < V.size(); i++) {
    int u = get_index(G, V[i]);
    if (u < n0) {
      continue;
    }
    auto Nu = get_neighbors(G, V[i]);
    for (size_t j = i + 1; j < V.size(); j++) {
      int v = get_index(G, V[j]);
      if (v < n0) {
        continue;
      }
      auto Nv = get_neighbors(G, V[j]);
      if (Nu.size() < Nv.size() && includes(Nv.begin(), Nv.end(), Nu.begin(), Nu.end())) {
        H.add_edge(u, v);
      } else if (Nv.size() < Nu.size() && includes(Nu.begin(), Nu.end(), Nv.begin(), Nv.end())) {
        H.add_edge(v, u);
      }
    }
  }
  return get_rank_from_DAG(H);
}

PairingScheme sort_by_probability(
    Graph &G, const int &n0, const Parameters &params, const double &threshold) {
  const auto permutations =
      get_permutation_probabilities_sampling(G, n0, params, SamplingMethod::UNIFORM, SIGMA_TRIES);
  const auto p_uv = get_p_uv_from_permutations(permutations, get_graph_size(G), n0);
  int n = get_graph_size(G);

  PairingScheme out;
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      if (i == j) {
        continue;
      }
      const auto &p_ij = p_uv.find(make_pair(i, j));
      if (p_ij != p_uv.end() && p_ij->second > threshold) {
        out.push_back(make_pair(i, j));
      }
    }
  }
  return out;
}

BinningScheme sort_by_probability_sum(Graph &G, const int &n0, const Parameters &params) {
  const auto permutations =
      get_permutation_probabilities_sampling(G, n0, params, SamplingMethod::UNIFORM, SIGMA_TRIES);
  const auto p_uv = get_p_uv_from_permutations(permutations, get_graph_size(G), n0);
  int n = get_graph_size(G);

  priority_queue<pair<long double, int>> p_v;
  for (int i = n0; i < n; i++) {
    long double p_i = 0;
    for (int j = n0; j < n; j++) {
      if (i == j) {
        continue;
      }
      const auto &p_ij = p_uv.find(make_pair(i, j));
      p_i += p_ij != p_uv.end() ? p_ij->second : 0.0L;
    }
    p_v.push(make_pair(p_i, i));
  }

  BinningScheme out;
  while (!p_v.empty()) {
    const auto &p_i = p_v.top();
    Bin current_bin{p_i.second};
    out.push_back(current_bin);
    p_v.pop();
  }
  return out;
}

DensityPrecision get_density_precision(const PairingScheme &solution, const int &count) {
  double total = solution.size(), correct = 0;
  for (const auto &uv : solution) {
    if (uv.first < uv.second) {
      correct++;
    }
  }
  return make_pair(total / (count * (count - 1) / 2), correct / total);
}

DensityPrecision get_density_precision(const BinningScheme &solution, const int &count) {
  double total = 0, correct = 0;
  for (size_t i = 0; i < solution.size(); i++) {
    for (size_t j = i + 1; j < solution.size(); j++) {
      total += solution[i].size() * solution[j].size();
      for (const auto &u : solution[i]) {
        for (const auto &v : solution[j]) {
          if (u < v) {
            correct++;
          }
        }
      }
    }
  }
  return make_pair(total / (count * (count - 1) / 2), correct / total);
}

DensityPrecision temporal_algorithm_single(
    const Graph &G0, const int &n, const Parameters &params, const TemporalAlgorithm &algorithm) {
  Graph G(G0);
  generate_graph(G, n, params);

  int n0 = get_graph_size(G0);
  switch (algorithm) {
    case DEGREE_SORT:
      return get_density_precision(sort_by_degree(G, n0), n - n0);
    case DEGREE_PEEL:
      return get_density_precision(peel_by_degree(G, n0), n - n0);
    case NEIGHBORHOOD_SORT:
      return get_density_precision(sort_by_neighborhood(G, n0), n - n0);
    case NEIGHBORHOOD_PEEL:
      return get_density_precision(peel_by_neighborhood(G, n0), n - n0);
    case NEIGHBORHOOD_RANK:
      return get_density_precision(rank_by_neighborhood(G, n0), n - n0);
    case PROBABILITY_SORT:
      // TODO(kturowski): parametrize by different values of params than used to generate G
      // TODO(kturowski): parametrize by different values of threshold than 0.5
      return get_density_precision(sort_by_probability(G, n0, params, 0.5), n - n0);
    case PROBABILITY_SUM_SORT:
      // TODO(kturowski): parametrize by different values of params than used to generate G
      return get_density_precision(sort_by_probability_sum(G, n0, params), n - n0);
    default:
      throw invalid_argument("Invalid algorithm: " + LONG_ALGORITHM_NAME.find(algorithm)->second);
  }
}

void print(
    const vector<DensityPrecision> solution, const int &n, const int &n0,
    const Parameters &params, const TemporalAlgorithm &algorithm, ostream &out_file,
    bool verbose = false) {
  cout << "Graph - n: " << n << ", n0: " << n0
      << ", parameters: " << params.to_string() << endl;
  cout << "Algorithm: " << LONG_ALGORITHM_NAME.find(algorithm)->second << endl;

  auto mean_density_precision =
      accumulate(
          solution.begin(), solution.end(), DensityPrecision(0.0, 0.0),
          [] (DensityPrecision &value, const DensityPrecision &density_precision) {
              value.first += density_precision.first, value.second += density_precision.second;
              return value;
          });
  cout << "Mean density: " << fixed << setw(6) << setprecision(3)
      << mean_density_precision.first / solution.size()
      << " mean precision: " << fixed << setw(6) << setprecision(3)
      << mean_density_precision.second / solution.size()
      << endl;
  if (verbose) {
    for (const auto &density_precision : solution) {
      cout << "density: " << fixed << setw(6) << setprecision(3) << density_precision.first
          << " precision: " << fixed << setw(6) << setprecision(3) << density_precision.second
          << endl;
    }
  }

  out_file << SHORT_ALGORITHM_NAME.find(algorithm)->second << " ";
  for (const auto &density_precision : solution) {
    out_file << density_precision.first << "," << density_precision.second << " ";
  }
  out_file << endl;
}

void temporal_algorithm(
    const int &n, const int &n0, const Parameters &params, const TemporalAlgorithm &algorithm,
    ostream &out_file) {
  Graph G0 = generate_seed(n0, 1.0);
  vector<DensityPrecision> solution(G_TRIES);
  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    solution[i] = temporal_algorithm_single(G0, n, params, algorithm);
    #pragma omp critical
    {
      if ((i + 1) % 1000 == 0) {
        cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
      }
    }
  }
  print(solution, n, n0, params, algorithm, out_file);
}

int main(int, char *argv[]) {
  try {
    string action(argv[1]), algorithm(argv[2]), mode(argv[3]);
    int n = stoi(argv[4]), n0 = stoi(argv[5]);
    Parameters params;
    params.initialize(mode, argv + 6);
    string name(TEMP_FOLDER + get_synthetic_filename(n, n0, params, "TA"));
    if (action != "synthetic") {
      throw invalid_argument("Invalid action: " + action);
    }
    if (algorithm == "all") {
      ofstream out_file(name, ios_base::app);
      for (const auto &alg : REVERSE_ALGORITHM_NAME) {
        temporal_algorithm(n, n0, params, alg.second, out_file);
      }
    } else if (REVERSE_ALGORITHM_NAME.count(algorithm)) {
      ofstream out_file(name, ios_base::app);
      temporal_algorithm(
          n, n0, params, REVERSE_ALGORITHM_NAME.find(algorithm)->second, out_file);
    } else {
      throw invalid_argument("Invalid algorithm: " + algorithm);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }

  return 0;
}
