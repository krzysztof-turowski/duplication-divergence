// Algorithms computating the temporal order for various duplication-divergence models.
// Compile: g++ dd_temporal_algorithms.cpp -O3 -o ./dd_temporal_algorithms
// Run: ./dd_temporal_algorithms synthetic all MODE_GEN n n0 PARAMETERS_GEN MODE PARAMETERS
//   or ./dd_temporal_algorithms synthetic ALGORITHM MODE_GEN n n0 PARAMETERS_GEN [MODE PARAMETERS]
//   or ./dd_temporal_algorithms real_data all FILE MODE PARAMETERS
//   or ./dd_temporal_algorithms real_data ALGORITHM FILE [MODE PARAMETERS]

#include "./dd_temporal.h"

#if defined(glpk)
  #include "./dd_glpk.h"
#elif defined(gurobi)
  #include "./dd_gurobi.h"
#endif

#include <deque>
#include <queue>

using namespace std;

typedef set<int> Bin;
typedef deque<Bin> BinningScheme;
typedef vector<pair<int, int>> PairingScheme;

const int G_TRIES = 20, SIGMA_TRIES = 100000;
const int AGE_ZERO = 0, AGE_TWO = 2, AGE_MAX = numeric_limits<int>::max();

enum TemporalAlgorithm {
  DEGREE_SORT, DEGREE_PEEL, NEIGHBORHOOD_SORT, NEIGHBORHOOD_PEEL_SL, NEIGHBORHOOD_PEEL_LF,
  NEIGHBORHOOD_RANK, PROBABILITY_SORT, PROBABILITY_SUM_SORT, LP_SOLUTION_SORT
};

const map<TemporalAlgorithm, string> SHORT_ALGORITHM_NAME = {
  { TemporalAlgorithm::DEGREE_SORT, "sort_by_degree" },
  { TemporalAlgorithm::DEGREE_PEEL, "peel_by_degree" },
  { TemporalAlgorithm::NEIGHBORHOOD_SORT, "sort_by_neighborhood" },
  { TemporalAlgorithm::NEIGHBORHOOD_PEEL_SL, "peel_by_neighborhood_sl" },
  { TemporalAlgorithm::NEIGHBORHOOD_PEEL_LF, "peel_by_neighborhood_lf" },
  { TemporalAlgorithm::NEIGHBORHOOD_RANK, "rank_by_neighborhood" },
  { TemporalAlgorithm::PROBABILITY_SORT, "sort_by_probability" },
  { TemporalAlgorithm::PROBABILITY_SUM_SORT, "sort_by_probability_sum" },
  { TemporalAlgorithm::LP_SOLUTION_SORT, "sort_by_lp_solution" },
};

const map<TemporalAlgorithm, string> LONG_ALGORITHM_NAME = {
  { TemporalAlgorithm::DEGREE_SORT, "Sort vertices by degree descending" },
  { TemporalAlgorithm::DEGREE_PEEL, "Peel vertices by degree (smallest last)" },
  { TemporalAlgorithm::NEIGHBORHOOD_SORT, "Sort vertices by neighborhood subset descending" },
  { TemporalAlgorithm::NEIGHBORHOOD_PEEL_SL,
      "Peel vertices by neighborhood subset relation (smallest last)" },
  { TemporalAlgorithm::NEIGHBORHOOD_PEEL_LF,
      "Peel vertices by neighborhood subset relation (largest first)" },
  { TemporalAlgorithm::NEIGHBORHOOD_RANK, "Rank vertices in DAG of neighborhood subset relation" },
  { TemporalAlgorithm::PROBABILITY_SORT, "Sort vertices if p_uv > threshold" },
  { TemporalAlgorithm::PROBABILITY_SUM_SORT, "Sort vertices by sum of p_uv" },
  { TemporalAlgorithm::LP_SOLUTION_SORT, "Sort vertices if x_uv > threshold" },
};

const map<string, TemporalAlgorithm> REVERSE_ALGORITHM_NAME = {
  { "sort_by_degree", TemporalAlgorithm::DEGREE_SORT },
  { "peel_by_degree", TemporalAlgorithm::DEGREE_PEEL },
  { "sort_by_neighborhood", TemporalAlgorithm::NEIGHBORHOOD_SORT },
  { "peel_by_neighborhood_sl", TemporalAlgorithm::NEIGHBORHOOD_PEEL_SL },
  { "peel_by_neighborhood_lf", TemporalAlgorithm::NEIGHBORHOOD_PEEL_LF },
  { "rank_by_neighborhood", TemporalAlgorithm::NEIGHBORHOOD_RANK },
  { "sort_by_probability", TemporalAlgorithm::PROBABILITY_SORT },
  { "sort_by_probability_sum", TemporalAlgorithm::PROBABILITY_SUM_SORT },
  { "sort_by_lp_solution", TemporalAlgorithm::LP_SOLUTION_SORT },
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
    for (size_t i = 0; i < sources.size(); i++) {
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

inline std::string get_age_name(const std::string &graph_name) {
  return std::regex_replace(graph_name, std::regex("^G"), "PH");
}

vector<int> read_age(const string &age_name) {
  ifstream age_file(age_name);
  if (age_file.fail()) {
    throw invalid_argument("Missing " + age_name + " file");
  }
  unsigned v;
  int age_v;
  vector<int> age;
  while (!age_file.eof()) {
    age_file >> v >> age_v;
    if (v >= age.size()) {
      age.resize(v + 1);
    }
    age[v] = age_v;
  }
  age_file.close();
  return age;
}

int count_age(const vector<int> &age, const int &age_value) {
  return count_if(age.begin(), age.end(), [&](const int &v){ return v == age_value; });
}

std::tuple<Graph, std::vector<int>, int> read_graph_with_age(
    const std::string &graph_name, const std::string &age_name,
    const int &min_age, const int &max_age) {
  std::ifstream graph_file(graph_name);
  if (graph_file.fail()) {
    throw std::invalid_argument("Missing " + graph_name + " file");
  }
  vector<int> all_nodes_age(read_age(age_name));
  int n0 = count_age(all_nodes_age, AGE_ZERO);

  Graph G;
  std::vector<int> age(all_nodes_age.size());
  std::map<int, Vertex> V;
  int u, v, first = 0, second = n0;
  while (graph_file >> u >> v) {
    if (all_nodes_age[u] >= min_age && all_nodes_age[u] <= max_age && !V.count(u)) {
      int &index = all_nodes_age[u] == AGE_ZERO ? first : second;
      V.insert(make_pair(u, add_vertex(G, index))), age[index] = all_nodes_age[u], ++index;
    }
    if (all_nodes_age[v] >= min_age && all_nodes_age[v] <= max_age && !V.count(v)) {
      int &index = all_nodes_age[v] == AGE_ZERO ? first : second;
      V.insert(make_pair(v, add_vertex(G, index))), age[index] = all_nodes_age[v], ++index;
    }
    if (u != v && V.count(u) && V.count(v)) {
      add_edge(G, V[u], V[v]);
    }
  }
  graph_file.close();
  age.resize(second);
  return make_tuple(G, age, n0);
}

void relabel_g0_first(Graph &G, const int &n0, const vector<int> &age) {
  int first = 0, second = n0;
  vector<Vertex> V(get_vertices(G));
  for (size_t i = 0; i < V.size(); i++) {
    if (age[get_index(G, V[i])] == AGE_ZERO) {
      set_index(G, V[i], first), first++;
    } else {
      set_index(G, V[i], second), second++;
    }
  }
}

vector<int> get_reverse_permutation(Graph &G) {
  vector<int> S(get_graph_size(G), 0);
  vector<Vertex> V(get_vertices(G));
  for (size_t i = 0; i < V.size(); i++) {
    S[get_index(G, V[i])] = i;
  }
  return S;
}

int get_total_pairs(const vector<int> &age) {
  int total = 0;
  for (const auto &u : age) {
    for (const auto &v : age) {
      if (u > AGE_ZERO && v > AGE_ZERO && u < v) {
        total++;
      }
    }
  }
  return total;
}

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

BinningScheme peel_by_neighborhood_sl(Graph &G, const int &n0) {
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

BinningScheme peel_by_neighborhood_lf(Graph &G, const int &n0) {
  BinningScheme out;
  while (get_graph_size(G) > n0) {
    set<Vertex> best_vertices;
    for (const auto &u : get_vertices(G)) {
      if (get_index(G, u) < n0) {
        continue;
      }
      auto Nu = get_neighbors(G, u);
      bool has_ancestors = false;
      for (const auto &v : get_vertices(G)) {
        if (u == v || get_index(G, v) < n0) {
          continue;
        }
        auto Nv = get_neighbors(G, v);
        if (Nu.size() < Nv.size() && includes(Nv.begin(), Nv.end(), Nu.begin(), Nu.end())) {
          has_ancestors = true;
          break;
        }
      }
      if (!has_ancestors) {
        best_vertices.insert(u);
      }
    }
    out.push_back(Bin());
    Bin &current_bin = out.back();
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
  vector<int> S = get_reverse_permutation(G);
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      if (i == j) {
        continue;
      }
      const auto &p_ij = p_uv.find(make_pair(i, j));
      if (p_ij != p_uv.end() && p_ij->second > threshold) {
        out.push_back(make_pair(S[i], S[j]));
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
  vector<int> S = get_reverse_permutation(G);
  for (int i = n0; i < n; i++) {
    long double p_i = 0;
    for (int j = n0; j < n; j++) {
      if (i == j) {
        continue;
      }
      const auto &p_ij = p_uv.find(make_pair(i, j));
      p_i += (p_ij != p_uv.end() ? p_ij->second : 0.0L);
    }
    p_v.push(make_pair(p_i, S[i]));
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

PairingScheme sort_by_lp_solution(
    Graph &G, const int &n0, const Parameters &params,
    const double &epsilon, const double &threshold) {
  const auto permutations =
      get_permutation_probabilities_sampling(G, n0, params, SamplingMethod::UNIFORM, SIGMA_TRIES);
  const auto p_uv = get_p_uv_from_permutations(permutations, get_graph_size(G), n0);
  int n = get_graph_size(G);
  map<pair<int, int>, double> x_uv;
  double solution;
  tie(solution, x_uv) = LP_solve(p_uv, n, n0, epsilon, true);

  PairingScheme out;
  vector<int> S = get_reverse_permutation(G);
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      if (i == j) {
        continue;
      }
      const auto &x_ij = x_uv.find(make_pair(i, j));
      if (x_ij != x_uv.end() && x_ij->second > threshold) {
        out.push_back(make_pair(S[i], S[j]));
      }
    }
  }
  return out;
}

DensityPrecision get_density_precision(
    const PairingScheme &solution, const vector<int> &node_age) {
  double count = get_total_pairs(node_age), total = 0, correct = 0;
  for (const auto &uv : solution) {
    if (node_age[uv.first] <= AGE_ZERO || node_age[uv.second] <= AGE_ZERO) {
      continue;
    }
    if (node_age[uv.first] < node_age[uv.second]) {
      correct++, total++;
    } else if (node_age[uv.first] > node_age[uv.second]) {
      total++;
    }
  }
  return make_pair(total / count, correct / total);
}

DensityPrecision get_density_precision(
    const BinningScheme &solution, const vector<int> &node_age) {
  double count = get_total_pairs(node_age), total = 0, correct = 0;
  for (size_t i = 0; i < solution.size(); i++) {
    for (size_t j = i + 1; j < solution.size(); j++) {
      for (const auto &u : solution[i]) {
        if (node_age[u] <= AGE_ZERO) {
          continue;
        }
        for (const auto &v : solution[j]) {
          if (node_age[v] <= AGE_ZERO) {
            continue;
          }
          if (node_age[u] < node_age[v]) {
            correct++, total++;
          } else if (node_age[u] > node_age[v]) {
            total++;
          }
        }
      }
    }
  }
  return make_pair(total / count, correct / total);
}

DensityPrecision temporal_algorithm_single(
    Graph &G, const int &n0, const vector<int> &node_age,
    const TemporalAlgorithm &algorithm, const Parameters &params) {
  switch (algorithm) {
    case DEGREE_SORT:
      return get_density_precision(sort_by_degree(G, n0), node_age);
    case DEGREE_PEEL:
      return get_density_precision(peel_by_degree(G, n0), node_age);
    case NEIGHBORHOOD_SORT:
      return get_density_precision(sort_by_neighborhood(G, n0), node_age);
    case NEIGHBORHOOD_PEEL_SL:
      return get_density_precision(peel_by_neighborhood_sl(G, n0), node_age);
    case NEIGHBORHOOD_PEEL_LF:
      return get_density_precision(peel_by_neighborhood_lf(G, n0), node_age);
    case NEIGHBORHOOD_RANK:
      return get_density_precision(rank_by_neighborhood(G, n0), node_age);
    case PROBABILITY_SORT:
      // TODO(kturowski): parametrize by different values of threshold than 0.5
      return get_density_precision(sort_by_probability(G, n0, params, 0.5), node_age);
    case PROBABILITY_SUM_SORT:
      return get_density_precision(sort_by_probability_sum(G, n0, params), node_age);
    case LP_SOLUTION_SORT:
      // TODO(kturowski): parametrize by different values of epsilon than 1.0
      // TODO(kturowski): parametrize by different values of threshold than 0.5
      return get_density_precision(sort_by_lp_solution(G, n0, params, 1.0, 0.5), node_age);
    default:
      throw invalid_argument("Invalid algorithm: " + LONG_ALGORITHM_NAME.find(algorithm)->second);
  }
}

void print(
    const string &name, const TemporalAlgorithm &algorithm,
    const vector<DensityPrecision> solution, ostream &out_file, bool verbose = false) {
  cout << name << endl;
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

void synthetic_data(
    const int &n, const int &n0, const Parameters &params, const TemporalAlgorithm &algorithm) {
  Graph G0(generate_seed(n0, 0.6));
  vector<DensityPrecision> density_precision_values(G_TRIES);
  vector<int> node_age(n, 0);
  for (int i = n0; i < n; i++) {
    node_age[i] = i;
  }
  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    Graph G(G0);
    generate_graph(G, n, params);
    // TODO(kturowski): parametrize by different values of params than used to generate G
    density_precision_values[i] = temporal_algorithm_single(G, n0, node_age, algorithm, params);
    #pragma omp critical
    {
      if ((i + 1) % 1000 == 0) {
        cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
      }
    }
  }
  ofstream out_file(TEMP_FOLDER + get_synthetic_filename(n, n0, params, "TA"), ios_base::app);
  print("Synthetic data: " + params.to_string(), algorithm, density_precision_values, out_file);
}

void real_world_data(
    const string &graph_name, const string &age_name,
    const TemporalAlgorithm &algorithm, const Parameters &params) {
  Graph G;
  vector<int> age;
  int n0;
  tie(G, age, n0) =
      read_graph_with_age(FILES_FOLDER + graph_name, FILES_FOLDER + age_name, AGE_ZERO, AGE_MAX);
  auto density_precision_value = temporal_algorithm_single(G, n0, age, algorithm, params);
  ofstream out_file(TEMP_FOLDER + get_real_filename(graph_name, params.mode, "TA"), ios_base::app);
  print(graph_name, algorithm, vector<DensityPrecision>{density_precision_value}, out_file);
}

int main(int, char *argv[]) {
  try {
    string action(argv[1]), algorithm_name(argv[2]);
    if (action == "synthetic") {
      string mode_0(argv[3]);
      int n = stoi(argv[4]), n0 = stoi(argv[5]);
      Parameters params_0;
      params_0.initialize(mode_0, argv + 6);
      if (algorithm_name == "all") {
        for (const auto &algorithm : REVERSE_ALGORITHM_NAME) {
          synthetic_data(n, n0, params_0, algorithm.second);
        }
      } else if (REVERSE_ALGORITHM_NAME.count(algorithm_name)) {
        synthetic_data(n, n0, params_0, REVERSE_ALGORITHM_NAME.find(algorithm_name)->second);
      } else {
        throw invalid_argument("Invalid algorithm: " + algorithm_name);
      }
    } else if (action == "real_data") {
      string graph_name(argv[3]), mode(argv[4]);
      Parameters params;
      params.initialize(mode, argv + 5);
      if (algorithm_name == "all") {
        for (const auto &algorithm : REVERSE_ALGORITHM_NAME) {
          real_world_data(
              graph_name, get_age_name(graph_name), algorithm.second, params);
        }
      } else if (REVERSE_ALGORITHM_NAME.count(algorithm_name)) {
        real_world_data(
            graph_name, get_age_name(graph_name),
            REVERSE_ALGORITHM_NAME.find(algorithm_name)->second, params);
      } else {
        throw invalid_argument("Invalid algorithm: " + algorithm_name);
      }
    } else {
      throw invalid_argument("Invalid action: " + action);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }

  return 0;
}
