// Algorithms computating the temporal order for various duplication-divergence models.
// Compile: g++ dd_temporal_algorithms.cpp -O3 -o ./dd_temporal_algorithms
// Run: ./dd_temporal_algorithms synthetic all MODE_GEN n n0 PARAMETERS_GEN MODE PARAMETERS
//   or ./dd_temporal_algorithms synthetic ALGORITHM MODE_GEN n n0 PARAMETERS_GEN [MODE PARAMETERS]
//   or ./dd_temporal_algorithms real_data all FILE MODE PARAMETERS
//   or ./dd_temporal_algorithms real_data ALGORITHM FILE [MODE PARAMETERS]

#include "./dd_input.h"
#include "./dd_temporal.h"
#include "./dd_perfect_pairs.h"

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
typedef vector<VertexPair> PairingScheme;

int G_TRIES, SIGMA_TRIES;

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

class PartialOrderScore {
 public:
  double density = 0.0, precision = 0.0, correlation = 0.0, gamma = 0.0;
};

class Settings {
 public:
  double threshold, epsilon, perfect_pairs_fraction;
  int bin_size;

  explicit Settings(TEnv &environment) {
    perfect_pairs_fraction =
        read_double(environment, "-perfect:", nan(""), "Fraction of perfect pairs");
    threshold = read_double(environment, "-threshold:", nan(""), "Threshold for uv");
    epsilon = read_double(environment, "-epsilon:", nan(""), "Solution density");
    bin_size = read_int(environment, "-binsize:", 1, "Maximum bin size");
  }
};

PartialOrderScore get_average(const vector<PartialOrderScore> &scores) {
  PartialOrderScore average;
  for (const auto &score : scores) {
    average.density += score.density, average.precision += score.precision;
    average.correlation += score.correlation, average.gamma += score.gamma;
  }
  average.density /= scores.size(), average.precision /= scores.size();
  average.correlation /= scores.size(), average.gamma /= scores.size();
  return average;
}

inline std::string get_age_name(const std::string &graph_name) {
  return std::regex_replace(graph_name, std::regex("^G"), "age");
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
  DAG H(get_graph_size(G));
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
    Graph &G, const int &n0, const Parameters &params, const double &p_uv_threshold,
    const std::set<VertexPair> &perfect_pairs) {
  int n = get_graph_size(G);
  const auto &DAG = get_DAG_from_perfect_pairs(perfect_pairs, n);
  auto permutations =
      get_log_permutation_probabilities_sampling(
          G, n0, params, DAG, SamplingMethod::UNIFORM, SIGMA_TRIES);
  normalize_log_probabilities(permutations);
  const auto &p_uv = get_p_uv_from_permutations(permutations, n, n0);

  PairingScheme out;
  vector<int> S = get_reverse_permutation(G);
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      if (i == j) {
        continue;
      }
      const auto &p_ij = p_uv.find(make_pair(i, j));
      if (p_ij != p_uv.end() && p_ij->second > p_uv_threshold) {
        out.push_back(make_pair(S[i], S[j]));
      }
    }
  }
  return out;
}

BinningScheme sort_by_probability_sum(
    Graph &G, const int &n0, const Parameters &params, const size_t &bin_size,
    const std::set<VertexPair> &perfect_pairs) {
  int n = get_graph_size(G);
  const auto &DAG = get_DAG_from_perfect_pairs(perfect_pairs, n);
  auto permutations =
      get_log_permutation_probabilities_sampling(
          G, n0, params, DAG, SamplingMethod::UNIFORM, SIGMA_TRIES);
  normalize_log_probabilities(permutations);
  const auto &p_uv = get_p_uv_from_permutations(permutations, n, n0);

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
  Bin current_bin;
  while (!p_v.empty()) {
    const auto &p_i = p_v.top();
    current_bin.insert(p_i.second);
    if (current_bin.size() == bin_size) {
      out.push_back(current_bin), current_bin.clear();
    }
    p_v.pop();
  }
  out.push_back(current_bin);
  return out;
}

PairingScheme sort_by_lp_solution(
    Graph &G, const int &n0, const Parameters &params, const double &epsilon,
    const double &x_uv_threshold, const std::set<VertexPair> &perfect_pairs) {
  int n = get_graph_size(G);
  const auto &DAG = get_DAG_from_perfect_pairs(perfect_pairs, n);
  auto permutations =
      get_log_permutation_probabilities_sampling(
          G, n0, params, DAG, SamplingMethod::UNIFORM, SIGMA_TRIES);
  normalize_log_probabilities(permutations);
  const auto p_uv = get_p_uv_from_permutations(permutations, n, n0);

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
      if (x_ij != x_uv.end() && x_ij->second > x_uv_threshold) {
        out.push_back(make_pair(S[i], S[j]));
      }
    }
  }
  return out;
}

void get_density_precision(
    const PairingScheme &solution, const vector<int> &node_age,
    const set<VertexPair> &perfect_pairs, PartialOrderScore &score) {
  double count = get_total_pairs(node_age) - perfect_pairs.size(), total = 0, correct = 0;
  for (const auto &uv : solution) {
    if (perfect_pairs.count(uv) || perfect_pairs.count(make_pair(uv.second, uv.first))) {
      continue;
    }
    if (node_age[uv.first] <= AGE_ZERO || node_age[uv.second] <= AGE_ZERO) {
      continue;
    }
    if (node_age[uv.first] < node_age[uv.second]) {
      correct++, total++;
    } else if (node_age[uv.first] > node_age[uv.second]) {
      total++;
    }
  }
  PartialOrderScore out;
  score.precision = correct / total, score.density = total / count;
}

void get_density_precision(
    const BinningScheme &solution, const vector<int> &node_age,
    const set<VertexPair> &perfect_pairs, PartialOrderScore &score) {
  double count = get_total_pairs(node_age) - perfect_pairs.size(), total = 0, correct = 0;
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
          if (perfect_pairs.count(make_pair(u, v)) || perfect_pairs.count(make_pair(v, u))) {
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
  score.precision = correct / total, score.density = total / count;
}

double get_goodman_gamma(
    const BinningScheme &solution, const vector<int> &node_age,
    const set<VertexPair> &perfect_pairs) {
  double data_size = 0, concordant_pairs = 0, discordant_pairs = 0, solution_ties = 0;
  map<int, int> values;
  for (size_t i = 0; i < solution.size(); i++) {
    int nodes = 0;
    for (const auto &u : solution[i]) {
      if (u > AGE_ZERO) {
        values[u]++, nodes++;
      }
    }
    data_size += nodes, solution_ties += nodes * (nodes - 1) / 2;
    for (size_t j = i + 1; j < solution.size(); j++) {
      for (const auto &u : solution[i]) {
        if (node_age[u] <= AGE_ZERO) {
          continue;
        }
        for (const auto &v : solution[j]) {
          if (node_age[v] <= AGE_ZERO) {
            continue;
          }
          if (perfect_pairs.count(make_pair(u, v)) || perfect_pairs.count(make_pair(v, u))) {
            continue;
          }
          if (node_age[u] < node_age[v]) {
            concordant_pairs++;
          } else if (node_age[u] > node_age[v]) {
            discordant_pairs++;
          }
        }
      }
    }
  }
  double original_ties = 0, all_pairs = data_size * (data_size - 1) / 2;
  for (const auto &u : values) {
    original_ties += u.second * (u.second - 1) / 2;
  }
  return (concordant_pairs - discordant_pairs)
      / (sqrt(all_pairs - solution_ties) * sqrt(all_pairs - original_ties));
}

double get_pearson_correlation(const BinningScheme &solution, const vector<int> &node_age) {
  double total_node_age = accumulate(node_age.begin(), node_age.end(), 0.0);

  vector<int> estimated_age(node_age.size());
  double total_estimated_age = 0, data_size = 0;
  for (size_t i = 0; i < solution.size(); i++) {
    for (const auto &v : solution[i]) {
      if (node_age[v] <= AGE_ZERO) {
        continue;
      }
      estimated_age[v] = i, total_estimated_age += i, data_size++;
    }
  }
  total_node_age /= data_size, total_estimated_age /= data_size;
  double rho_num = 0, rho_den = 0, rho_den_estimated = 0;
  for (size_t i = 0; i < node_age.size(); i++) {
    if (node_age[i] <= AGE_ZERO) {
      continue;
    }
    double diff = node_age[i] - total_node_age;
    double estimated_diff = estimated_age[i] - total_estimated_age;
    rho_num += diff * estimated_diff;
    rho_den += diff * diff, rho_den_estimated += estimated_diff * estimated_diff;
  }
  return rho_num / sqrt(rho_den * rho_den_estimated);
}

PartialOrderScore get_score(
    const PairingScheme &solution, const vector<int> &node_age,
    const set<VertexPair> &perfect_pairs = set<VertexPair>()) {
  PartialOrderScore score;
  get_density_precision(solution, node_age, perfect_pairs, score);
  score.correlation = nan(""), score.gamma = nan("");
  return score;
}

PartialOrderScore get_score(
    const BinningScheme &solution, const vector<int> &node_age,
    const set<VertexPair> &perfect_pairs = set<VertexPair>()) {
  PartialOrderScore score;
  score.correlation = get_pearson_correlation(solution, node_age);
  score.gamma = get_goodman_gamma(solution, node_age, perfect_pairs);
  get_density_precision(solution, node_age, perfect_pairs, score);
  return score;
}

PartialOrderScore temporal_algorithm_single(
    Graph &G, const int &n0, const vector<int> &node_age,
    const TemporalAlgorithm &algorithm, const Parameters &params,
    const Settings &settings) {
  switch (algorithm) {
    case DEGREE_SORT:
      return get_score(sort_by_degree(G, n0), node_age);
    case DEGREE_PEEL:
      return get_score(peel_by_degree(G, n0), node_age);
    case NEIGHBORHOOD_SORT:
      return get_score(sort_by_neighborhood(G, n0), node_age);
    case NEIGHBORHOOD_PEEL_SL:
      return get_score(peel_by_neighborhood_sl(G, n0), node_age);
    case NEIGHBORHOOD_PEEL_LF:
      return get_score(peel_by_neighborhood_lf(G, n0), node_age);
    case NEIGHBORHOOD_RANK:
      return get_score(rank_by_neighborhood(G, n0), node_age);
    case PROBABILITY_SORT: {
      auto perfect_pairs = get_perfect_pairs(node_age, settings.perfect_pairs_fraction);
      return get_score(
          sort_by_probability(G, n0, params, settings.threshold, perfect_pairs),
          node_age, perfect_pairs);
    }
    case PROBABILITY_SUM_SORT: {
      auto perfect_pairs = get_perfect_pairs(node_age, settings.perfect_pairs_fraction);
      return get_score(
          sort_by_probability_sum(G, n0, params, settings.bin_size, perfect_pairs),
          node_age, perfect_pairs);
    }
    case LP_SOLUTION_SORT: {
      auto perfect_pairs = get_perfect_pairs(node_age, settings.perfect_pairs_fraction);
      return get_score(
          sort_by_lp_solution(G, n0, params, settings.epsilon, settings.threshold, perfect_pairs),
          node_age, perfect_pairs);
    }
    default:
      throw invalid_argument("Invalid algorithm: " + LONG_ALGORITHM_NAME.find(algorithm)->second);
  }
}

void print(
    const string &name, const TemporalAlgorithm &algorithm, const Settings &settings,
    const vector<PartialOrderScore> &scores, ostream &out_file, bool verbose = false) {
  cout << name << endl;
  cout << "Algorithm: " << LONG_ALGORITHM_NAME.find(algorithm)->second << endl;
  cout << "Perfect pairs fraction: " << settings.perfect_pairs_fraction << endl;
  if (!isnan(settings.threshold)) {
    cout << "Threshold: " << fixed << setw(6) << setprecision(3) << settings.threshold << endl;
  }
  if (!isnan(settings.epsilon)) {
    cout << "Epsilon: " << fixed << setw(6) << setprecision(3) << settings.epsilon << endl;
  }

  auto mean_score = get_average(scores);
  cout << "Mean density: " << fixed << setw(6) << setprecision(3)
      << mean_score.density
      << " mean precision: " << fixed << setw(6) << setprecision(3)
      << mean_score.precision
      << endl;
  if (!isnan(mean_score.correlation)) {
    cout << "Mean Pearson correlation: " << fixed << setw(6) << setprecision(3)
        << mean_score.gamma << endl;
  }
  if (!isnan(mean_score.gamma)) {
    cout << "Mean Goodman gamma: " << fixed << setw(6) << setprecision(3)
        << mean_score.gamma << endl;
  }
  if (verbose) {
    for (const auto &score : scores) {
      cout << "density: " << fixed << setw(6) << setprecision(3) << score.density
          << " precision: " << fixed << setw(6) << setprecision(3) << score.precision
          << endl;
    }
  }

  out_file << SHORT_ALGORITHM_NAME.find(algorithm)->second;
  if (!isnan(settings.perfect_pairs_fraction)) {
    out_file << "-pp:" << fixed << setprecision(3) << settings.perfect_pairs_fraction;
  }
  if (settings.bin_size > 1) {
    out_file << "-b:" << fixed << setprecision(3) << settings.bin_size;
  }
  if (!isnan(settings.threshold)) {
    out_file << "-th:" << fixed << setprecision(3) << settings.threshold;
  }
  if (!isnan(settings.epsilon)) {
    out_file << "-eps:" << fixed << setprecision(3) << settings.epsilon;
  }
  for (const auto &score : scores) {
    out_file << " " << score.density << "," << score.precision;
  }
  out_file << endl;
}

void synthetic_data(
    const int &n, const int &n0, const Parameters &params, const double &p0,
    const TemporalAlgorithm &algorithm, const Settings &settings) {
  Graph G0(generate_seed(n0, p0));
  vector<PartialOrderScore> scores(G_TRIES);
  vector<int> node_age(n, 0);
  for (int i = n0; i < n; i++) {
    node_age[i] = i;
  }
  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    Graph G(G0);
    generate_graph(G, n, params);
    // TODO(kturowski): parametrize by different values of params than used to generate G
    scores[i] = temporal_algorithm_single(G, n0, node_age, algorithm, params, settings);
    #pragma omp critical
    {
      if ((i + 1) % 1000 == 0) {
        cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
      }
    }
  }
  ofstream out_file(TEMP_FOLDER + get_synthetic_filename(n, n0, params, "TA"), ios_base::app);
  print("Synthetic data: " + params.to_string(), algorithm, settings, scores, out_file);
}

void real_world_data(
    const string &graph_name, const string &age_name, const TemporalAlgorithm &algorithm,
    const Parameters &params, const Settings &settings) {
  Graph G;
  vector<int> age;
  int n0;
  tie(G, age, n0) =
      read_graph_with_age(FILES_FOLDER + graph_name, FILES_FOLDER + age_name, AGE_ZERO, AGE_MAX);
  auto score = temporal_algorithm_single(G, n0, age, algorithm, params, settings);
  ofstream out_file(TEMP_FOLDER + get_real_filename(graph_name, params.mode, "TA"), ios_base::app);
  print(graph_name, algorithm, settings, vector<PartialOrderScore>{score}, out_file);
}

int main(int argc, char **argv) {
  try {
    Env = prepare_environment(argc, argv);
    G_TRIES = read_int(Env, "-gt:", 1, "G_TRIES");
    SIGMA_TRIES = read_int(Env, "-st:", 1, "SIGMA_TRIES");
    string action = read_action(Env);
    string algorithm_name = read_string(
        Env, "-algorithm:", "all", "Temporal algorithm to run");
    Settings settings(Env);
    if (action == "synthetic") {
      const int n = read_n(Env), n0 = read_n0(Env);
      const double p0 = read_p0(Env);
      Parameters params_0 = read_parameters(Env);

      if (algorithm_name == "all") {
        for (const auto &algorithm : REVERSE_ALGORITHM_NAME) {
          synthetic_data(n, n0, params_0, p0, algorithm.second, settings);
        }
      } else if (REVERSE_ALGORITHM_NAME.count(algorithm_name)) {
        synthetic_data(
            n, n0, params_0, p0, REVERSE_ALGORITHM_NAME.find(algorithm_name)->second, settings);
      } else {
        throw invalid_argument("Invalid algorithm: " + algorithm_name);
      }
    } else if (action == "real_data") {
      string graph_name = read_graph_name(Env);
      Parameters params = read_parameters(Env);
      if (algorithm_name == "all") {
        for (const auto &algorithm : REVERSE_ALGORITHM_NAME) {
          real_world_data(
              graph_name, get_age_name(graph_name), algorithm.second, params, settings);
        }
      } else if (REVERSE_ALGORITHM_NAME.count(algorithm_name)) {
        real_world_data(
            graph_name, get_age_name(graph_name),
            REVERSE_ALGORITHM_NAME.find(algorithm_name)->second, params, settings);
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
