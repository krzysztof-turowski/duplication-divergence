// Tools for computation the temporal order bound for various duplication-divergence models.
// Compile: g++ dd_temporal_bound.cpp -O3 -lgmpxx -lgmp -lglpk -o ./dd_temporal_bound
// Run: ./dd_temporal_bound [exact|ALGORITHM] MODE n n0 PARAMETERS

// TODO(kturowski): deal gurobi output suppression and output to cout instead of cerr

#include "./dd_input.h"
#include "./dd_temporal.h"

#if defined(glpk)
  #include "./dd_glpk.h"
#elif defined(gurobi)
  #include "./dd_gurobi.h"
#endif

#include <queue>

using namespace std;

int G_TRIES, SIGMA_TRIES;
const double EPS_MIN = 0.05, EPS_STEP = 0.05;
const long double RANDOM_WALK_THRESHOLD = 1.0;

void print_density_precision(
    const string &name, const vector<double> &density, const vector<double> &precision,
    const int &n, const int &n0, const Parameters &params, ostream &out_file) {
  cerr << "Graph - n: " << n << ", n0: " << n0
      << ", parameters: " << params.to_string() << endl;
  cerr << "Method: " << name << endl;
  for (size_t i = 0; i < precision.size(); i++) {
    cerr << "density: " << fixed << setw(6) << setprecision(3) << density[i]
        << " precision: " << fixed << setw(6) << setprecision(3) << precision[i] << endl;
  }
  out_file << name << " ";
  for (size_t i = 0; i < precision.size(); i++) {
    out_file << density[i] << "," << precision[i] << " ";
  }
  out_file << endl;
}

vector<double> LP_bound_exact_single(
    const Graph &G0, const int &n, const Parameters &params, const vector<double> &epsilon,
    const set<VertexPair> &perfect_pairs) {
  Graph G(G0);
  generate_graph(G, n, params);

  vector<int> S = generate_permutation(n, get_graph_size(G0));
  apply_permutation(G, S);

  auto permutations = get_log_permutation_probabilities(G, get_graph_size(G0), params);
  normalize_log_probabilities(permutations);
  print_best_permutations(
      permutations, G, params, get_graph_size(G0), "exact", PERMUTATION_COUNT_LIMIT);
  auto p_uv = get_p_uv_from_permutations(permutations, n, get_graph_size(G0));
  set_perfect_pairs(p_uv, perfect_pairs);
  vector<double> solutions;
  for (const double &eps : epsilon) {
    double solution;
    tie(solution, ignore) = LP_solve(p_uv, n, get_graph_size(G0), eps);
    solutions.push_back(solution);
  }
  return solutions;
}

vector<double> LP_bound_random_walk_single(
    const Graph &G0, const int &n, const Parameters &params, const vector<double> &epsilon,
    const set<VertexPair> &perfect_pairs) {
  const int SIGMA_RANDOM_WALK_LENGTH = 10000, SIGMA_BEST = 100;
  Graph G(G0);
  generate_graph(G, n, params);

  int n0 = get_graph_size(G0);

  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> pick_distribution(0.0, 1.0);
  std::map<mpz_class, long double> permutations;
  auto sample_permutations = get_log_permutation_probabilities_sampling(
      G, get_graph_size(G0), params, get_DAG_from_perfect_pairs(perfect_pairs, n),
      SamplingMethod::UNIFORM, SIGMA_TRIES);
  std::priority_queue<std::pair<long double, mpz_class>> Q;
  for (const auto &permutation : sample_permutations) {
    Q.push(std::make_pair(permutation.second, permutation.first));
  }
  for (int i = 0; i < SIGMA_BEST; i++) {
    vector<int> S = decode_permutation(Q.top().second, n);
    long double score = get_log_permutation_probability(G, n0, params, reverse_permutation(S));
    mpz_class best_sigma = Q.top().second;
    long double best_score = score;
    for (int j = 0; j < SIGMA_RANDOM_WALK_LENGTH; j++) {
      std::uniform_int_distribution<int> swap_distribution(n0, n - 1);
      int u = swap_distribution(generator), v = swap_distribution(generator);
      if (u == v) {
        continue;
      }
      swap(S[u], S[v]);
      score = get_log_permutation_probability(G, n0, params, reverse_permutation(S));
      if (best_score < score) {
        best_score = score, best_sigma = encode_permutation(S);
      } else {
        swap(S[u], S[v]);
      }
    }
    Q.pop();
  }
  normalize_log_probabilities(permutations);
  print_best_permutations(permutations, G, params, n0, "random_walk", PERMUTATION_COUNT_LIMIT);
  auto p_uv = get_p_uv_from_permutations(permutations, n, n0);
  set_perfect_pairs(p_uv, perfect_pairs);
  vector<double> solutions;
  for (const double &eps : epsilon) {
    double solution;
    tie(solution, ignore) = LP_solve(p_uv, n, get_graph_size(G0), eps);
    solutions.push_back(solution);
  }
  return solutions;
}

vector<double> LP_bound_approximate_single(
    const Graph &G0, const int &n, const Parameters &params, const SamplingMethod &algorithm,
    const vector<double> &epsilon, const set<VertexPair> &perfect_pairs) {
  Graph G(G0);
  generate_graph(G, n, params);

  vector<int> S = generate_permutation(n, get_graph_size(G0));
  apply_permutation(G, S);

  auto permutations =
      get_log_permutation_probabilities_sampling(
          G, get_graph_size(G0), params, get_DAG_from_perfect_pairs(perfect_pairs, n),
          algorithm, SIGMA_TRIES);
  normalize_log_probabilities(permutations);
  print_best_permutations(
      permutations, G, params, get_graph_size(G0),
      SAMPLING_METHOD_NAME.find(algorithm)->second, PERMUTATION_COUNT_LIMIT);
  auto p_uv = get_p_uv_from_permutations(permutations, n, get_graph_size(G0));
  set_perfect_pairs(p_uv, perfect_pairs);
  vector<double> solutions;
  for (const double &eps : epsilon) {
    double solution;
    tie(solution, ignore) = LP_solve(p_uv, n, get_graph_size(G0), eps);
    solutions.push_back(solution);
  }
  return solutions;
}

void LP_bound_exact(
    const int &n, const int &n0, const Parameters &params, const double &p0, ostream &out_file) {
  Graph G0 = generate_seed(n0, p0);
  vector<double> epsilon;
  for (double eps = EPS_MIN; eps <= 1.0 + 10e-9; eps += EPS_STEP) {
    epsilon.push_back(eps);
  }
  vector<double> solution(epsilon.size(), 0.0);
  vector<vector<double>> solutions(G_TRIES);
  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    solutions[i] = LP_bound_exact_single(G0, n, params, epsilon, set<VertexPair>());
    #pragma omp critical
    {
      cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
    }
  }
  for (int i = 0; i < G_TRIES; i++) {
    transform(
        solution.begin(), solution.end(), solutions[i].begin(), solution.begin(),
        std::plus<double>());
  }
  for (auto &s : solution) {
    s /= G_TRIES;
  }
  print_density_precision("exact", epsilon, solution, n, n0, params, out_file);
}

void LP_bound_random_walk(
    const int &n, const int &n0, const Parameters &params, const double &p0, ostream &out_file) {
  Graph G0 = generate_seed(n0, p0);
  vector<double> epsilon;
  for (double eps = EPS_MIN; eps <= 1.0 + 10e-9; eps += EPS_STEP) {
    epsilon.push_back(eps);
  }
  vector<double> solution(epsilon.size(), 0.0);
  vector<vector<double>> solutions(G_TRIES);
  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    solutions[i] = LP_bound_random_walk_single(G0, n, params, epsilon, set<VertexPair>());
    #pragma omp critical
    {
      cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
    }
  }
  for (int i = 0; i < G_TRIES; i++) {
    transform(
        solution.begin(), solution.end(), solutions[i].begin(), solution.begin(),
        std::plus<double>());
  }
  for (auto &sol : solution) {
    sol /= G_TRIES;
  }
  print_density_precision(
      "random_walk-" + to_string(G_TRIES) + "-" + to_string(SIGMA_TRIES),
      epsilon, solution, n, n0, params, out_file);
}

void LP_bound_approximate(
    const int &n, const int &n0, const Parameters &params, const double &p0,
    const SamplingMethod &algorithm, ostream &out_file) {
  Graph G0 = generate_seed(n0, p0);
  vector<double> epsilon;
  for (double eps = EPS_MIN; eps <= 1.0 + 10e-9; eps += EPS_STEP) {
    epsilon.push_back(eps);
  }
  vector<double> solution(epsilon.size(), 0.0);
  vector<vector<double>> solutions(G_TRIES);
  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    solutions[i] =
        LP_bound_approximate_single(G0, n, params, algorithm, epsilon, set<VertexPair>());
    #pragma omp critical
    {
      cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
    }
  }
  for (int i = 0; i < G_TRIES; i++) {
    transform(
        solution.begin(), solution.end(), solutions[i].begin(), solution.begin(),
        std::plus<double>());
  }
  for (auto &sol : solution) {
    sol /= G_TRIES;
  }
  print_density_precision(
      SAMPLING_METHOD_NAME.find(algorithm)->second + "-" + to_string(G_TRIES)
          + "-" + to_string(SIGMA_TRIES),
      epsilon, solution, n, n0, params, out_file);
}

int main(int argc, char **argv) {
  try {
    TEnv environment = prepare_environment(argc, argv);
    G_TRIES = read_int(environment, "-gt:", 100, "G_TRIES");
    SIGMA_TRIES = read_int(environment, "-st:", 100000, "SIGMA_TRIES");

    string algorithm = read_string(
        environment, "-algorithm:", "all", "Sampling algorithm to run");
    string mode = read_mode(environment);
    const int n = read_n(environment), n0 = read_n0(environment);
    const double p0 = read_p0(environment);
    Parameters params = read_parameters(environment);
    string name(TEMP_FOLDER + get_synthetic_filename(n, n0, params, "TC"));
    if (algorithm == "exact") {
      if (!validate_problem_size(n, n0)) {
        throw out_of_range(
            "Graph too large for exact mode: n = " + to_string(n) + ", n0 = " + to_string(n0));
      }
      ofstream out_file(name, ios_base::app);
      LP_bound_exact(n, n0, params, p0, out_file);
    } else if (algorithm == "random_walk") {
      ofstream out_file(name, ios_base::app);
      LP_bound_random_walk(n, n0, params, p0, out_file);
    } else if (SAMPLING_METHOD_REVERSE_NAME.count(algorithm)) {
      ofstream out_file(name, ios_base::app);
      LP_bound_approximate(
          n, n0, params, p0, SAMPLING_METHOD_REVERSE_NAME.find(algorithm)->second, out_file);
    } else {
      throw invalid_argument("Invalid algorithm: " + algorithm);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }

  return 0;
}
