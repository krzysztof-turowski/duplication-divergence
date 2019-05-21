/*
Tool for computing the upper bound for temporal order inference for various duplication-divergence models.

Compile: make dd_temporal_bound

Syntax: ./dd_temporal_bound <options>
<options> are
-algorithm:
  exact: generate all admissible permutations for a given graph under a given model
  local-unif-sampling: sample admissible permutations according to the local uniform sampling rule
  high-prob-sampling: sample admissible permutations according to the high probability sampling rule
-st: Number of independent admissible permutations sampled by the algorithm `local-unif-sampling` or `high-prob-sampling`.
-gt: Number of independently generated graphs to test.
-mode: {pure_duplication, pastor_satorras, chung_lu}. In case of `synthetic` action, the mode (type) of the duplication-divergence graph model.
<parameters>: Depending on `mode`, the parameters `p,q,r` of the duplication-divergence graph model.
-n: The size of a graph.
-n0, -p0: The parameters for generating a seed graph.

Example run:
  ./dd_temporal_bound -algorithm:high-prob-sampling -n:100 -n0:10 -mode:pastor_satorras -p:0.5 -r:2.0 -p0:0.6 -gt:100 -st:1000
*/ 

#include "./dd_input.h"
#include "./dd_perfect_pairs.h"
#include "./dd_temporal.h"

#if defined(glpk)
  #include "./dd_glpk.h"
#elif defined(gurobi)
  #include "./dd_gurobi.h"
#endif

#include <queue>

int G_TRIES, SIGMA_TRIES;
const double EPS_MIN = 0.05, EPS_STEP = 0.05;

void print_density_precision(
    const std::string &name, const std::vector<double> &density,
    const std::vector<double> &precision, const int &n, const int &n0, const Parameters &params,
    std::ostream &out_file) {
  std::cerr << "Graph - n: " << n << ", n0: " << n0
      << ", parameters: " << params.to_string() << std::endl;
  std::cerr << "Method: " << name << std::endl;
  for (size_t i = 0; i < precision.size(); i++) {
    std::cerr << "density: " << std::fixed << std::setw(6) << std::setprecision(3) << density[i]
        << " precision: " << std::fixed << std::setw(6) << std::setprecision(3) << precision[i]
        << std::endl;
  }
  out_file << name << " ";
  for (size_t i = 0; i < precision.size(); i++) {
    out_file << density[i] << "," << precision[i] << " ";
  }
  out_file << std::endl;
}

std::vector<double> LP_bound_exact_single(
    const Graph &G0, const int &n, const Parameters &params, const std::vector<double> &epsilon,
    const std::set<VertexPair> &perfect_pairs) {
  Graph G(G0);
  generate_graph(G, n, params);

  std::vector<int> S = generate_permutation(n, get_graph_size(G0));
  apply_permutation(G, S);

  auto permutations = get_log_permutation_probabilities(G, get_graph_size(G0), params);
  normalize_log_probabilities(permutations);
  auto p_uv = get_p_uv_from_permutations(permutations, n, get_graph_size(G0));
  set_perfect_pairs(p_uv, perfect_pairs);
  std::vector<double> solutions;
  for (const double &eps : epsilon) {
    double solution;
    std::tie(solution, std::ignore) = LP_solve(p_uv, n, get_graph_size(G0), eps);
    solutions.push_back(solution);
  }
  return solutions;
}

std::vector<double> LP_bound_approximate_single(
    const Graph &G0, const int &n, const Parameters &params, const SamplingMethod &algorithm,
    const std::vector<double> &epsilon, const std::set<VertexPair> &perfect_pairs) {
  Graph G(G0);
  generate_graph(G, n, params);

  std::vector<int> S = generate_permutation(n, get_graph_size(G0));
  apply_permutation(G, S);

  auto permutations =
      get_log_permutation_probabilities_sampling(
          G, get_graph_size(G0), params, get_DAG_from_perfect_pairs(perfect_pairs, n),
          algorithm, SIGMA_TRIES);
  normalize_log_probabilities(permutations);
  auto p_uv = get_p_uv_from_permutations(permutations, n, get_graph_size(G0));
  set_perfect_pairs(p_uv, perfect_pairs);
  std::vector<double> solutions;
  for (const double &eps : epsilon) {
    double solution;
    std::tie(solution, std::ignore) = LP_solve(p_uv, n, get_graph_size(G0), eps);
    solutions.push_back(solution);
  }
  return solutions;
}

void LP_bound_exact(
    const int &n, const int &n0, const Parameters &params, const double &p0,
    std::ostream &out_file) {
  Graph G0 = generate_seed(n0, p0);
  std::vector<double> epsilon;
  for (double eps = EPS_MIN; eps <= 1.0 + 10e-9; eps += EPS_STEP) {
    epsilon.push_back(eps);
  }
  std::vector<double> solution(epsilon.size(), 0.0);
  std::vector<std::vector<double>> solutions(G_TRIES);
  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    solutions[i] = LP_bound_exact_single(G0, n, params, epsilon, std::set<VertexPair>());
    #pragma omp critical
    {
      std::cerr << "Finished run " << i + 1 << "/" << G_TRIES << std::endl;
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

void LP_bound_approximate(
    const int &n, const int &n0, const Parameters &params, const double &p0,
    const SamplingMethod &algorithm, std::ostream &out_file) {
  Graph G0 = generate_seed(n0, p0);
  std::vector<double> epsilon;
  for (double eps = EPS_MIN; eps <= 1.0 + 10e-9; eps += EPS_STEP) {
    epsilon.push_back(eps);
  }
  std::vector<double> solution(epsilon.size(), 0.0);
  std::vector<std::vector<double>> solutions(G_TRIES);
  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    solutions[i] =
        LP_bound_approximate_single(G0, n, params, algorithm, epsilon, std::set<VertexPair>());
    #pragma omp critical
    {
      std::cerr << "Finished run " << i + 1 << "/" << G_TRIES << std::endl;
    }
  }
  for (int i = 0; i < G_TRIES; i++) {
    std::transform(
        solution.begin(), solution.end(), solutions[i].begin(), solution.begin(),
        std::plus<double>());
  }
  for (auto &sol : solution) {
    sol /= G_TRIES;
  }
  print_density_precision(
      SAMPLING_METHOD_NAME.find(algorithm)->second + "-" + std::to_string(G_TRIES)
          + "-" + std::to_string(SIGMA_TRIES),
      epsilon, solution, n, n0, params, out_file);
}

int main(int argc, char **argv) {
  try {
    Env = prepare_environment(argc, argv);
    G_TRIES = read_int(Env, "-gt:", 1, "G_TRIES");
    SIGMA_TRIES = read_int(Env, "-st:", 1, "SIGMA_TRIES");

    std::string algorithm = read_string(Env, "-algorithm:", "all", "Sampling algorithm to run");
    const int n = read_n(Env), n0 = read_n0(Env);
    const double p0 = read_p0(Env);
    Parameters params = read_parameters(Env);
    std::string name(TEMP_FOLDER + get_synthetic_filename(n, n0, params, "TC"));
    if (algorithm == "exact") {
      if (!validate_problem_size(n, n0)) {
        throw std::out_of_range(
            "Graph too large for exact mode: n = " + std::to_string(n)
                + ", n0 = " + std::to_string(n0));
      }
      std::ofstream out_file(name, std::ios_base::app);
      LP_bound_exact(n, n0, params, p0, out_file);
    } else if (SAMPLING_METHOD_REVERSE_NAME.count(algorithm)) {
      std::ofstream out_file(name, std::ios_base::app);
      LP_bound_approximate(
          n, n0, params, p0, SAMPLING_METHOD_REVERSE_NAME.find(algorithm)->second, out_file);
    } else {
      throw std::invalid_argument("Invalid algorithm: " + algorithm);
    }
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
  }

  return 0;
}
