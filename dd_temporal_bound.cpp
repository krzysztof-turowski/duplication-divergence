// Tools for computation the temporal order bound for various duplication-divergence models.
// Compile: g++ dd_temporal_bound.cpp -O3 -lgmpxx -lgmp -lglpk -o ./dd_temporal_bound
// Run: ./dd_temporal_bound [exact|ALGORITHM] MODE n n0 PARAMETERS

// TODO(kturowski): deal gurobi output suppression and output to cout instead of cerr

#include "./dd_temporal.h"

#if defined(glpk)
  #include "./dd_glpk.h"
#elif defined(gurobi)
  #include "./dd_gurobi.h"
#endif

using namespace std;

const int G_TRIES = 20, SIGMA_TRIES = 100000;
const double EPS_MIN = 0.05, EPS_STEP = 0.05;

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
    const Graph &G0, const int &n, const Parameters &params, const vector<double> &epsilon) {
  Graph G(G0);
  generate_graph(G, n, params);

  vector<int> S = generate_permutation(n, get_graph_size(G0));
  apply_permutation(G, S);

  auto permutations = get_permutation_probabilities(G, get_graph_size(G0), params);
  auto p_uv = get_p_uv_from_permutations(permutations, n, get_graph_size(G0));
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
    const vector<double> &epsilon) {
  Graph G(G0);
  generate_graph(G, n, params);

  vector<int> S = generate_permutation(n, get_graph_size(G0));
  apply_permutation(G, S);

  auto permutations =
      get_permutation_probabilities_sampling(
          G, get_graph_size(G0), params, algorithm, SIGMA_TRIES);
  auto p_uv = get_p_uv_from_permutations(permutations, n, get_graph_size(G0));
  vector<double> solutions;
  for (const double &eps : epsilon) {
    double solution;
    tie(solution, ignore) = LP_solve(p_uv, n, get_graph_size(G0), eps);
    solutions.push_back(solution);
  }
  return solutions;
}

void LP_bound_exact(
    const int &n, const int &n0, const Parameters &params, ostream &out_file) {
  Graph G0 = generate_seed(n0, 1.0);
  vector<double> epsilon;
  for (double eps = EPS_MIN; eps <= 1.0 + 10e-9; eps += EPS_STEP) {
    epsilon.push_back(eps);
  }
  vector<double> solution(epsilon.size(), 0.0);
  vector<vector<double>> solutions(G_TRIES);
  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    solutions[i] = LP_bound_exact_single(G0, n, params, epsilon);
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

void LP_bound_approximate(
    const int &n, const int &n0, const Parameters &params, const SamplingMethod &algorithm,
    ostream &out_file) {
  Graph G0 = generate_seed(n0, 0.6);
  vector<double> epsilon;
  for (double eps = EPS_MIN; eps <= 1.0 + 10e-9; eps += EPS_STEP) {
    epsilon.push_back(eps);
  }
  vector<double> solution(epsilon.size(), 0.0);
  vector<vector<double>> solutions(G_TRIES);
  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    solutions[i] = LP_bound_approximate_single(G0, n, params, algorithm, epsilon);
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

int main(int, char *argv[]) {
  try {
    string action(argv[1]), mode(argv[2]);
    int n = stoi(argv[3]), n0 = stoi(argv[4]);
    Parameters params;
    params.initialize(mode, argv + 5);
    string name(TEMP_FOLDER + get_synthetic_filename(n, n0, params, "TC"));
    if (action == "exact") {
      if (!validate_problem_size(n, n0)) {
        throw out_of_range(
            "Graph too large for exact mode: n = " + to_string(n) + ", n0 = " + to_string(n0));
      }
      ofstream out_file(name, ios_base::app);
      LP_bound_exact(n, n0, params, out_file);
    } else if (SAMPLING_METHOD_REVERSE_NAME.count(action)) {
      ofstream out_file(name, ios_base::app);
      LP_bound_approximate(
          n, n0, params, SAMPLING_METHOD_REVERSE_NAME.find(action)->second, out_file);
    } else {
      throw invalid_argument("Invalid action: " + action);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }

  return 0;
}
