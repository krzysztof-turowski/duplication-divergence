// Tools for computation the temporal order bound for various duplication-divergence models.
// Compile: g++ dd_temporal_bound.cpp -O3 -lgmpxx -lgmp -lglpk -o ./dd_temporal_bound
// Run: ./dd_temporal_bound [exact|ALGORITHM|check_convergence] MODE n n0 PARAMETERS

// TODO(kturowski): deal gurobi output suppression and output to cout instead of cerr

#include "./dd_temporal.h"

#if defined(glpk)
  #include "./dd_glpk.h"
#elif defined(gurobi)
  #include "./dd_gurobi.h"
#endif

using namespace std;

const int G_TRIES = 20, SIGMA_TRIES = 200000;
const double EPS_MIN = 0.2, EPS_STEP = 0.1;
const int MIN_TRIES_TEST = 10, MAX_TRIES_TEST = 20000;

void print_density_precision(
    const string &name, const vector<double> &density, const vector<double> &precision,
    const int &n, const int &n0, const Parameters &params, ostream &out_file) {
  cout << "Graph - n: " << n << ", n0: " << n0
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
  Graph G0 = generate_seed(n0, 1.0);
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

class ErrorStruct {
 private:
  int p_uv_counter = 0;
  map<int, long double> p_uv_mse, p_uv_lambda;

  template <typename T>
  long double mean_square_error(
      const map<T, long double> &opt, const map<T, long double> &apx,
      const int &n0, const int &n) {
    long double mse = 0;
    for (int i = n0; i < n; i++) {
      for (int j = n0; j < n; j++) {
        if (i == j) {
          continue;
        }
        auto uv(make_pair(i, j));
        long double opt_uv = opt.count(uv) ? opt.find(uv)->second : 0.0L;
        long double apx_uv = apx.count(uv) ? apx.find(uv)->second : 0.0L;
        mse += powl(opt_uv - apx_uv, 2);
      }
    }
    return mse / opt.size();
  }

  template <typename T>
  long double max_relative_error(
      const map<T, long double> &opt, const map<T, long double> &apx,
      const int &n0, const int &n) {
    long double mre = 0;
    for (int i = n0; i < n; i++) {
      for (int j = n0; j < n; j++) {
        if (i == j) {
          continue;
        }
        auto uv(make_pair(i, j));
        long double opt_uv = opt.count(uv) ? opt.find(uv)->second : 0.0L;
        long double apx_uv = apx.count(uv) ? apx.find(uv)->second : 0.0L;
        if (opt_uv > 0.0L) {
          mre = max(mre, fabsl(apx_uv / opt_uv - 1.0L));
        }
      }
    }
    return mre;
  }

 public:
  void add_p_uv(
      const int &tries,
      const map<VertexPair, long double> &p_uv_opt,
      const map<VertexPair, long double> &p_uv_apx,
      const int &n0, const int &n) {
    this->p_uv_mse[tries] += mean_square_error(p_uv_opt, p_uv_apx, n0, n);
    this->p_uv_lambda[tries] += max_relative_error(p_uv_opt, p_uv_apx, n0, n);
    this->p_uv_counter++;
  }

  long double get_p_uv_mse(const int &tries) const {
    return this->p_uv_mse.find(tries)->second / p_uv_counter;
  }

  long double get_p_uv_lambda(const int &tries) const {
    return this->p_uv_lambda.find(tries)->second / p_uv_counter;
  }
};

void print_errors(
    const vector<int> &sigma_tries, const map<SamplingMethod, ErrorStruct> &errors,
    function<double(map<SamplingMethod, ErrorStruct>, SamplingMethod, int)> get_value) {
  cout << setw(6) << "n" << " ";
  for (const auto &algorithm : SAMPLING_METHOD_NAME) {
    cout << setw(20) << algorithm.second << " ";
  }
  cout << endl;
  for (const int &tries : sigma_tries) {
    cout << setw(6) << tries << " ";
    for (const auto &algorithm : SAMPLING_METHOD_NAME) {
      cout << fixed << setw(20) << setprecision(9)
          << get_value(errors, algorithm.first, tries) << " ";
    }
    cout << endl;
  }
}

void print_errors(
    const vector<int> &sigma_tries, const map<SamplingMethod, ErrorStruct> &errors) {
  auto p_uv_mse = [](
      const map<SamplingMethod, ErrorStruct> &e,
      const SamplingMethod &method, const int &tries) -> long double {
        return e.find(method)->second.get_p_uv_mse(tries);
      };
  cout << "Mean square errors for p_uv: " << endl;
  print_errors(sigma_tries, errors, p_uv_mse);

  auto p_uv_lambda = [](
      const map<SamplingMethod, ErrorStruct> &e,
      const SamplingMethod &method, const int &tries) -> long double {
        return e.find(method)->second.get_p_uv_lambda(tries);
      };
  cout << "Max relative errors for p_uv: " << endl;
  print_errors(sigma_tries, errors, p_uv_lambda);
}

inline bool validate_problem_size(const int &n, const int &n0) {
  return exp(lgamma(n) - lgamma(n0)) <= 10e8;
}

void compare_probabilities(const int &n, const int &n0, const Parameters &params) {
  bool exact_mode = validate_problem_size(n, n0);

  Graph G0 = generate_seed(n0, 1.0);
  vector<int> sigma_tries;
  for (int tries = MIN_TRIES_TEST; tries <= MAX_TRIES_TEST; tries *= 2) {
    sigma_tries.push_back(tries);
  }

  map<SamplingMethod, ErrorStruct> errors;

  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    Graph G(G0);
    generate_graph(G, n, params);

    vector<int> S = generate_permutation(n, get_graph_size(G0));
    apply_permutation(G, S);

    map<VertexPair, long double> p_uv_opt;
    if (exact_mode) {
      auto permutations_opt = get_permutation_probabilities(G, get_graph_size(G0), params);
      p_uv_opt = get_p_uv_from_permutations(permutations_opt, n, get_graph_size(G0));
    }
    for (const auto &algorithm : SAMPLING_METHOD_NAME) {
      for (const int &tries : sigma_tries) {
        auto permutations_apx =
            get_permutation_probabilities_sampling(
                G, get_graph_size(G0), params, algorithm.first, tries);
        auto p_uv_apx = get_p_uv_from_permutations(permutations_apx, n, get_graph_size(G0));

        if (!exact_mode) {
          auto permutations_opt =
                get_permutation_probabilities_sampling(
                    G, get_graph_size(G0), params, algorithm.first, tries);
          p_uv_opt =
              get_p_uv_from_permutations(permutations_opt, n, get_graph_size(G0));
        }
        #pragma omp critical
        {
          ErrorStruct &error = errors[algorithm.first];
          error.add_p_uv(tries, p_uv_opt, p_uv_apx, n0, n);
        }
      }
    }
    #pragma omp critical
    {
      cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
    }
  }
  print_errors(sigma_tries, errors);
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
    } else if (action == "check_convergence") {
      compare_probabilities(n, n0, params);
    } else {
      throw invalid_argument("Invalid action: " + action);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }

  return 0;
}
