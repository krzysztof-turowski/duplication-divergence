// Tools for computation the temporal order bound for various duplication-divergence models.
// Compile: g++ dd_temporal_bound.cpp -O3 -lgmpxx -lgmp -lglpk -o ./dd_temporal_bound
// Run: ./dd_temporal_bound [exact|ALGORITHM|check_convergence] MODE n n0 PARAMETERS

// TODO(kturowski): deal gurobi output suppression and output to cout instead of cerr

#include "./dd_graph.h"

#if defined(glpk)
  #include "./dd_glpk.h"
#elif defined(gurobi)
  #include "./dd_gurobi.h"
#endif

#include <gmpxx.h>

#include <random>

using namespace std;

const int G_TRIES = 20, SIGMA_TRIES = 200000;
const double EPS_MIN = 0.2, EPS_STEP = 0.1;
const int MIN_TRIES_TEST = 10, MAX_TRIES_TEST = 20000;

enum SamplingMethod { WIUF, UNIFORM };

const map<SamplingMethod, string> SAMPLING_METHOD_NAME = {
  { SamplingMethod::WIUF, "wiuf" },
  { SamplingMethod::UNIFORM, "uniform" },
};

const map<string, SamplingMethod> SAMPLING_METHOD_REVERSE_NAME = {
  { "wiuf", SamplingMethod::WIUF },
  { "uniform", SamplingMethod::UNIFORM },
};

vector<int> generate_permutation(const int &n, const int &n0) {
  random_device device;
  mt19937 generator(device());
  vector<int> S(n);
  for (int i = 0; i < n; i++) {
    S[i] = i;
  }
  for (int i = n - 1; i > n0; i--) {
    uniform_int_distribution<int> swap_distribution(n0, i);
    int index = swap_distribution(generator);
    swap(S[i], S[index]);
  }
  return S;
}

vector<int> decode_permutation(const mpz_class &sigma, const int &n) {
  vector<int> S(n);
  mpz_class value(sigma);
  for (int i = 1; i <= n; i++) {
    mpz_class base(i + 1), remainder = value % base;
    value /= base;
    S[n - i] = remainder.get_si();
  }
  for (int i = n - 1; i >= 0; i--) {
    for (int j = i + 1; j < n; j++) {
      if (S[i] <= S[j]) {
        S[j]++;
      }
    }
  }
  return S;
}

mpz_class encode_permutation(const vector<int> &S) {
  int n = S.size();
  mpz_class sigma(0);
  vector<int> V(S);
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      if (V[i] < V[j]) {
        V[j]--;
      }
    }
  }
  mpz_class base(1);
  for (int i = 1; i <= n; i++) {
    base *= mpz_class(i);
    sigma += mpz_class(V[n - i]) * base;
  }
  return sigma;
}

void apply_permutation(Graph &G, const vector<int> &S) {
  vector<Vertex> V(get_vertices(G));
  for (auto v : V) {
    set_index(v, S[get_index(v)]);
  }
}

map<mpz_class, double> get_permutation_probabilities(
    Graph &G, const int &n0, const Parameters &params,
    NeighborhoodStructure &aux, vector<int> &S, const double &p_sigma) {
  map<mpz_class, double> permutations;
  if (get_graph_size(G) == n0) {
    mpz_class sigma = encode_permutation(S);
    permutations.insert(make_pair(sigma, p_sigma));
    return permutations;
  }

  vector<Vertex> V(get_vertices(G));
  Graph H;
  for (auto v : V) {
    if (get_index(v) < n0) {
      continue;
    }
    double p_v = get_transition_probability(G, params, v, aux);
    if (p_v > 0.0) {
      set<Vertex> neighbors_v(get_neighbors(G, v));
      aux.remove_vertex(neighbors_v), S[get_graph_size(G) - 1] = get_index(v);
      move_vertex(H, G, v);
      assert(aux.verify(G));

      auto permutations_v = get_permutation_probabilities(G, n0, params, aux, S, p_sigma * p_v);
      permutations.insert(permutations_v.begin(), permutations_v.end());

      move_vertex(G, H, v), aux.restore_vertex(neighbors_v), S[get_graph_size(G) - 1] = -1;
      for (auto &u : neighbors_v) {
        add_edge(G, v, u);
      }
      assert(aux.verify(G));
    }
  }
  return permutations;
}

map<mpz_class, double> get_permutation_probabilities(
    const Graph &G, const int &n0, const Parameters &params) {
  Graph H(G);
  NeighborhoodStructure aux(H);
  vector<int> S(get_graph_size(H), -1);
  for (int i = 0; i < n0; i++) {
    S[i] = i;
  }
  auto permutations = get_permutation_probabilities(H, n0, params, aux, S, 1.0);
  double total_probability = accumulate(
      permutations.begin(), permutations.end(), 0.0,
      [] (double value, const map<mpz_class, double>::value_type &permutation) {
          return value + permutation.second;
      });
  for (auto &permutation : permutations) {
    permutation.second /= total_probability;
  }
  return permutations;
}

tuple<Vertex, double> sample_vertex(
    const vector<Vertex> &V, const vector<double> &P,
    const SamplingMethod &algorithm, mt19937 &generator) {
  switch (algorithm) {
    case WIUF: {
      double P_sum = accumulate(P.begin(), P.end(), 0.0);
      discrete_distribution<int> choose_vertex(P.begin(), P.end());
      int index = choose_vertex(generator);
      return make_tuple(V[index], P_sum);
    }
    case UNIFORM: {
      vector<int> C(P.size());
      transform(
          P.begin(), P.end(), C.begin(), [](const double &value) -> int { return value != 0.0; });
      int C_sum = accumulate(C.begin(), C.end(), 0.0);
      discrete_distribution<int> choose_vertex(C.begin(), C.end());
      int index = choose_vertex(generator);
      return make_tuple(V[index], C_sum * P[index]);
    }
    default:
      throw invalid_argument("Invalid algorithm: " + SAMPLING_METHOD_NAME.find(algorithm)->second);
  }
}

pair<mpz_class, double> get_permutation_sample(
    const Graph &G, const int &n0, const Parameters &params, const SamplingMethod &algorithm) {
  random_device device;
  mt19937 generator(device());
  Graph H(G);
  NeighborhoodStructure aux(H);

  vector<int> S(get_graph_size(H), -1);
  for (int i = 0; i < n0; i++) {
    S[i] = i;
  }
  double p_sigma = 1.0, pv;
  while (get_graph_size(H) > n0) {
    vector<Vertex> V, U(get_vertices(G));
    vector<double> P;
    for (const auto &v : U) {
      if (get_index(v) < n0) {
        continue;
      }
      V.push_back(v);
    }
    for (auto v : V) {
      P.push_back(get_transition_probability(H, params, v, aux));
    }
    Vertex v;
    tie(v, pv) = sample_vertex(V, P, algorithm, generator);
    S[get_graph_size(H) - 1] = get_index(v), p_sigma *= pv;
    assert(aux.verify(H));
    aux.remove_vertex(get_neighbors(H, v)), delete_vertex(H, v);
    assert(aux.verify(H));
  }
  return make_pair(encode_permutation(S), p_sigma);
}

map<mpz_class, double> get_permutation_probabilities_sampling(
    const Graph &G, const int &n0, const Parameters &params, const SamplingMethod &algorithm,
    const int &tries) {
  map<mpz_class, double> permutations;
  for (int i = 0; i < tries; i++) {
    auto sigma_with_probability = get_permutation_sample(G, n0, params, algorithm);
    permutations[sigma_with_probability.first] += sigma_with_probability.second;
    if ((i + 1) % 10000 == 0) {
      #pragma omp critical
      {
        cerr << "Finished tries " << i + 1 << "/" << tries << endl;
      }
    }
  }
  double total_probability = accumulate(
      permutations.begin(), permutations.end(), 0.0,
      [] (double value, const map<mpz_class, double>::value_type &permutation) {
          return value + permutation.second;
      });
  for (auto &permutation : permutations) {
    permutation.second /= total_probability;
  }
  return permutations;
}

map<pair<int, int>, double> get_p_uv_from_permutations(
    const map<mpz_class, double> &permutations, const int &n, const int &n0) {
  map<pair<int, int>, double> p_uv;
  for (auto &permutation : permutations) {
    vector<int> S = decode_permutation(permutation.first, n);
    for (int i = n0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
          p_uv[make_pair(S[i], S[j])] += permutation.second;
      }
    }
  }
  return p_uv;
}

void print_density_precision(
    const string &name, const vector<double> &density, const vector<double> &precision,
    const int &n, const int &n0, const Parameters &params, ostream &out_file) {
  cout << "Graph - n: " << n << ", n0: " << n0
      << ", parameters: " << params.to_string() << endl;
  cerr << "Method: " << name << endl;
  for (int i = 0; i < static_cast<int>(precision.size()); i++) {
    cerr << "density: " << fixed << setw(6) << setprecision(3) << density[i]
        << " precision: " << fixed << setw(6) << setprecision(3) << precision[i] << endl;
  }
  out_file << name << " ";
  for (int i = 0; i < static_cast<int>(precision.size()); i++) {
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
    solutions.push_back(LP_solve(p_uv, n, get_graph_size(G0), eps));
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
    solutions.push_back(LP_solve(p_uv, n, get_graph_size(G0), eps));
  }
  return solutions;
}

void LP_bound_exact(
    const int &n, const int &n0, const Parameters &params, ostream &out_file) {
  Graph G0 = generate_seed_graph(n0, 1.0);
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
  Graph G0 = generate_seed_graph(n0, 1.0);
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

template <typename T>
double mean_square_error(const map<T, double> &opt, const map<T, double> &apx) {
  double mse = 0;
  for (auto &sigma : opt) {
    mse += pow(opt.find(sigma.first)->second - apx.find(sigma.first)->second, 2);
  }
  return mse / opt.size();
}

template <typename T>
double mean_square_error_uniform(const map<T, double> &apx, const double &p, const double &N) {
  double mse = 0;
  for (auto &sigma : apx) {
    mse += pow(p - sigma.second, 2) / p;
  }
  mse += p * p * (N - apx.size());
  return mse / N;
}

template <typename T>
double max_relative_error(const map<T, double> &opt, const map<T, double> &apx) {
  double mre = 0;
  for (auto &uv : opt) {
    double opt_uv = opt.find(uv.first)->second;
    double apx_uv = apx.find(uv.first)->second;
    mre = max(mre, fabs(apx_uv / opt_uv - 1));
  }
  return mre;
}

template <typename T>
double max_relative_error_uniform(const map<T, double> &apx, const double &p) {
  double mre = 0;
  for (auto &uv : apx) {
    mre = max(mre, fabs(uv.second / p - 1));
  }
  return mre;
}

class ErrorStruct {
 private:
  int permutations_counter = 0, p_uv_counter = 0;
  map<int, double> permutations_mse, p_uv_mse, permutations_lambda, p_uv_lambda;

 public:
  void add_permutations(
      const int &tries,
      const map<mpz_class, double> &permutations_opt,
      const map<mpz_class, double> &permutations_apx) {
    this->permutations_mse[tries] += mean_square_error(permutations_opt, permutations_apx);
    this->permutations_lambda[tries] += max_relative_error(permutations_opt, permutations_apx);
    this->permutations_counter++;
  }

  void add_permutations_uniform(
      const int &tries, const map<mpz_class, double> &permutations_apx,
      const int &n, const int &n0) {
    double p = exp(-lgamma(n - n0 + 1));
    this->permutations_mse[tries] += mean_square_error_uniform(permutations_apx, p, 1 / p);
    this->permutations_lambda[tries] += max_relative_error_uniform(permutations_apx, p);
    this->permutations_counter++;
  }

  void add_p_uv(
      const int &tries,
      const map<pair<int, int>, double> &p_uv_opt,
      const map<pair<int, int>, double> &p_uv_apx) {
    this->p_uv_mse[tries] += mean_square_error(p_uv_opt, p_uv_apx);
    this->p_uv_lambda[tries] += max_relative_error(p_uv_opt, p_uv_apx);
    this->p_uv_counter++;
  }

  void add_p_uv_uniform(
      const int &tries, const map<pair<int, int>, double> &p_uv_apx, const int &n, const int &n0) {
    this->p_uv_mse[tries] +=
        mean_square_error_uniform(p_uv_apx, 0.5, (n - n0) * (n - n0 - 1));
    this->p_uv_lambda[tries] += max_relative_error_uniform(p_uv_apx, 0.5);
    this->p_uv_counter++;
  }

  double get_permutations_mse(const int &tries) const {
    return this->permutations_mse.find(tries)->second / permutations_counter;
  }

  double get_p_uv_mse(const int &tries) const {
    return this->p_uv_mse.find(tries)->second / p_uv_counter;
  }

  double get_permutations_lambda(const int &tries) const {
    return this->permutations_lambda.find(tries)->second / permutations_counter;
  }

  double get_p_uv_lambda(const int &tries) const {
    return this->p_uv_lambda.find(tries)->second / p_uv_counter;
  }
};

void print_errors(
    const vector<int> &sigma_tries, const map<SamplingMethod, ErrorStruct> &errors,
    function<double(map<SamplingMethod, ErrorStruct>, SamplingMethod, int)> get_value) {
  cout << setw(6) << "n" << " ";
  for (const auto &algorithm : SAMPLING_METHOD_NAME) {
    cout << setw(12) << algorithm.second << " ";
  }
  cout << endl;
  for (const int &tries : sigma_tries) {
    cout << setw(6) << tries << " ";
    for (const auto &algorithm : SAMPLING_METHOD_NAME) {
      cout << fixed << setw(12) << setprecision(9)
          << get_value(errors, algorithm.first, tries) << " ";
    }
    cout << endl;
  }
}

void print_errors(
    const vector<int> &sigma_tries, const map<SamplingMethod, ErrorStruct> &errors) {
  auto permutations_mse = [](
      const map<SamplingMethod, ErrorStruct> &e,
      const SamplingMethod &method, const int &tries) -> double {
        return e.find(method)->second.get_permutations_mse(tries);
      };
  cout << "Mean square errors for permutations: " << endl;
  print_errors(sigma_tries, errors, permutations_mse);

  auto p_uv_mse = [](
      const map<SamplingMethod, ErrorStruct> &e,
      const SamplingMethod &method, const int &tries) -> double {
        return e.find(method)->second.get_p_uv_mse(tries);
      };
  cout << "Mean square errors for p_uv: " << endl;
  print_errors(sigma_tries, errors, p_uv_mse);

  auto permutations_lambda = [](
      const map<SamplingMethod, ErrorStruct> &e,
      const SamplingMethod &method, const int &tries) -> double {
        return e.find(method)->second.get_permutations_lambda(tries);
      };
  cout << "Max relative errors for permutations: " << endl;
  print_errors(sigma_tries, errors, permutations_lambda);

  auto p_uv_lambda = [](
      const map<SamplingMethod, ErrorStruct> &e,
      const SamplingMethod &method, const int &tries) -> double {
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

  Graph G0 = generate_seed_graph(n0, 1.0);
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

    map<mpz_class, double> permutations_opt;
    map<pair<int, int>, double> p_uv_opt;
    if (exact_mode) {
      permutations_opt = get_permutation_probabilities(G, get_graph_size(G0), params);
      p_uv_opt = get_p_uv_from_permutations(permutations_opt, n, get_graph_size(G0));
    }
    for (const auto &algorithm : SAMPLING_METHOD_NAME) {
      auto permutations_apx_prev =
            get_permutation_probabilities_sampling(
                G, get_graph_size(G0), params, algorithm.first, MIN_TRIES_TEST);
      auto p_uv_apx_prev =
          get_p_uv_from_permutations(permutations_apx_prev, n, get_graph_size(G0));
      for (const int &tries : sigma_tries) {
        auto permutations_apx =
            get_permutation_probabilities_sampling(
                G, get_graph_size(G0), params, algorithm.first, tries);
        auto p_uv_apx = get_p_uv_from_permutations(permutations_apx, n, get_graph_size(G0));

        #pragma omp critical
        {
          ErrorStruct &error = errors[algorithm.first];
          if (exact_mode) {
            error.add_permutations(tries, permutations_opt, permutations_apx);
            error.add_p_uv(tries, p_uv_opt, p_uv_apx);
          } else {
            error.add_permutations(tries, permutations_apx_prev, permutations_apx);
            error.add_p_uv(tries, p_uv_apx_prev, p_uv_apx);
            permutations_apx_prev = permutations_apx, p_uv_apx_prev = p_uv_apx;
          }
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
