// Tools for checking p_uv estimation for various duplication-divergence models.

#include "./dd_temporal.h"

#include <queue>

using namespace std;

const int G_TRIES = 1, SIGMA_TRIES = 10000;
const int MIN_TRIES_TEST = 1000, MAX_TRIES_TEST = 1000;
const int PERMUTATION_SIZE_LIMIT = 10, PERMUTATION_COUNT_LIMIT = 10;

class ErrorStruct {
 private:
  int p_uv_counter = 0;
  map<int, long double> p_uv_mse, p_uv_mae, p_uv_mre;

  long double mean_square_error(
      const map<VertexPair, long double> &opt, const map<VertexPair, long double> &apx,
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

  long double max_absolute_error(
      const map<VertexPair, long double> &opt, const map<VertexPair, long double> &apx,
      const int &n0, const int &n) {
    long double mae = 0;
    for (int i = n0; i < n; i++) {
      for (int j = n0; j < n; j++) {
        if (i == j) {
          continue;
        }
        auto uv(make_pair(i, j));
        long double opt_uv = opt.count(uv) ? opt.find(uv)->second : 0.0L;
        long double apx_uv = apx.count(uv) ? apx.find(uv)->second : 0.0L;
        mae = max(mae, fabsl(apx_uv - opt_uv));
      }
    }
    return mae;
  }

  long double max_relative_error(
      const map<VertexPair, long double> &opt, const map<VertexPair, long double> &apx,
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
    this->p_uv_mae[tries] += max_absolute_error(p_uv_opt, p_uv_apx, n0, n);
    this->p_uv_mre[tries] += max_relative_error(p_uv_opt, p_uv_apx, n0, n);
    this->p_uv_counter++;
  }

  long double get_p_uv_mse(const int &tries) const {
    return this->p_uv_mse.find(tries)->second / p_uv_counter;
  }

  long double get_p_uv_mae(const int &tries) const {
    return this->p_uv_mae.find(tries)->second / p_uv_counter;
  }

  long double get_p_uv_mre(const int &tries) const {
    return this->p_uv_mre.find(tries)->second / p_uv_counter;
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

  auto p_uv_mae = [](
      const map<SamplingMethod, ErrorStruct> &e,
      const SamplingMethod &method, const int &tries) -> long double {
        return e.find(method)->second.get_p_uv_mae(tries);
      };
  cout << "Max absolute errors for p_uv: " << endl;
  print_errors(sigma_tries, errors, p_uv_mae);

  auto p_uv_mre = [](
      const map<SamplingMethod, ErrorStruct> &e,
      const SamplingMethod &method, const int &tries) -> long double {
        return e.find(method)->second.get_p_uv_mre(tries);
      };
  cout << "Max relative errors for p_uv: " << endl;
  print_errors(sigma_tries, errors, p_uv_mre);
}

void print_best_permutations(
    const map<mpz_class, long double> &permutations, const Graph &G,
    const Parameters &params, const int &n0, const string &algorithm_name, const int &limit) {
  priority_queue<pair<long double, mpz_class>> Q;
  int n = get_graph_size(G);
  for (auto &permutation : permutations) {
    Q.push(make_pair(permutation.second, permutation.first));
  }
  cout << "Best " << min(limit, static_cast<int>(Q.size())) << " permutations for "
      << algorithm_name << " method:" << endl;
  for (int i = 0; i < limit && !Q.empty(); i++) {
    const auto &permutation = Q.top();
    const auto V(decode_permutation(permutation.second, n));
    if (n <= PERMUTATION_SIZE_LIMIT) {
      for (const auto &v : V) {
        cout << v << " ";
      }
    } else {
      cout << "Permutation " << i << " ";
    }
    cout << permutation.first << " "
        << get_log_permutation_probability(G, n0, params, reverse_permutation(V)) << endl;
    Q.pop();
  }
}

void print_compare(
    const map<VertexPair, long double> &opt, const map<VertexPair, long double> &apx,
    const int &n0, const int &n) {
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      auto uv(make_pair(i, j));
      cout << setw(4) << i << setw(4) << j
          << setw(20) << setprecision(9) << (opt.count(uv) ? opt.find(uv)->second : 0.0L) << " "
          << setw(20) << setprecision(9) << (apx.count(uv) ? apx.find(uv)->second : 0.0L) << endl;
    }
  }
}

void check_convergence(const int &n, const int &n0, const Parameters &params) {
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

    map<mpz_class, long double> permutations_opt;
    map<VertexPair, long double> p_uv_opt;
    if (exact_mode) {
      permutations_opt =
          get_log_permutation_probabilities(G, get_graph_size(G0), params);
      normalize_log_probabilities(permutations_opt);
      p_uv_opt = get_p_uv_from_permutations(permutations_opt, n, get_graph_size(G0));
    }
    for (const auto &algorithm : SAMPLING_METHOD_NAME) {
      for (const int &tries : sigma_tries) {
        auto permutations_apx =
            get_log_permutation_probabilities_sampling(
                G, get_graph_size(G0), params, algorithm.first, tries);
        if (!exact_mode) {
          permutations_opt =
                get_log_permutation_probabilities_sampling(
                    G, get_graph_size(G0), params, algorithm.first, tries);
          #pragma omp critical
          {
            print_best_permutations(
                permutations_opt, G, params, n0, algorithm.second, PERMUTATION_COUNT_LIMIT);
            print_best_permutations(
                permutations_apx, G, params, n0, algorithm.second, PERMUTATION_COUNT_LIMIT);
          }
          normalize_log_probabilities(permutations_opt);
          p_uv_opt =
              get_p_uv_from_permutations(permutations_opt, n, get_graph_size(G0));
        }

        normalize_log_probabilities(permutations_apx);
        auto p_uv_apx = get_p_uv_from_permutations(permutations_apx, n, get_graph_size(G0));
        #pragma omp critical
        {
          ErrorStruct &error = errors[algorithm.first];
          error.add_p_uv(tries, p_uv_opt, p_uv_apx, n0, n);
          // print_compare(p_uv_opt, p_uv_apx, n0, n);
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

void check_permutations(const int &n, const int &n0, const Parameters &params) {
  Graph G0 = generate_seed(n0, 1.0);
  bool exact_mode = validate_problem_size(n, n0);

  for (int i = 0; i < G_TRIES; i++) {
    Graph G(G0);
    generate_graph(G, n, params);

    vector<int> S = generate_permutation(n, get_graph_size(G0));
    apply_permutation(G, S);

    if (exact_mode) {
      auto permutations_opt = get_log_permutation_probabilities(G, get_graph_size(G0), params);
      print_best_permutations(permutations_opt, G, params, n0, "exact", PERMUTATION_COUNT_LIMIT);
    }
    for (const auto &algorithm : SAMPLING_METHOD_NAME) {
      auto permutations_apx =
          get_log_permutation_probabilities_sampling(
              G, get_graph_size(G0), params, algorithm.first, SIGMA_TRIES);
      print_best_permutations(
          permutations_apx, G, params, n0, algorithm.second, PERMUTATION_COUNT_LIMIT);
    }
  }
}

int main(int, char *argv[]) {
  try {
    string action(argv[1]), mode(argv[2]);
    int n = stoi(argv[3]), n0 = stoi(argv[4]);
    Parameters params;
    params.initialize(mode, argv + 5);
    if (action == "convergence") {
      check_convergence(n, n0, params);
    } else if (action == "permutations") {
      check_permutations(n, n0, params);
    } else {
      throw invalid_argument("Invalid action: " + action);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }

  return 0;
}
