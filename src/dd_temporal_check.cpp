// Tools for checking p_uv estimation for various duplication-divergence models.
// Compile: g++ dd_temporal_check.cpp -O3 -o ./dd_temporal_check
// Example run:
//  ./dd_temporal_check -action:permutations -n:100 -n0:10
//      -mode:pastor_satorras -p:0.5 -r:2.0 -p0:0.6

#include "./dd_input.h"
#include "./dd_perfect_pairs.h"
#include "./dd_temporal.h"

#include <queue>

int G_TRIES, SIGMA_TRIES, MIN_SIGMA_TRIES, MAX_SIGMA_TRIES;

class ErrorStruct {
 private:
  int p_uv_counter = 0;
  std::map<int, long double> p_uv_mse, p_uv_mae, p_uv_mre;

  long double mean_square_error(
      const std::map<VertexPair, long double> &opt, const std::map<VertexPair, long double> &apx,
      const int &n0, const int &n) {
    long double mse = 0;
    for (int i = n0; i < n; i++) {
      for (int j = n0; j < n; j++) {
        if (i == j) {
          continue;
        }
        auto uv(std::make_pair(i, j));
        long double opt_uv = opt.count(uv) ? opt.find(uv)->second : 0.0L;
        long double apx_uv = apx.count(uv) ? apx.find(uv)->second : 0.0L;
        mse += powl(opt_uv - apx_uv, 2);
      }
    }
    return mse / opt.size();
  }

  long double max_absolute_error(
      const std::map<VertexPair, long double> &opt, const std::map<VertexPair, long double> &apx,
      const int &n0, const int &n) {
    long double mae = 0;
    for (int i = n0; i < n; i++) {
      for (int j = n0; j < n; j++) {
        if (i == j) {
          continue;
        }
        auto uv(std::make_pair(i, j));
        long double opt_uv = opt.count(uv) ? opt.find(uv)->second : 0.0L;
        long double apx_uv = apx.count(uv) ? apx.find(uv)->second : 0.0L;
        mae = std::max(mae, fabsl(apx_uv - opt_uv));
      }
    }
    return mae;
  }

  long double max_relative_error(
      const std::map<VertexPair, long double> &opt, const std::map<VertexPair, long double> &apx,
      const int &n0, const int &n) {
    long double mre = 0;
    for (int i = n0; i < n; i++) {
      for (int j = n0; j < n; j++) {
        if (i == j) {
          continue;
        }
        auto uv(std::make_pair(i, j));
        long double opt_uv = opt.count(uv) ? opt.find(uv)->second : 0.0L;
        long double apx_uv = apx.count(uv) ? apx.find(uv)->second : 0.0L;
        if (opt_uv > 0.0L) {
          mre = std::max(mre, fabsl(apx_uv / opt_uv - 1.0L));
        }
      }
    }
    return mre;
  }

 public:
  void add_p_uv(
      const int &tries,
      const std::map<VertexPair, long double> &p_uv_opt,
      const std::map<VertexPair, long double> &p_uv_apx,
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
    const std::vector<int> &sigma_tries, const std::map<SamplingMethod, ErrorStruct> &errors,
    std::function<double(std::map<SamplingMethod, ErrorStruct>, SamplingMethod, int)> get_value) {
  std::cout << std::setw(6) << "n" << " ";
  for (const auto &algorithm : SAMPLING_METHOD_NAME) {
    std::cout << std::setw(20) << algorithm.second << " ";
  }
  std::cout << std::endl;
  for (const int &tries : sigma_tries) {
    std::cout << std::setw(6) << tries << " ";
    for (const auto &algorithm : SAMPLING_METHOD_NAME) {
      std::cout << std::fixed << std::setw(20) << std::setprecision(9)
          << get_value(errors, algorithm.first, tries) << " ";
    }
    std::cout << std::endl;
  }
}

void print_errors(
    const std::vector<int> &sigma_tries, const std::map<SamplingMethod, ErrorStruct> &errors) {
  auto p_uv_mse = [](
      const std::map<SamplingMethod, ErrorStruct> &e,
      const SamplingMethod &method, const int &tries) -> long double {
        return e.find(method)->second.get_p_uv_mse(tries);
      };
  std::cout << "Mean square errors for p_uv: " << std::endl;
  print_errors(sigma_tries, errors, p_uv_mse);

  auto p_uv_mae = [](
      const std::map<SamplingMethod, ErrorStruct> &e,
      const SamplingMethod &method, const int &tries) -> long double {
        return e.find(method)->second.get_p_uv_mae(tries);
      };
  std::cout << "Max absolute errors for p_uv: " << std::endl;
  print_errors(sigma_tries, errors, p_uv_mae);

  auto p_uv_mre = [](
      const std::map<SamplingMethod, ErrorStruct> &e,
      const SamplingMethod &method, const int &tries) -> long double {
        return e.find(method)->second.get_p_uv_mre(tries);
      };
  std::cout << "Max relative errors for p_uv: " << std::endl;
  print_errors(sigma_tries, errors, p_uv_mre);
}

void print_compare(
    const std::map<VertexPair, long double> &opt, const std::map<VertexPair, long double> &apx,
    const int &n0, const int &n) {
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      auto uv(std::make_pair(i, j));
      std::cout << std::setw(4) << i << std::setw(4) << j
          << std::setw(20) << std::setprecision(9)
              << (opt.count(uv) ? opt.find(uv)->second : 0.0L) << " "
          << std::setw(20) << std::setprecision(9)
              << (apx.count(uv) ? apx.find(uv)->second : 0.0L) << std::endl;
    }
  }
}

void check_convergence(const int &n, const int &n0, const double &p0, const Parameters &params) {
  bool exact_mode = validate_problem_size(n, n0);

  Graph G0 = generate_seed(n0, p0);
  std::vector<int> sigma_tries;
  for (int tries = MIN_SIGMA_TRIES; tries <= MAX_SIGMA_TRIES; tries *= 2) {
    sigma_tries.push_back(tries);
  }

  std::map<SamplingMethod, ErrorStruct> errors;
  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    Graph G(G0);
    generate_graph(G, n, params);

    std::vector<int> S = generate_permutation(n, get_graph_size(G0));
    apply_permutation(G, S);

    DAG perfect_pairs(get_graph_size(G));

    std::map<mpz_class, long double> permutations_opt;
    std::map<VertexPair, long double> p_uv_opt;
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
                G, get_graph_size(G0), params, perfect_pairs, algorithm.first, tries);
        if (!exact_mode) {
          permutations_opt =
                get_log_permutation_probabilities_sampling(
                    G, get_graph_size(G0), params, perfect_pairs, algorithm.first, tries);
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
      std::cerr << "Finished run " << i + 1 << "/" << G_TRIES << std::endl;
    }
  }
  print_errors(sigma_tries, errors);
}

void check_permutations(const int &n, const int &n0, const double &p0, const Parameters &params) {
  Graph G0 = generate_seed(n0, p0);
  bool exact_mode = validate_problem_size(n, n0);

  for (int i = 0; i < G_TRIES; i++) {
    Graph G(G0);
    generate_graph(G, n, params);

    std::vector<int> S = generate_permutation(n, get_graph_size(G0));
    apply_permutation(G, S);

    DAG perfect_pairs(get_graph_size(G));

    if (exact_mode) {
      auto permutations_opt = get_log_permutation_probabilities(G, get_graph_size(G0), params);
      print_best_permutations(permutations_opt, G, params, n0, "exact", PERMUTATION_COUNT_LIMIT);
    }
    for (const auto &algorithm : SAMPLING_METHOD_NAME) {
      auto permutations_apx =
          get_log_permutation_probabilities_sampling(
              G, get_graph_size(G0), params, perfect_pairs, algorithm.first, SIGMA_TRIES);
      print_best_permutations(
          permutations_apx, G, params, n0, algorithm.second, PERMUTATION_COUNT_LIMIT);
    }
  }
}

void check_random_walk(const int &n, const int &n0, const double &p0, const Parameters &params) {
  Graph G0 = generate_seed(n0, p0);

  for (int i = 0; i < G_TRIES; i++) {
    Graph G(G0);
    generate_graph(G, n, params);

    std::vector<int> S = generate_permutation(n, n0);
    apply_permutation(G, S);

    std::map<mpz_class, long double> permutations;
    auto sigma = encode_permutation(S);
    long double score = get_log_permutation_probability(G, n0, params, reverse_permutation(S));
    permutations.insert(std::make_pair(sigma, score));

    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_real_distribution<double> pick_distribution(0.0, 1.0);
    long double best_score = score;
    while (static_cast<int>(permutations.size()) < SIGMA_TRIES) {
      std::uniform_int_distribution<int> swap_distribution(n0, n - 1);
      int u = swap_distribution(generator), v = swap_distribution(generator);
      if (u == v) {
        continue;
      }
      std::swap(S[u], S[v]), sigma = encode_permutation(S);
      score = get_log_permutation_probability(G, n0, params, reverse_permutation(S));
      if (pick_distribution(generator) <= exp2(score - best_score)) {
        permutations.insert(std::make_pair(sigma, score));
        best_score = std::max(best_score, score);
        #pragma omp critical
        {
          if (permutations.size() % 10 == 0) {
            std::cerr << "Finished run " << permutations.size() << "/" << SIGMA_TRIES << std::endl;
          }
        }
      } else {
        std::swap(S[u], S[v]);
      }
    }
    print_best_permutations(permutations, G, params, n0, "random_walk", PERMUTATION_COUNT_LIMIT);
    normalize_log_probabilities(permutations);
    auto p_uv = get_p_uv_from_permutations(permutations, n, n0);
    for (int u = n0; u < n; u++) {
      for (int v = n0; v < n; v++) {
        auto uv(std::make_pair(u, v));
        std::cout << std::setw(4) << u << std::setw(4) << v
            << std::setw(20) << std::setprecision(9)
            << (p_uv.count(uv) ? p_uv.find(uv)->second : 0.0L) << std::endl;
      }
    }
  }
}

int main(int argc, char **argv) {
  try {
    Env = prepare_environment(argc, argv);
    G_TRIES = read_int(Env, "-gt:", 1, "G_TRIES");
    std::string action = read_action(Env);
    const int n = read_n(Env), n0 = read_n0(Env);
    const double p0 = read_p0(Env);
    const Parameters params = read_parameters(Env);
    if (action == "convergence") {
      MIN_SIGMA_TRIES = read_int(Env, "-gt:", 1, "MIN_SIGMA_TRIES");
      MAX_SIGMA_TRIES = read_int(Env, "-gt:", 1, "MAX_SIGMA_TRIES");
      check_convergence(n, n0, p0, params);
    } else if (action == "permutations") {
      SIGMA_TRIES = read_int(Env, "-st:", 1, "SIGMA_TRIES");
      check_permutations(n, n0, p0, params);
    } else if (action == "random_walk") {
      SIGMA_TRIES = read_int(Env, "-st:", 1, "SIGMA_TRIES");
      check_random_walk(n, n0, p0, params);
    } else {
      throw std::invalid_argument("Invalid action: " + action);
    }
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
  }

  return 0;
}
