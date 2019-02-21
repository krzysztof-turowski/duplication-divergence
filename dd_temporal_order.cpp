// Tools for computation the temporal order for various duplication-divergence models.
// Compile: g++ dd_temporal_order.cpp -O3 -lgmpxx -lgmp -lglpk -o ./dd_temporal_order
// Run: ./dd_temporal_order exact_bound MODE n n0 PARAMETERS

// TODO(kturowski): deal gurobi output suppression and output to cout instead of cerr

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-default"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"
#if __GNUC__ >= 7
  #pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif
#if __GNUC__ >= 6
  #pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif
#ifndef __clang__
  #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include "./dd_koala.h"
#pragma GCC diagnostic pop

#if defined(glpk)
  #include "./dd_glpk.h"
#elif defined(gurobi)
  #include "./dd_gurobi.h"
#endif

#include <gmpxx.h>

#include <future>
#include <random>

using namespace std;

typedef Koala::Graph<int, int> Graph;
typedef Koala::Graph<int, int>::PVertex Vertex;

const int G_TRIES = 100, SIGMA_TRIES = 100;
const bool G_PARALLEL = false, SIGMA_PARALLEL = false;
const double EPS_MIN = 0.2, EPS_STEP = 0.025;

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
  vector<Vertex> V(G.getVertNo());
  G.getVerts(V.begin());
  for (auto v : V) {
    v->setInfo(S[v->getInfo()]);
  }
}

map<mpz_class, double> get_permutation_probabilities(
    Graph &G, const int &n0, const Parameters &params,
    NeighborhoodStructure &aux, vector<int> &S, const double &p_sigma) {
  map<mpz_class, double> permutations;
  if (G.getVertNo() == n0) {
    mpz_class sigma = encode_permutation(S);
    permutations.insert(make_pair(sigma, p_sigma));
    return permutations;
  }

  vector<Vertex> V(G.getVertNo());
  G.getVerts(V.begin());
  Graph H;
  for (auto v : V) {
    if (v->getInfo() < n0) {
      continue;
    }
    double p_v = get_transition_probability(G, params, v, aux);
    if (p_v > 0.0) {
      set<Vertex> neighbors_v = G.getNeighSet(v);
      aux.remove_vertex(v, neighbors_v), S[G.getVertNo() - 1] = v->getInfo(), H.move(G, v);
      assert(aux.verify(G) == true);

      map<mpz_class, double> permutations_v =
        get_permutation_probabilities(G, n0, params, aux, S, p_sigma * p_v);
      permutations.insert(permutations_v.begin(), permutations_v.end());

      G.move(H, v), aux.restore_vertex(v, neighbors_v), S[G.getVertNo() - 1] = -1;
      for (auto &u : neighbors_v) {
        G.addEdge(v, u);
      }
      assert(aux.verify(G) == true);
    }
  }
  return permutations;
}

map<mpz_class, double> get_permutation_probabilities(
    const Graph &G, const int &n0, const Parameters &params) {
  Graph H(G);
  NeighborhoodStructure aux(H);
  vector<int> S(H.getVertNo(), -1);
  for (int i = 0; i < n0; i++) {
    S[i] = i;
  }
  map<mpz_class, double> permutations = get_permutation_probabilities(H, n0, params, aux, S, 1.0);
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

pair<mpz_class, double> get_permutation_sample(
    const Graph &G, const int &n0, const Parameters &params) {
  random_device device;
  mt19937 generator(device());
  Graph H(G);
  NeighborhoodStructure aux(H);

  vector<int> S(H.getVertNo(), -1);
  for (int i = 0; i < n0; i++) {
    S[i] = i;
  }
  double p_sigma = 1.0;
  while (H.getVertNo() > n0) {
    vector<Vertex> V;
    vector<double> P;
    for (auto v = H.getVert(); v; v = H.getVertNext(v)) {
      if (v->getInfo() < n0) {
        continue;
      }
      V.push_back(v);
      P.push_back(get_transition_probability(H, params, v, aux));
    }
    double P_sum = accumulate(P.begin(), P.end(), 0.0);
    discrete_distribution<int> choose_vertex(P.begin(), P.end());
    int index = choose_vertex(generator);
    S[H.getVertNo() - 1] = V[index]->getInfo(), p_sigma *= P_sum;
    assert(aux.verify(H) == true);
    aux.remove_vertex(V[index], H.getNeighSet(V[index])), H.delVert(V[index]);
    assert(aux.verify(H) == true);
  }
  return make_pair(encode_permutation(S), p_sigma);
}

map<mpz_class, double> get_permutation_probabilities_sampling(
    const Graph &G, const int &n0, const Parameters &params, const int &tries) {
  map<mpz_class, double> permutations;
  // TODO(unknown): parallelize
  for (int i = 0; i < tries; i++) {
    pair<mpz_class, double> sigma_with_probability = get_permutation_sample(G, n0, params);
    permutations[sigma_with_probability.first] += sigma_with_probability.second;
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

void print(const vector<double> &epsilon, const vector<double> &solution) {
  for (int i = 0; i < static_cast<int>(solution.size()); i++) {
    cerr << fixed << setw(6) << setprecision(3) << epsilon[i] << " "
        << fixed << setw(6) << setprecision(3) << solution[i] << endl;
  }
  // export to file
}

vector<double> LP_bound_exact_single(
    const Graph &G0, const int &n, const Parameters &params, const vector<double> &epsilon) {
  Graph G(G0);
  generate_graph_koala(G, n, params);

  vector<int> S = generate_permutation(n, G0.getVertNo());
  apply_permutation(G, S);

  map<mpz_class, double> permutations = get_permutation_probabilities(G, G0.getVertNo(), params);
  map<pair<int, int>, double> p_uv = get_p_uv_from_permutations(permutations, n, G0.getVertNo());
  vector<double> solutions;
  for (const double &eps : epsilon) {
    solutions.push_back(LP_solve(p_uv, n, G0.getVertNo(), eps));
  }
  return solutions;
}

vector<double> LP_bound_approximate_single(
    const Graph &G0, const int &n, const Parameters &params, const vector<double> &epsilon) {
  Graph G(G0);
  generate_graph_koala(G, n, params);

  vector<int> S = generate_permutation(n, G0.getVertNo());
  apply_permutation(G, S);

  map<mpz_class, double> permutations =
      get_permutation_probabilities_sampling(G, G0.getVertNo(), params, SIGMA_TRIES);
  map<pair<int, int>, double> p_uv = get_p_uv_from_permutations(permutations, n, G0.getVertNo());
  vector<double> solutions;
  for (const double &eps : epsilon) {
    solutions.push_back(LP_solve(p_uv, n, G0.getVertNo(), eps));
  }
  return solutions;
}

void LP_bound_exact(const int &n, const int &n0, const Parameters &params) {
  Graph G0 = generate_seed_koala(n0, 1.0);
  vector<double> epsilon;
  for (double eps = EPS_STEP; eps <= 1.0 + 10e-9; eps += EPS_STEP) {
    epsilon.push_back(eps);
  }
  vector<double> solution(epsilon.size(), 0.0);
  if (G_PARALLEL) {
    vector<future<vector<double>>> futures(G_TRIES);
    for (int i = 0; i < G_TRIES; i++) {
      futures[i] = async(
          launch::async, &LP_bound_exact_single,
          cref(G0), cref(n), cref(params), cref(epsilon));
    }
    for (int i = 0; i < G_TRIES; i++) {
      vector<double> solution_single = futures[i].get();
      transform(
          solution.begin(), solution.end(), solution_single.begin(), solution.begin(),
          std::plus<double>());
      cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
    }
  } else {
    for (int i = 0; i < G_TRIES; i++) {
      vector<double> solution_single = LP_bound_exact_single(G0, n, params, epsilon);
      transform(
          solution.begin(), solution.end(), solution_single.begin(), solution.begin(),
          std::plus<double>());
      cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
    }
  }
  for (auto &s : solution) {
    s /= G_TRIES;
  }
  print(epsilon, solution);
}

void LP_bound_approximate(const int &n, const int &n0, const Parameters &params) {
  Graph G0 = generate_seed_koala(n0, 1.0);
  vector<double> epsilon;
  for (double eps = EPS_MIN; eps <= 1.0 + 10e-9; eps += EPS_STEP) {
    epsilon.push_back(eps);
  }
  vector<double> solution(epsilon.size(), 0.0);
  if (G_PARALLEL) {
    vector<future<vector<double>>> futures(G_TRIES);
    for (int i = 0; i < G_TRIES; i++) {
      futures[i] = async(
          launch::async, &LP_bound_approximate_single,
          cref(G0), cref(n), cref(params), cref(epsilon));
    }
    for (int i = 0; i < G_TRIES; i++) {
      vector<double> solution_single = futures[i].get();
      transform(
          solution.begin(), solution.end(), solution_single.begin(), solution.begin(),
          std::plus<double>());
      cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
    }
  } else {
    for (int i = 0; i < G_TRIES; i++) {
      vector<double> solution_single = LP_bound_approximate_single(G0, n, params, epsilon);
      transform(
          solution.begin(), solution.end(), solution_single.begin(), solution.begin(),
          std::plus<double>());
      cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
    }
  }
  for (auto &sol : solution) {
    sol /= G_TRIES;
  }
  print(epsilon, solution);
}

double mean_square_error(const map<mpz_class, double> &opt, const map<mpz_class, double> &apx) {
  double mse = 0;
  for (auto &sigma : opt) {
    mse += pow(opt.find(sigma.first)->second - apx.find(sigma.first)->second, 2);
  }
  return mse / opt.size();
}

void compare_probabilities(const int &n, const int &n0, const Parameters &params) {
  Graph G0 = generate_seed_koala(n0, 1.0);
  map<int, double> mse;
  for (int sigma_tries = 10; sigma_tries < 12000; sigma_tries *= 2) {
    mse.insert(make_pair(sigma_tries, 0.0));
  }
  for (int i = 0; i < G_TRIES; i++) {
    Graph G(G0);
    generate_graph_koala(G, n, params);

    vector<int> S = generate_permutation(n, G0.getVertNo());
    apply_permutation(G, S);

    map<mpz_class, double> permutations_opt =
        get_permutation_probabilities(G, G0.getVertNo(), params);
    for (auto &it : mse) {
      map<mpz_class, double> permutations_apx =
          get_permutation_probabilities_sampling(G, G0.getVertNo(), params, it.first);
      it.second += mean_square_error(permutations_opt, permutations_apx);
    }
    cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
  }
  for (auto &it : mse) {
    cerr << setw(6) << it.first << " "
        << fixed << setw(13) << setprecision(10) << it.second / G_TRIES << endl;
  }
}

void validate_problem_size(const int &n, const int &n0) {
  if (exp(lgamma(n) - lgamma(n0)) > 10e8) {
    throw out_of_range(
        "Graph too large for exact_bound mode: n = " + to_string(n) + ", n0 = " + to_string(n0));
  }
}

int main(int, char *argv[]) {
  try {
    string action(argv[1]), mode(argv[2]);
    int n = stoi(argv[3]), n0 = stoi(argv[4]);
    Parameters params;
    params.initialize(mode, argv + 5);
    if (action == "exact_bound") {
      validate_problem_size(n, n0);
      LP_bound_exact(n, n0, params);
    } else if (action == "approximate_bound") {
      LP_bound_approximate(n, n0, params);
    } else if (action == "compare_probabilities") {
      validate_problem_size(n, n0);
      compare_probabilities(n, n0, params);
    } else {
      throw invalid_argument("Invalid action: " + action);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }

  return 0;
}
