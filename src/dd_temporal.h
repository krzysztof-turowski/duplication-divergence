#pragma once

#include "./dd_dag.h"
#include "./dd_graph.h"

#include <gmpxx.h>

#include <algorithm>
#include <cassert>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

typedef std::pair<double, double> DensityPrecision;
typedef std::pair<int, int> VertexPair;
typedef CompleteNeighborhoodStructure NeighborhoodStructure;

const int PERMUTATION_SIZE_LIMIT = 100, PERMUTATION_COUNT_LIMIT = 10;

enum SamplingMethod { WIUF, UNIFORM, MIN_DISCARD };

const std::map<SamplingMethod, std::string> SAMPLING_METHOD_NAME = {
  { SamplingMethod::WIUF, "high-prob-sampling" },
  { SamplingMethod::UNIFORM, "local-unif-sampling" },
  { SamplingMethod::MIN_DISCARD, "discard-minimum-sampling" },
};

const std::map<std::string, SamplingMethod> SAMPLING_METHOD_REVERSE_NAME = {
  { "high-prob-sampling", SamplingMethod::WIUF },
  { "local-unif-sampling", SamplingMethod::UNIFORM },
  { "discard-minimum-sampling", SamplingMethod::MIN_DISCARD },
};

inline bool validate_problem_size(const int &n, const int &n0) {
  return exp(lgamma(n) - lgamma(n0)) <= 10e8;
}

std::vector<int> generate_permutation(const int &n, const int &n0) {
  std::random_device device;
  std::mt19937 generator(device());
  std::vector<int> S(n);
  for (int i = 0; i < n; i++) {
    S[i] = i;
  }
  for (int i = n - 1; i > n0; i--) {
    std::uniform_int_distribution<int> swap_distribution(n0, i);
    int index = swap_distribution(generator);
    std::swap(S[i], S[index]);
  }
  return S;
}

std::vector<int> decode_permutation(const mpz_class &sigma, const int &n) {
  std::vector<int> S(n);
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

mpz_class encode_permutation(const std::vector<int> &S) {
  int n = S.size();
  mpz_class sigma(0);
  std::vector<int> V(S);
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

void apply_permutation(Graph &G, const std::vector<int> &S) {
  std::vector<Vertex> V(get_vertices(G));
  for (auto &v : V) {
    set_index(G, v, S[get_index(G, v)]);
  }
}

std::vector<int> reverse_permutation(const std::vector<int> &S) {
  std::vector<int> V(S.size());
  for (size_t i = 0; i < S.size(); i++) {
    V[S[i]] = i;
  }
  return V;
}

void normalize_log_probabilities(std::map<mpz_class, long double> &permutations) {
  long double max_log_value = -std::numeric_limits<long double>::infinity();
  for (const auto &permutation : permutations) {
    max_log_value = std::max(max_log_value, permutation.second);
  }
  long double total_probability = 0.0L;
  for (auto &permutation : permutations) {
    permutation.second = exp2l(permutation.second - max_log_value);
    total_probability += permutation.second;
  }
  for (auto &permutation : permutations) {
    permutation.second /= total_probability;
  }
}

std::map<mpz_class, long double> get_log_permutation_probabilities(
    Graph &G, const int &n0, const Parameters &params,
    NeighborhoodStructure &aux, std::vector<int> &S, const long double &p_sigma) {
  std::map<mpz_class, long double> permutations;
  if (get_graph_size(G) == n0) {
    mpz_class sigma = encode_permutation(S);
    permutations.insert(std::make_pair(sigma, p_sigma));
    return permutations;
  }

  std::vector<Vertex> V(get_vertices(G));
  Graph H;
  for (auto &v : V) {
    if (get_index(G, v) < n0) {
      continue;
    }
    long double p_v = get_log_transition_probability(G, params, v, aux);
    if (std::isfinite(p_v)) {
      std::set<Vertex> neighbors_v(get_neighbors(G, v));
      aux.remove_vertex(neighbors_v), S[get_graph_size(G) - 1] = get_index(G, v);
      move_vertex(H, G, v);
      assert(aux.verify(G));

      auto permutations_v =
          get_log_permutation_probabilities(G, n0, params, aux, S, p_sigma + p_v);
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

std::map<mpz_class, long double> get_log_permutation_probabilities(
    const Graph &G, const int &n0, const Parameters &params) {
  Graph H(G);
  NeighborhoodStructure aux(H);
  std::vector<int> S(get_graph_size(H), -1);
  for (int i = 0; i < n0; i++) {
    S[i] = i;
  }
  return get_log_permutation_probabilities(H, n0, params, aux, S, 0.0L);
}

long double get_log_permutation_probability(
    const Graph &G, const int &n0, const Parameters &params, const std::vector<int> &rev_S) {
  Graph H(G);
  CompleteNeighborhoodStructure aux(H);
  long double p_sigma = 0.0L;
  std::map<int, Vertex> V;
  for (const auto &v : get_vertices(H)) {
    V.insert(std::make_pair(rev_S[get_index(H, v)], v));
  }
  for (int i = V.size() - 1; i >= n0; i--) {
    long double p_v = get_log_transition_probability(H, params, V[i], aux);
    aux.remove_vertex(get_neighbors(H, V[i])), delete_vertex(H, V[i]);
    assert(aux.verify(H));
    p_sigma += p_v;
  }
  return p_sigma;
}

std::tuple<Vertex, double> sample_vertex(
    const Graph &G, const int &n0, const Parameters &params, const DAG &perfect_pairs,
    const NeighborhoodStructure &aux, const SamplingMethod &algorithm, std::mt19937 &generator) {
  std::vector<Vertex> V;
  std::vector<long double> omega;
  for (const auto &v : get_vertices(G)) {
    if (get_index(G, v) < n0) {
      continue;
    }
    if (!perfect_pairs.is_source(get_index(G, v))) {
      continue;
    }
    V.push_back(v);
  }
  switch (algorithm) {
    case WIUF: {
      for (const auto &v : V) {
        omega.push_back(exp2l(get_log_transition_probability(G, params, v, aux)));
      }
      long double omega_sum = accumulate(omega.begin(), omega.end(), 0.0L);
      std::discrete_distribution<int> choose_vertex(omega.begin(), omega.end());
      int index = choose_vertex(generator);
      return std::make_tuple(V[index], log2l(omega_sum));
    }
    case UNIFORM: {
      std::vector<int> C(V.size());
      std::transform(
          V.begin(), V.end(), C.begin(),
          [&](const Vertex &v) -> int { return is_feasible(G, params, v, aux); });
      int C_sum = accumulate(C.begin(), C.end(), 0);
      std::discrete_distribution<int> choose_vertex(C.begin(), C.end());
      int index = choose_vertex(generator);
      return std::make_tuple(
          V[index], log2l(C_sum) + get_log_transition_probability(G, params, V[index], aux));
    }
    case MIN_DISCARD: {
      std::vector<long double> C(V.size());
      std::transform(
          V.begin(), V.end(), C.begin(),
          [&](const Vertex &v) -> long double { return get_discard_score(G, params, v, aux); });
      long double C_sum = accumulate(C.begin(), C.end(), 0.0L);
      std::discrete_distribution<int> choose_vertex(C.begin(), C.end());
      int index = choose_vertex(generator);
      return std::make_tuple(
          V[index], log2l(C_sum) + get_log_transition_probability(G, params, V[index], aux));
    }
    default:
      throw std::invalid_argument(
          "Invalid algorithm: " + SAMPLING_METHOD_NAME.find(algorithm)->second);
  }
}

std::pair<mpz_class, long double> get_log_permutation_sample(
    const Graph &G, const int &n0, const Parameters &params, const DAG &perfect_pairs_G,
    const NeighborhoodStructure &aux_G, const SamplingMethod &algorithm) {
  std::random_device device;
  std::mt19937 generator(device());
  Graph H(G);
  NeighborhoodStructure aux(aux_G);
  DAG perfect_pairs(perfect_pairs_G);

  std::vector<int> S(get_graph_size(H), -1);
  for (int i = 0; i < n0; i++) {
    S[i] = i;
  }
  long double p_sigma = 0.0L, pv;
  while (get_graph_size(H) > n0) {
    Vertex v;
    std::tie(v, pv) =
        sample_vertex(H, n0, params, perfect_pairs, aux, algorithm, generator);
    S[get_graph_size(H) - 1] = get_index(H, v), p_sigma += pv;
    assert(aux.verify(H));
    perfect_pairs.remove_vertex(get_index(H, v));
    aux.remove_vertex(get_neighbors(H, v)), delete_vertex(H, v);
    assert(aux.verify(H));
  }
  return std::make_pair(encode_permutation(S), p_sigma);
}

std::map<mpz_class, long double> get_log_permutation_probabilities_sampling(
    const Graph &G, const int &n0, const Parameters &params, const DAG &perfect_pairs,
    const SamplingMethod &algorithm, const int &tries) {
  std::map<mpz_class, long double> permutations;
  NeighborhoodStructure aux(G);
  #pragma omp parallel for shared(permutations)
  for (int i = 0; i < tries; i++) {
    auto sigma_with_probability =
        get_log_permutation_sample(G, n0, params, perfect_pairs, aux, algorithm);
    auto permutation = permutations.find(sigma_with_probability.first);
    if (permutation != permutations.end()) {
      permutation->second = add_exp_log(permutation->second, sigma_with_probability.second);
    } else {
      permutations.insert(sigma_with_probability);
    }
    #pragma omp critical
    {
      if ((i + 1) % 10000 == 0) {
        std::cerr << "Finished tries " << i + 1 << "/" << tries << std::endl;
      }
    }
  }
  return permutations;
}

std::map<VertexPair, long double> get_p_uv_from_permutations(
    const std::map<mpz_class, long double> &permutations, const int &n, const int &n0) {
  std::map<VertexPair, long double> p_uv;
  const long double P_SIGMA_THRESHOLD = 10e-20L;
  for (const auto &permutation : permutations) {
    if (permutation.second < P_SIGMA_THRESHOLD) {
      continue;
    }
    const auto &S = decode_permutation(permutation.first, n);
    for (int j = n0; j < n; j++) {
      for (int k = j + 1; k < n; k++) {
        p_uv[std::make_pair(S[j], S[k])] += permutation.second;
      }
    }
  }
  return p_uv;
}

void print_best_permutations(
    const std::map<mpz_class, long double> &permutations, const Graph &G,
    const Parameters &params, const int &n0, const std::string &algorithm_name, const int &limit) {
  std::priority_queue<std::pair<long double, mpz_class>> Q;
  int n = get_graph_size(G);
  for (auto &permutation : permutations) {
    Q.push(std::make_pair(permutation.second, permutation.first));
  }
  std::cout << "Best " << std::min(limit, static_cast<int>(Q.size())) << "/" << permutations.size()
      << " permutations for " << algorithm_name << " method:" << std::endl;
  for (int i = 0; i < limit && !Q.empty(); i++) {
    const auto &permutation = Q.top();
    const auto V(decode_permutation(permutation.second, n));
    if (n <= PERMUTATION_SIZE_LIMIT) {
      for (const auto &v : V) {
        std::cout << v << " ";
      }
    } else {
      std::cout << "Permutation " << i << " ";
    }
    std::cout << std::setw(20) << std::setprecision(9) << permutation.first << " "
        << std::setw(20) << std::setprecision(9)
        << get_log_permutation_probability(G, n0, params, reverse_permutation(V)) << std::endl;
    Q.pop();
  }
}
