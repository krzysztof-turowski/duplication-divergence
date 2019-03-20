#pragma once

#include "./dd_graph.h"

#include <gmpxx.h>

#include <cassert>
#include <map>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

typedef std::pair<double, double> DensityPrecision;
typedef std::pair<int, int> VertexPair;

enum SamplingMethod { WIUF, UNIFORM };

const std::map<SamplingMethod, std::string> SAMPLING_METHOD_NAME = {
  { SamplingMethod::WIUF, "wiuf" },
  { SamplingMethod::UNIFORM, "uniform" },
};

const std::map<std::string, SamplingMethod> SAMPLING_METHOD_REVERSE_NAME = {
  { "wiuf", SamplingMethod::WIUF },
  { "uniform", SamplingMethod::UNIFORM },
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

std::map<mpz_class, long double> get_permutation_probabilities(
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
          get_permutation_probabilities(G, n0, params, aux, S, p_sigma + p_v);
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

std::map<mpz_class, long double> get_permutation_probabilities(
    const Graph &G, const int &n0, const Parameters &params) {
  Graph H(G);
  CompleteNeighborhoodStructure aux(H);
  std::vector<int> S(get_graph_size(H), -1);
  for (int i = 0; i < n0; i++) {
    S[i] = i;
  }
  auto permutations = get_permutation_probabilities(H, n0, params, aux, S, 0.0L);
  long double total_probability = accumulate(
      permutations.begin(), permutations.end(), 0.0L,
      [] (long double value, const std::map<mpz_class, long double>::value_type &permutation) {
          return value + exp2l(permutation.second);
      });
  for (auto &permutation : permutations) {
    permutation.second = exp2l(permutation.second) / total_probability;
  }
  return permutations;
}

std::tuple<Vertex, double> sample_vertex(
    const Graph &G, const int &n0, const Parameters &params, const NeighborhoodStructure &aux,
    const SamplingMethod &algorithm, std::mt19937 &generator) {
  std::vector<Vertex> V;
  std::vector<long double> omega;
  for (const auto &v : get_vertices(G)) {
    if (get_index(G, v) < n0) {
      continue;
    }
    V.push_back(v);
  }
  switch (algorithm) {
    case WIUF: {
      for (const auto &v : V) {
        omega.push_back(exp2l(get_log_transition_probability(G, params, v, aux)));
      }
      long double omega_sum = accumulate(omega.begin(), omega.end(), 0.0);
      std::discrete_distribution<int> choose_vertex(omega.begin(), omega.end());
      int index = choose_vertex(generator);
      return std::make_tuple(V[index], log2l(omega_sum));
    }
    case UNIFORM: {
      std::vector<int> C(V.size());
      std::transform(
          V.begin(), V.end(), C.begin(),
          [&](const Vertex &v) -> int { return is_feasible(G, params, v, aux); });
      int C_sum = accumulate(C.begin(), C.end(), 0.0);
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

std::pair<mpz_class, long double> get_permutation_sample(
    const Graph &G, const int &n0, const Parameters &params,
    const CompleteNeighborhoodStructure &aux_G, const SamplingMethod &algorithm) {
  std::random_device device;
  std::mt19937 generator(device());
  Graph H(G);
  CompleteNeighborhoodStructure aux(aux_G);

  std::vector<int> S(get_graph_size(H), -1);
  for (int i = 0; i < n0; i++) {
    S[i] = i;
  }
  long double p_sigma = 0.0L, pv;
  while (get_graph_size(H) > n0) {
    Vertex v;
    std::tie(v, pv) = sample_vertex(H, n0, params, aux, algorithm, generator);
    S[get_graph_size(H) - 1] = get_index(H, v), p_sigma += pv;
    assert(aux.verify(H));
    aux.remove_vertex(get_neighbors(H, v)), delete_vertex(H, v);
    assert(aux.verify(H));
  }
  return std::make_pair(encode_permutation(S), p_sigma);
}

std::map<mpz_class, long double> get_permutation_probabilities_sampling(
    const Graph &G, const int &n0, const Parameters &params, const SamplingMethod &algorithm,
    const int &tries) {
  std::map<mpz_class, long double> permutations;
  CompleteNeighborhoodStructure aux(G);
  #pragma omp parallel for
  for (int i = 0; i < tries; i++) {
    auto sigma_with_probability = get_permutation_sample(G, n0, params, aux, algorithm);
    auto permutation = permutations.find(sigma_with_probability.first);
    if (permutation != permutations.end()) {
      permutation->second = add_exp_log(permutation->second, sigma_with_probability.second);
    } else {
      permutations.insert(sigma_with_probability);
    }
    if ((i + 1) % 10000 == 0) {
      #pragma omp critical
      {
        std::cerr << "Finished tries " << i + 1 << "/" << tries << std::endl;
      }
    }
  }
  auto max_log_value = std::max_element(
      permutations.begin(), permutations.end(),
      [] (const std::map<mpz_class, long double>::value_type &first,
          const std::map<mpz_class, long double>::value_type &second) {
          return first.second < second.second;
      })->second;
  for (auto &permutation : permutations) {
    permutation.second = exp2l(permutation.second - max_log_value);
  }
  auto total_probability = std::accumulate(
      permutations.begin(), permutations.end(), 0.0L,
      [] (long double &value, const std::map<mpz_class, long double>::value_type &permutation) {
          return value + permutation.second;
      });
  for (auto &permutation : permutations) {
    permutation.second /= total_probability;
  }
  return permutations;
}

std::map<VertexPair, long double> get_p_uv_from_permutations(
    const std::map<mpz_class, long double> &permutations, const int &n, const int &n0) {
  std::map<VertexPair, long double> p_uv;
  for (auto &permutation : permutations) {
    std::vector<int> S = decode_permutation(permutation.first, n);
    for (int i = n0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
          p_uv[std::make_pair(S[i], S[j])] += permutation.second;
      }
    }
  }
  return p_uv;
}
