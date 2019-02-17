// Tools for computation the temporal order for various duplication-divergence models.
// Compile: g++ dd_temporal_order.cpp -O3 -o ./dd_temporal_order
// Run: ./dd_temporal_order exact_bound MODE n n0 PARAMETERS - e.g. ./dd_temporal_order exact_bound pastor_satorras 100 20 0.5 2.0

#include "dd_koala.h"

#include <random>

#include <gmpxx.h>

using namespace std;

typedef Koala::Graph<int, int> Graph;
typedef Koala::Graph<int, int>::PVertex Vertex;

vector<int> generate_permutation(const int &n, const int &n0) {
  random_device device;
  mt19937 generator(device());
  vector<int> S(n);
  for (int i = 0; i < n; i++) {
    S[i] = i;
  }
  for (int i = n - 1; i > n0; i--) {
    std::uniform_int_distribution<int> swap_distribution(n0, i);
    int index = swap_distribution(generator);
    swap(S[i], S[index]);
  }
  return S;
}

vector<int> invert_permutation(const vector<int> &S) {
  vector<int> V(S.size());
  for (int i = 0; i < static_cast<int>(S.size()); i++) {
    V[S[i]] = i;
  }
  return V;
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

mpz_class encode_permutation(const vector<int> &S, const int &n) {
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

void LP_bound_exact(const int &n, const int &n0, const Parameters &params) {
  Graph G0 = generate_seed_koala(n0, 1.0);
  Graph G(G0);
  generate_graph_koala(G, n, params);

  Graph H(G);
  vector<int> S = generate_permutation(n, n0);
  apply_permutation(H, S);
  // compute P(Pi_n = sigma_n | Pi_n(G_n) = h_n, G_n0 = g_n0)
  // compute p_uv coefficients
  // construct and solve LP
  // export to file and to stdout
}

int main(int argc, char *argv[]) {
  string action(argv[1]), mode(argv[2]);
  int n = stoi(argv[3]), n0 = stoi(argv[4]);
  Parameters params;
  params.initialize(mode, argv + 5);
  if (action == "exact_bound") {
    LP_bound_exact(n, n0, params);
  }
  else {
    assert(0);
  }

  return 0;
}
