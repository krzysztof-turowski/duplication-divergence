// Tools for computation the temporal order for various duplication-divergence models.
// Compile: g++ dd_temporal_order.cpp -O3 -lgmpxx -lgmp -lglpk -o ./dd_temporal_order
// Run: ./dd_temporal_order exact_bound MODE n n0 PARAMETERS - e.g. ./dd_temporal_order exact_bound pastor_satorras 100 20 0.5 2.0

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-default"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "dd_koala.h"
#pragma GCC diagnostic pop

#include <random>

#include <glpk.h>
#include <gmpxx.h>

using namespace std;

typedef Koala::Graph<int, int> Graph;
typedef Koala::Graph<int, int>::PVertex Vertex;

const int TRIES = 100;
const double EPS_STEP = 0.05;

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
    Graph &G, const int &n0, const Parameters &params, NeighborhoodStructure &aux, vector<int> &S, const double &p) {
  map<mpz_class, double> permutations;
  if (G.getVertNo() == n0) {
    mpz_class sigma = encode_permutation(S);
    permutations.insert(make_pair(sigma, p));
    return permutations;
  }

  multimap<Vertex, Vertex> candidates;
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

      map<mpz_class, double> permutations_v = get_permutation_probabilities(G, n0, params, aux, S, p * p_v);
      permutations.insert(permutations_v.begin(), permutations_v.end());

      G.move(H, v), aux.restore_vertex(v, neighbors_v), S[G.getVertNo() - 1] = -1;
      for (auto &u : neighbors_v) {
        G.addEdge(v, u);
      }
    }
  }
  return permutations;
}

map<mpz_class, double> get_permutation_probabilities(const Graph &G, const int &n0, const Parameters &params) {
  Graph H(G);
  NeighborhoodStructure aux(H);
  vector<int> S(H.getVertNo(), -1);
  for (int i = 0; i < n0; i++) {
    S[i] = i;
  }
  map<mpz_class, double> permutations = get_permutation_probabilities(H, n0, params, aux, S, 1.0);

  double total_probability = accumulate(
    permutations.begin(), permutations.end(), 0.0,
    [] (double value, const map<mpz_class, double>::value_type &permutation) { return value + permutation.second; });
  for (auto &permutation : permutations) {
    permutation.second /= total_probability;
  }
  return permutations;
}

map<pair<int, int>, double> get_p_uv_from_permutations(const map<mpz_class, double> &permutations, const int &n, const int &n0) {
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

string LP_row_name(const string &prefix, const initializer_list<int> &vertices) {
  ostringstream out;
  out << prefix;
  if (vertices.size() > 0) {
    out << "_{";
    for (auto v : vertices) {
      out << to_string(v) << ",";
    }
    out.seekp(-1, ios_base::end), out << "}";
  }
  return out.str();
}

double LP_solve(const map<pair<int, int>, double> &p_uv, const int &n, const int &n0, const double &epsilon) {
  glp_prob *LP = glp_create_prob();
  glp_set_prob_name(LP, ("Solve " + to_string(epsilon)).c_str());
  glp_set_obj_dir(LP, GLP_MAX);
  double variables = (n - n0) * (n - n0 - 1), density = epsilon * variables / 2;

  auto get_variable_index = [&](const int &u, const int &v) {
    return (u - n0) * (n - n0) + (v - n0);
  };
  auto name_variable = [&](const int &u, const int &v) {
    return "x_{" + to_string(u) + "," + to_string(v) + "}";
  };

  // Objective function
  glp_add_cols(LP, (n - n0) * (n - n0));
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      auto index = get_variable_index(i, j);
      glp_set_col_name(LP, index + 1, name_variable(i, j).c_str());
      if (i != j) {
        glp_set_col_bnds(LP, index + 1, GLP_DB, 0.0, 1.0);
      }
      if (i < j) {
        glp_set_obj_coef(LP, index + 1, p_uv.find(make_pair(i, j))->second / density);
      }
    }
  }

  vector<int> X, Y;
  vector<double> A;
  int row = 1;
  glp_add_rows(LP, (n - n0) * (n - n0 - 1) * (n - n0 - 2) + (n - n0) * (n - n0 - 1) / 2 + 1);
  // Antisymmetry
  for (int i = n0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      glp_set_row_name(LP, row, LP_row_name("A", {i, j}).c_str());
      X.push_back(row), Y.push_back(get_variable_index(i, j) + 1), A.push_back(1);
      X.push_back(row), Y.push_back(get_variable_index(j, i) + 1), A.push_back(1);
      glp_set_row_bnds(LP, row, GLP_UP, 1.0, 1.0);
      row++;
    }
  }
  // Transitivity
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      for (int k = n0; k < n; k++) {
        if (i != j && j != k && i != k) {
          glp_set_row_name(LP, row, LP_row_name("T", {i, j, k}).c_str());
          X.push_back(row), Y.push_back(get_variable_index(i, j) + 1), A.push_back(1);
          X.push_back(row), Y.push_back(get_variable_index(j, k) + 1), A.push_back(1);
          X.push_back(row), Y.push_back(get_variable_index(i, k) + 1), A.push_back(-1);
          glp_set_row_bnds(LP, row, GLP_UP, 1.0, 1.0), row++;
        }
      }
    }
  }
  // Density
  glp_set_row_name(LP, row, LP_row_name("D", {}).c_str());
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      if (i != j) {
        X.push_back(row), Y.push_back(get_variable_index(i, j) + 1), A.push_back(1);
      }
    }
  }
  glp_set_row_bnds(LP, row, GLP_FX, density, density), row++;

  glp_load_matrix(LP, A.size(), &X[0] - 1, &Y[0] - 1, &A[0] - 1);
  // glp_write_lp(LP, NULL, "/dev/stdout");
  glp_term_out(0);
  glp_simplex(LP, NULL);

  double solution = glp_get_obj_val(LP);
  glp_delete_prob(LP);
  glp_free_env();
  return solution;
}

vector<double> LP_bound_exact_single(const Graph &G0, const int &n, const Parameters &params, const vector<double> &epsilon) {
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

void LP_bound_exact(const int &n, const int &n0, const Parameters &params) {
  Graph G0 = generate_seed_koala(n0, 1.0);
  vector<double> epsilon;
  for (double eps = EPS_STEP; eps <= 1.0 + 10e-9; eps += EPS_STEP) {
    epsilon.push_back(eps);
  }

  // TODO: parallelize
  vector<double> solution(epsilon.size(), 0.0);
  for (int i = 0; i < TRIES; i++) {
    vector<double> solution_single = LP_bound_exact_single(G0, n, params, epsilon);
    transform(solution.begin(), solution.end(), solution_single.begin(), solution.begin(), std::plus<double>());
  }
  
  for (int i = 0; i < static_cast<int>(solution.size()); i++) {
    cout << fixed << setw(6) << setprecision(3) << epsilon[i] << " "
        << fixed << setw(6) << setprecision(3) << solution[i] / TRIES << endl;
  }
  // export to file
}

int main(int, char *argv[]) {
  try {
    string action(argv[1]), mode(argv[2]);
    int n = stoi(argv[3]), n0 = stoi(argv[4]);
    Parameters params;
    params.initialize(mode, argv + 5);
    if (action == "exact_bound") {
      LP_bound_exact(n, n0, params);
    }
    else {
      throw invalid_argument("Invalid action: " + action);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }

  return 0;
}
