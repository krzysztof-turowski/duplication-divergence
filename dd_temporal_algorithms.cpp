// Algorithms computating the temporal order for various duplication-divergence models.
// Compile: g++ dd_temporal_algorithms.cpp -O3 -o ./dd_temporal_algorithms
// Run: ./dd_temporal_algorithms synthetic all MODE n n0 PARAMETERS
//   or ./dd_temporal_algorithms synthetic ALGORITHM MODE n n0 PARAMETERS

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

#include <vector>

using namespace std;

typedef Koala::Graph<int, int> Graph;
typedef Koala::Graph<int, int>::PVertex Vertex;

const int G_TRIES = 10000;

enum TemporalAlgorithm { DEGREE_SORT };

const std::map<TemporalAlgorithm, std::string> SHORT_ALGORITHM_NAME = {
  { TemporalAlgorithm::DEGREE_SORT, "sort_by_degree" },
};

const std::map<TemporalAlgorithm, std::string> LONG_ALGORITHM_NAME = {
  { TemporalAlgorithm::DEGREE_SORT, "Sort vertices by degree" },
};

const std::map<std::string, TemporalAlgorithm> REVERSE_ALGORITHM_NAME = {
  { "sort_by_degree", TemporalAlgorithm::DEGREE_SORT },
};

void apply_permutation(Graph &G, const vector<int> &S) {
  vector<Vertex> V(G.getVertNo());
  G.getVerts(V.begin());
  for (auto v : V) {
    v->setInfo(S[v->getInfo()]);
  }
}

vector<set<int>> sort_by_degree(const Graph &G, const int &n0) {
  vector<set<int>> out(G.Delta() + 1);
  for (Vertex v = G.getVert(); v; v = G.getVertNext(v)) {
    if (v->getInfo() < n0) {
      out[G.Delta() - G.deg(v)].insert(v->getInfo());
    }
  }
  return out;
}

pair<double, double> score(const vector<set<int>> &solution) {
  double total = 0, correct = 0, n = 0;
  for (int i = 0; i < static_cast<int>(solution.size()); i++) {
    n += solution[i].size();
    for (int j = i + 1; j < static_cast<int>(solution.size()); j++) {
      total += solution[i].size() * solution[j].size();
      for (const auto &u : solution[i]) {
        for (const auto &v : solution[j]) {
          if (u < v) {
            correct++;
          }
        }
      }
    }
  }
  return make_pair(total / (n * (n - 1) / 2), correct / total);
}

pair<double, double> temporal_algorithm_single(
    const Graph &G0, const int &n, const Parameters &params, const TemporalAlgorithm &algorithm) {
  Graph G(G0);
  generate_graph_koala(G, n, params);

  switch (algorithm) {
    case DEGREE_SORT:
      return score(sort_by_degree(G, G0.getVertNo()));
    default:
      throw invalid_argument("Invalid algorithm: " + LONG_ALGORITHM_NAME.find(algorithm)->second);
  }
}

void print(
    const vector<pair<double, double>> solution, const int &n, const int &n0,
    const Parameters &params, const TemporalAlgorithm &algorithm, ostream &out_file) {
  cout << "Graph - n: " << n << ", n0: " << n0
      << ", parameters: " << params.to_string() << endl;
  cout << "Algorithm: " << LONG_ALGORITHM_NAME.find(algorithm)->second << endl;
  // TODO(unknown): Add printing mean
  for (const auto &density_precision : solution) {
    cout << "density: " << fixed << setw(6) << setprecision(3) << density_precision.first
        << " precision: " << fixed << setw(6) << setprecision(3) << density_precision.second
        << endl;
  }
  out_file << SHORT_ALGORITHM_NAME.find(algorithm)->second << " ";
  for (const auto &density_precision : solution) {
    out_file << density_precision.first << "," << density_precision.second << " ";
  }
  out_file << endl;
}

void temporal_algorithm(
    const int &n, const int &n0, const Parameters &params, const TemporalAlgorithm &algorithm,
    ostream &out_file) {
  Graph G0 = generate_seed_koala(n0, 1.0);
  // TODO(kturowski): parallelize
  vector<pair<double, double>> solution(G_TRIES);
  for (int i = 0; i < G_TRIES; i++) {
    solution[i] = temporal_algorithm_single(G0, n, params, algorithm);
  }
  print(solution, n, n0, params, algorithm, out_file);
}

int main(int, char *argv[]) {
  try {
    string action(argv[1]), algorithm(argv[2]), mode(argv[3]);
    int n = stoi(argv[4]), n0 = stoi(argv[5]);
    Parameters params;
    params.initialize(mode, argv + 6);
    string name(TEMP_FOLDER + get_synthetic_filename(n, n0, params, "TA"));
    if (algorithm == "all") {
      ofstream out_file(name, ios_base::app);
      for (const auto &alg : REVERSE_ALGORITHM_NAME) {
        temporal_algorithm(n, n0, params, alg.second, out_file);
      }
    } else if (REVERSE_ALGORITHM_NAME.count(algorithm)) {
      ofstream out_file(name, ios_base::app);
      temporal_algorithm(
          n, n0, params, REVERSE_ALGORITHM_NAME.find(algorithm)->second, out_file);
    } else {
      throw invalid_argument("Invalid algorithm: " + algorithm);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }

  return 0;
}
