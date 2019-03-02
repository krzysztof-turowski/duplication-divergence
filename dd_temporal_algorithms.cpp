// Algorithms computating the temporal order for various duplication-divergence models.
// Compile: g++ dd_temporal_algorithms.cpp -O3 -o ./dd_temporal_algorithms
// Run: ./dd_temporal_algorithms synthetic all MODE n n0 PARAMETERS
//   or ./dd_temporal_algorithms synthetic ALGORITHM MODE n n0 PARAMETERS

#include "./dd_graph.h"

#include <vector>

using namespace std;

typedef pair<double, double> DensityPrecision;

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
  vector<Vertex> V(get_vertices(G));
  for (auto &v : V) {
    set_index(G, v, S[get_index(G, v)]);
  }
}

vector<set<int>> sort_by_degree(const Graph &G, const int &n0) {
  int n = get_graph_size(G);
  vector<set<int>> out(n);
  for (const auto &v : get_vertices(G)) {
    int index = get_index(G, v);
    if (index >= n0) {
      out[n - 1 - get_degree(G, v)].insert(index);
    }
  }
  return out;
}

DensityPrecision get_density_precision(const vector<set<int>> &solution) {
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

DensityPrecision temporal_algorithm_single(
    const Graph &G0, const int &n, const Parameters &params, const TemporalAlgorithm &algorithm) {
  Graph G(G0);
  generate_graph(G, n, params);

  switch (algorithm) {
    case DEGREE_SORT: {
      return get_density_precision(sort_by_degree(G, G0.getVertNo()));
    }
    default:
      throw invalid_argument("Invalid algorithm: " + LONG_ALGORITHM_NAME.find(algorithm)->second);
  }
}

void print(
    const vector<DensityPrecision> solution, const int &n, const int &n0,
    const Parameters &params, const TemporalAlgorithm &algorithm, ostream &out_file,
    bool verbose = false) {
  cout << "Graph - n: " << n << ", n0: " << n0
      << ", parameters: " << params.to_string() << endl;
  cout << "Algorithm: " << LONG_ALGORITHM_NAME.find(algorithm)->second << endl;

  auto mean_density_precision =
      accumulate(
          solution.begin(), solution.end(), DensityPrecision(0.0, 0.0),
          [] (DensityPrecision &value, const DensityPrecision &density_precision) {
              value.first += density_precision.first, value.second += density_precision.second;
              return value;
          });
  cout << "Mean density: " << fixed << setw(6) << setprecision(3)
      << mean_density_precision.first / solution.size()
      << " mean precision: " << fixed << setw(6) << setprecision(3)
      << mean_density_precision.second / solution.size()
      << endl;
  if (verbose) {
    for (const auto &density_precision : solution) {
      cout << "density: " << fixed << setw(6) << setprecision(3) << density_precision.first
          << " precision: " << fixed << setw(6) << setprecision(3) << density_precision.second
          << endl;
    }
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
  Graph G0 = generate_seed(n0, 1.0);
  vector<DensityPrecision> solution(G_TRIES);
  #pragma omp parallel for
  for (int i = 0; i < G_TRIES; i++) {
    solution[i] = temporal_algorithm_single(G0, n, params, algorithm);
    #pragma omp critical
    {
      if ((i + 1) % 1000 == 0) {
        cerr << "Finished run " << i + 1 << "/" << G_TRIES << endl;
      }
    }
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
