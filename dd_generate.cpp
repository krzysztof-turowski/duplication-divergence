// Tool for graph generation for various duplication-divergence models.
// Compile: g++ dd_generate.cpp -O3 -o ./dd_generate
// Run: ./dd_generate MODE n n0 PARAMETERS - e.g. ./dd_generate pastor_satorras 100 20 0.5 2.0

#include "./dd_header.h"

#include <exception>

using namespace std;

inline string name(const int &n, const int &n0, const Parameters &params) {
  return to_string(n) + "-" + to_string(n0) + "-" + params.to_filename();
}

void export_graph(const string &name, const vector<set<int>> &G) {
  ofstream G_out_file(name);
  for (size_t i = 0; i < G.size(); i++) {
    G_out_file << i << " " << i << endl;
    for (auto j : G[i]) {
      G_out_file << i << " " << j << endl;
    }
  }
  G_out_file.close();
}

void generate_graph(const int &n, const int &n0, const Parameters &params) {
  vector<set<int>> G = generate_seed_simple(n0, 1.0);
  export_graph(FILES_FOLDER + "G0-" + name(n, n0, params) + ".txt", G);
  generate_graph_simple(G, n, params);
  export_graph(FILES_FOLDER + "G-" + name(n, n0, params) + ".txt", G);
}

int main(int, char *argv[]) {
  try {
    string mode(argv[1]);
    int n = stoi(argv[2]), n0 = stoi(argv[3]);
    Parameters params;
    params.initialize(mode, argv + 4);
    generate_graph(n, n0, params);
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }
  return 0;
}
