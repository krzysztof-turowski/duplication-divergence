// Tool for graph generation for various duplication-divergence models.
// Compile: g++ dd_generate.cpp -O3 -o ./dd_generate
// Example run: ./dd_generate -n:100 -n0:10 -mode:pastor_satorras -p:0.5 -r:2.0 -p0:0.6

#include "./dd_input.h"
#include "./dd_header.h"

#include <exception>

typedef std::vector<std::set<unsigned>> Graph;

inline std::string name(const int &n, const int &n0, const Parameters &params) {
  return std::to_string(n) + "-" + std::to_string(n0) + "-" + params.to_filename();
}

void export_graph(const std::string &name, const Graph &G) {
  std::ofstream G_out_file(name);
  for (std::size_t i = 0; i < G.size(); i++) {
    G_out_file << i << " " << i << std::endl;
    for (auto j : G[i]) {
      if (i < j) {
        G_out_file << i << " " << j << std::endl;
      }
    }
  }
  G_out_file.close();
}

void generate_graph(const int &n, const int &n0, const double &p0, const Parameters &params) {
  Graph G = generate_seed_simple(n0, p0);
  export_graph(FILES_FOLDER + "G0-" + name(n, n0, params) + ".txt", G);
  generate_graph_simple(G, n, params);
  export_graph(FILES_FOLDER + "G-" + name(n, n0, params) + ".txt", G);
}

int main(int argc, char **argv) {
  try {
    Env = prepare_environment(argc, argv);
    const int n = read_n(Env), n0 = read_n0(Env);
    const double p0 = read_p0(Env);
    const Parameters params = read_parameters(Env);
    generate_graph(n, n0, p0, params);
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
  }
  return 0;
}
