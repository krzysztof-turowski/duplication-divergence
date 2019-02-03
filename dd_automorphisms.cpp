#include "dd_header.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "dd_automorphisms.h"
#pragma GCC diagnostic pop

#include <chrono>

const int PVAL_TRIES = 100;

enum Algorithm { NAUTY_DENSE, NAUTY_SPARSE, ALL };

double log_automorphisms(const std::vector<std::set<int>> &G, const Algorithm &algorithm = Algorithm::NAUTY_DENSE) {
  switch(algorithm) {
    default:
      break;
    case NAUTY_DENSE:
      return log_automorphisms_dense(G);
    case NAUTY_SPARSE:
      return log_automorphisms_sparse(G);
    case ALL:
      auto first = std::chrono::system_clock::now();
      double dense = log_automorphisms_dense(G);
      auto second = std::chrono::system_clock::now();
      double sparse = log_automorphisms_sparse(G);
      auto third = std::chrono::system_clock::now();
      double run_dense = std::chrono::duration_cast<std::chrono::milliseconds>(second - first).count() / 1000.0;
      double run_sparse = std::chrono::duration_cast<std::chrono::milliseconds>(third - second).count() / 1000.0;
      assert(sparse == dense);
      std::cerr << "Time: " << run_sparse << " " << run_dense << std::endl;
      return dense;
  }
  assert(0);
}

std::vector<double> log_automorphisms(const std::vector<std::set<int>> &G0, const int &n, const Parameters &params, const int &tries) {
  // TODO: provide parallelized version
  std::vector<double> log_aut_H(tries);
  for (int i = 0; i < tries; i++) {
    std::cerr << "Run " << i << ": ";
    std::vector<std::set<int>> H = G0;
    generate_graph(H, n, params);
    log_aut_H[i] = log_automorphisms(H);
    std::cerr << log_aut_H[i] << std::endl;
  }
  return log_aut_H;
}

void log_automorphisms_p_value(const std::string &graph_name, const std::string &seed_name, const Parameters &params) {
  std::vector<std::set<int>> G = read_graph(FILES_FOLDER + graph_name);
  std::vector<std::set<int>> G0 = read_graph(FILES_FOLDER + seed_name);
  double log_aut_G = log_automorphisms(G);
  std::vector<double> log_aut_H = log_automorphisms(G0, G.size(), params, PVAL_TRIES);
  double p_lower = std::count_if(log_aut_H.begin(), log_aut_H.end(), [&](const double &value){ return value < log_aut_G; }) / PVAL_TRIES;
  double p_upper = std::count_if(log_aut_H.begin(), log_aut_H.end(), [&](const double &value){ return value > log_aut_G; }) / PVAL_TRIES;
  double p_value = 2 * std::max(p_lower, p_upper);
  std::cout << graph_name << " " << std::accumulate(log_aut_H.begin(), log_aut_H.end(), 0.0) / PVAL_TRIES << " " << p_value << std::endl;
}

void log_automorphisms(const std::string &graph_name) {
  std::vector<std::set<int>> G = read_graph(FILES_FOLDER + graph_name);
  std::cout << graph_name << " " << log_automorphisms(G) << std::endl;
}

int main(int argc, char *argv[]) {
  try {
    std::string action(argv[1]), mode(argv[2]);
    if (action == "real_seed") {
      // TODO: make datasets parameter-variable
      // TODO: make parameters ranges
      Parameters params;
      params.initialize(mode, argv + 3);
      log_automorphisms_p_value("G-100-20-PS-0.70-2.00.txt", "G0-100-20-PS-0.70-2.00.txt", params);
      log_automorphisms_p_value("G-c-elegans.txt", "G0-c-elegans.txt", params);
      log_automorphisms_p_value("G-a-thaliana.txt", "G0-a-thaliana.txt", params);
      log_automorphisms_p_value("G-c-elegans.txt", "G0-c-elegans.txt", params);
      log_automorphisms_p_value("G-d-melanogaster.txt", "G0-d-melanogaster.txt", params);
      log_automorphisms_p_value("G-homo-sapiens.txt", "G0-homo-sapiens.txt", params);
      log_automorphisms_p_value("G-mus-musculus.txt", "G0-mus-musculus.txt", params);
      log_automorphisms_p_value("G-s-cerevisiae.txt", "G0-s-cerevisiae.txt", params);
      log_automorphisms_p_value("G-s-pombe.txt", "G0-s-pombe.txt", params);
    }
    else if (action == "real_graph") {
      std::string graph_name(argv[3]);
      log_automorphisms(graph_name);
    }
    else {
      assert(0);
    }
  } catch (std::exception &e) {
    std::cout << "ERROR: " << e.what() << std::endl;
  }
  return 0;
}
