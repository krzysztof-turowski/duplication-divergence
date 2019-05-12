/*
Tool for computation automorphisms for various duplication-divergence models.

Compile: make dd_automorphisms

Syntax: ./dd_automorphisms <options>
<options> are
-action:
  real_graph: Find the logarithm of the number of symmetries in the given real_graph file
  real_seed: Find the expected logarithm of the number of symmetries of a synthetic graph with given parameters and given seed file.
  synthetic: Find the average logarithm of the number of symmetries of a synthetic graph with given size and given parameters.
-graph: Graph file name and file should be in edge list format. In case of `real_seed` option, seed graph name is extracted from the graph file name
-mode: {pure_duplication, pastor_satorras, chung_lu}. In case of `real_seed` or `synthetic` action, the mode (type) of the duplication-divergence graph model.
-st: When action is `real_seed` or `synthetic`, number of independent tries to calculate empirical average and/or p-value.
<parameters>: Depending on `mode`, the parameters `p,q,r` of the duplication-divergence model.
-n: The size of a graph in the case of `synthetic` action.
-n0, -p0: The parameters for generating a seed graph in the case of `synthetic` action.

Example runs:
  ./dd_automorphisms -action:real_graph -graph:G-s-cerevisiae.txt
  ./dd_automorphisms -action:real_seed -graph:G-s-cerevisiae.txt -mode:pastor_satorras -p:0.5 -r:2.0 -p0:0.6 -st:1000
  ./dd_automorphisms -action:synthetic -n:200 -n0:20 -st:100 -mode:pastor_satorras
*/

#include "./dd_input.h"
#include "./dd_header.h"
#include "./dd_automorphisms.h"

#include <algorithm>
#include <chrono>
#include <string>

typedef std::tuple<double, double, double, double> AutomorphismsInfo;
typedef std::tuple<double, double, double, double> PValuesInfo;

int PVAL_TRIES, AVG_TRIES;

enum AutomorphismsDetection { NAUTY_DENSE, NAUTY_SPARSE, TRACES, ALL };

const AutomorphismsDetection ALGORITHM = AutomorphismsDetection::NAUTY_DENSE;

double log_automorphisms(const Graph &G) {
  switch (ALGORITHM) {
    case NAUTY_DENSE:
      return log_automorphisms_dense(G);
    case NAUTY_SPARSE:
      return log_automorphisms_sparse(G);
    case TRACES:
      return log_automorphisms_traces(G);
    case ALL: {
        auto first = std::chrono::system_clock::now();
        double dense = log_automorphisms_dense(G);
        auto second = std::chrono::system_clock::now();
        double sparse = log_automorphisms_sparse(G);
        auto third = std::chrono::system_clock::now();
        double traces = log_automorphisms_traces(G);
        auto fourth = std::chrono::system_clock::now();
        double run_dense =
            std::chrono::duration_cast<std::chrono::milliseconds>(second - first).count() / 1000.0;
        double run_sparse =
            std::chrono::duration_cast<std::chrono::milliseconds>(third - second).count() / 1000.0;
        double run_traces =
            std::chrono::duration_cast<std::chrono::milliseconds>(fourth - third).count() / 1000.0;
        if (dense != sparse || dense != traces) {
          throw std::logic_error(
              "Automorphism number do not match: dense " + std::to_string(dense)
              + ", sparse " + std::to_string(sparse) + ", traces " + std::to_string(traces));
        }
        std::cerr << "Time: " << run_dense << " " << run_sparse << " " << run_traces << std::endl;
        return dense;
      }
    default:
      throw std::invalid_argument("Invalid automorphisms detection algorithm");
  }
}

double log_automorphisms_from_isolated_nodes(const Graph &G) {
  return lgamma(isolated_nodes(G) + 1);
}

double log_automorphisms_from_cherries(const Graph &G) {
  std::vector<int> V(G.size());
  for (size_t i = 0; i < G.size(); i++) {
    if (G[i].size() == 1) {
      V[*(G[i].begin())] += 1;
    }
  }
  double out =
      std::accumulate(
          V.begin(), V.end(), 0.0,
          [] (double &value, const double &element) { return value + lgamma(element + 1); });
  return out;
}

double log_automorphisms_from_copies(const Graph &G) {
  std::vector<bool> V(G.size(), false);
  double out = 0, count;
  for (size_t i = 0; i < G.size(); i++) {
    if (V[i]) {
      continue;
    }
    count = 1, V[i] = true;
    for (size_t j = i + 1; j < G.size(); j++) {
      if (G[i] == G[j]) {
        count++, V[j] = true;
      }
    }
    out += lgamma(count);
  }
  return out;
}

AutomorphismsInfo log_automorphisms_single(const Graph &G) {
  return AutomorphismsInfo(
      log_automorphisms(G), log_automorphisms_from_isolated_nodes(G),
      log_automorphisms_from_cherries(G), log_automorphisms_from_copies(G));
}

AutomorphismsInfo log_automorphisms_single(
    const Graph &G0, const int &n, const Parameters &params) {
  Graph H(G0);
  generate_graph_simple(H, n, params);
  return log_automorphisms_single(H);
}

std::vector<AutomorphismsInfo> log_automorphisms(
    const Graph &G0, const int &n, const Parameters &params, const int &tries) {
  std::vector<AutomorphismsInfo> log_aut_H(tries);
  #pragma omp parallel for num_threads(THREADS)
  for (int i = 0; i < tries; i++) {
    log_aut_H[i] = log_automorphisms_single(G0, n, params);
    # pragma omp critical
    {
      std::cerr << "Run " << i + 1 << "/" << tries << ": "
          << std::get<0>(log_aut_H[i]) << std::endl;
    }
  }
  return log_aut_H;
}

template <typename T>
void print(const std::string &graph_name, const T &info, const bool &all_types) {
  std::cout << std::left << std::setw(25) << graph_name << " " << std::right;
  std::cout << std::fixed << std::setw(8) << std::setprecision(3) << std::get<0>(info) << " ";
  if (all_types) {
    std::cout << std::fixed << std::setw(8) << std::setprecision(3) << std::get<1>(info) << " ";
    std::cout << std::fixed << std::setw(8) << std::setprecision(3) << std::get<2>(info) << " ";
    std::cout << std::fixed << std::setw(8) << std::setprecision(3) << std::get<3>(info) << " ";
  }
  std::cout << std::endl;
}

double get_average(
    const std::vector<AutomorphismsInfo> &log_aut_H,
    std::function<double(const AutomorphismsInfo&)> get_value) {
  return std::accumulate(
      log_aut_H.begin(), log_aut_H.end(), 0.0,
      [&](double &value, const AutomorphismsInfo &info){
          return value + get_value(info);
      }) / log_aut_H.size();
}

AutomorphismsInfo get_average(const std::vector<AutomorphismsInfo> &log_aut_H) {
  auto get_log_aut = [](const AutomorphismsInfo &info) { return std::get<0>(info); };
  auto get_log_aut_isolated = [](const AutomorphismsInfo &info) { return std::get<1>(info); };
  auto get_log_aut_cherries = [](const AutomorphismsInfo &info) { return std::get<2>(info); };
  auto get_log_aut_copies = [](const AutomorphismsInfo &info) { return std::get<3>(info); };

  return AutomorphismsInfo(
      get_average(log_aut_H, get_log_aut),
      get_average(log_aut_H, get_log_aut_isolated),
      get_average(log_aut_H, get_log_aut_cherries),
      get_average(log_aut_H, get_log_aut_copies));
}

double get_p_value(
    const AutomorphismsInfo log_aut_G, const std::vector<AutomorphismsInfo> &log_aut_H,
    std::function<double(const AutomorphismsInfo&)> get_value) {
  double log_aut_G_value = get_value(log_aut_G);
  double p_lower =
      std::count_if(
          log_aut_H.begin(), log_aut_H.end(),
          [&](const AutomorphismsInfo &info){ return get_value(info) < log_aut_G_value; });
  double p_upper =
      std::count_if(
          log_aut_H.begin(), log_aut_H.end(),
          [&](const AutomorphismsInfo &info){ return get_value(info) > log_aut_G_value; });
  return 2 * std::min(p_lower, p_upper) / log_aut_H.size();
}

PValuesInfo get_p_value(
    const AutomorphismsInfo log_aut_G, const std::vector<AutomorphismsInfo> &log_aut_H) {
  auto get_log_aut = [](const AutomorphismsInfo &info) { return std::get<0>(info); };
  auto get_log_aut_isolated = [](const AutomorphismsInfo &info) { return std::get<1>(info); };
  auto get_log_aut_cherries = [](const AutomorphismsInfo &info) { return std::get<2>(info); };
  auto get_log_aut_copies = [](const AutomorphismsInfo &info) { return std::get<3>(info); };

  return PValuesInfo(
      get_p_value(log_aut_G, log_aut_H, get_log_aut),
      get_p_value(log_aut_G, log_aut_H, get_log_aut_isolated),
      get_p_value(log_aut_G, log_aut_H, get_log_aut_cherries),
      get_p_value(log_aut_G, log_aut_H, get_log_aut_copies));
}

void log_automorphisms_p_value(
    const std::string &graph_name, const std::string &seed_name, const Parameters &params) {
  Graph G = read_graph_simple(FILES_FOLDER + graph_name);
  Graph G0 = read_graph_simple(FILES_FOLDER + seed_name);
  AutomorphismsInfo log_aut_G = log_automorphisms_single(G);
  std::vector<AutomorphismsInfo> log_aut_H = log_automorphisms(G0, G.size(), params, PVAL_TRIES);

  AutomorphismsInfo log_aut_avg_values = get_average(log_aut_H);
  print(graph_name, log_aut_avg_values, false);

  PValuesInfo p_values = get_p_value(log_aut_G, log_aut_H);
  print("", p_values, false);
}

void log_automorphisms(const std::string &graph_name) {
  Graph G(read_graph_simple(FILES_FOLDER + graph_name));
  print(graph_name, log_automorphisms_single(G), false);
}

void synthetic_data(const int &n, const int &n0, const double &p0, const Parameters &params) {
  Graph G0(generate_seed_simple(n0, p0));
  switch (mode) {
    case PASTOR_SATORRAS:
      std::vector<AutomorphismsInfo> log_aut_G = log_automorphisms(G0, n, params, AVG_TRIES);
      AutomorphismsInfo log_aut_avg_values = get_average(log_aut_G);
      print("Synthetic data: " + params.to_string(), log_aut_avg_values, false);
      break;
    default:
      throw std::invalid_argument("Invalid mode: " + LONG_NAME.find(mode)->second);
  }
}

int main(int argc, char **argv) {
  try {
    Env = prepare_environment(argc, argv);
    std::string action =
        read_string(Env, "-action:", "", "Action: real_graph, real_seed, synthetic");
    // TODO(unknown): test parameters ranges
    if (action == "synthetic") {
      AVG_TRIES = read_int(Env, "-st:", 1, "Number of tries");
      const int n = read_n(Env), n0 = read_n0(Env);
      const double p0 = read_p0(Env);
      Parameters params = read_parameters(Env);
      synthetic_data(n, n0, p0, params);
    } else if (action == "real_seed") {
      PVAL_TRIES = read_int(Env, "-st:", 1, "p-value tries");
      std::string graph_name = read_graph_name(Env);
      Parameters params = read_parameters(Env);
      log_automorphisms_p_value(graph_name, get_seed_name(graph_name), params);
    } else if (action == "real_graph") {
      std::string graph_name = read_graph_name(Env);
      log_automorphisms(graph_name);
    } else {
      throw std::invalid_argument("Invalid action: " + action);
    }
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
  }
  return 0;
}
