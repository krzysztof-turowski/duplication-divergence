/*
Generate information about the basic parameters of the graph: mean degree, number of open triangles, number of triangles for various duplication-divergence models.

Compile: make dd_graphlets

Syntax: ./dd_graphlets <options>
<options> are
-action:
  synthetic: synthetic random graph generations from duplication-divergence model
  real_seed: synthetic graphs with given parameters and given seed file.
-graph: If action is `real_seed`, then provide graph file name (located in `files/` folder). File should be in edge list format.
-gt: Number of independent graphs generated from a given seed for a given parameters.
-mode: {pure_duplication, pastor_satorras, chung_lu}. In case of `synthetic` action, the mode (type) of the duplication-divergence graph model.
<parameters>: Depending on `mode`, the parameters `p,q,r` of the DD model.
-n: The size of a graph in the case of `synthetic` action.
-n0, -p0: The parameters for generating a seed graph in the case of `synthetic` action.

Example run:
  ./dd_graphlets -action:real_seed -graph:G-a-thaliana.txt -gt:100 -mode:pastor_satorras -p:0.98 -r:0.49
*/ 

#include "./dd_input.h"
#include "./dd_graphlets.h"

typedef std::tuple<double, double, double> PValuesInfo;

int TRIES;

double get_p_value(
    const DataObject &graphlets_G, const std::vector<DataObject> &graphlets_H,
    std::function<double(const DataObject&)> get_value) {
  double graphlet_value = get_value(graphlets_G);
  double p_lower =
      std::count_if(
          graphlets_H.begin(), graphlets_H.end(),
          [&](const DataObject &info){ return get_value(info) < graphlet_value; });
  double p_upper =
      std::count_if(
          graphlets_H.begin(), graphlets_H.end(),
          [&](const DataObject &info){ return get_value(info) > graphlet_value; });
  return 2 * std::min(p_lower, p_upper) / graphlets_H.size();
}

PValuesInfo get_p_value(
    const DataObject &graphlets_G, const std::vector<DataObject> &graphlets_H) {
  auto get_average_degree = [](const DataObject &info) { return info.average_degree; };
  auto get_open_triangles = [](const DataObject &info) { return info.open_triangles; };
  auto get_triangles = [](const DataObject &info) { return info.triangles; };

  return PValuesInfo(
      get_p_value(graphlets_G, graphlets_H, get_average_degree),
      get_p_value(graphlets_G, graphlets_H, get_open_triangles),
      get_p_value(graphlets_G, graphlets_H, get_triangles));
}

void print(const std::string &name, const std::vector<DataObject> &graphlets) {
  std::cout << std::fixed << std::setprecision(3) << name << std::endl;
  for (const auto &v : graphlets) {
    std::cout << v.to_string() << std::endl;
  }
}

void print(const PValuesInfo &pvalues) {
  std::cout << "p-value for average degree: " << std::fixed << std::setw(8) << std::setprecision(3)
      << std::get<0>(pvalues) << std::endl;
  std::cout << "p-value for open triangles: " << std::fixed << std::setw(8) << std::setprecision(3)
      << std::get<1>(pvalues) << std::endl;
  std::cout << "p-value for triangles: " << std::fixed << std::setw(8) << std::setprecision(3)
      << std::get<2>(pvalues) << std::endl;
}

void synthetic_data(const int &n, const int &n0, const double &p0, const Parameters &params) {
  Graph G0 = generate_seed_simple(n0, p0);
  std::vector<DataObject> graphlets(TRIES);
  #pragma omp parallel for
  for (int i = 0; i < TRIES; i++) {
    graphlets[i] = get_params_for_synthetic_graph(G0, n, params);
    #pragma omp critical
    {
      std::cerr << "Run " << i + 1 << "/" << TRIES << std::endl;
    }
  }
  print("Synthetic data: " + params.to_string(), graphlets);
}

void real_seed_data(
    const std::string &graph_name, const std::string &seed_name, const Parameters &params) {
  Graph G = read_graph_simple(FILES_FOLDER + graph_name);
  Graph G0 = read_graph_simple(FILES_FOLDER + seed_name);
  DataObject graphlets_G = get_params_for_graph(G);
  std::vector<DataObject> graphlets_H(TRIES);
  #pragma omp parallel for
  for (int i = 0; i < TRIES; i++) {
    graphlets_H[i] = get_params_for_synthetic_graph(G0, G.size(), params);
    #pragma omp critical
    {
      std::cerr << "Run " << i + 1 << "/" << TRIES << std::endl;
    }
  }
  print("File: " + graph_name, std::vector<DataObject>{graphlets_G});
  print("", graphlets_H);
  PValuesInfo pvalues = get_p_value(graphlets_G, graphlets_H);
  print(pvalues);
}

int main(int argc, char **argv) {
  try {
    Env = prepare_environment(argc, argv);
    TRIES = read_int(Env, "-gt:", 1, "Number of tries");
    std::string action = read_action(Env);
    if (action == "synthetic") {
      const int n = read_n(Env), n0 = read_n0(Env);
      const double p0 = read_p0(Env);
      std::unique_ptr<Parameters> params = read_parameters(Env);
      synthetic_data(n, n0, p0, *params);
    } else if (action == "real_seed") {
      std::string graph_name = read_graph_name(Env);
      std::unique_ptr<Parameters> params = read_parameters(Env);
      real_seed_data(graph_name, get_seed_name(graph_name), *params);
    } else {
      throw std::invalid_argument("Invalid action: " + action);
    }
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
  }
  return 0;
}
