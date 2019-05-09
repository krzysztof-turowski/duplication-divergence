/*
Tool for computation the Maximum Likelihood Estimator for various duplication-divergence models.
This is an implementation and our extension of the idea proposed in the following paper.
Wiuf, Carsten, Markus Brameier, Oskar Hagberg, and Michael PH Stumpf. "A likelihood approach to analysis of network data." Proceedings of the National Academy of Sciences 103, no. 20 (2006): 7566-7570.

Compile: make dd_ml_estimation

Syntax: ./dd_ml_estimation <options>
<options> are
-action:
   synthetic: synthetic random graph generations from DD-model
   real_data: real world graph
-graph: If action is `real_data`, give graph file name and file should be in edge list format.
-st: Number of independent tries for estimating the log-likelihood function. We use the idea proposed by Wiuf et al, and hence these tries correspond to importance sampling tries.
-mode: {pure_duplication, pastor_satorras, chung_lu}. In case of `synthetic` action, the mode (type) of the DD-graph model.
<parameters>: Depending on `mode`, the parameters `n,p,q,r`of the DD model.
-n0, -p0: These are the parameters for generating seed graph in the case of `synthetic` action.

Example runs:
 ./dd_ml_estimation -action:synthetic -n:100 -n0:10
     -mode:pastor_satorras -p:0.5 -r:2.0 -p0:0.6 -st:1000
 ./dd_ml_estimation -action:real_data -graph:G-s-cerevisiae.txt
     -mode:pastor_satorras -st:1000
 */

#include "./dd_input.h"
#include "./dd_graph.h"

const double STEP_P = 0.1, STEP_R = 1.0;
int IS_TRIES;

class LikelihoodValue {
 public:
  Parameters params;
  long double likelihood;

  LikelihoodValue(const Parameters& params_v, long double const &likelihood_v)
      : params(params_v), likelihood(likelihood_v) { }
};

long double likelihood_value(
    const Graph &G, const int &n0, const Parameters &params, const Parameters &params_0) {
  std::random_device device;
  std::mt19937 generator(device());
  Graph H(G);
  CompleteNeighborhoodStructure aux(H);

  long double ML_value = 0;
  while (get_graph_size(H) > n0) {
    std::vector<Vertex> V = get_vertices(H);
    std::vector<long double> P = get_transition_probability(H, params_0, aux);
    long double P_sum = std::accumulate(P.begin(), P.end(), 0.0L);
    if (P_sum == 0.0L) {
      return -std::numeric_limits<long double>::infinity();
    }
    std::discrete_distribution<int> choose_vertex(P.begin(), P.end());
    int index = choose_vertex(generator);
    long double log_p_transition = get_log_transition_probability(H, params, V[index], aux);
    if (!std::isfinite(log_p_transition)) {
      return -std::numeric_limits<long double>::infinity();
    }
    ML_value += log2l(P_sum) + log_p_transition - log2l(P[index]);

    aux.remove_vertex(get_neighbors(H, V[index]));
    delete_vertex(H, V[index]);
    aux.verify(H);
  }
  return ML_value;
}

LikelihoodValue importance_sampling(
    const Graph &G, const int &n0, const Parameters &params, const Parameters &params_0) {
  std::vector<long double> likelihood_values(IS_TRIES);
  #pragma omp parallel for
  for (int i = 0; i < IS_TRIES; i++) {
    likelihood_values[i] = likelihood_value(G, n0, params, params_0);
    #pragma omp critical
    {
      std::cerr << "Run " << i + 1 << "/" << IS_TRIES << std::endl;
    }
  }
  long double likelihood_score = -std::numeric_limits<long double>::infinity();
  for (const auto &value : likelihood_values) {
    likelihood_score = add_exp_log(likelihood_score, value);
  }
  return LikelihoodValue(params, likelihood_score - log2l(IS_TRIES));
}

// TODO(unknown): split into gradient search and sampling likelihood over some parameters interval
std::vector<LikelihoodValue> find_likelihood_values(
    const Graph &G, const int &n0, const Mode &mode) {
  Parameters params_0;
  std::vector<LikelihoodValue> likelihood_values;
  switch (mode) {
    case PASTOR_SATORRAS: {
        // TODO(unknown): parametrize parameter ranges and driving values
        const double MAX_R = n0 - 1;
        params_0.initialize_pastor_satorras(0.5, MAX_R / 2);
        for (double p = 0; p <= 1.0 + EPS; p += STEP_P) {
          for (double r = 0; r <= MAX_R + EPS; r += STEP_R) {
            Parameters params;
            params.initialize_pastor_satorras(p, r);
            LikelihoodValue value = importance_sampling(G, n0, params, params_0);
            std::cerr << "Finished run for " << value.params.to_string() << std::endl;
            likelihood_values.push_back(value);
          }
        }
      }
      return likelihood_values;
    default:
      throw std::invalid_argument("Invalid mode: " + LONG_NAME.find(mode)->second);
  }
}

void print(
    const std::string &name, const std::vector<LikelihoodValue> &likelihood_values,
    std::ostream &out_file) {
  std::cout << name << std::endl;
  auto ML = *max_element(
      likelihood_values.begin(), likelihood_values.end(),
      [&] (const LikelihoodValue& lhs, const LikelihoodValue& rhs) {
          return lhs.likelihood < rhs.likelihood;
      });
  std::cout << "Best found: " << ML.params.to_string()
      << " ML score: " << ML.likelihood << std::endl;
  for (auto const& value : likelihood_values) {
    std::cout << value.params.to_string() << " ML score: " << value.likelihood << std::endl;
  }
  for (auto const& value : likelihood_values) {
    out_file << value.params.to_csv() << "," << value.likelihood << " ";
  }
}

void synthetic_data(const int &n, const int &n0, const double &p0, const Parameters &params) {
  Graph G = generate_seed(n0, p0);
  generate_graph(G, n, params);
  auto likelihood_values = find_likelihood_values(G, n0, params.mode);
  std::ofstream out_file(TEMP_FOLDER + get_synthetic_filename(n, n0, params, "ML"));
  print("Synthetic data: " + params.to_string(), likelihood_values, out_file);
}

void real_world_data(
    const std::string &graph_name, const std::string &seed_name, const Mode &mode) {
  Graph G = read_graph(FILES_FOLDER + graph_name);
  int n0 = read_graph_size(FILES_FOLDER + seed_name);
  auto likelihood_values = find_likelihood_values(G, n0, mode);
  std::ofstream out_file(TEMP_FOLDER + get_real_filename(graph_name, mode, "ML"));
  print(graph_name, likelihood_values, out_file);
}

int main(int argc, char **argv) {
  try {
    Env = prepare_environment(argc, argv);
    IS_TRIES = read_int(Env, "-st:", 1, "Important sampling tries");
    std::string action = read_action(Env);
    if (action == "synthetic") {
      const int n = read_n(Env), n0 = read_n0(Env);
      const double p0 = read_p0(Env);
      Parameters params = read_parameters(Env);
      synthetic_data(n, n0, p0, params);
    } else if (action == "real_data") {
      std::string graph_name = read_graph_name(Env);
      std::string mode = read_mode(Env);
      real_world_data(graph_name, get_seed_name(graph_name), REVERSE_NAME.find(mode)->second);
    } else {
      throw std::invalid_argument("Invalid action: " + action);
    }
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
  }
  return 0;
}
