// Tool for computation the Maximum Likelihood Estimator for various duplication-divergence models.
// Compile: g++ dd_ml_estimation.cpp -O3 -o ./dd_ml_estimation
// Run: ./dd_ml_estimation synthetic MODE n n0 PARAMETERS or ./dd_ml_estimation real_data FILE MODE

#include "./dd_graph.h"

using namespace std;

const double STEP_P = 0.1, STEP_R = 1.0;
const int IS_TRIES = 100;

class LikelihoodValue {
 public:
  Parameters params;
  long double likelihood;

  LikelihoodValue(const Parameters& params_v, long double const &likelihood_v)
      : params(params_v), likelihood(likelihood_v) { }
};

long double likelihood_value(
    const Graph &G, const int &n0, const Parameters &params, const Parameters &params_0) {
  random_device device;
  mt19937 generator(device());
  Graph H(G);
  CompleteNeighborhoodStructure aux(H);

  long double ML_value = 0;
  while (get_graph_size(H) > n0) {
    vector<Vertex> V = get_vertices(H);
    vector<long double> P = get_transition_probability(H, params_0, aux);
    long double P_sum = accumulate(P.begin(), P.end(), 0.0L);
    if (P_sum == 0.0L) {
      return 0.0L;
    }
    discrete_distribution<int> choose_vertex(P.begin(), P.end());
    int index = choose_vertex(generator);
    long double p_transition = get_transition_probability(H, params, V[index], aux);
    if (p_transition == 0.0L) {
      return 0.0L;
    }
    long double log_p_transition = log(p_transition);
    ML_value += log(P_sum) + log_p_transition - log(P[index]);

    aux.remove_vertex(get_neighbors(H, V[index]));
    delete_vertex(H, V[index]);
    aux.verify(H);
  }
  return exp(ML_value);
}

LikelihoodValue importance_sampling(
    const Graph &G, const int &n0, const Parameters &params, const Parameters &params_0) {
  vector<long double> likelihood_values(IS_TRIES);
  #pragma omp parallel for
  for (int i = 0; i < IS_TRIES; i++) {
    likelihood_values[i] = likelihood_value(G, n0, params, params_0);
    #pragma omp critical
    {
      std::cerr << "Run " << i + 1 << "/" << IS_TRIES << std::endl;
    }
  }
  long double likelihood_score =
      accumulate(likelihood_values.begin(), likelihood_values.end(), 0.0L) / IS_TRIES;
  return LikelihoodValue(params, likelihood_score);
}

// TODO(unknown): split into gradient search and sampling likelihood over some parameters interval
vector<LikelihoodValue> find_likelihood_values(const Graph &G, const int &n0, const Mode &mode) {
  Parameters params_0;
  vector<LikelihoodValue> likelihood_values;
  switch (mode) {
    case PASTOR_SATORRAS:
      {
        // TODO(unknown): parametrize parameter ranges and driving values
        const double MAX_R = n0 - 1;
        params_0.initialize_pastor_satorras(0.5, MAX_R / 2);
        for (double p = 0; p <= 1.0 + EPS; p += STEP_P) {
          for (double r = 0; r <= MAX_R + EPS; r += STEP_R) {
            Parameters params;
            params.initialize_pastor_satorras(p, r);
            LikelihoodValue value = importance_sampling(G, n0, params, params_0);
            cerr << "Finished run for " << value.params.to_string() << endl;
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
    const string &name, const vector<LikelihoodValue> &likelihood_values, ostream &out_file) {
  cout << name << endl;
  auto ML = *max_element(
      likelihood_values.begin(), likelihood_values.end(),
      [&] (const LikelihoodValue& lhs, const LikelihoodValue& rhs) {
          return lhs.likelihood < rhs.likelihood;
      });
  cout << "Best found: " << ML.params.to_string() << " ML score: " << ML.likelihood << endl;
  for (auto const& value : likelihood_values) {
    cout << value.params.to_string() << " ML score: " << value.likelihood << endl;
  }
  for (auto const& value : likelihood_values) {
    out_file << value.params.to_csv() << "," << log(value.likelihood) << " ";
  }
}

void synthetic_data(const int &n, const int &n0, const Parameters &params) {
  Graph G = generate_seed(n0, 1.0);
  generate_graph(G, n, params);
  auto likelihood_values = find_likelihood_values(G, n0, params.mode);
  ofstream out_file(TEMP_FOLDER + get_synthetic_filename(n, n0, params, "ML"));
  print("Synthetic data: " + params.to_string(), likelihood_values, out_file);
}

void real_world_data(const string &graph_name, const string &seed_name, const Mode &mode) {
  Graph G = read_graph(FILES_FOLDER + graph_name);
  int n0 = read_graph_size(FILES_FOLDER + seed_name);
  auto likelihood_values = find_likelihood_values(G, n0, mode);
  ofstream out_file(TEMP_FOLDER + get_real_filename(graph_name, mode, "ML"));
  print(graph_name, likelihood_values, out_file);
}

int main(int, char *argv[]) {
  try {
    string action(argv[1]);
    if (action == "synthetic") {
      string mode(argv[2]);
      int n = stoi(argv[3]), n0 = stoi(argv[4]);
      Parameters params;
      params.initialize(mode, argv + 5);
      synthetic_data(n, n0, params);
    } else if (action == "real_data") {
      string graph_name(argv[2]), mode(argv[3]);
      real_world_data(graph_name, get_seed_name(graph_name), REVERSE_NAME.find(mode)->second);
    } else {
      throw std::invalid_argument("Invalid action: " + action);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }
  return 0;
}
