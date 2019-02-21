// Tool for computation the Maximum Likelihood Estimator for various duplication-divergence models.
// Compile: g++ dd_ml_estimation.cpp -O3 -o ./dd_ml_estimation
// Run: ./dd_ml_estimation synthetic MODE n n0 PARAMETERS or ./dd_ml_estimation real_data MODE

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

#include <future>

using namespace std;

typedef Koala::Graph<int, int> Graph;
typedef Koala::Graph<int, int>::PVertex Vertex;

const double STEP_P = 0.1, STEP_R = 1.0;
const int IS_TRIES = 100;
const bool ML_PARALLEL = false;

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
  NeighborhoodStructure aux(H);

  long double ML_value = 0;
  while (H.getVertNo() > n0) {
    vector<long double> P = get_transition_probability(H, params_0, aux);
    long double P_sum = accumulate(P.begin(), P.end(), 0.0L);
    if (P_sum == 0.0L) {
      return 0.0L;
    }
    discrete_distribution<int> choose_vertex(P.begin(), P.end());
    int index = choose_vertex(generator);
    Vertex v = H.vertByNo(index);
    long double p_transition = get_transition_probability(H, params, v, aux);
    if (p_transition == 0.0L) {
      return 0.0L;
    }
    long double log_p_transition = log(p_transition);
    ML_value += log(P_sum) + log_p_transition - log(P[index]);

    aux.remove_vertex(v, H.getNeighSet(v));
    H.delVert(v);
    aux.verify(H);
  }
  return exp(ML_value);
}

LikelihoodValue importance_sampling(
    const Graph &G, const int &n0, const Parameters &params, const Parameters &params_0) {
  vector<long double> likelihood_values(IS_TRIES);
  if (ML_PARALLEL) {
    vector<future<long double>> futures(IS_TRIES);
    for (int i = 0; i < IS_TRIES; i++) {
      futures[i] =
          async(launch::async, &likelihood_value, cref(G), cref(n0), cref(params), cref(params_0));
    }
    for (int i = 0; i < IS_TRIES; i++) {
      likelihood_values[i] = futures[i].get();
    }
  } else {
    for (int i = 0; i < IS_TRIES; i++) {
      likelihood_values[i] = likelihood_value(G, n0, params, params_0);
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
  Graph G = generate_seed_koala(n0, 1.0);
  generate_graph_koala(G, n, params);
  auto likelihood_values = find_likelihood_values(G, n0, params.mode);
  ofstream out_file(TEMP_FOLDER + "synthetic_" + SHORT_NAME.find(params.mode)->second + "-ML.txt");
  print("Synthetic data: " + params.to_string(), likelihood_values, out_file);
}

void real_world_data(const string &graph_name, const string &seed_name, const Mode &mode) {
  Graph G = read_graph_koala(FILES_FOLDER + graph_name);
  int n0 = read_graph_size(FILES_FOLDER + seed_name);
  auto likelihood_values = find_likelihood_values(G, n0, mode);
  ofstream out_file(TEMP_FOLDER + graph_name.substr(0, graph_name.find_last_of("."))
      + "_" + SHORT_NAME.find(mode)->second + "-ML.txt");
  print(graph_name, likelihood_values, out_file);
}

int main(int, char *argv[]) {
  try {
    string action(argv[1]), mode(argv[2]);
    if (action == "synthetic") {
      int n = stoi(argv[3]), n0 = stoi(argv[4]);
      Parameters params;
      params.initialize(mode, argv + 5);
      synthetic_data(n, n0, params);
    } else if (action == "real_data") {
      // TODO(unknown): make datasets parameter-variable
      Mode mode_v = REVERSE_NAME.find(mode)->second;
      real_world_data("G-100-20-PS-0.1-0.3.txt", "G0-100-20-PS-0.1-0.3.txt", mode_v);
      real_world_data("G-100-20-PS-0.7-2.txt", "G0-100-20-PS-0.7-2.txt", mode_v);
      real_world_data("G-100-20-PS-0.99-3.txt", "G0-100-20-PS-0.99-3.txt", mode_v);
      real_world_data("G-a-thaliana.txt", "G0-a-thaliana.txt", mode_v);
      real_world_data("G-c-elegans.txt", "G0-c-elegans.txt", mode_v);
      real_world_data("G-d-melanogaster.txt", "G0-d-melanogaster.txt", mode_v);
      real_world_data("G-homo-sapiens.txt", "G0-homo-sapiens.txt", mode_v);
      real_world_data("G-mus-musculus.txt", "G0-mus-musculus.txt", mode_v);
      real_world_data("G-s-cerevisiae.txt", "G0-s-cerevisiae.txt", mode_v);
      real_world_data("G-s-pombe.txt", "G0-s-pombe.txt", mode_v);
    } else {
      throw std::invalid_argument("Invalid action: " + action);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }
  return 0;
}
