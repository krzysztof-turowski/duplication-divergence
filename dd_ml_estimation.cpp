// Tool for computation the Maximum Likelihood Estimator for various duplication-divergence models.
// Compile: g++ dd_ml_estimation.cpp -O3 -o ./dd_ml_estimation
// Run: ./dd_ml_estimation synthetic MODE n n0 PARAMETERS - e.g. ./dd_ml_estimation synthetic pastor_satorras 100 20 0.5 2.0
//      ./dd_ml_estimation real_data MODE                 - e.g. ./dd_ml_estimation real_data pastor_satorras

#include "dd_koala.h"

#include <future>

using namespace std;

typedef Koala::Graph<int, int> Graph;
typedef Koala::Graph<int, int>::PVertex Vertex;

const double STEP_P = 0.1, STEP_R = 1.0, EPS = 10e-9;
const int IS_TRIES = 100;
const bool ML_PARALLEL = false;

class LikelihoodValue {
public:
  Parameters params;
  double likelihood;

  LikelihoodValue(const Parameters& params_v, const double likelihood_v) : params(params_v), likelihood(likelihood_v) { }
};

double likelihood_value(const Graph &G, const int &n0, const Parameters &params, const Parameters &params_0) {
  random_device device;
  mt19937 generator(device());
  int n = G.getVertNo();
  Graph H(G);
  NeighborhoodStructure aux(H);

  double ML_value = 0;
  while (H.getVertNo() > n0) {
    vector<double> P = get_transition_probability(H, params_0, aux);
    double P_sum = accumulate(P.begin(), P.end(), 0.0);
    if (P_sum == 0.0L) {
      return -numeric_limits<double>::infinity();
    }
    discrete_distribution<int> choose_vertex(P.begin(), P.end());
    int index = choose_vertex(generator);
    Vertex v = H.vertByNo(index);
    ML_value += log(P_sum) + log(get_transition_probability(H, params, v, aux)) - log(H.getVertNo()) - log(P[index]);

    aux.remove_vertex(v, G.getNeighSet(v));
    H.delVert(v);
  }
  return ML_value;
}

LikelihoodValue importance_sampling(const Graph &G, const int &n0, const Parameters &params, const Parameters &params_0) {
  vector<double> likelihood_values(IS_TRIES);
  if (ML_PARALLEL) {
    vector<future<double>> futures(IS_TRIES);
    for(int i = 0; i < IS_TRIES; i++) {
      futures[i] = async(launch::async, &likelihood_value, cref(G), cref(n0), cref(params), cref(params_0));
    }
    for(int i = 0; i < IS_TRIES; i++) {
      likelihood_values[i] = futures[i].get();
    }
  } else {
    for(int i = 0; i < IS_TRIES; i++) {
      likelihood_values[i] = likelihood_value(G, n0, params, params_0);
    }
  }
  double likelihood_score = accumulate(likelihood_values.begin(), likelihood_values.end(), 0.0) / IS_TRIES;
  return LikelihoodValue(params, likelihood_score);
}

// TODO: split into 2 functionalities: finding ML (gradient search etc.) and sampling likelihood over some parameters interval
vector<LikelihoodValue> find_likelihood_values(Graph &G, const int &n0, const Mode &mode) {
  Parameters params_0;
  vector<LikelihoodValue> likelihood_values;
  switch (mode) {
    case PASTOR_SATORRAS:
      {
        // TODO: parametrize parameter ranges and driving values
        const double MAX_R = n0;
        params_0.initialize_pastor_satorras(0.5, MAX_R / 2);
        for (double p = 0.0; p <= 1.0 + EPS; p += STEP_P) {
          for (double r = 0.0; r <= MAX_R + EPS; r += STEP_R) {
            Parameters params;
            params.initialize_pastor_satorras(p, r);
            LikelihoodValue value = importance_sampling(G, n0, params, params_0);
            cerr << "Finished run for " << value.params.to_string() << endl;
            // TODO: consider filtering out infeasible values
            likelihood_values.push_back(value);
          }
        }
      }
      return likelihood_values;
    default:
      throw std::invalid_argument("Invalid mode: " + mode);
  }
}

void print(const string &name, const vector<LikelihoodValue> &likelihood_values, ostream &out_file) {
  cout << name << endl;
  auto ML = *max_element(
      likelihood_values.begin(), likelihood_values.end(),
      [&] (const LikelihoodValue& lhs, const LikelihoodValue& rhs) { return lhs.likelihood < rhs.likelihood; });
  cout << "Best found: " << ML.params.to_string() << " ML score: " << ML.likelihood << endl;
  for (auto const& value : likelihood_values) {
    cout << value.params.to_string() << " ML score: " << value.likelihood << endl;
  }
  for (auto const& value : likelihood_values) {
    out_file << value.params.to_csv() << "," << value.likelihood << " ";
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
  ofstream out_file(TEMP_FOLDER + graph_name.substr(0, graph_name.find_last_of(".")) + "_" + SHORT_NAME.find(mode)->second + "-ML.txt");
  print(graph_name, likelihood_values, out_file);
}

int main(int argc, char *argv[]) {
  try {
    string action(argv[1]), mode(argv[2]);
    if (action == "synthetic") {
      int n = stoi(argv[3]), n0 = stoi(argv[4]);
      Parameters params;
      params.initialize(mode, argv + 5);
      synthetic_data(n, n0, params);
    }
    else if (action == "real_data") {
      // TODO: make datasets parameter-variable
      real_world_data("G-100-20-PS-0.1-0.3.txt", "G0-100-20-PS-0.1-0.3.txt", REVERSE_NAME.find(mode)->second);
      real_world_data("G-100-20-PS-0.7-2.txt", "G0-100-20-PS-0.7-2.txt", REVERSE_NAME.find(mode)->second);
      real_world_data("G-100-20-PS-0.99-3.txt", "G0-100-20-PS-0.99-3.txt", REVERSE_NAME.find(mode)->second);
      real_world_data("G-a-thaliana.txt", "G0-a-thaliana.txt", REVERSE_NAME.find(mode)->second);
      real_world_data("G-c-elegans.txt", "G0-c-elegans.txt", REVERSE_NAME.find(mode)->second);
      real_world_data("G-d-melanogaster.txt", "G0-d-melanogaster.txt", REVERSE_NAME.find(mode)->second);
      real_world_data("G-homo-sapiens.txt", "G0-homo-sapiens.txt", REVERSE_NAME.find(mode)->second);
      real_world_data("G-mus-musculus.txt", "G0-mus-musculus.txt", REVERSE_NAME.find(mode)->second);
      real_world_data("G-s-cerevisiae.txt", "G0-s-cerevisiae.txt", REVERSE_NAME.find(mode)->second);
      real_world_data("G-s-pombe.txt", "G0-s-pombe.txt", REVERSE_NAME.find(mode)->second);
    }
    else {
      throw std::invalid_argument("Invalid action: " + action);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }
  return 0;
}
