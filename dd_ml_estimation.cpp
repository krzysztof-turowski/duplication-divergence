// Tool for computation the Maximum Likelihood Estimator for various duplication-divergence models.
// Compile: g++ dd_ml_estimation.cpp -O3 -o ./dd_ml_estimation
// Run: ./dd_ml_estimation synthetic MODE n n0 PARAMETERS - e.g. ./dd_ml_estimation synthetic pastor_satorras 100 20 0.5 2.0
//      ./dd_ml_estimation real_data MODE                 - e.g. ./dd_ml_estimation real_data pastor_satorras

#include "dd_header.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-default"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "lib/koala/graph/graph.h"
#pragma GCC diagnostic pop

#include <future>

using namespace std;

typedef Koala::Graph<int, int> Graph;
typedef Koala::Graph<int, int>::PVertex Vertex;

const double STEP_P = 0.1, STEP_R = 1.0, EPS = 10e-9;
const int IS_TRIES = 100;
const bool ML_PARALLEL = false;

template <typename T>
struct counting_iterator
{
    size_t count;
    T dummy;

    counting_iterator() : count(0) { }
    counting_iterator& operator++() { ++count; return *this; }
    T& operator*() { return dummy; }
};

class LikelihoodValue {
public:
  Parameters params;
  double likelihood;

  LikelihoodValue(const Parameters& params_v, const double likelihood_v) : params(params_v), likelihood(likelihood_v) { }
};

inline int get_index(const int &n, const Vertex &v, const Vertex &u) {
  return min(v->getInfo(), u->getInfo()) * n + max(v->getInfo(), u->getInfo());
}

Graph generate_graph_koala(const int &n, const int &n0, const Parameters &params) {
  random_device device;
  mt19937 generator(device());
  uniform_real_distribution<double> edge_distribution(0.0, 1.0);
  double p0 = 1.0;

  Graph G;
  vector<Vertex> V;
  for (int i = 0; i < n; i++) {
    V.push_back(G.addVert(i));
  }
  for (int i = 0; i < n0; i++) {
    for (int j = i + 1; j < n0; j++) {
      if (edge_distribution(generator) <= p0) {
        G.addEdge(V[i], V[j]);
      }
    }
  }
  for (int i = n0; i < n; i++) {
    uniform_int_distribution<int> parent_distribution(0, i - 1);
    int parent = parent_distribution(generator);
    set<Vertex> neighbors = G.getNeighSet(V[parent]);
    if (params.mode == Mode::PURE_DUPLICATION) {
      for (auto v : neighbors) {
        if (edge_distribution(generator) <= params.p) {
          G.addEdge(V[i], v);
        }
      }
    }
    else if (params.mode == Mode::PURE_DUPLICATION_CONNECTED) {
      while(true) {
        for (auto v : neighbors) {
          if (edge_distribution(generator) <= params.p) {
            G.addEdge(V[i], v);
          }
        }
        if (G.getNeighNo(V[i]) > 0) {
          break;
        }
        neighbors = G.getNeighSet(V[parent_distribution(generator)]);
      }
    }
    else if (params.mode == Mode::CHUNG_LU) {
      for (auto v : neighbors) {
        if (edge_distribution(generator) <= params.p) {
          G.addEdge(V[i], v);
        }
      }
      if (edge_distribution(generator) <= params.q) {
        G.addEdge(V[i], V[parent]);
      }
    }
    else if (params.mode == Mode::PASTOR_SATORRAS) {
      for (int j = 0; j < i; j++) {
        if (neighbors.count(V[j])) {
          if (edge_distribution(generator) <= params.p) {
            G.addEdge(V[i], V[j]);
          }
        }
        else {
          if (edge_distribution(generator) <= params.r / i) {
            G.addEdge(V[i], V[j]);
          }
        }
      }
    }
    else {
      assert(0);
    }
  }
  return G;
}

Graph read_graph_koala(const string &graph_name) {
  ifstream graph_file(graph_name);
  if (graph_file.fail()) {
    throw invalid_argument("Missing " + graph_name + " file");
  }
  Graph G;
  vector<Vertex> V;
  int u, v;
  while (!graph_file.eof())
  {
    graph_file >> u >> v;
    if (v >= G.getVertNo()) {
      for (int i = G.getVertNo(); i <= v; i++) {
        V.push_back(G.addVert(i));
      }
    }
    if (u != v) {
      G.addEdge(V[u], V[v]);
    }
  }
  graph_file.close();
  return G;
}

double omega(const Graph &G, vector<int> &V, const int &n, const Parameters &params, const Vertex &v) {
  double i = G.getVertNo(), out = 0;
  for (Vertex u = G.getVert(); u; u = G.getVertNext(u)) {
    if (u != v) {
      if (params.mode == Mode::PASTOR_SATORRAS) {
        double a = V[get_index(n, v, u)], b = G.getNeighNo(v) - a, c = G.getNeighNo(u) - a, d = i + a - b - c;
        out += pow(params.p, a) * pow(params.r / i, b) * pow(1 - params.p, c) * pow(1 - (params.r / i), d);
      }
      else {
        assert(0);
      }
    }
  }
  return out;
}

vector<double> omega(const Graph &G, vector<int> &V, const int &n, const Parameters &params) {
  vector<double> out;
  for (Vertex v = G.getVert(); v; v = G.getVertNext(v)) {
    out.push_back(omega(G, V, n, params, v));
  }
  return out;
}

double likelihood_value(const Graph &G, const int &n0, const Parameters &params, const Parameters &params_0) {
  random_device device;
  mt19937 generator(device());
  int n = G.getVertNo();
  Graph H;
  H.copy(G);
  vector<int> V(n * n);
  for (Vertex v = H.getVert(); v; v = H.getVertNext(v)) {
    set<Vertex> N_v = H.getNeighSet(v);
    for (Vertex u = H.getVert(); u; u = H.getVertNext(u)) {
      if (v < u) {
        set<Vertex> N_u = H.getNeighSet(u);
        V[get_index(n, v, u)] = set_intersection(N_v.begin(), N_v.end(), N_u.begin(), N_u.end(), counting_iterator<Vertex>()).count;
      }
    }
  }
  double ML_value = 0;
  while (H.getVertNo() > n0) {
    vector<double> P = omega(H, V, n, params_0);
    double P_sum = accumulate(P.begin(), P.end(), 0.0);
    if (P_sum == 0.0L) {
      return -numeric_limits<double>::infinity();
    }
    discrete_distribution<int> choose_vertex(P.begin(), P.end());
    int index = choose_vertex(generator);
    Vertex v = H.vertByNo(index);
    ML_value += log(P_sum) + log(omega(H, V, n, params, v)) - log(H.getVertNo()) - log(P[index]);

    set<Vertex> neighbors = G.getNeighSet(v);
    for (Vertex w : neighbors) {
      for (Vertex u : neighbors) {
        if (w < u) {
          --V[get_index(n, w, u)];
        }
      }
    }
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
      break;
    default:
      assert(0);
      break;
  }
  return likelihood_values;
}

void print(const string &name, const vector<LikelihoodValue> &likelihood_values, ostream &out_file) {
  cout << name << endl;
  auto ML = *max_element(
      likelihood_values.begin(), likelihood_values.end(),
      [] (const LikelihoodValue& lhs, const LikelihoodValue& rhs) { return lhs.likelihood > rhs.likelihood; });
  cout << "Best found: " << ML.params.to_string() << " ML score: " << ML.likelihood << endl;
  for (auto const& value : likelihood_values) {
    cout << value.params.to_string() << " ML score: " << value.likelihood << endl;
  }
  for (auto const& value : likelihood_values) {
    out_file << value.params.to_csv() << "," << value.likelihood << " ";
  }
}

void synthetic_data(const int &n, const int &n0, const Parameters &params) {
  Graph G = generate_graph_koala(n, n0, params);
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

int main(int, char *argv[]) {
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
      assert(0);
    }
  } catch (exception &e) {
    cout << "ERROR: " << e.what() << endl;
  }
  return 0;
}
