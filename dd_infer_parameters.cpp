// Tool for inference the values of parameters for various duplication-divergence models.
// Compile: g++ dd_infer_parameters.cpp -O3 -o ./dd_infer_parameters
// Run: ./dd_infer_parameters synthetic MODE n n0 PARAMETERS - e.g. ./dd_infer_parameters synthetic pastor_satorras 100 20 0.5 2.0
//      ./dd_infer_parameters real_data MODE                 - e.g. ./dd_infer_parameters real_data pastor_satorras

#include "dd_header.h"

#include <functional>
#include <iomanip>

using namespace std;

const double R_STEP = 1.0;
const double Q_STEP = 0.1;
const double EPS = 10e-9;
const double P_DISTANCE = 10e-3;

class DataObject {
public:
  double open_triangles = 0, triangles = 0, average_degree = 0, average_degree_squared = 0;
  int no_vertices = 0, low_degree = 0, high_degree = 0;

  void print() const {
    printf("%5d vertices, %5d min degree, %5d max degree\n", no_vertices, low_degree, high_degree);
    printf("%10.0lf open triangles, %10.0lf triangles\n", open_triangles, triangles);
    printf("%10.3lf average degree, %10.3lf average degree squared\n", average_degree, average_degree_squared);
  }
};

void degree_distribution(const vector<set<int>> &G, DataObject &data) {
  data.no_vertices = G.size(), data.low_degree = G.size(), data.high_degree = 0;
  for (auto v : G) {
    data.low_degree = min(data.low_degree, static_cast<int>(v.size()));
    data.high_degree = max(data.high_degree, static_cast<int>(v.size()));
    data.average_degree += v.size();
    data.average_degree_squared += (v.size() * static_cast<double>(v.size()));
  }
  data.average_degree /= data.no_vertices;
  data.average_degree_squared /= data.no_vertices;
}

void count_triangles(const vector<set<int>> &G, DataObject &data) {
  double triangles = 0;
  for (auto v : G) {
    data.open_triangles += v.size() * (v.size() - 1) / 2;

    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = next(it, 1); jt != v.end(); ++jt) {
        if (G[*it].count(*jt)) {
          data.triangles++;
        }
      }
    }
  }
  data.triangles /= 3;
}

DataObject get_params_for_graph(const vector<set<int>> &G) {
  DataObject data;
  degree_distribution(G, data);
  count_triangles(G, data);

  data.print();
  return data;
}

// TODO: replace Score by Parameters
class Score {
public:
  Parameters params;
  double apx, score;
  Score(const Parameters &params_v, const double &opt_v, const double &apx_v) : params(params_v), apx(apx_v), score(abs(opt_v - apx_v)) { }

  bool operator <(const Score &other) const {
    return score < other.score;
  }
};

inline bool contains(const double &low, const double &high, const double &value) {
  return min(low, high) <= value && max(low, high) >= value;
}

bool check_close_solutions(const vector<Score> &V, const Score &value) {
  for (int i = 0; V[i].score < value.score; i++) {
    if (abs(V[i].params.p - value.params.p) < P_DISTANCE) {
      return true;
    }
  }
  return false;
}

void print(const string &name, const double &g0_value, const double &g_value, vector<Score> &V, ostream &out_file, bool erase = false) {
  if (erase) {
    sort(V.begin(), V.end());
    V.erase(remove_if(V.begin(), V.end(), [&](Score &value) { 
      return check_close_solutions(V, value); }), V.end());
  }

  cout << fixed << setprecision(3) << name << " - G0: " << g0_value  << ", G: " << g_value << endl;
  if (!V.empty()) {
    for(auto v : V) {
      cout << v.params.to_string() << endl;
    }
  }
  else {
    cout << "There are no suitable parameter values" << endl;
  }
  for(auto v : V) {
    out_file << v.params.to_csv();
  }
  out_file << endl;
}

void chung_lu_estimate_iterative(DataObject &apx_data, const double &p, const double &q, const int &n0, const int &n) {
  for (int i = n0; i < n; i++) {
    double id = i, D = apx_data.average_degree, D2 = apx_data.average_degree_squared, S2 = apx_data.open_triangles, C3 = apx_data.triangles;
    apx_data.triangles = C3 * (1 + (3 * p * p) / id) + D * p * q;
    apx_data.average_degree = D * (1 + (2 * p - 1) / (id + 1)) + 2 * q / (id + 1);
    apx_data.average_degree_squared =
        D2 * (1 + (2 * p + p * p - 1) / (id + 1))
        + D * (2 * q + 2 * p + 2 * p * q - p * p) / (id + 1) + 2 * q / (id + 1);
    apx_data.open_triangles = S2 * (1 + (2 * p + p * p) / id) + D * (p * q + p + q);
    // apx_data.open_triangles = (id + 1) * (apx_data.average_degree_squared - apx_data.average_degree) / 2;
  }
}

Score chung_lu_binary_search_p(
    const DataObject &g0_data, const DataObject &g_data, const double &high_value, const double &q,
    function<double(DataObject const&)> get_value) {
  double p_low = 0.0, p_high = 1.0, error = EPS;
  DataObject apx_data;
  while (p_high - p_low > error) {
    double p_mid = (p_high + p_low) / 2;
    apx_data = g0_data;
    chung_lu_estimate_iterative(apx_data, p_mid, q, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_data), high_value, get_value(g_data))) {
      p_low = p_mid;
    }
    else {
      p_high = p_mid;
    }
  }
  Parameters params;
  params.initialize_chung_lu(p_low, q);
  return Score(params, get_value(g_data), get_value(apx_data));
}

vector<Score> chung_lu_get_scores(const DataObject &g0_data, const DataObject &g_data, function<double (DataObject const&)> get_value) {
  vector<Score> params;
  // TODO: variable step, increased if points too close
  for (double q = 0.0; q <= 1.0 + EPS; q += Q_STEP) {
    DataObject apx_min(g0_data), apx_max(g0_data);
    chung_lu_estimate_iterative(apx_min, 0.0, q, g0_data.no_vertices, g_data.no_vertices);
    chung_lu_estimate_iterative(apx_max, 1.0, q, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_min), get_value(apx_max), get_value(g_data))) {
      Score value = chung_lu_binary_search_p(g0_data, g_data, get_value(apx_max), q, get_value);
      params.push_back(value);
    }
  }
  return params;
}

void chung_lu_estimate_parameter(
    const string &name, const DataObject &g0_data, const DataObject &g_data,
    function<double (DataObject const&)> get_value, ofstream &out_file) {
  vector<Score> S = chung_lu_get_scores(g0_data, g_data, get_value);
  print(name, get_value(g0_data), get_value(g_data), S, out_file);
}

void chung_lu_estimate(const DataObject &g0_data, const DataObject &g_data, const vector<set<int>> &G0, ofstream &out_file) {
  auto D_lambda = [](const DataObject &data) { return data.average_degree; };
  chung_lu_estimate_parameter("Average degree", g_data, g0_data, D_lambda, out_file);
  
  auto D2_lambda = [](const DataObject &data) { return data.average_degree_squared; };
  chung_lu_estimate_parameter("Average degree squared", g_data, g0_data, D2_lambda, out_file);

  auto C3_lambda = [](const DataObject &data) { return data.triangles; };
  chung_lu_estimate_parameter("Triangles", g_data, g0_data, C3_lambda, out_file);
}

void pastor_satorras_estimate_iterative(DataObject &apx_data, const double &p, const double &r, const int &n0, const int &n) {
  for (int i = n0; i < n; i++) {
    double id = i, D = apx_data.average_degree, D2 = apx_data.average_degree_squared, S2 = apx_data.open_triangles, C3 = apx_data.triangles;
    apx_data.triangles =
        C3 * (1 + 3 * p * p / id - 6 * p * r / (id * id) + 3 * r * r / (id * id * id))
        + D2 * (p * r / id - r * r / (id * id)) + D * r * r / (2 * id);
    apx_data.open_triangles =
        S2 * (1 + (2 * p + p * p) / id - 2 * (p + 1) * r / (id * id) + r * r / (id * id * id))
        + D * (p * r + p + r - (p * r + r + r * r) / id + r * r / (id * id))
        + r * r / 2 - r * r / (2 * id);
    apx_data.average_degree_squared =
        D2 * (1 + (2 * p + p * p - 1) / (id + 1) - 2 * r * (1 + p) / (id * (id + 1)) + r * r / (id * id * (id + 1)))
        + D * ((2 * p - p * p + 2 * p * r + 2 * r) / (id + 1) - (2 * r + 2 * r * r) / (id * (id + 1))
        + r * r / (id * id * (id + 1))) + (2 * r + 2 * r * r) / (id + 1) - r * r / (id * (id + 1));
    apx_data.average_degree = D * (1 + (2 * p - 1) / (id + 1) - (2 * r) / (id * (id + 1))) + 2 * r / (id + 1);
    // apx_data.open_triangles = (id + 1) * (apx_data.average_degree_squared - apx_data.average_degree) / 2;
  }
}

Score pastor_satorras_binary_search_p(
    const DataObject &g0_data, const DataObject &g_data, const double &high_value, const double &r,
    function<double(DataObject const&)> get_value) {
  double p_low = 0.0, p_high = 1.0, error = EPS;
  DataObject apx_data;
  while (p_high - p_low > error) {
    double p_mid = (p_high + p_low) / 2;
    apx_data = g0_data;
    pastor_satorras_estimate_iterative(apx_data, p_mid, r, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_data), high_value, get_value(g_data))) {
      p_low = p_mid;
    }
    else {
      p_high = p_mid;
    }
  }
  Parameters params;
  params.initialize_pastor_satorras(p_low, r);
  return Score(params, get_value(g_data), get_value(apx_data));
}

Score pastor_satorras_estimate_binary_search_r(
    const DataObject &g0_data, const DataObject &g_data, const double &high_value, const double &p,
    function<double(DataObject const&)> get_value) {
  double r_low = 0.0, r_high = g0_data.no_vertices, error = EPS;
  DataObject apx_data;
  while (r_high - r_low > error) {
    double r_mid = (r_high + r_low) / 2;
    apx_data = g0_data;
    pastor_satorras_estimate_iterative(apx_data, p, r_mid, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_data), high_value, get_value(g_data))) {
      r_low = r_mid;
    }
    else {
      r_high = r_mid;
    }
  }
  Parameters params;
  params.initialize_pastor_satorras(p, r_low);
  return Score(params, get_value(g_data), get_value(apx_data));
}

vector<Score> pastor_satorras_get_scores(const DataObject &g0_data, const DataObject &g_data, function<double (DataObject const&)> get_value) {
  vector<Score> params;
  // TODO: variable step, increased if points too close
  double r_max = g0_data.no_vertices;
  for (double r = 0.0; r <= r_max + EPS; r += R_STEP) {
    DataObject apx_min(g0_data), apx_max(g0_data);
    pastor_satorras_estimate_iterative(apx_min, 0.0, r, g0_data.no_vertices, g_data.no_vertices);
    pastor_satorras_estimate_iterative(apx_max, 1.0, r, g0_data.no_vertices, g_data.no_vertices);
    if (g0_data.triangles > apx_min.triangles || g0_data.open_triangles > apx_min.open_triangles) {
      continue;
    }
    if (contains(get_value(apx_min), get_value(apx_max), get_value(g_data))) {
      Score value = pastor_satorras_binary_search_p(g0_data, g_data, get_value(apx_max), r, get_value);
      // TODO: add tolerance interval
      params.push_back(value);
    }
  }
  return params;
}

void pastor_satorras_estimate_parameter(
    const string &name, const DataObject &g0_data, const DataObject &g_data,
    function<double (DataObject const&)> get_value, ofstream &out_file) {
  vector<Score> S = pastor_satorras_get_scores(g0_data, g_data, get_value);
  print(name, get_value(g0_data), get_value(g_data), S, out_file);
}

void pastor_satorras_estimate(const DataObject &g0_data, const DataObject &g_data, const vector<set<int>> &G0, ofstream &out_file) {
  auto D_lambda = [](const DataObject &data) { return data.average_degree; };
  pastor_satorras_estimate_parameter("Average degree", g_data, g0_data, D_lambda, out_file);
  
  auto D2_lambda = [](const DataObject &data) { return data.average_degree_squared; };
  pastor_satorras_estimate_parameter("Average degree squared", g_data, g0_data, D2_lambda, out_file);

  auto C3_lambda = [](const DataObject &data) { return data.triangles; };
  pastor_satorras_estimate_parameter("Triangles", g_data, g0_data, C3_lambda, out_file);
}

void process_graph(const vector<set<int>> &G, const vector<set<int>> &G0, const Mode &mode, ofstream &out_file) {
  DataObject g0_data = get_params_for_graph(G0), g_data = get_params_for_graph(G);
  switch (mode) {
    case CHUNG_LU:
      chung_lu_estimate(g_data, g0_data, G0, out_file);
      break;
    case PASTOR_SATORRAS:
      pastor_satorras_estimate(g_data, g0_data, G0, out_file);
      break;
    default:
      assert(0);
  }
}

void synthetic_data(const int &n, const int &n0, const Parameters &params) {
  ofstream out_file(TEMP_FOLDER + "synthetic_" + SHORT_NAME.find(params.mode)->second + ".txt");
  vector<set<int>> G0 = generate_seed(n0, 1.0);
  vector<set<int>> G = G0;
  generate_graph(G, n, params);
  string name = "Synthetic data: " + params.to_string();
  cout << name << endl, out_file << name << endl;
  process_graph(G, G0, params.mode, out_file);
}

void real_world_data(const string &graph_name, const string &seed_name, const Mode &mode) {
  vector<set<int>> G = read_graph(FILES_FOLDER + graph_name);
  vector<set<int>> G0 = read_graph(FILES_FOLDER + seed_name);
  ofstream out_file(TEMP_FOLDER + graph_name.substr(0, graph_name.find_last_of(".")) + "_" + SHORT_NAME.find(mode)->second + ".txt");
  cout << "File: " << graph_name << endl, out_file << graph_name << endl;
  process_graph(G, G0, mode, out_file);
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
      real_world_data("G-100-20-PS-0.70-2.00.txt", "G0-100-20-PS-0.70-2.00.txt", REVERSE_NAME.find(mode)->second);
      real_world_data("G-c-elegans.txt", "G0-c-elegans.txt", REVERSE_NAME.find(mode)->second);
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
