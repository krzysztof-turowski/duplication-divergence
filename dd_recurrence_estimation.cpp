// Tool for inference the values of parameters for various duplication-divergence models.
// Compile: g++ dd_recurrence_estimation.cpp -O3 -o ./dd_recurrence_estimation
// Run: ./dd_recurrence_estimation synthetic MODE n n0 PARAMETERS
//   or ./dd_recurrence_estimation real_data FILE MODE

#include "./dd_header.h"

#include <functional>
#include <tuple>

using namespace std;

typedef vector<set<int>> Graph;

const double R_STEP = 1.0;
const double Q_STEP = 0.1;
const double EPS = 10e-9;
const double P_DISTANCE = 10e-3;
const double TI_ALPHA = 0.05;
const int TI_TRIES = 100;
const double PERCENTILE_95 = 1.96, PERCENTILE_99 = 2.575;

enum ToleranceInterval { DIRECT, EMPIRICAL_VARIANCE };

const ToleranceInterval TI_ALGORITHM = ToleranceInterval::DIRECT;

class DataObject {
 public:
  double open_triangles = 0, triangles = 0, average_degree = 0, average_degree_squared = 0;
  int no_vertices = 0, low_degree = 0, high_degree = 0;

  void print() const {
    printf("%5d vertices, %5d min degree, %5d max degree\n", no_vertices, low_degree, high_degree);
    printf("%10.0lf open triangles, %10.0lf triangles\n", open_triangles, triangles);
    printf(
        "%10.3lf average degree, %10.3lf average degree squared\n",
        average_degree, average_degree_squared);
  }
};

void degree_distribution(const Graph &G, DataObject &data) {
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

void count_triangles(const Graph &G, DataObject &data) {
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

DataObject get_params_for_graph(const Graph &G, bool verbose = false) {
  DataObject data;
  degree_distribution(G, data);
  count_triangles(G, data);

  if (verbose) {
    data.print();
  }
  return data;
}

DataObject get_params_for_synthetic_graph(
    const Graph &G0, const int &n, const Parameters &params) {
  Graph G = G0;
  generate_graph(G, n, params);
  return get_params_for_graph(G);
}

tuple<DataObject, DataObject> get_empirical_interval(
    DataObject &g_data, const Graph &G0, const Parameters &params,
    function<double& (DataObject&)> get_value) {
  if (!(TI_ALPHA < 0.5 && TI_ALPHA * TI_TRIES >= 1 && (1 - TI_ALPHA) * TI_TRIES >= 1)) {
    throw invalid_argument(
        "Invalid tolerance interval constants: TI_ALPHA = " + to_string(TI_ALPHA)
            + ", TI_TRIES = " + to_string(TI_TRIES));
  }
  vector<DataObject> values(TI_TRIES);
  #pragma omp parallel for
  for (int i = 0; i < TI_TRIES; i++) {
    values[i] = get_params_for_synthetic_graph(G0, g_data.no_vertices, params);
    #pragma omp critical
    {
      std::cerr << "Run " << i + 1 << "/" << TI_TRIES << std::endl;
    }
  }
  switch (TI_ALGORITHM) {
    case ToleranceInterval::DIRECT: {
        sort(
            values.begin(), values.end(), [&](DataObject &a, DataObject &b) {
                return get_value(a) < get_value(b);
            });
        int low = floor(TI_ALPHA * TI_TRIES), high = (TI_TRIES - 1) - floor(TI_ALPHA * TI_TRIES);
        return make_tuple(
            get_value(values[low]) < get_value(g_data) ? values[low] : g_data,
            get_value(values[high]) > get_value(g_data) ? values[high] : g_data);
      }
    case ToleranceInterval::EMPIRICAL_VARIANCE: {
        double empirical_first_moment =
            accumulate(
                values.begin(), values.end(), 0.0,
                [&](double value, DataObject &a) { return value + get_value(a); }) / values.size();
        double empirical_second_moment =
            accumulate(
                values.begin(), values.end(), 0.0,
                [&](double value, DataObject &a) {
                    return value + get_value(a) * get_value(a);
                }) / values.size();
        double empirical_std_deviation =
            sqrt(empirical_second_moment - empirical_first_moment * empirical_first_moment);
        DataObject g_data_low(g_data), g_data_high(g_data);
        get_value(g_data_low) -= PERCENTILE_95 * empirical_std_deviation;
        get_value(g_data_high) += PERCENTILE_95 * empirical_std_deviation;
        return make_tuple(g_data_low, g_data_high);
      }
    default:
      throw invalid_argument("Invalid tolerance interval computation algorithm");
  }
}

inline bool contains(const double &low, const double &high, const double &value) {
  return min(low, high) <= value && max(low, high) >= value;
}

void print(
    const string &name, const double &g0_value, const double &g_value,
    const vector<Parameters> &V, ostream &out_file) {
  cout << fixed << setprecision(3) << name << " - G0: " << g0_value  << ", G: " << g_value << endl;
  if (!V.empty()) {
    for (const auto &v : V) {
      cout << v.to_string() << endl;
    }
  } else {
    cout << "There are no suitable parameter values" << endl;
  }
  for (const auto &v : V) {
    out_file << v.to_csv() << " ";
  }
  out_file << endl;
}

void print(
    const string &name, const double &g0_value, const double &g_value,
    const vector<Parameters> &V_low, const vector<Parameters> &V, const vector<Parameters> &V_high,
    ostream &out_file) {
  cout << fixed << setprecision(3) << name << " - G0: " << g0_value  << ", G: " << g_value << endl;
  if (!V.empty()) {
    for (int i = 0; i < static_cast<int>(V.size()); i++) {
      cout << V[i].to_string(V_low[i], V_high[i]) << endl;
    }
  } else {
    cout << "There are no suitable parameter values" << endl;
  }
  for (int i = 0; i < static_cast<int>(V.size()); i++) {
    out_file << V_low[i].to_csv() << ";" << V[i].to_csv() << ";" << V_high[i].to_csv() << " ";
  }
  out_file << endl;
}

void chung_lu_estimate_iterative(
    DataObject &apx_data, const double &p, const double &q, const int &n0, const int &n) {
  for (int i = n0; i < n; i++) {
    double id = i, D = apx_data.average_degree, D2 = apx_data.average_degree_squared,
        S2 = apx_data.open_triangles, C3 = apx_data.triangles;
    apx_data.triangles = C3 * (1 + (3 * p * p) / id) + D * p * q;
    apx_data.average_degree = D * (1 + (2 * p - 1) / (id + 1)) + 2 * q / (id + 1);
    apx_data.average_degree_squared =
        D2 * (1 + (2 * p + p * p - 1) / (id + 1))
        + D * (2 * q + 2 * p + 2 * p * q - p * p) / (id + 1) + 2 * q / (id + 1);
    apx_data.open_triangles = S2 * (1 + (2 * p + p * p) / id) + D * (p * q + p + q);
    // apx_data.open_triangles =
        // (id + 1) * (apx_data.average_degree_squared - apx_data.average_degree) / 2;
  }
}

Parameters chung_lu_binary_search_p(
    DataObject &g0_data, DataObject &g_data, const double &q,
    function<double& (DataObject&)> get_value) {
  double p_low = 0.0, p_high = 1.0, error = EPS;
  DataObject apx_data;
  while (p_high - p_low > error) {
    double p_mid = (p_high + p_low) / 2;
    apx_data = g0_data;
    chung_lu_estimate_iterative(apx_data, p_mid, q, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_data), numeric_limits<double>::infinity(), get_value(g_data))) {
      p_low = p_mid;
    } else {
      p_high = p_mid;
    }
  }
  Parameters params;
  params.initialize_chung_lu(p_low, q);
  return params;
}

Parameters chung_lu_binary_search_q(
    DataObject &g0_data, DataObject &g_data, const double &p,
    function<double& (DataObject&)> get_value) {
  double q_low = 0.0, q_high = 1.0, error = EPS;
  DataObject apx_data;
  while (q_high - q_low > error) {
    double q_mid = (q_high + q_low) / 2;
    apx_data = g0_data;
    chung_lu_estimate_iterative(apx_data, p, q_mid, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_data), numeric_limits<double>::infinity(), get_value(g_data))) {
      q_low = q_mid;
    } else {
      q_high = q_mid;
    }
  }
  Parameters params;
  params.initialize_chung_lu(p, q_low);
  return params;
}

vector<Parameters> chung_lu_get_parameters(
    DataObject &g0_data, DataObject &g_data, function<double& (DataObject&)> get_value) {
  vector<Parameters> S;
  // TODO(krzysztof-turowski): variable step, increased if points too close
  for (double q = 0.0; q <= 1.0 + EPS; q += Q_STEP) {
    DataObject apx_min(g0_data), apx_max(g0_data);
    chung_lu_estimate_iterative(apx_min, 0.0, q, g0_data.no_vertices, g_data.no_vertices);
    chung_lu_estimate_iterative(apx_max, 1.0, q, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_min), get_value(apx_max), get_value(g_data))) {
      S.push_back(chung_lu_binary_search_p(g0_data, g_data, q, get_value));
    }
  }
  return S;
}

void chung_lu_estimate_parameter(
    const string &name, DataObject &g_data, DataObject &g0_data, const Graph &G0,
    function<double& (DataObject&)> get_value, ofstream &out_file, bool tolerance_interval) {
  vector<Parameters> S = chung_lu_get_parameters(g0_data, g_data, get_value);
  if (tolerance_interval) {
    vector<Parameters> S_low, S_high;
    for (auto params : S) {
      DataObject g_data_low, g_data_high;
      tie(g_data_low, g_data_high) = get_empirical_interval(g_data, G0, params, get_value);
      S_low.push_back(chung_lu_binary_search_q(g0_data, g_data_low, params.p, get_value));
      S_high.push_back(chung_lu_binary_search_q(g0_data, g_data_high, params.p, get_value));
    }
    print(name, get_value(g0_data), get_value(g_data), S_low, S, S_high, out_file);
  } else {
    print(name, get_value(g0_data), get_value(g_data), S, out_file);
  }
}

void chung_lu_estimate(
    DataObject &g_data, DataObject &g0_data, const Graph &G0, ofstream &out_file) {
  auto D_lambda = [](DataObject &data) -> double& { return data.average_degree; };
  chung_lu_estimate_parameter("Average degree", g_data, g0_data, G0, D_lambda, out_file, true);

  auto D2_lambda = [](DataObject &data) -> double& { return data.average_degree_squared; };
  chung_lu_estimate_parameter(
      "Average degree squared", g_data, g0_data, G0, D2_lambda, out_file, false);

  auto S2_lambda = [](DataObject &data) -> double& { return data.open_triangles; };
  chung_lu_estimate_parameter("Open triangles", g_data, g0_data, G0, S2_lambda, out_file, true);

  auto C3_lambda = [](DataObject &data) -> double& { return data.triangles; };
  chung_lu_estimate_parameter("Triangles", g_data, g0_data, G0, C3_lambda, out_file, true);
}

void pastor_satorras_estimate_iterative(
    DataObject &apx_data, const double &p, const double &r, const int &n0, const int &n) {
  for (int i = n0; i < n; i++) {
    double id = i, D = apx_data.average_degree, D2 = apx_data.average_degree_squared,
        S2 = apx_data.open_triangles, C3 = apx_data.triangles;
    apx_data.triangles =
        C3 * (1 + 3 * p * p / id - 6 * p * r / (id * id) + 3 * r * r / (id * id * id))
        + D2 * (p * r / id - r * r / (id * id)) + D * r * r / (2 * id);
    apx_data.open_triangles =
        S2 * (1 + (2 * p + p * p) / id - 2 * (p + 1) * r / (id * id)
            + r * r / (id * id * id))
        + D * (p * r + p + r - (p * r + r + r * r) / id + r * r / (id * id))
        + r * r / 2 - r * r / (2 * id);
    apx_data.average_degree_squared =
        D2 * (1 + (2 * p + p * p - 1) / (id + 1) - 2 * r * (1 + p) / (id * (id + 1))
            + r * r / (id * id * (id + 1)))
        + D * ((2 * p - p * p + 2 * p * r + 2 * r) / (id + 1)
            - (2 * r + 2 * r * r) / (id * (id + 1))
        + r * r / (id * id * (id + 1))) + (2 * r + 2 * r * r) / (id + 1) - r * r / (id * (id + 1));
    apx_data.average_degree =
        D * (1 + (2 * p - 1) / (id + 1) - (2 * r) / (id * (id + 1))) + 2 * r / (id + 1);
    // apx_data.open_triangles =
        // (id + 1) * (apx_data.average_degree_squared - apx_data.average_degree) / 2;
  }
}

Parameters pastor_satorras_binary_search_p(
    DataObject &g0_data, DataObject &g_data, const double &r,
    function<double& (DataObject&)> get_value) {
  double p_low = 0.0, p_high = 1.0, error = EPS;
  DataObject apx_data;
  while (p_high - p_low > error) {
    double p_mid = (p_high + p_low) / 2;
    apx_data = g0_data;
    pastor_satorras_estimate_iterative(
        apx_data, p_mid, r, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_data), numeric_limits<double>::infinity(), get_value(g_data))) {
      p_low = p_mid;
    } else {
      p_high = p_mid;
    }
  }
  Parameters params;
  params.initialize_pastor_satorras(p_low, r);
  return params;
}

Parameters pastor_satorras_binary_search_r(
    DataObject &g0_data, DataObject &g_data, const double &p,
    function<double& (DataObject&)> get_value) {
  double r_low = 0.0, r_high = g0_data.no_vertices, error = EPS;
  DataObject apx_data;
  while (r_high - r_low > error) {
    double r_mid = (r_high + r_low) / 2;
    apx_data = g0_data;
    pastor_satorras_estimate_iterative(
        apx_data, p, r_mid, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_data), numeric_limits<double>::infinity(), get_value(g_data))) {
      r_low = r_mid;
    } else {
      r_high = r_mid;
    }
  }
  Parameters params;
  params.initialize_pastor_satorras(p, r_low);
  return params;
}

vector<Parameters> pastor_satorras_get_parameters(
    DataObject &g0_data, DataObject &g_data, function<double& (DataObject&)> get_value) {
  vector<Parameters> S;
  // TODO(krzysztof-turowski): variable step, increased if points too close
  double r_max = g0_data.no_vertices;
  for (double r = 0.0; r <= r_max + EPS; r += R_STEP) {
    DataObject apx_min(g0_data), apx_max(g0_data);
    pastor_satorras_estimate_iterative(apx_min, 0.0, r, g0_data.no_vertices, g_data.no_vertices);
    pastor_satorras_estimate_iterative(apx_max, 1.0, r, g0_data.no_vertices, g_data.no_vertices);
    if (g0_data.triangles > apx_min.triangles || g0_data.open_triangles > apx_min.open_triangles) {
      continue;
    }
    if (contains(get_value(apx_min), get_value(apx_max), get_value(g_data))) {
      S.push_back(pastor_satorras_binary_search_p(g0_data, g_data, r, get_value));
    }
  }
  return S;
}

void pastor_satorras_estimate_parameter(
    const string &name, DataObject &g_data, DataObject &g0_data, const Graph &G0,
    function<double& (DataObject&)> get_value, ofstream &out_file, bool tolerance_interval) {
  vector<Parameters> S = pastor_satorras_get_parameters(g0_data, g_data, get_value);
  if (tolerance_interval) {
    vector<Parameters> S_low, S_high;
    for (auto params : S) {
      DataObject g_data_low, g_data_high;
      tie(g_data_low, g_data_high) = get_empirical_interval(g_data, G0, params, get_value);
      S_low.push_back(pastor_satorras_binary_search_r(g0_data, g_data_low, params.p, get_value));
      S_high.push_back(pastor_satorras_binary_search_r(g0_data, g_data_high, params.p, get_value));
    }
    print(name, get_value(g0_data), get_value(g_data), S_low, S, S_high, out_file);
  } else {
    print(name, get_value(g0_data), get_value(g_data), S, out_file);
  }
}

void pastor_satorras_estimate(
    DataObject &g_data, DataObject &g0_data, const Graph &G0, ofstream &out_file) {
  auto D_lambda = [](DataObject &data) -> double& { return data.average_degree; };
  pastor_satorras_estimate_parameter(
      "Average degree", g_data, g0_data, G0, D_lambda, out_file, true);

  auto D2_lambda = [](DataObject &data) -> double& { return data.average_degree_squared; };
  pastor_satorras_estimate_parameter(
      "Average degree squared", g_data, g0_data, G0, D2_lambda, out_file, false);

  auto S2_lambda = [](DataObject &data) -> double& { return data.open_triangles; };
  pastor_satorras_estimate_parameter(
      "Open triangles", g_data, g0_data, G0, S2_lambda, out_file, true);

  auto C3_lambda = [](DataObject &data) -> double& { return data.triangles; };
  pastor_satorras_estimate_parameter("Triangles", g_data, g0_data, G0, C3_lambda, out_file, true);
}

void process_graph(
    const Graph &G, const Graph &G0, const Mode &mode, ofstream &out_file) {
  DataObject g0_data = get_params_for_graph(G0, true), g_data = get_params_for_graph(G, true);
  switch (mode) {
    case CHUNG_LU:
      chung_lu_estimate(g_data, g0_data, G0, out_file);
      break;
    case PASTOR_SATORRAS:
      pastor_satorras_estimate(g_data, g0_data, G0, out_file);
      break;
    default:
      throw invalid_argument("Invalid mode: " + LONG_NAME.find(mode)->second);
  }
}

void synthetic_data(const int &n, const int &n0, const Parameters &params) {
  Graph G0 = generate_seed(n0, 1.0);
  Graph G = G0;
  generate_graph(G, n, params);
  ofstream out_file(TEMP_FOLDER + get_synthetic_filename(n, n0, params, ""));
  cout << "Synthetic data: " + params.to_string() << endl;
  process_graph(G, G0, params.mode, out_file);
}

void real_world_data(const string &graph_name, const string &seed_name, const Mode &mode) {
  Graph G = read_graph(FILES_FOLDER + graph_name);
  Graph G0 = read_graph(FILES_FOLDER + seed_name);
  ofstream out_file(TEMP_FOLDER + get_real_filename(graph_name, mode, ""));
  cout << "File: " << graph_name << endl;
  process_graph(G, G0, mode, out_file);
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
      throw invalid_argument("Invalid action: " + action);
    }
  } catch (const exception &e) {
    cerr << "ERROR: " << e.what() << endl;
  }
  return 0;
}
