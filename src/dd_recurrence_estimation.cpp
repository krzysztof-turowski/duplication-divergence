/*
Tool for inference the values of parameters for various duplication-divergence models.

Compile: make dd_recurrence_estimation

Syntax: ./dd_recurrence_estimation <options>
<options> are
-action:
  synthetic: synthetic random graph generations from duplication-divergence model
  real_data: real-world network
-graph: If action is `real_data`, then provide graph file name (located in `files/` folder). File should be in edge list format.
-st: Number of independent tries for estimating the empirical tolerance interval for parameters.
-mode: {pure_duplication, pastor_satorras, chung_lu}. In case of `synthetic` action, the mode (type) of the duplication-divergence graph model.
<parameters>: Depending on `mode`, the parameters `p,q,r` of the DD model.
-n: The size of a graph in the case of `synthetic` action.
-n0, -p0: The parameters for generating a seed graph in the case of `synthetic` action.

Example runs:
  ./dd_recurrence_estimation -action:synthetic -n:100 -n0:10 -mode:pastor_satorras -p:0.5 -r:2.0 -p0:0.6 -st:1000
  ./dd_recurrence_estimation -action:real_data -graph:G-test.txt -mode:pastor_satorras -st:1000
*/

#include "./dd_input.h"
#include "./dd_graphlets.h"

#include <functional>
#include <tuple>

const double R_STEP = 0.01, R_EXP = 1.5;
const double Q_STEP = 0.001, Q_EXP = 1.5;
const double EPS = 10e-9;
const double P_DISTANCE = 10e-2;
const double TI_ALPHA = 0.05;
const double PERCENTILE_95 = 1.96, PERCENTILE_99 = 2.575;
int TI_TRIES;

enum ToleranceInterval { DIRECT, EMPIRICAL_VARIANCE };

const ToleranceInterval TI_ALGORITHM = ToleranceInterval::DIRECT;

std::tuple<DataObject, DataObject> get_empirical_interval(
    DataObject &g_data, const Graph &G0, const Parameters &params,
    std::function<double& (DataObject&)> get_value) {
  if (!(TI_ALPHA < 0.5 && TI_ALPHA * TI_TRIES >= 1 && (1 - TI_ALPHA) * TI_TRIES >= 1)) {
    throw std::invalid_argument(
        "Invalid tolerance interval constants: TI_ALPHA = " + std::to_string(TI_ALPHA)
            + ", TI_TRIES = " + std::to_string(TI_TRIES));
  }
  std::vector<DataObject> values(TI_TRIES);
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
        std::sort(
            values.begin(), values.end(), [&](DataObject &a, DataObject &b) {
                return get_value(a) < get_value(b);
            });
        int low = floor(TI_ALPHA * TI_TRIES), high = (TI_TRIES - 1) - floor(TI_ALPHA * TI_TRIES);
        return std::make_tuple(
            get_value(values[low]) < get_value(g_data) ? values[low] : g_data,
            get_value(values[high]) > get_value(g_data) ? values[high] : g_data);
      }
    case ToleranceInterval::EMPIRICAL_VARIANCE: {
        double empirical_first_moment =
            std::accumulate(
                values.begin(), values.end(), 0.0,
                [&](double value, DataObject &a) { return value + get_value(a); }) / values.size();
        double empirical_second_moment =
            std::accumulate(
                values.begin(), values.end(), 0.0,
                [&](double value, DataObject &a) {
                    return value + get_value(a) * get_value(a);
                }) / values.size();
        double empirical_std_deviation =
            sqrt(empirical_second_moment - empirical_first_moment * empirical_first_moment);
        DataObject g_data_low(g_data), g_data_high(g_data);
        get_value(g_data_low) -= PERCENTILE_95 * empirical_std_deviation;
        get_value(g_data_high) += PERCENTILE_95 * empirical_std_deviation;
        return std::make_tuple(g_data_low, g_data_high);
      }
    default:
      throw std::invalid_argument("Invalid tolerance interval computation algorithm");
  }
}

inline bool contains(const double &low, const double &high, const double &value) {
  return std::min(low, high) <= value && std::max(low, high) >= value;
}

void print(
    const std::string &name, const double &g0_value, const double &g_value,
    const std::vector<Parameters> &V, std::ostream &out_file) {
  std::cout << std::fixed << std::setprecision(3) << name
      << " - G0: " << g0_value << ", G: " << g_value << std::endl;
  if (!V.empty()) {
    for (const auto &v : V) {
      std::cout << v.to_string() << std::endl;
    }
  } else {
    std::cout << "There are no suitable parameter values" << std::endl;
  }
  for (const auto &v : V) {
    out_file << v.to_csv() << " ";
  }
  out_file << std::endl;
}

void print(
    const std::string &name, const double &g0_value, const double &g_value,
    const std::vector<Parameters> &V_low, const std::vector<Parameters> &V,
    const std::vector<Parameters> &V_high, std::ostream &out_file) {
  std::cout << std::fixed << std::setprecision(3) << name
      << " - G0: " << g0_value << ", G: " << g_value << std::endl;
  if (!V.empty()) {
    for (size_t i = 0; i < V.size(); i++) {
      std::cout << V[i].to_string(V_low[i], V_high[i]) << std::endl;
    }
  } else {
    std::cout << "There are no suitable parameter values" << std::endl;
  }
  for (size_t i = 0; i < V.size(); i++) {
    out_file << V_low[i].to_csv() << ";" << V[i].to_csv() << ";" << V_high[i].to_csv() << " ";
  }
  out_file << std::endl;
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
    std::function<double& (DataObject&)> get_value) {
  double p_low = 0.0, p_high = 1.0, error = EPS;
  DataObject apx_data;
  while (p_high - p_low > error) {
    double p_mid = (p_high + p_low) / 2;
    apx_data = g0_data;
    chung_lu_estimate_iterative(apx_data, p_mid, q, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_data), std::numeric_limits<double>::infinity(), get_value(g_data))) {
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
    std::function<double& (DataObject&)> get_value) {
  double q_low = 0.0, q_high = 1.0, error = EPS;
  DataObject apx_data;
  while (q_high - q_low > error) {
    double q_mid = (q_high + q_low) / 2;
    apx_data = g0_data;
    chung_lu_estimate_iterative(apx_data, p, q_mid, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_data), std::numeric_limits<double>::infinity(), get_value(g_data))) {
      q_low = q_mid;
    } else {
      q_high = q_mid;
    }
  }
  Parameters params;
  params.initialize_chung_lu(p, q_low);
  return params;
}

std::vector<Parameters> chung_lu_get_parameters(
    DataObject &g0_data, DataObject &g_data, std::function<double& (DataObject&)> get_value) {
  std::vector<Parameters> S;
  double q_step = Q_STEP;
  for (double q = 0.0; q <= 1.0 + EPS; q += q_step) {
    DataObject apx_min(g0_data), apx_max(g0_data);
    chung_lu_estimate_iterative(apx_min, 0.0, q, g0_data.no_vertices, g_data.no_vertices);
    chung_lu_estimate_iterative(apx_max, 1.0, q, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_min), get_value(apx_max), get_value(g_data))) {
      S.push_back(chung_lu_binary_search_p(g0_data, g_data, q, get_value));
      if (S.size() >= 2 && fabs(S[S.size() - 1].p - S[S.size() - 2].p) < P_DISTANCE) {
        q_step *= Q_EXP;
      }
    }
  }
  return S;
}

void chung_lu_estimate_parameter(
    const std::string &name, DataObject &g_data, DataObject &g0_data, const Graph &G0,
    std::function<double& (DataObject&)> get_value, std::ofstream &out_file,
    bool tolerance_interval) {
  std::vector<Parameters> S = chung_lu_get_parameters(g0_data, g_data, get_value);
  if (tolerance_interval) {
    std::vector<Parameters> S_low, S_high;
    for (auto params : S) {
      DataObject g_data_low, g_data_high;
      std::tie(g_data_low, g_data_high) = get_empirical_interval(g_data, G0, params, get_value);
      S_low.push_back(chung_lu_binary_search_q(g0_data, g_data_low, params.p, get_value));
      S_high.push_back(chung_lu_binary_search_q(g0_data, g_data_high, params.p, get_value));
    }
    print(name, get_value(g0_data), get_value(g_data), S_low, S, S_high, out_file);
  } else {
    print(name, get_value(g0_data), get_value(g_data), S, out_file);
  }
}

void chung_lu_estimate(
    DataObject &g_data, DataObject &g0_data, const Graph &G0, std::ofstream &out_file) {
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
    std::function<double& (DataObject&)> get_value) {
  double p_low = 0.0, p_high = 1.0, error = EPS;
  DataObject apx_data;
  while (p_high - p_low > error) {
    double p_mid = (p_high + p_low) / 2;
    apx_data = g0_data;
    pastor_satorras_estimate_iterative(
        apx_data, p_mid, r, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_data), std::numeric_limits<double>::infinity(), get_value(g_data))) {
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
    std::function<double& (DataObject&)> get_value) {
  double r_low = 0.0, r_high = g0_data.no_vertices, error = EPS;
  DataObject apx_data;
  while (r_high - r_low > error) {
    double r_mid = (r_high + r_low) / 2;
    apx_data = g0_data;
    pastor_satorras_estimate_iterative(
        apx_data, p, r_mid, g0_data.no_vertices, g_data.no_vertices);
    if (contains(get_value(apx_data), std::numeric_limits<double>::infinity(), get_value(g_data))) {
      r_low = r_mid;
    } else {
      r_high = r_mid;
    }
  }
  Parameters params;
  params.initialize_pastor_satorras(p, r_low);
  return params;
}

std::vector<Parameters> pastor_satorras_get_parameters(
    DataObject &g0_data, DataObject &g_data, std::function<double& (DataObject&)> get_value) {
  std::vector<Parameters> S;
  double r_max = g0_data.no_vertices, r_step = R_STEP;
  for (double r = 0.0; r <= r_max + EPS; r += r_step) {
    DataObject apx_min(g0_data), apx_max(g0_data);
    pastor_satorras_estimate_iterative(apx_min, 0.0, r, g0_data.no_vertices, g_data.no_vertices);
    pastor_satorras_estimate_iterative(apx_max, 1.0, r, g0_data.no_vertices, g_data.no_vertices);
    if (g0_data.triangles > apx_min.triangles || g0_data.open_triangles > apx_min.open_triangles) {
      continue;
    }
    if (contains(get_value(apx_min), get_value(apx_max), get_value(g_data))) {
      S.push_back(pastor_satorras_binary_search_p(g0_data, g_data, r, get_value));
      if (S.size() >= 2 && fabs(S[S.size() - 1].p - S[S.size() - 2].p) < P_DISTANCE) {
        r_step *= R_EXP;
      }
    }
  }
  return S;
}

void pastor_satorras_estimate_parameter(
    const std::string &name, DataObject &g_data, DataObject &g0_data, const Graph &G0,
    std::function<double& (DataObject&)> get_value, std::ofstream &out_file,
    bool tolerance_interval) {
  std::vector<Parameters> S = pastor_satorras_get_parameters(g0_data, g_data, get_value);
  if (tolerance_interval) {
    std::vector<Parameters> S_low, S_high;
    for (auto params : S) {
      DataObject g_data_low, g_data_high;
      std::tie(g_data_low, g_data_high) = get_empirical_interval(g_data, G0, params, get_value);
      S_low.push_back(pastor_satorras_binary_search_r(g0_data, g_data_low, params.p, get_value));
      S_high.push_back(pastor_satorras_binary_search_r(g0_data, g_data_high, params.p, get_value));
    }
    print(name, get_value(g0_data), get_value(g_data), S_low, S, S_high, out_file);
  } else {
    print(name, get_value(g0_data), get_value(g_data), S, out_file);
  }
}

void pastor_satorras_estimate(
    DataObject &g_data, DataObject &g0_data, const Graph &G0, std::ofstream &out_file) {
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
    const Graph &G, const Graph &G0, const Mode &mode, std::ofstream &out_file) {
  DataObject g0_data = get_params_for_graph(G0, true), g_data = get_params_for_graph(G, true);
  switch (mode) {
    case CHUNG_LU:
      chung_lu_estimate(g_data, g0_data, G0, out_file);
      break;
    case PASTOR_SATORRAS:
      pastor_satorras_estimate(g_data, g0_data, G0, out_file);
      break;
    default:
      throw std::invalid_argument("Invalid mode: " + LONG_NAME.find(mode)->second);
  }
}

void synthetic_data(const int &n, const int &n0, const double &p0, const Parameters &params) {
  Graph G0 = generate_seed_simple(n0, p0);
  Graph G = G0;
  generate_graph_simple(G, n, params);
  std::ofstream out_file(TEMP_FOLDER + get_synthetic_filename(n, n0, params, "RE"));
  std::cout << "Synthetic data: " + params.to_string() << std::endl;
  process_graph(G, G0, params.mode, out_file);
}

void real_world_data(
    const std::string &graph_name, const std::string &seed_name, const Mode &mode) {
  Graph G = read_graph_simple(FILES_FOLDER + graph_name);
  Graph G0 = read_graph_simple(FILES_FOLDER + seed_name);
  std::ofstream out_file(TEMP_FOLDER + get_real_filename(graph_name, mode, "RE"));
  std::cout << "File: " << graph_name << std::endl;
  process_graph(G, G0, mode, out_file);
}

int main(int argc, char **argv) {
  try {
    Env = prepare_environment(argc, argv);
    TI_TRIES = read_int(Env, "-st:", 1, "Tolerance interval tries");
    std::string action = read_action(Env);
    if (action == "synthetic") {
      const int n = read_n(Env), n0 = read_n0(Env);
      const double p0 = read_p0(Env);
      std::unique_ptr<Parameters> params = read_parameters(Env);
      synthetic_data(n, n0, p0, *params);
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
