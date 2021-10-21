#pragma once

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <regex>
#include <set>
#include <string>
#include <sstream>
#include <vector>

const std::string FILES_FOLDER = "files/", TEMP_FOLDER = "temp/";
const int PRECISION_P = 3, PRECISION_Q = 3, PRECISION_R = 2, WIDTH_R = 6;

enum Mode { INVALID, PURE_DUPLICATION, PURE_DUPLICATION_CONNECTED, CHUNG_LU, PASTOR_SATORRAS, STICKY };

const std::map<Mode, std::string> SHORT_NAME = {
  { Mode::PURE_DUPLICATION, "PD" },
  { Mode::PURE_DUPLICATION_CONNECTED, "PDC" },
  { Mode::CHUNG_LU, "CL" },
  { Mode::PASTOR_SATORRAS, "PS" },
  { Mode::STICKY, "STICKY" },
};

const std::map<Mode, std::string> LONG_NAME = {
  { Mode::PURE_DUPLICATION, "Pure duplication" },
  { Mode::PURE_DUPLICATION_CONNECTED, "Pure duplication without isolated vertices" },
  { Mode::CHUNG_LU, "Chung-Lu" },
  { Mode::PASTOR_SATORRAS, "Pastor-Satorras" },
  { Mode::STICKY, "STICKY" },
};

const std::map<std::string, Mode> REVERSE_NAME = {
  { "pure_duplication", Mode::PURE_DUPLICATION },
  { "pure_duplication_connected", Mode::PURE_DUPLICATION_CONNECTED },
  { "chung_lu", Mode::CHUNG_LU },
  { "pastor_satorras", Mode::PASTOR_SATORRAS },
  { "sticky", Mode::STICKY },
};

class Parameters {
 public:
  Mode mode;
  double p, q, r;
  std::vector<int> degrees;

  Parameters() : mode(Mode::INVALID), p(nan("")), q(nan("")), r(nan("")) { }

  void initialize(const std::string &mode_v, char *argv[]) {
    if (mode_v == "pure_duplication") {
      initialize_pure_duplication(std::stod(argv[0]));
    } else if (mode_v == "chung_lu") {
      initialize_chung_lu(std::stod(argv[0]), std::stod(argv[1]));
    } else if (mode_v == "pastor_satorras") {
      initialize_pastor_satorras(std::stod(argv[0]), std::stod(argv[1]));
    } else {
      throw std::invalid_argument("Invalid mode: " + mode_v);
    }
  }

  void initialize_pure_duplication(const double &p_v) {
    this->mode = Mode::PURE_DUPLICATION;
    this->p = p_v;
    this->q = this->r = nan("");
  }

  void initialize_pure_duplication_connected(const double &p_v) {
    this->mode = Mode::PURE_DUPLICATION_CONNECTED;
    this->p = p_v;
    this->q = this->r = nan("");
  }

  void initialize_chung_lu(const double &p_v, const double &q_v) {
    this->mode = Mode::CHUNG_LU;
    this->p = p_v, this->q = q_v;
    this->r = nan("");
  }

  void initialize_pastor_satorras(const double &p_v, const double &r_v) {
    this->mode = Mode::PASTOR_SATORRAS;
    this->p = p_v, this->r = r_v;
    this->q = nan("");
  }

  void initialize_sticky(std::vector<int> &&_degrees) {
    this->mode = Mode::STICKY;
    degrees = _degrees;

    this->p = this->r = this->q = nan("");
  }

  std::string to_string() const {
    std::stringstream out;
    out << LONG_NAME.find(this->mode)->second << " ";
    out << "p = " << std::fixed << std::setprecision(PRECISION_P) << this->p << " ";
    if (!std::isnan(this->q)) {
      out << "q = " << std::fixed << std::setprecision(PRECISION_Q) << this->q << " ";
    }
    if (!std::isnan(this->r)) {
      out << "r = " << std::fixed << std::setw(WIDTH_R) << std::setprecision(PRECISION_R)
          << this->r << " ";
    }
    return out.str();
  }

  std::string to_string(const Parameters &low, const Parameters &high) const {
    std::stringstream out;
    out << LONG_NAME.find(this->mode)->second << " ";
    out << "p_min = " << std::fixed << std::setprecision(PRECISION_P) << low.p << " "
        << "p = " << std::fixed << std::setprecision(PRECISION_P) << this->p << " "
        << "p_max = " << std::fixed << std::setprecision(PRECISION_P) << high.p << " ";
    if (!std::isnan(this->q)) {
      out << "q = " << std::fixed << std::setprecision(PRECISION_Q) << this->q << " ";
    }
    if (!std::isnan(this->r)) {
      out << "r = " << std::fixed << std::setw(WIDTH_R) << std::setprecision(PRECISION_R)
          << this->r << " ";
    }
    return out.str();
  }

  std::string to_filename() const {
    std::stringstream out;
    out << SHORT_NAME.find(this->mode)->second << "-";
    out << std::fixed << std::setprecision(PRECISION_P) << this->p;
    if (!std::isnan(this->q)) {
      out << "-" << std::fixed << std::setprecision(PRECISION_Q) << this->q;
    }
    if (!std::isnan(this->r)) {
      out << "-" << std::fixed << std::setprecision(PRECISION_R) << this->r;
    }
    return out.str();
  }

  std::string to_csv() const {
    std::stringstream out;
    out << this->p << ",";
    if (!std::isnan(this->q)) {
      out << this->q;
    }
    out << ",";
    if (!std::isnan(this->r)) {
      out << this->r;
    }
    return out.str();
  }
};

inline std::string get_synthetic_filename(
    const int &n, const int &n0, const Parameters &params, const std::string &suffix) {
  return "synthetic-" + std::to_string(n) + "-"  + std::to_string(n0)
      + "-" + params.to_filename() + (suffix.length() > 0 ? "-" + suffix : "") + ".txt";
}

inline std::string get_real_filename(
    const std::string &graph_name, const Mode &mode, const std::string &suffix) {
  return graph_name.substr(0, graph_name.find_last_of(".")) + "-"
      + SHORT_NAME.find(mode)->second + (suffix.length() > 0 ? "-" + suffix : "") + ".txt";
}

inline std::string get_seed_name(const std::string &graph_name) {
  return std::regex_replace(graph_name, std::regex("^G"), "G0");
}

std::vector<std::set<unsigned>> generate_sticky_graph(const Parameters &params, std::uniform_real_distribution<double> distribution) {
  return std::vector<std::set<unsigned>>();
}

std::vector<std::set<unsigned>> generate_seed_simple(const int &n0, const double &p0) {
  std::vector<std::set<unsigned>> G(n0);
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  for (int i = 0; i < n0; i++) {
    for (int j = i + 1; j < n0; j++) {
      if (edge_distribution(generator) <= p0) {
        G[i].insert(j), G[j].insert(i);
      }
    }
  }
  return G;
}

std::vector<std::set<unsigned>> generate_graph_simple(
    std::vector<std::set<unsigned>> &G, const int &n, const Parameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  if (params.mode == Mode::STICKY) {
      return generate_sticky_graph(params, edge_distribution);
  }

  for (int i = G.size(); i < n; i++) {
    std::uniform_int_distribution<int> parent_distribution(0, i - 1);
    int parent = parent_distribution(generator);
    G.resize(i + 1);
    if (params.mode == Mode::PURE_DUPLICATION) {
      for (auto j : G[parent]) {
        if (edge_distribution(generator) <= params.p) {
          G[i].insert(j), G[j].insert(i);
        }
      }
    } else if (params.mode == Mode::PURE_DUPLICATION_CONNECTED) {
      while (true) {
        for (auto j : G[parent]) {
          if (edge_distribution(generator) <= params.p) {
            G[i].insert(j), G[j].insert(i);
          }
        }
        if (!G[i].empty()) {
          break;
        }
        parent = parent_distribution(generator);
      }
    } else if (params.mode == Mode::CHUNG_LU) {
      for (auto j : G[parent]) {
        if (edge_distribution(generator) <= params.p) {
          G[i].insert(j), G[j].insert(i);
        }
      }
      if (edge_distribution(generator) <= params.q) {
        G[i].insert(parent), G[parent].insert(i);
      }
    } else if (params.mode == Mode::PASTOR_SATORRAS) {
      for (int j = 0; j < i; j++) {
        if (G[parent].count(j)) {
          if (edge_distribution(generator) <= params.p) {
            G[i].insert(j), G[j].insert(i);
          }
        } else {
          if (edge_distribution(generator) <= params.r / i) {
            G[i].insert(j), G[j].insert(i);
          }
        }
      }
    } else {
      throw std::invalid_argument("Invalid mode: " + LONG_NAME.find(params.mode)->second);
    }
  }
  return G;
}

std::vector<std::set<unsigned>> read_graph_simple(const std::string &graph_name) {
  std::ifstream graph_file(graph_name);
  if (graph_file.fail()) {
    throw std::invalid_argument("Missing " + graph_name + " file");
  }
  std::vector<std::set<unsigned>> G;
  unsigned u, v;
  while (graph_file >> u >> v) {
    if (v >= G.size()) {
      G.resize(v + 1);
    }
    if (u != v) {
      G[u].insert(v), G[v].insert(u);
    }
  }
  graph_file.close();
  return G;
}

int read_graph_size(const std::string &graph_name) {
  std::ifstream graph_file(graph_name);
  if (graph_file.fail()) {
    throw std::invalid_argument("Missing " + graph_name + " file");
  }
  int n = 0, u, v;
  while (graph_file >> u >> v) {
    n = std::max(n, v + 1);
  }
  graph_file.close();
  return n;
}
