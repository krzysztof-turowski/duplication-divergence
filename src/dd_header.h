#pragma once

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using Graph = std::vector<std::set<unsigned>>;

const std::string FILES_FOLDER = "files/", TEMP_FOLDER = "temp/";
const int PRECISION_P = 3, PRECISION_Q = 3, PRECISION_R = 2, WIDTH_R = 6;

enum Mode {
  INVALID,
  PURE_DUPLICATION,
  PURE_DUPLICATION_CONNECTED,
  CHUNG_LU,
  PASTOR_SATORRAS,
  STICKY,
  BA,
  COPY_GRAPH
};

const std::map<Mode, std::string> SHORT_NAME = {
  { Mode::PURE_DUPLICATION, "PD" },
  { Mode::PURE_DUPLICATION_CONNECTED, "PDC" },
  { Mode::CHUNG_LU, "CL" },
  { Mode::PASTOR_SATORRAS, "PS" },
  { Mode::STICKY, "STICKY" },
  { Mode::BA, "BA" },
  { Mode::COPY_GRAPH, "COPY" },
};

const std::map<Mode, std::string> LONG_NAME = {
  { Mode::PURE_DUPLICATION, "Pure duplication" },
  { Mode::PURE_DUPLICATION_CONNECTED, "Pure duplication without isolated vertices" },
  { Mode::CHUNG_LU, "Chung-Lu" },
  { Mode::PASTOR_SATORRAS, "Pastor-Satorras" },
  { Mode::STICKY, "STICKY" },
  { Mode::BA, "Barabasi-Albert" },
  { Mode::COPY_GRAPH, "Copy graph" },
};

const std::map<std::string, Mode> REVERSE_NAME = {
  { "pure_duplication", Mode::PURE_DUPLICATION },
  { "pure_duplication_connected", Mode::PURE_DUPLICATION_CONNECTED },
  { "chung_lu", Mode::CHUNG_LU },
  { "pastor_satorras", Mode::PASTOR_SATORRAS },
  { "sticky", Mode::STICKY },
  { "ba", Mode::BA },
  { "copy_graph", Mode::COPY_GRAPH },
};

class Parameters {
 public:
  Mode mode;
  double p, q, r;
  std::vector<int> degrees;
  int m;

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
    this->m = -1;
  }

  void initialize_pure_duplication_connected(const double &p_v) {
    this->mode = Mode::PURE_DUPLICATION_CONNECTED;
    this->p = p_v;
    this->q = this->r = nan("");
    this->m = -1;
  }

  void initialize_chung_lu(const double &p_v, const double &q_v) {
    this->mode = Mode::CHUNG_LU;
    this->p = p_v, this->q = q_v;
    this->r = nan("");
    this->m = -1;
  }

  void initialize_pastor_satorras(const double &p_v, const double &r_v) {
    this->mode = Mode::PASTOR_SATORRAS;
    this->p = p_v, this->r = r_v;
    this->q = nan("");
    this->m = -1;
  }

  void initialize_sticky(std::vector<int> &&_degrees) {
    this->mode = Mode::STICKY;
    degrees = std::move(_degrees);

    this->p = this->r = this->q = nan("");
    this->m = -1;
  }

  void initialize_ba(int const &_m) {
    this->mode = Mode::BA;
    this->p = this->r = this->q = nan("");
    this->m = _m;
  }

  virtual std::string to_string() const {
    std::stringstream out;
    out << LONG_NAME.find(this->mode)->second << " ";
    if (!std::isnan(this->p)) {
      out << "p = " << std::fixed << std::setprecision(PRECISION_P) << this->p << " ";
    }
    if (!std::isnan(this->q)) {
      out << "q = " << std::fixed << std::setprecision(PRECISION_Q) << this->q << " ";
    }
    if (!std::isnan(this->r)) {
      out << "r = " << std::fixed << std::setw(WIDTH_R) << std::setprecision(PRECISION_R) << this->r
          << " ";
    }
    if (this->m != -1) {
      out << "m = " << this->m << " ";
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
      out << "r = " << std::fixed << std::setw(WIDTH_R) << std::setprecision(PRECISION_R) << this->r
          << " ";
    }
    return out.str();
  }

  virtual std::string to_filename() const {
    std::stringstream out;
    out << SHORT_NAME.find(this->mode)->second;
    if (!std::isnan(this->p)) {
      out << "-" << std::fixed << std::setprecision(PRECISION_P) << this->p;
    }
    if (!std::isnan(this->q)) {
      out << "-" << std::fixed << std::setprecision(PRECISION_Q) << this->q;
    }
    if (!std::isnan(this->r)) {
      out << "-" << std::fixed << std::setprecision(PRECISION_R) << this->r;
    }
    if (this->m != -1) {
      out << "-" << this->m << " ";
    }
    return out.str();
  }

  virtual std::string to_csv() const {
    std::stringstream out;
    if (!std::isnan(this->p)) {
      out << this->p;
    }
    out << ",";
    if (!std::isnan(this->q)) {
      out << this->q;
    }
    out << ",";
    if (!std::isnan(this->r)) {
      out << this->r;
    }
    out << ",";
    if (!std::isnan(this->m)) {
      out << this->m;
    }
    return out.str();
  }
};

inline std::string get_synthetic_filename(
    const int &n, const int &n0, const Parameters &params, const std::string &suffix) {
  return "synthetic-" + std::to_string(n) + "-" + std::to_string(n0) + "-" + params.to_filename()
         + (suffix.length() > 0 ? "-" + suffix : "") + ".txt";
}

inline std::string get_real_filename(
    const std::string &graph_name, const Mode &mode, const std::string &suffix) {
  return graph_name.substr(0, graph_name.find_last_of(".")) + "-" + SHORT_NAME.find(mode)->second
         + (suffix.length() > 0 ? "-" + suffix : "") + ".txt";
}

inline std::string get_seed_name(const std::string &graph_name) {
  return std::regex_replace(graph_name, std::regex("^G"), "G0");
}

Graph read_graph_simple(const std::string &graph_name) {
  std::ifstream graph_file(graph_name);
  if (graph_file.fail()) {
    throw std::invalid_argument("Missing " + graph_name + " file");
  }
  Graph G;
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
