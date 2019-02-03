#pragma once

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <sstream>
#include <vector>

const std::string FILES_FOLDER = "files/", TEMP_FOLDER = "temp/";

enum Mode { PURE_DUPLICATION, PURE_DUPLICATION_CONNECTED, CHUNG_LU, PASTOR_SATORRAS };

const std::map<Mode, std::string> SHORT_NAME = {
  { Mode::PURE_DUPLICATION, "PD" },
  { Mode::PURE_DUPLICATION_CONNECTED, "PDC" },
  { Mode::CHUNG_LU, "CL" },
  { Mode::PASTOR_SATORRAS, "PS" },
};

const std::map<Mode, std::string> LONG_NAME = {
  { Mode::PURE_DUPLICATION, "Pure duplication" },
  { Mode::PURE_DUPLICATION_CONNECTED, "Pure duplication without isolated vertices" },
  { Mode::CHUNG_LU, "Chung-Lu" },
  { Mode::PASTOR_SATORRAS, "Pastor-Satorras" },
};

const std::map<std::string, Mode> REVERSE_NAME = {
  { "pure_duplication", Mode::PURE_DUPLICATION },
  { "pure_duplication_connected", Mode::PURE_DUPLICATION_CONNECTED },
  { "chung_lu", Mode::CHUNG_LU },
  { "pastor_satorras", Mode::PASTOR_SATORRAS },
};

class Parameters {
public:
  Mode mode;
  double p, q, r;

  Parameters() : p(nan("")), q(nan("")), r(nan("")) { }

  void initialize(const std::string &mode_v, char *argv[]) {
    if (mode_v == "chung_lu") {
      initialize_chung_lu(std::stod(argv[0]), std::stod(argv[1]));
    }
    else if (mode_v == "pastor_satorras") {
      initialize_pastor_satorras(std::stod(argv[0]), std::stod(argv[1]));
    }
    else {
      assert(0);
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

  std::string to_string() const {
    std::stringstream out;
    out << LONG_NAME.find(this->mode)->second << " ";
    out << "p = " << this->p << " ";
    if (!std::isnan(this->q)) {
      out << "q = " << this->q << " ";
    }
    if (!std::isnan(this->r)) {
      out << "r = " << this->r << " ";
    }
    return out.str();
  }

  std::string to_string(const Parameters &low, const Parameters &high) const {
    std::stringstream out;
    out << LONG_NAME.find(this->mode)->second << " ";
    out << "p_min = " << low.p << " " << "p = " << this->p << " " << "p_max = " << high.p << " ";
    if (!std::isnan(this->q)) {
      out << "q = " << this->q << " ";
    }
    if (!std::isnan(this->r)) {
      out << "r = " << this->r << " ";
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
      out  << this->r;
    }
    return out.str();
  }
};

std::vector<std::set<int>> generate_seed(const int &n0, const double &p0) {
  std::vector<std::set<int>> G(n0);
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
 
std::vector<std::set<int>> generate_graph(std::vector<std::set<int>> &G, const int &n, const Parameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

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
    }
    else if (params.mode == Mode::PURE_DUPLICATION_CONNECTED) {
      while(true) {
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
    }
    else if (params.mode == Mode::CHUNG_LU) {
      for (auto j : G[parent]) {
        if (edge_distribution(generator) <= params.p) {
          G[i].insert(j), G[j].insert(i);
        }
      }
      if (edge_distribution(generator) <= params.q) {
        G[i].insert(parent), G[parent].insert(i);
      }
    }
    else if (params.mode == Mode::PASTOR_SATORRAS) {
      for (int j = 0; j < i; j++) {
        if (G[parent].count(j)) {
          if (edge_distribution(generator) <= params.p) {
            G[i].insert(j), G[j].insert(i);
          }
        }
        else {
          if (edge_distribution(generator) <= params.r / i) {
            G[i].insert(j), G[j].insert(i);
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

std::vector<std::set<int>> read_graph(const std::string &graph_name) {
  std::ifstream graph_file(graph_name);
  if (graph_file.fail()) {
    throw std::invalid_argument("Missing " + graph_name + " file");
  }
  std::vector<std::set<int>> G;
  int u, v;
  while (!graph_file.eof())
  {
    graph_file >> u >> v;
    if (v >= static_cast<int>(G.size())) {
      G.resize(v + 1);
    }
    if (u != v) {
      G[u].insert(v), G[v].insert(u);
    }
  }
  graph_file.close();
  return G;
}