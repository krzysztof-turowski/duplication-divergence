#pragma once

#include "generators/parameters.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>

const std::string FILES_FOLDER = "files/", TEMP_FOLDER = "temp/";

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

SimpleGraph read_graph_simple(const std::string &graph_name) {
  std::ifstream graph_file(graph_name);
  if (graph_file.fail()) {
    throw std::invalid_argument("Missing " + graph_name + " file");
  }
  SimpleGraph G;
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
