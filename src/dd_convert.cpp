// Tool for graph conversion from popular formats.
// Compile: g++ dd_convert.cpp -O3 -o ./dd_convert
// Run: ./dd_convert -graph:FILENAME -age:FILENAME -out:FILE_SUFFIX

#include "./dd_input.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

std::set<int> read_all_graph_labels(std::ifstream &in_G_file) {
  std::set<int> vertices;
  int u, v;
  while (in_G_file >> u >> v) {
    vertices.insert(u), vertices.insert(v);
  }
  return vertices;
}

std::multimap<unsigned, int> read_age_to_id_map(
    std::ifstream &in_age_file, const std::set<int> &vertices) {
  std::multimap<unsigned, int> age_to_id;
  int v, age;
  while (in_age_file >> v >> age) {
    if (vertices.count(v)) {
      age_to_id.insert(std::make_pair(age, v));
    }
  }
  return age_to_id;
}

std::tuple<std::vector<int>, std::vector<unsigned>> get_vertices_age_order(
    std::multimap<unsigned, int> &age_to_id, const size_t &n_max) {
  std::vector<int> V;
  std::vector<unsigned> A;
  unsigned age_max = std::numeric_limits<unsigned>::max();
  for (const auto &age_v : age_to_id) {
    if (age_max < age_v.first) {
      break;
    }
    V.push_back(age_v.second), A.push_back(age_v.first);
    if (V.size() == n_max) {
      age_max = A.back();
    }
  }
  return std::make_tuple(V, A);
}

void print_graph(
    const std::string &out_name, std::ifstream &in_G_file, const std::vector<int> &V,
    const std::vector<unsigned> &A, const unsigned &age_zero) {
  std::ofstream out_G_file("G-" + out_name), out_G0_file("G0-" + out_name);
  int n = V.size(), u, v;
  for (int i = 0; i < n; i++) {
    out_G_file << i << " " << i << std::endl;
  }
  for (size_t i = 0; A[i] <= age_zero; i++) {
    out_G0_file << i << " " << i << std::endl;
  }
  while (in_G_file >> u >> v) {
    int u_new = std::find(V.begin(), V.end(), u) - V.begin();
    int v_new = std::find(V.begin(), V.end(), v) - V.begin();
    if (u_new >= n || v_new >= n) {
      continue;
    }
    if (u_new > v_new) {
      std::swap(u_new, v_new);
    }
    if (A[u_new] <= age_zero && A[v_new] <= age_zero) {
      out_G0_file << u_new << " " << v_new << std::endl;
    }
    out_G_file << u_new << " " << v_new << std::endl;
  }
  out_G_file.close(), out_G0_file.close();
}

void print_age(
    const std::string &out_name, const std::vector<unsigned> &A, const int &age_zero) {
  std::ofstream out_age_file("PH-" + out_name);
  int n = A.size();
  for (int i = 0; i < n; i++) {
    int age = A[i];
    out_age_file << i << " " << (age != -1 ? std::max(age - age_zero, 0) : -1) << std::endl;
  }
  out_age_file.close();
}

void convert_edge_list(
    const std::string &graph_name, const std::string &age_name, const std::string &out_name,
    const size_t &n_max = std::numeric_limits<int>::max(), const double &n0_fraction = 0.0) {
  std::ifstream in_G_file(graph_name), in_age_file(age_name);

  std::set<int> vertices = read_all_graph_labels(in_G_file);
  std::multimap<unsigned, int> age_to_id = read_age_to_id_map(in_age_file, vertices);
  std::vector<int> V;
  std::vector<unsigned> A;
  tie(V, A) = get_vertices_age_order(age_to_id, n_max);

  int age_zero = A[n0_fraction * V.size()];
  in_G_file.clear(), in_G_file.seekg(0);
  print_graph(out_name, in_G_file, V, A, age_zero);
  in_G_file.close();

  print_age(out_name, A, age_zero);
  in_age_file.close();
}

int main(int argc, char **argv) {
  Env = prepare_environment(argc, argv);
  std::string graph_name = read_string(Env, "-graph:", "", "Graph data filename");
  std::string age_name = read_string(Env, "-age:", "", "Age data filename");
  std::string out_name = read_string(Env, "-out:", "", "Output filename suffix");
  convert_edge_list(graph_name, age_name, out_name);

  return 0;
}
