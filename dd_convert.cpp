#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

set<int> read_all_graph_labels(ifstream &in_G_file) {
  set<int> vertices;
  int u, v;
  while (in_G_file >> u >> v) {
    vertices.insert(u), vertices.insert(v);
  }
  return vertices;
}

multimap<unsigned, int> read_age_to_id_map(
    ifstream &in_age_file, const set<int> &vertices) {
  multimap<unsigned, int> age_to_id;
  int v, age;
  while (in_age_file >> v >> age) {
    if (vertices.count(v)) {
      age_to_id.insert(make_pair(age, v));
    }
  }
  return age_to_id;
}

tuple<vector<int>, vector<unsigned>> get_vertices_age_order(
    multimap<unsigned, int> &age_to_id, const int &n_max) {
  vector<int> V;
  vector<unsigned> A;
  unsigned age_max = numeric_limits<unsigned>::max();
  for (const auto &age_v : age_to_id) {
    if (age_max < age_v.first) {
      break;
    }
    V.push_back(age_v.second), A.push_back(age_v.first);
    if (V.size() == n_max) {
      age_max = A.back();
    }
  }
  return make_tuple(V, A);
}

void print_graph(
    const string &out_name, ifstream &in_G_file, const vector<int> &V,
    const vector<unsigned> &A, const int &age_zero) {
  ofstream out_G_file("G-" + out_name), out_G0_file("G0-" + out_name);
  int n = V.size(), u, v;
  for (int i = 0; i < n; i++) {
    out_G_file << i << " " << i << endl;
  }
  for (int i = 0; A[i] <= age_zero; i++) {
    out_G0_file << i << " " << i << endl;
  }
  while (in_G_file >> u >> v) {
    int u_new = find(V.begin(), V.end(), u) - V.begin();
    int v_new = find(V.begin(), V.end(), v) - V.begin();
    if (u_new >= n || v_new >= n) {
      continue;
    }
    if (u_new > v_new) {
      swap(u_new, v_new);
    }
    if (A[u_new] <= age_zero && A[v_new] <= age_zero) {
      out_G0_file << u_new << " " << v_new << endl;
    }
    out_G_file << u_new << " " << v_new << endl;
  }
  out_G_file.close(), out_G0_file.close();
}

void print_age(
    const string &out_name, const vector<unsigned> &A, const int &age_zero) {
  ofstream out_age_file("PH-" + out_name);
  int n = A.size();
  for (int i = 0; i < n; i++) {
    int age = A[i];
    out_age_file << i << " " << (age != -1 ? max(age - age_zero, 0) : -1) << endl;
  }
  out_age_file.close();
}

void convert_edge_list(
    const string &graph_name, const string &age_name, const string &out_name,
    const int &n_max = numeric_limits<int>::max(), const double &n0_fraction = 0.0) {
  ifstream in_G_file(graph_name), in_age_file(age_name);

  set<int> vertices = read_all_graph_labels(in_G_file);
  multimap<unsigned, int> age_to_id = read_age_to_id_map(in_age_file, vertices);
  vector<int> V;
  vector<unsigned> A;
  tie(V, A) = get_vertices_age_order(age_to_id, n_max);

  int age_zero = A[n0_fraction * V.size()];
  in_G_file.clear(), in_G_file.seekg(0);
  print_graph(out_name, in_G_file, V, A, age_zero);
  in_G_file.close();

  print_age(out_name, A, age_zero);
  in_age_file.close();
}

int main(int, char **argv) {
  convert_edge_list(argv[1], argv[2], argv[3], 10000, 0.1);

  return 0;
}
