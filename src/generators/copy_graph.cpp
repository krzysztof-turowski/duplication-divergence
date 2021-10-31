#include "copy_graph.h"

const int PRECISION_C = 3;

CopyGraphParameters::CopyGraphParameters(const int &_a, const int &_b, const double &_c)
    : a(_a), b(_b), c(_c) {
  this->mode = Mode::COPY_GRAPH;
}

std::string CopyGraphParameters::to_string() const {
  std::stringstream out;
  out << LONG_NAME.find(this->mode)->second << " ";
  out << "a = " << this->a << " ";
  out << "b = " << this->b << " ";
  out << "c = " << std::fixed << std::setprecision(PRECISION_C) << this->c << " ";
  return out.str();
}

std::string CopyGraphParameters::to_filename() const {
  std::stringstream out;
  out << SHORT_NAME.find(this->mode)->second;
  out << "-" << this->a;
  out << "-" << this->b;
  out << "-" << std::fixed << std::setprecision(PRECISION_C) << this->c;
  return out.str();
}

std::string CopyGraphParameters::to_csv() const {
  std::stringstream out;
  out << this->a;
  out << ",";
  out << this->b;
  out << ",";
  out << this->c;
  return out.str();
}

void generate_copy_graph(Graph &G, const int &n, const CopyGraphParameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  const size_t seed_size = G.size();
  G.resize(n);

  for (int i = seed_size; i < n; i++) {
    std::uniform_int_distribution<int> parent_distribution(0, i - 1);

    const auto &a_param = std::min(params.a, i);
    while (G[i].size() < static_cast<size_t>(a_param)) {
      const int new_edge = parent_distribution(generator);
      if (G[i].find(new_edge) != G[i].end()) {
        continue;
      }
      G[i].insert(new_edge);
      G[new_edge].insert(i);
    }

    const auto &b_param = std::min(params.b, i - 1);
    std::set<int> copy_sources;
    while (copy_sources.size() < static_cast<size_t>(b_param)) {
      const int new_copy_source = parent_distribution(generator);
      if (copy_sources.find(new_copy_source) != copy_sources.end()) {
        continue;
      }
      copy_sources.insert(new_copy_source);

      const auto &c_param = params.c;
      for (auto &&new_edge : G[new_copy_source]) {
        if (edge_distribution(generator) <= c_param && G[i].find(new_edge) == G[i].end()) {
          G[i].insert(new_edge);
          G[new_edge].insert(i);
        }
      }
    }
  }
}
