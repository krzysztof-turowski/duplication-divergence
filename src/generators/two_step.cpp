#include "two_step.h"

TwoStepParameters::TwoStepParameters(const double &_alpha, const size_t &_target_edges)
    : alpha(_alpha), target_edges(_target_edges) {
  this->mode = Mode::TWO_STEP;
}

std::string TwoStepParameters::to_string() const {
  std::stringstream out;
  out << LONG_NAME.find(this->mode)->second << " ";
  out << "α = " << this->alpha << " ";
  out << "m = " << this->target_edges << " ";
  return out.str();
}

std::string TwoStepParameters::to_filename() const {
  std::stringstream out;
  out << SHORT_NAME.find(this->mode)->second;
  out << "-" << this->alpha;
  out << "-" << this->target_edges;
  return out.str();
}

std::string TwoStepParameters::to_csv() const {
  std::stringstream out;
  out << this->alpha;
  out << ",";
  out << this->target_edges;
  return out.str();
}

static int count_common_neighbors(const Graph &G, const int &i, const int &j) {
  int result = 0;
  for (auto &&neighbor : G[i]) {
    if (G[j].find(neighbor) != G[j].end()) {
      result++;
    }
  }
  return result;
}

Graph generate_two_step_graph(const int &n, const TwoStepParameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  std::uniform_int_distribution<int> vertex_distribution(0, n - 1);

  auto G = Graph(n);
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      G[i].insert(j);
      G[j].insert(i);
    }
  }

  int number_of_edges = n * (n - 1) / 2;

  while (number_of_edges > n) {
    const int i = vertex_distribution(generator);
    if (G[i].size() <= 0) {
      continue;
    }
    std::uniform_int_distribution<int> edge_distribution(0, G[i].size() - 1);
    const int j_index = edge_distribution(generator);
    auto it = G[i].begin();
    std::advance(it, j_index);
    const int j = *it;
    if (G[j].size() <= 1) {
      continue;
    }
    double p_j = std::pow(G[j].size(), -params.alpha);
    double N_i = 0;
    for (auto &&edge : G[i]) {
      N_i += std::pow(G[edge].size(), -params.alpha);
    }

    if (distribution(generator) <= p_j / N_i) {
      G[i].erase(j);
      G[j].erase(i);
      number_of_edges--;
    }
  }

  while (number_of_edges < static_cast<int>(params.target_edges)) {
    int i = vertex_distribution(generator);
    int j = vertex_distribution(generator);
    while (i == j || G[i].find(j) != G[i].end()) {
      i = vertex_distribution(generator);
      j = vertex_distribution(generator);
    }

    int64_t nc = count_common_neighbors(G, i, j);
    double p_ij = G[i].size() * G[j].size() == 0
                      ? 1
                      : static_cast<double>(nc * nc) / (G[i].size() * G[j].size());

    if (distribution(generator) <= p_ij) {
      G[i].insert(j);
      G[j].insert(i);
      number_of_edges++;
    }
  }

  return G;
}
