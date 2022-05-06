#include "kumar_linear.h"

const int PRECISION_COPY_FACTOR = 3;

KumarLinearParameters::KumarLinearParameters(const int &d, const double &cf)
    : degree(d), copy_factor(cf) {
  this->mode = Mode::KUMAR_LINEAR;
}

std::string KumarLinearParameters::to_string() const {
  std::stringstream out;
  out << this->long_name() << " ";
  out << "d = " << this->degree << " ";
  out << "α = " << std::fixed << std::setprecision(PRECISION_COPY_FACTOR) << this->copy_factor
      << " ";
  return out.str();
}

std::string KumarLinearParameters::to_filename() const {
  std::stringstream out;
  out << this->short_name();
  out << "-" << this->degree;
  out << "-" << std::fixed << std::setprecision(PRECISION_COPY_FACTOR) << this->copy_factor;
  return out.str();
}

std::string KumarLinearParameters::to_csv() const {
  std::stringstream out;
  out << this->degree;
  out << ",";
  out << this->copy_factor;
  return out.str();
}

std::string KumarLinearParameters::short_name() const {
  return "KUMARL";
}

std::string KumarLinearParameters::long_name() const {
  return "Kumar linear copy";
}

void generate_kumar_linear_graph(
    SimpleGraph &G, const int &n, const KumarLinearParameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  const size_t seed_size = G.size();
  G.resize(n);

  for (int i = seed_size; i < n; i++) {
    std::uniform_int_distribution<int> parent_distribution(0, i - 1);

    const auto &d_param = std::min(params.degree, i);
    const int parent = parent_distribution(generator);
    for (auto it = G[parent].begin(); it != G[parent].end(); it++) {
      int new_edge;
      if (edge_distribution(generator) <= params.copy_factor) {
        new_edge = parent_distribution(generator);
      } else {
        new_edge = *it;
      }

      G[i].insert(new_edge);
      G[new_edge].insert(i);
    }
  }
}
