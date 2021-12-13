#include "sticky.h"

StickyParameters::StickyParameters(std::vector<int> &&_degrees) {
  this->mode = Mode::STICKY;
  degrees = std::move(_degrees);
}

StickyParameters::StickyParameters(const size_t &n, const double &gamma) {
  this->mode = Mode::STICKY;

  std::vector<double> weights(n);
  for (size_t i = 0; i < weights.size(); i++) {
    weights[i] = i == 0 ? 0 : std::pow(static_cast<double>(i), -gamma);
  }

  std::random_device device;
  std::mt19937 generator(device());
  std::discrete_distribution<int> degree_distribution(weights.begin(), weights.end());

  degrees.resize(n);
  for (size_t i = 0; i < degrees.size(); i++) {
    degrees[i] = degree_distribution(generator);
  }
}

std::string StickyParameters::to_string() const {
  std::stringstream out;
  return this->long_name();
}

std::string StickyParameters::to_filename() const {
  std::stringstream out;
  return this->short_name();
}

std::string StickyParameters::to_csv() const {
  return "";
}

Graph generate_sticky_graph(const StickyParameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  const auto &degrees = params.degrees;
  Graph G(degrees.size());
  int64_t degree_sum = 0;
  for (auto &&deg : params.degrees) {
    degree_sum += deg;
  }

  for (size_t i = 0; i < G.size(); i++) {
    for (size_t j = 0; j < G.size(); j++) {
      if (distribution(generator) * degree_sum <= static_cast<double>(degrees[i]) * degrees[j]) {
        G[i].insert(j), G[j].insert(i);
      }
    }
  }

  return G;
}
