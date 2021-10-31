#include "sticky.h"

StickyParameters::StickyParameters(std::vector<int> &&_degrees) {
  this->mode = Mode::STICKY;
  degrees = std::move(_degrees);
}

std::string StickyParameters::to_string() const {
  std::stringstream out;
  return LONG_NAME.find(this->mode)->second;
}

std::string StickyParameters::to_filename() const {
  std::stringstream out;
  return SHORT_NAME.find(this->mode)->second;
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
      if (distribution(generator) * degree_sum * degree_sum
          <= static_cast<double>(degrees[i]) * degrees[j]) {
        G[i].insert(j), G[j].insert(i);
      }
    }
  }

  return G;
}
