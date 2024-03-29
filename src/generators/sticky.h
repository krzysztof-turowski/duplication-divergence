#pragma once

#include "parameters.h"
#include <random>

class StickyParameters : public Parameters {
 public:
  std::vector<int> degrees;
  double gamma;

  explicit StickyParameters(std::vector<int> &&_degrees);
  StickyParameters(const size_t &n, const double &_gamma);

  std::string to_string() const;

  std::string to_filename() const;

  std::string to_csv() const;
};

SimpleGraph generate_sticky_graph(const StickyParameters &params);
