#pragma once

#include "parameters.h"
#include <random>

class TwoStepParameters : public Parameters {
 public:
  double alpha;
  size_t target_edges;

  TwoStepParameters(const double &_alpha, const size_t &target_edges);

  std::string to_string() const;

  std::string to_filename() const;

  std::string to_csv() const;
};

SimpleGraph generate_two_step_graph(const int &n, const TwoStepParameters &params);
