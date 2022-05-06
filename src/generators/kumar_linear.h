#pragma once

#include "parameters.h"
#include <random>

extern const int PRECISION_C;

class KumarLinearParameters : public Parameters {
 public:
  int degree;
  double copy_factor;

  KumarLinearParameters(const int &d, const double &cf);

  std::string to_string() const;

  std::string to_filename() const;

  std::string to_csv() const;

  std::string short_name() const;
  std::string long_name() const;
};

void generate_kumar_linear_graph(SimpleGraph &G, const int &n, const KumarLinearParameters &params);
