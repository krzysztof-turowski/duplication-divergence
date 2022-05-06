#pragma once

#include "parameters.h"
#include <random>

extern const int PRECISION_C;

class CopyGraphParameters : public Parameters {
 public:
  int a;
  int b;
  double c;

  CopyGraphParameters(const int &_a, const int &_b, const double &_c);

  std::string to_string() const;

  std::string to_filename() const;

  std::string to_csv() const;
};

void generate_copy_graph(SimpleGraph &G, const int &n, const CopyGraphParameters &params);
