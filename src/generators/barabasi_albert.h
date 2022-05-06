#pragma once

#include "parameters.h"
#include <random>

class BarabasiAlbertParameters : public Parameters {
 public:
  int m;

  explicit BarabasiAlbertParameters(const int &_m);

  std::string to_string() const;

  std::string to_filename() const;

  std::string to_csv() const;
};

void generate_ba_graph(SimpleGraph &G, const int &n, const BarabasiAlbertParameters &params);
