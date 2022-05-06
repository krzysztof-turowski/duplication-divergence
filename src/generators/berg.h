#pragma once

#include "parameters.h"
#include <random>

class BergParameters : public Parameters {
 public:
  double time_unit;
  double evolution_duration;
  double average_connectivity;
  double duplication_rate;
  double link_addition_rate;
  double link_detachment_rate;

  explicit BergParameters(const double &tu, const double &ed, const double &ac, const double &dr,
      const double &lar, const double &ldr);

  std::string to_string() const;

  std::string to_filename() const;

  std::string to_csv() const;

  std::string short_name() const;
  std::string long_name() const;
};

void generate_berg_graph(SimpleGraph &G, const BergParameters &params);
