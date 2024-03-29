#pragma once

#include "../graph/graph.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

extern const int PRECISION_P, PRECISION_Q, PRECISION_R, WIDTH_R;

enum Mode {
  INVALID,
  PURE_DUPLICATION,
  PURE_DUPLICATION_CONNECTED,
  CHUNG_LU,
  PASTOR_SATORRAS,
  STICKY,
  BA,
  COPY_GRAPH,
  TWO_STEP,
  BERG,
  KUMAR_LINEAR,
};

extern const std::map<Mode, std::string> SHORT_NAME;

extern const std::map<Mode, std::string> LONG_NAME;

extern const std::map<std::string, Mode> REVERSE_NAME;

class Parameters {
 public:
  Mode mode;
  double p, q, r;

  Parameters();

  void initialize(const std::string &mode_v, char *argv[]);

  void initialize_pure_duplication(const double &p_v);

  void initialize_pure_duplication_connected(const double &p_v);

  void initialize_chung_lu(const double &p_v, const double &q_v);

  void initialize_pastor_satorras(const double &p_v, const double &r_v);

  virtual std::string to_string() const;

  std::string to_string(const Parameters &low, const Parameters &high) const;

  virtual std::string to_filename() const;

  virtual std::string to_csv() const;

  virtual std::string short_name() const;

  virtual std::string long_name() const;
};
