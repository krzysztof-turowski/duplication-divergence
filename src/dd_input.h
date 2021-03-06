#pragma once

#include "./dd_header.h"

#include <Snap.h>

#include <string>

TEnv prepare_environment(int argc, char **argv) {
  return TEnv(argc, argv, TNotify::StdNotify);
}

inline int read_int(
    TEnv &environment, const std::string &prefix, const int &default_value,
    const std::string &description) {
  return environment.GetIfArgPrefixInt(
      TStr(prefix.c_str()), default_value, TStr(description.c_str()));
}

inline double read_double(
    TEnv &environment, const std::string &prefix, const double &default_value,
    const std::string &description) {
  return environment.GetIfArgPrefixFlt(
      TStr(prefix.c_str()), default_value, TStr(description.c_str()));
}

inline std::string read_string(
    TEnv &environment, const std::string &prefix, const std::string &default_value,
    const std::string &description) {
  return environment.GetIfArgPrefixStr(
      TStr(prefix.c_str()), TStr(default_value.c_str()), TStr(description.c_str())).CStr();
}

inline std::string read_action(TEnv &environment) {
  return read_string(
      environment, "-action:", "synthetic", "Actions: {synthetic, real_data}");
}

inline int read_n(TEnv &environment) {
  return read_int(environment, "-n:", 0, "Number of nodes");
}

inline int read_n0(TEnv &environment) {
  return read_int(environment, "-n0:", 0, "Number of nodes in seed graph");
}

inline double read_p0(TEnv &environment) {
  return read_double(environment, "-p0:", 1.0, "Probability p_0 for seed graph");
}

inline std::string read_mode(TEnv &environment) {
  return read_string(
      environment, "-mode:", "pure_duplication",
      "Mode: {pure_duplication, pastor_satorras, chung_lu}");
}

inline Parameters read_parameters(TEnv &environment) {
  Parameters params;
  std::string mode = read_mode(environment);
  double p = read_double(environment, "-p:", 0.0, "Probability p for duplication");
  if (mode == "pure_duplication") {
    params.initialize_pure_duplication(p);
  } else if (mode == "chung_lu") {
    double q = read_double(environment, "-q:", 0.0, "Probability q for connection to parent");
    params.initialize_chung_lu(p, q);
  } else if (mode == "pastor_satorras") {
    double r = read_double(
        environment, "-r:", 0.0, "Probability r/n for connection to all other vertices");
    params.initialize_pastor_satorras(p, r);
  } else {
    throw std::invalid_argument("Invalid mode: " + mode);
  }
  return params;
}

inline std::string read_graph_name(TEnv &environment) {
  return read_string(environment, "-graph:", "", "Real-world graph name");
}
