// Header for LP computation of the temporal upper bound using Gurobi.

#pragma once

#include <gurobi_c++.h>

#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

std::string LP_row_name(const std::string &prefix, const std::initializer_list<int> &vertices) {
  std::ostringstream out;
  out << prefix;
  if (vertices.size() > 0) {
    out << "_{";
    for (auto v : vertices) {
      out << std::to_string(v) << ",";
    }
    out.seekp(-1, std::ios_base::end), out << "}";
  }
  return out.str();
}

inline std::string LP_name_variable(const int &u, const int &v) {
  return "x_{" + std::to_string(u) + "," + std::to_string(v) + "}";
}

inline int LP_get_variable_index(const int &u, const int &v, const int &n, const int &n0) {
  return (u - n0) * (n - n0) + (v - n0);
}

double LP_solve(
    const std::map<std::pair<int, int>, double> &p_uv, const int &n, const int &n0,
    const double &epsilon) {
  GRBEnv* environment = new GRBEnv();
  GRBModel *LP = new GRBModel(*environment);
  LP->set(GRB_StringAttr_ModelName, "Solve " + std::to_string(epsilon));
  LP->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
  double density = epsilon * (n - n0) * (n - n0 - 1) / 2;

  // Objective function
  std::vector<GRBVar> vars((n - n0) * (n - n0) + 1);
  int s_index = (n - n0) * (n - n0);
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      auto index = LP_get_variable_index(i, j, n, n0);
      if (i != j) {
        vars[index] =
            LP->addVar(
                0.0, 1.0, p_uv.find(std::make_pair(i, j))->second,
                GRB_CONTINUOUS, LP_name_variable(i, j));
      } else {
        vars[index] = LP->addVar(0.0, 0.0, 0.0, GRB_CONTINUOUS, LP_name_variable(i, j));
      }
    }
  }
  vars[s_index] = LP->addVar(0.0, 1 / density, 0.0, GRB_CONTINUOUS, "s");

  // Antisymmetry
  for (int i = n0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      GRBLinExpr row =
          vars[LP_get_variable_index(i, j, n, n0)] + vars[LP_get_variable_index(j, i, n, n0)];
      LP->addConstr(row, GRB_LESS_EQUAL, vars[s_index], LP_row_name("A", { i, j }));
    }
  }
  // Transitivity
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      for (int k = n0; k < n; k++) {
        if (i != j && j != k && i != k) {
          GRBLinExpr row =
              vars[LP_get_variable_index(i, j, n, n0)] + vars[LP_get_variable_index(j, k, n, n0)]
                  - vars[LP_get_variable_index(i, k, n, n0)];
          LP->addConstr(row, GRB_LESS_EQUAL, vars[s_index], LP_row_name("T", { i, j, k }));
        }
      }
    }
  }
  // Density
  GRBLinExpr row = 0;
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      if (i != j) {
        row += vars[LP_get_variable_index(i, j, n, n0)];
      }
    }
  }
  LP->addConstr(row, GRB_EQUAL, 1.0, LP_row_name("D", {}));

  LP->set(GRB_IntParam_OutputFlag, 0);
  LP->optimize();
  int status = LP->get(GRB_IntAttr_Status);
  if (status == GRB_OPTIMAL) {
    double solution = LP->get(GRB_DoubleAttr_ObjVal);
    delete LP, delete environment;
    return solution;
  } else {
    delete LP, delete environment;
    throw std::domain_error("Invalid LP status: " + status);
  }
}
