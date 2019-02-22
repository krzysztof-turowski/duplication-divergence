// Header for LP computation of the temporal upper bound using GLPK.

#pragma once

#include <glpk.h>

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
  glp_prob *LP = glp_create_prob();
  glp_set_prob_name(LP, ("Solve " + std::to_string(epsilon)).c_str());
  glp_set_obj_dir(LP, GLP_MAX);
  double density = epsilon * (n - n0) * (n - n0 - 1) / 2;

  // Objective function
  glp_add_cols(LP, (n - n0) * (n - n0) + 1);
  int s_index = (n - n0) * (n - n0);
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      auto index = LP_get_variable_index(i, j, n, n0);
      glp_set_col_name(LP, index + 1, LP_name_variable(i, j).c_str());
      if (i != j) {
        glp_set_col_bnds(LP, index + 1, GLP_DB, 0.0, 1.0);
        glp_set_obj_coef(LP, index + 1, p_uv.find(std::make_pair(i, j))->second);
      }
    }
  }
  glp_set_col_name(LP, s_index + 1, "s");
  glp_set_col_bnds(LP, s_index + 1, GLP_DB, 0.0, 1 / density);

  std::vector<int> X, Y;
  std::vector<double> A;
  int row = 1;
  glp_add_rows(LP, (n - n0) * (n - n0 - 1) * (n - n0 - 2) + (n - n0) * (n - n0 - 1) / 2 + 1);
  // Antisymmetry
  for (int i = n0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      glp_set_row_name(LP, row, LP_row_name("A", {i, j}).c_str());
      X.push_back(row), Y.push_back(LP_get_variable_index(i, j, n, n0) + 1), A.push_back(1);
      X.push_back(row), Y.push_back(LP_get_variable_index(j, i, n, n0) + 1), A.push_back(1);
      X.push_back(row), Y.push_back(s_index + 1), A.push_back(-1);
      glp_set_row_bnds(LP, row, GLP_UP, 0.0, 0.0), row++;
    }
  }
  // Transitivity
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      for (int k = n0; k < n; k++) {
        if (i != j && j != k && i != k) {
          glp_set_row_name(LP, row, LP_row_name("T", {i, j, k}).c_str());
          X.push_back(row), Y.push_back(LP_get_variable_index(i, j, n, n0) + 1), A.push_back(1);
          X.push_back(row), Y.push_back(LP_get_variable_index(j, k, n, n0) + 1), A.push_back(1);
          X.push_back(row), Y.push_back(LP_get_variable_index(i, k, n, n0) + 1), A.push_back(-1);
          X.push_back(row), Y.push_back(s_index + 1), A.push_back(-1);
          glp_set_row_bnds(LP, row, GLP_UP, 0.0, 0.0), row++;
        }
      }
    }
  }
  // Density
  glp_set_row_name(LP, row, LP_row_name("D", {}).c_str());
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      if (i != j) {
        X.push_back(row), Y.push_back(LP_get_variable_index(i, j, n, n0) + 1), A.push_back(1);
      }
    }
  }
  glp_set_row_bnds(LP, row, GLP_FX, 1.0, 1.0), row++;

  glp_load_matrix(LP, A.size(), &X[0] - 1, &Y[0] - 1, &A[0] - 1);
  glp_term_out(0);
  glp_simplex(LP, NULL);

  double solution = glp_get_obj_val(LP);
  glp_delete_prob(LP);
  glp_free_env();
  return solution;
}
