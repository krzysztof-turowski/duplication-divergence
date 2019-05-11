// Header for LP computation of the temporal upper bound using GLPK.

#pragma once

#include <glpk.h>

#include <map>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <tuple>
#include <vector>

const int64_t MAX_CONSTRAINTS = 5000000000L;

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

inline void add_transitivity_constraint(
    glp_prob *LP, std::vector<int> &X, std::vector<int> &Y, std::vector<double> &A, int &row,
    const int &n, const int &n0, const int &s_index, const int &i, const int &j, const int &k) {
  #pragma omp critical
  {
    glp_add_rows(LP, 1);
    glp_set_row_name(LP, row, LP_row_name("T", {i, j, k}).c_str());
    X.push_back(row), Y.push_back(LP_get_variable_index(i, j, n, n0) + 1);
    X.push_back(row), Y.push_back(LP_get_variable_index(j, k, n, n0) + 1);
    X.push_back(row), Y.push_back(LP_get_variable_index(i, k, n, n0) + 1);
    A.push_back(1), A.push_back(1), A.push_back(-1);
    X.push_back(row), Y.push_back(s_index + 1), A.push_back(-1);
    glp_set_row_bnds(LP, row, GLP_UP, 0.0, 0.0), row++;
  }
}

std::tuple<double, std::map<std::pair<int, int>, double>> LP_solve(
    const std::map<std::pair<int, int>, long double> &p_uv, const int &n, const int &n0,
    const double &epsilon, const bool get_solution = false) {
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
        const auto &p_ij = p_uv.find(std::make_pair(i, j));
        glp_set_col_bnds(LP, index + 1, GLP_DB, 0.0, 1.0);
        glp_set_obj_coef(
            LP, index + 1, (p_ij != p_uv.end()) ? static_cast<double>(p_ij->second) : 0.0);
      }
    }
  }
  glp_set_col_name(LP, s_index + 1, "s");
  glp_set_col_bnds(LP, s_index + 1, GLP_DB, 0.0, 1 / density);

  std::vector<int> X, Y;
  std::vector<double> A;
  int row = 1;
  glp_add_rows(LP, (n - n0) * (n - n0 - 1) / 2);
  // Antisymmetry
  #pragma omp parallel for
  for (int i = n0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      #pragma omp critical
      {
        glp_set_row_name(LP, row, LP_row_name("A", {i, j}).c_str());
        X.push_back(row), Y.push_back(LP_get_variable_index(i, j, n, n0) + 1), A.push_back(1);
        X.push_back(row), Y.push_back(LP_get_variable_index(j, i, n, n0) + 1), A.push_back(1);
        X.push_back(row), Y.push_back(s_index + 1), A.push_back(-1);
        glp_set_row_bnds(LP, row, GLP_UP, 0.0, 0.0), row++;
      }
    }
  }
  // Transitivity
  if (MAX_CONSTRAINTS >= pow(n - n0, 3.0)) {
    #pragma omp parallel for
    for (int i = n0; i < n; i++) {
      for (int j = n0; j < n; j++) {
        for (int k = n0; k < n; k++) {
          if (i != j && j != k && i != k) {
            add_transitivity_constraint(LP, X, Y, A, row, n, n0, s_index, i, j, k);
          }
        }
      }
    }
  } else {
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_int_distribution<int> index_distribution(n0, n - 1);
    #pragma omp parallel for
    for (int64_t constraint = 0; constraint < MAX_CONSTRAINTS; constraint++) {
      int i = index_distribution(generator), j = index_distribution(generator),
          k = index_distribution(generator);
      if (i == j || j == k || i == k) {
        continue;
      }
      add_transitivity_constraint(LP, X, Y, A, row, n, n0, s_index, i, j, k);
    }
  }
  // Density
  glp_add_rows(LP, 1);
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

  double objective = glp_get_obj_val(LP);
  std::map<std::pair<int, int>, double> solution;
  if (get_solution) {
    double s = glp_get_col_prim(LP, s_index + 1);
    for (int i = n0; i < n; i++) {
      for (int j = n0; j < n; j++) {
        if (i == j) {
          continue;
        }
        double y_ij = glp_get_col_prim(LP, LP_get_variable_index(i, j, n, n0) + 1);
        solution.insert(std::make_pair(std::make_pair(i, j), y_ij / s));
      }
    }
  }
  glp_delete_prob(LP);
  glp_free_env();
  return std::make_tuple(objective, solution);
}