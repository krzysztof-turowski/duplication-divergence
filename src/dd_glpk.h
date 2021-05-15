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

std::string LP_name(const std::string &prefix, const std::initializer_list<int> &vertices) {
  std::ostringstream out;
  out << prefix;
  if (vertices.size() > 0) {
    out << "_{";
    for (const auto &v : vertices) {
      out << std::to_string(v) << ",";
    }
    out.seekp(-1, std::ios_base::end), out << "}";
  }
  return out.str();
}

inline int LP_get_variable_index(const int &u, const int &v, const int &n, const int &n0) {
  return (u - n0) * (n - n0) + (v - n0);
}

inline int LP_get_variable_index(
    const int &u, const int &i, const int &v, const int &j, const int &n, const int &n0) {
  int x = LP_get_variable_index(u, i, n, n0), y = LP_get_variable_index(v, j, n, n0);
  return LP_get_variable_index(n, n0, n, n0) + x * (n - n0) * (n - n0) + y;
}

inline void add_transitivity_constraint(
    glp_prob *LP, std::vector<int> &X, std::vector<int> &Y, std::vector<double> &A, int &row,
    const int &n, const int &n0, const int &i, const int &j, const int &k) {
  glp_add_rows(LP, 1);
  glp_set_row_name(LP, row, LP_name("T", {i, j, k}).c_str());
  X.push_back(row), Y.push_back(LP_get_variable_index(i, j, n, n0) + 1);
  X.push_back(row), Y.push_back(LP_get_variable_index(j, k, n, n0) + 1);
  X.push_back(row), Y.push_back(LP_get_variable_index(i, k, n, n0) + 1);
  A.push_back(1), A.push_back(1), A.push_back(-1);
}

inline void add_density_constraint(
    glp_prob *LP, std::vector<int> &X, std::vector<int> &Y, std::vector<double> &A, int &row,
    const double &density, const int &n, const int &n0) {
  glp_add_rows(LP, 1);
  glp_set_row_name(LP, row, LP_name("D", {}).c_str());
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      if (i != j) {
        X.push_back(row), Y.push_back(LP_get_variable_index(i, j, n, n0) + 1), A.push_back(1);
      }
    }
  }
  glp_set_row_bnds(LP, row, GLP_UP, density, density), row++;
}

std::map<std::pair<int, int>, double> retrieve_solution(
    glp_prob *LP, const int &n, const int &n0, const double &s) {
  std::map<std::pair<int, int>, double> solution;
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      if (i == j) {
        continue;
      }
      double y_ij = glp_get_col_prim(LP, LP_get_variable_index(i, j, n, n0) + 1);
      solution.insert(std::make_pair(std::make_pair(i, j), y_ij / s));
    }
  }
  return solution;
}

std::tuple<double, std::map<std::pair<int, int>, double>> LP_ordering_solve(
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
      glp_set_col_name(LP, index + 1, LP_name("y", {i, j}).c_str());
      if (i != j) {
        const auto &p_ij = p_uv.find(std::make_pair(i, j));
        glp_set_col_bnds(LP, index + 1, GLP_LO, 0.0, 1.0);
        glp_set_obj_coef(
            LP, index + 1, p_ij != p_uv.end() ? static_cast<double>(p_ij->second) : 0.0);
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
        glp_set_row_name(LP, row, LP_name("A", {i, j}).c_str());
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
            #pragma omp critical
            {
              add_transitivity_constraint(LP, X, Y, A, row, n, n0, i, j, k);
              X.push_back(row), Y.push_back(s_index + 1), A.push_back(-1);
              glp_set_row_bnds(LP, row, GLP_UP, 0.0, 0.0), row++;
            }
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
      #pragma omp critical
      {
        add_transitivity_constraint(LP, X, Y, A, row, n, n0, i, j, k);
        X.push_back(row), Y.push_back(s_index + 1), A.push_back(-1);
        glp_set_row_bnds(LP, row, GLP_UP, 0.0, 0.0), row++;
      }
    }
  }
  // Density
  add_density_constraint(LP, X, Y, A, row, 1.0, n, n0);

  glp_load_matrix(LP, A.size(), &X[0] - 1, &Y[0] - 1, &A[0] - 1);
  glp_term_out(0);
  glp_simplex(LP, NULL);

  double objective = glp_get_obj_val(LP);
  std::map<std::pair<int, int>, double> solution;
  if (get_solution) {
    solution = retrieve_solution(LP, n, n0, glp_get_col_prim(LP, s_index + 1));
  }
  glp_delete_prob(LP);
  glp_free_env();
  return std::make_tuple(objective, solution);
}

std::tuple<double, std::map<std::pair<int, int>, double>> LP_binning_solve(
    const std::map<std::pair<int, int>, long double> &p_uv, const int &n, const int &n0,
    const double &epsilon, const bool get_solution = false) {
  (void)get_solution;
  glp_prob *LP = glp_create_prob();
  glp_set_prob_name(LP, ("Solve " + std::to_string(epsilon)).c_str());
  glp_set_obj_dir(LP, GLP_MAX);
  double density = epsilon * (n - n0) * (n - n0 - 1) / 2;

  // Objective function
  int var_count = pow(n - n0, 4.0) + pow(n - n0, 2.0) + 1;
  glp_add_cols(LP, var_count);
  for (int u = n0; u < n; u++) {
    for (int i = n0; i < n; i++) {
      auto index = LP_get_variable_index(u, i, n, n0);
      glp_set_col_name(LP, index + 1, LP_name("y", {u, i}).c_str());
      glp_set_col_bnds(LP, index + 1, GLP_LO, 0.0, 1.0);
    }
  }
  for (int u = n0; u < n; u++) {
    for (int i = n0; i < n; i++) {
      for (int v = n0; v < n; v++) {
        for (int j = n0; j < n; j++) {
          auto index = LP_get_variable_index(u, i, v, j, n, n0);
          glp_set_col_name(LP, index + 1, LP_name("w", {u, i, v, j}).c_str());
          glp_set_col_bnds(LP, index + 1, GLP_LO, 0.0, 1.0);
          if (u != v && i < j) {
            const auto &p_ij = p_uv.find(std::make_pair(u, v));
            glp_set_obj_coef(
                LP, index + 1, (p_ij != p_uv.end()) ? static_cast<double>(p_ij->second) : 0.0);
          }
        }
      }
    }
  }
  int s_index = var_count - 1;
  glp_set_col_name(LP, s_index + 1, "s");
  glp_set_col_bnds(LP, s_index + 1, GLP_DB, 0.0, 1 / density);

  std::vector<int> X, Y;
  std::vector<double> A;
  int row = 1;

  // Identity
  glp_add_rows(LP, (n - n0) * (n - n0));
  for (int u = n0; u < n; u++) {
    for (int i = n0; i < n; i++) {
      glp_set_row_name(LP, row, LP_name("I", {u, i}).c_str());
      X.push_back(row), Y.push_back(LP_get_variable_index(u, i, n, n0) + 1), A.push_back(1);
      X.push_back(row), Y.push_back(LP_get_variable_index(u, i, u, i, n, n0) + 1), A.push_back(-1);
      glp_set_row_bnds(LP, row, GLP_FX, 0.0, 0.0), row++;
    }
  }
  // Symmetry
  glp_add_rows(LP, pow(n - n0, 4.0));
  for (int u = n0; u < n; u++) {
    for (int i = n0; i < n; i++) {
      for (int v = n0; v < n; v++) {
        for (int j = i + 1; j < n; j++) {
          glp_set_row_name(LP, row, LP_name("S", {u, i, v, j}).c_str());
          X.push_back(row), Y.push_back(LP_get_variable_index(u, i, v, j, n, n0) + 1);
          X.push_back(row), Y.push_back(LP_get_variable_index(v, j, u, i, n, n0) + 1);
          A.push_back(1), A.push_back(-1);
          glp_set_row_bnds(LP, row, GLP_FX, 0.0, 0.0), row++;
        }
      }
    }
  }
  // y-density
  glp_add_rows(LP, n - n0);
  for (int u = n0; u < n; u++) {
    glp_set_row_name(LP, row, LP_name("yD", {u}).c_str());
    for (int i = n0; i < n; i++) {
      X.push_back(row), Y.push_back(LP_get_variable_index(u, i, n, n0) + 1), A.push_back(1);
    }
    X.push_back(row), Y.push_back(s_index + 1), A.push_back(-1);
    glp_set_row_bnds(LP, row, GLP_FX, 0.0, 0.0), row++;
  }
  // w-density
  glp_add_rows(LP, pow(n - n0, 3.0));
  for (int u = n0; u < n; u++) {
    for (int i = n0; i < n; i++) {
      for (int v = n0; v < n; v++) {
        glp_set_row_name(LP, row, LP_name("wD", {u, i, v}).c_str());
        for (int j = n0; j < n; j++) {
            X.push_back(row), Y.push_back(LP_get_variable_index(u, i, v, j, n, n0) + 1);
            A.push_back(1);
        }
        X.push_back(row), Y.push_back(LP_get_variable_index(u, i, n, n0) + 1), A.push_back(-1);
        glp_set_row_bnds(LP, row, GLP_FX, 0.0, 0.0), row++;
      }
    }
  }
  // Density
  glp_add_rows(LP, 1);
  glp_set_row_name(LP, row, LP_name("D", {}).c_str());
  for (int u = n0; u < n; u++) {
    for (int i = n0; i < n; i++) {
      for (int v = n0; v < n; v++) {
        for (int j = i + 1; j < n; j++) {
          if (u != v) {
            X.push_back(row), Y.push_back(LP_get_variable_index(u, i, v, j, n, n0) + 1);
            A.push_back(1);
          }
        }
      }
    }
  }
  glp_set_row_bnds(LP, row, GLP_FX, 1.0, 1.0), row++;

  glp_load_matrix(LP, A.size(), &X[0] - 1, &Y[0] - 1, &A[0] - 1);
  glp_term_out(0);
  glp_simplex(LP, NULL);

  double objective = glp_get_obj_val(LP);
  std::map<std::pair<int, int>, double> solution;
  glp_delete_prob(LP);
  glp_free_env();
  return std::make_tuple(objective, solution);
}

std::tuple<double, std::map<std::pair<int, int>, double>> IP_ordering_solve(
    const std::map<std::pair<int, int>, long double> &p_uv, const int &n, const int &n0,
    const double &epsilon, const bool get_solution = false) {
  glp_prob *IP = glp_create_prob();
  glp_set_prob_name(IP, ("Solve " + std::to_string(epsilon)).c_str());
  glp_set_obj_dir(IP, GLP_MAX);
  int density = epsilon * (n - n0) * (n - n0 - 1) / 2;

  // Objective function
  glp_add_cols(IP, (n - n0) * (n - n0));
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      auto index = LP_get_variable_index(i, j, n, n0);
      glp_set_col_name(IP, index + 1, LP_name("y", {i, j}).c_str());
      if (i != j) {
        const auto &p_ij = p_uv.find(std::make_pair(i, j));
        glp_set_col_kind(IP, index + 1, GLP_BV);
        glp_set_obj_coef(
            IP, index + 1, p_ij != p_uv.end() ? static_cast<double>(p_ij->second) : 0.0);
      }
    }
  }

  std::vector<int> X, Y;
  std::vector<double> A;
  int row = 1;
  glp_add_rows(IP, (n - n0) * (n - n0 - 1) / 2);
  // Antisymmetry
  #pragma omp parallel for
  for (int i = n0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      #pragma omp critical
      {
        glp_set_row_name(IP, row, LP_name("A", {i, j}).c_str());
        X.push_back(row), Y.push_back(LP_get_variable_index(i, j, n, n0) + 1), A.push_back(1);
        X.push_back(row), Y.push_back(LP_get_variable_index(j, i, n, n0) + 1), A.push_back(1);
        glp_set_row_bnds(IP, row, GLP_UP, 1, 1), row++;
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
            #pragma omp critical
            {
              add_transitivity_constraint(IP, X, Y, A, row, n, n0, i, j, k);
              glp_set_row_bnds(IP, row, GLP_UP, 1, 1), row++;
            }
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
      #pragma omp critical
      {
        add_transitivity_constraint(IP, X, Y, A, row, n, n0, i, j, k);
        glp_set_row_bnds(IP, row, GLP_UP, 1, 1), row++;
      }
    }
  }
  // Density
  add_density_constraint(IP, X, Y, A, row, density, n, n0);

  glp_load_matrix(IP, A.size(), &X[0] - 1, &Y[0] - 1, &A[0] - 1);
  glp_term_out(0);
  glp_simplex(IP, NULL);
  glp_intopt(IP, NULL);

  double objective = glp_mip_obj_val(IP) / density;
  std::map<std::pair<int, int>, double> solution;
  if (get_solution) {
    solution = retrieve_solution(IP, n, n0, 1);
  }
  glp_delete_prob(IP);
  glp_free_env();
  return std::make_tuple(objective, solution);
}
