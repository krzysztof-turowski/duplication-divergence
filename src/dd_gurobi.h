// Header for LP computation of the temporal upper bound using Gurobi.

#pragma once

#include <gurobi_c++.h>

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
    for (auto v : vertices) {
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

inline void add_asymmetry_constraint(
    GRBModel *LP, const std::vector<GRBVar> &vars, const int &n, const int &n0,
    const GRBLinExpr &s) {
  #pragma omp parallel for
  for (int i = n0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      #pragma omp critical
      {
        GRBLinExpr row =
            vars[LP_get_variable_index(i, j, n, n0)] + vars[LP_get_variable_index(j, i, n, n0)];
        LP->addConstr(row, GRB_LESS_EQUAL, s, LP_name("A", { i, j }));
      }
    }
  }
}

inline void add_transitivity_constraint(
    GRBModel *LP, const std::vector<GRBVar> &vars, const int &n, const int &n0,
    const GRBLinExpr &s, const int &i, const int &j, const int &k) {
  #pragma omp critical
  {
    GRBLinExpr row =
        vars[LP_get_variable_index(i, j, n, n0)]
            + vars[LP_get_variable_index(j, k, n, n0)]
            - vars[LP_get_variable_index(i, k, n, n0)];
    LP->addConstr(row, GRB_LESS_EQUAL, s, LP_name("T", { i, j, k }));
  }
}

inline void add_density_constraint(
    GRBModel *LP, const std::vector<GRBVar> &vars, const int &n, const int &n0,
    const double &density) {
  GRBLinExpr row = 0;
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      if (i != j) {
        row += vars[LP_get_variable_index(i, j, n, n0)];
      }
    }
  }
  LP->addConstr(row, GRB_LESS_EQUAL, density, LP_name("D", {}));
}

std::map<std::pair<int, int>, double> retrieve_solution(
    GRBModel *LP, const int &n, const int &n0, const double &s) {
  std::map<std::pair<int, int>, double> solution;
  for (int i = n0; i < n; i++) {
    for (int j = n0; j < n; j++) {
      if (i == j) {
        continue;
      }
      double y_ij = vars[LP_get_variable_index(i, j, n, n0)].get(GRB_DoubleAttr_X);
      solution.insert(std::make_pair(std::make_pair(i, j), y_ij / s));
    }
  }
  return solution;
}

std::tuple<double, std::map<std::pair<int, int>, double>> LP_ordering_solve(
    const std::map<std::pair<int, int>, long double> &p_uv, const int &n, const int &n0,
    const double &epsilon, const bool get_solution = false) {
  try {
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
          const auto &p_ij = p_uv.find(std::make_pair(i, j));
          vars[index] =
              LP->addVar(
                  0.0, 1.0, (p_ij != p_uv.end()) ? static_cast<double>(p_ij->second) : 0.0,
                  GRB_CONTINUOUS, LP_name("y", {i, j}));
        } else {
          vars[index] = LP->addVar(0.0, 0.0, 0.0, GRB_CONTINUOUS, LP_name("y", {i, j}));
        }
      }
    }
    vars[s_index] = LP->addVar(0.0, 1 / density, 0.0, GRB_CONTINUOUS, "s");

    // Antisymmetry
    add_asymmetry_constraint(LP, vars, n, n0, vars[s_index]);
    // Transitivity
    if (MAX_CONSTRAINTS >= pow(n - n0, 3.0)) {
      #pragma omp parallel for
      for (int i = n0; i < n; i++) {
        for (int j = n0; j < n; j++) {
          for (int k = n0; k < n; k++) {
            if (i != j && j != k && i != k) {
              add_transitivity_constraint(LP, vars, n, n0, vars[s_index], i, j, k);
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
        add_transitivity_constraint(LP, vars, n, n0, vars[s_index], i, j, k);
      }
    }
    // Density
    add_density_constraint(LP, vars, n, n0, 1.0);

    LP->set(GRB_IntParam_OutputFlag, 0);
    LP->optimize();
    int status = LP->get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL) {
      double objective = LP->get(GRB_DoubleAttr_ObjVal);
      std::map<std::pair<int, int>, double> solution;
      if (get_solution) {
        double s = vars[s_index].get(GRB_DoubleAttr_X);
        solution = retrieve_solution(LP, n, n0, s);
      }
      delete LP, delete environment;
      return std::make_tuple(objective, solution);
    } else {
      delete LP, delete environment;
      throw std::domain_error("Invalid LP status: " + std::to_string(status));
    }
  } catch (const GRBException &e) {
    throw std::domain_error(
        "LP solver exception code: " + std::to_string(e.getErrorCode())
            + ", message: " + e.getMessage());
  }
}

std::tuple<double, std::map<std::pair<int, int>, double>> LP_binning_solve(
    const std::map<std::pair<int, int>, long double> &p_uv, const int &n, const int &n0,
    const double &epsilon, const bool get_solution = false) {
  try {
    (void)get_solution;
    GRBEnv* environment = new GRBEnv();
    GRBModel *LP = new GRBModel(*environment);
    LP->set(GRB_StringAttr_ModelName, "Solve " + std::to_string(epsilon));
    LP->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    double density = epsilon * (n - n0) * (n - n0 - 1) / 2;

    // Objective function
    int var_count = pow(n - n0, 4.0) + pow(n - n0, 2.0) + 1;
    std::vector<GRBVar> vars(var_count);
    for (int u = n0; u < n; u++) {
      for (int i = n0; i < n; i++) {
        auto index = LP_get_variable_index(u, i, n, n0);
        vars[index] = LP->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, LP_name("y", {u, i}).c_str());
      }
    }
    for (int u = n0; u < n; u++) {
      for (int i = n0; i < n; i++) {
        for (int v = n0; v < n; v++) {
          for (int j = n0; j < n; j++) {
            auto index = LP_get_variable_index(u, i, v, j, n, n0);
            if (u != v && i < j) {
              const auto &p_ij = p_uv.find(std::make_pair(u, v));
              vars[index] =
                  LP->addVar(
                      0.0, 1.0, (p_ij != p_uv.end()) ? static_cast<double>(p_ij->second) : 0.0,
                      GRB_CONTINUOUS, LP_name("w", {u, i, v, j}).c_str());
            } else {
              vars[index] =
                  LP->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, LP_name("w", {u, i, v, j}).c_str());
            }
          }
        }
      }
    }
    int s_index = var_count - 1;
    vars[s_index] = LP->addVar(0.0, 1 / density, 0.0, GRB_CONTINUOUS, "s");

    // Identity
    for (int u = n0; u < n; u++) {
      for (int i = n0; i < n; i++) {
        GRBLinExpr row =
            vars[LP_get_variable_index(u, i, n, n0)]
                - vars[LP_get_variable_index(u, i, u, i, n, n0)];
        LP->addConstr(row, GRB_EQUAL, 0.0, LP_name("I", {u, i}));
      }
    }
    // Symmetry
    for (int u = n0; u < n; u++) {
      for (int i = n0; i < n; i++) {
        for (int v = n0; v < n; v++) {
          for (int j = i + 1; j < n; j++) {
            GRBLinExpr row =
                vars[LP_get_variable_index(u, i, v, j, n, n0)]
                    - vars[LP_get_variable_index(v, j, u, i, n, n0)];
            LP->addConstr(row, GRB_EQUAL, 0.0, LP_name("S", {u, i, v, j}));
          }
        }
      }
    }
    // y-density
    for (int u = n0; u < n; u++) {
      GRBLinExpr row = 0;
      for (int i = n0; i < n; i++) {
        row += vars[LP_get_variable_index(u, i, n, n0)];
      }
      row -= vars[s_index];
      LP->addConstr(row, GRB_EQUAL, 0.0, LP_name("yD", {u}));
    }
    // w-density
    for (int u = n0; u < n; u++) {
      for (int i = n0; i < n; i++) {
        for (int v = n0; v < n; v++) {
          GRBLinExpr row = 0;
          for (int j = n0; j < n; j++) {
            row += vars[LP_get_variable_index(u, i, v, j, n, n0)];
          }
          row -= vars[LP_get_variable_index(u, i, n, n0)];
          LP->addConstr(row, GRB_EQUAL, 0.0, LP_name("wD", {u, i, v}));
        }
      }
    }
    // Density
    GRBLinExpr row = 0;
    for (int u = n0; u < n; u++) {
      for (int i = n0; i < n; i++) {
        for (int v = n0; v < n; v++) {
          for (int j = i + 1; j < n; j++) {
            if (u != v) {
              row += vars[LP_get_variable_index(u, i, v, j, n, n0)];
            }
          }
        }
      }
    }
    LP->addConstr(row, GRB_EQUAL, 1.0, LP_name("D", {}));

    LP->set(GRB_IntParam_OutputFlag, 0);
    LP->optimize();
    int status = LP->get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL) {
      double objective = LP->get(GRB_DoubleAttr_ObjVal);
      std::map<std::pair<int, int>, double> solution;
      delete LP, delete environment;
      return std::make_tuple(objective, solution);
    } else {
      delete LP, delete environment;
      throw std::domain_error("Invalid LP status: " + std::to_string(status));
    }
  } catch (const GRBException &e) {
    throw std::domain_error(
        "LP solver exception code: " + std::to_string(e.getErrorCode())
            + ", message: " + e.getMessage());
  }
}

std::tuple<double, std::map<std::pair<int, int>, double>> IP_ordering_solve(
    const std::map<std::pair<int, int>, long double> &p_uv, const int &n, const int &n0,
    const double &epsilon, const bool get_solution = false) {
  try {
    GRBEnv* environment = new GRBEnv();
    GRBModel *IP = new GRBModel(*environment);
    IP->set(GRB_StringAttr_ModelName, "Solve " + std::to_string(epsilon));
    IP->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    int density = epsilon * (n - n0) * (n - n0 - 1) / 2;

    // Objective function
    std::vector<GRBVar> vars((n - n0) * (n - n0));
    for (int i = n0; i < n; i++) {
      for (int j = n0; j < n; j++) {
        auto index = LP_get_variable_index(i, j, n, n0);
        if (i != j) {
          const auto &p_ij = p_uv.find(std::make_pair(i, j));
          vars[index] =
              IP->addVar(
                  0.0, 1.0, (p_ij != p_uv.end()) ? static_cast<double>(p_ij->second) : 0.0,
                  GRB_BINARY, LP_name("y", {i, j}));
        } else {
          vars[index] = IP->addVar(0.0, 0.0, 0.0, GRB_INTEGER, LP_name("y", {i, j}));
        }
      }
    }

    // Antisymmetry
    add_asymmetry_constraint(IP, vars, n, n0, 1.0);
    // Transitivity
    #pragma omp parallel for
    for (int i = n0; i < n; i++) {
      for (int j = n0; j < n; j++) {
        for (int k = n0; k < n; k++) {
          if (i != j && j != k && i != k) {
            add_transitivity_constraint(IP, vars, n, n0, 1.0, i, j, k);
          }
        }
      }
    }
    // Density
    add_density_constraint(IP, vars, n, n0, density);

    IP->set(GRB_IntParam_OutputFlag, 0);
    IP->optimize();
    int status = IP->get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL) {
      double objective = IP->get(GRB_DoubleAttr_ObjVal) / density;
      std::map<std::pair<int, int>, double> solution;
      if (get_solution) {
        solution = retrieve_solution(IP, n, n0, 1);
      }
      delete IP, delete environment;
      return std::make_tuple(objective, solution);
    } else {
      delete IP, delete environment;
      throw std::domain_error("Invalid IP status: " + std::to_string(status));
    }
  } catch (const GRBException &e) {
    throw std::domain_error(
        "IP solver exception code: " + std::to_string(e.getErrorCode())
            + ", message: " + e.getMessage());
  }
}
