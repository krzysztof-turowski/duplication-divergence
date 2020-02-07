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

inline void add_transitivity_constraint(
    GRBModel *LP, const std::vector<GRBVar> &vars, const int &n, const int &n0,
    const int &s_index, const int &i, const int &j, const int &k) {
  #pragma omp critical
  {
    GRBLinExpr row =
        vars[LP_get_variable_index(i, j, n, n0)]
            + vars[LP_get_variable_index(j, k, n, n0)]
            - vars[LP_get_variable_index(i, k, n, n0)];
    LP->addConstr(row, GRB_LESS_EQUAL, vars[s_index], LP_name("T", { i, j, k }));
  }
}

std::tuple<double, std::map<std::pair<int, int>, double>> LP_solve(
    GRBModel *LP, const std::vector<GRBVar> &vars, const int &n, const int &n0,
    const int &s_index, const bool get_solution = false) {
  LP->set(GRB_IntParam_OutputFlag, 0);
  LP->optimize();
  int status = LP->get(GRB_IntAttr_Status);
  if (status == GRB_OPTIMAL) {
    double objective = LP->get(GRB_DoubleAttr_ObjVal);
    std::map<std::pair<int, int>, double> solution;
    if (get_solution) {
      double s = vars[s_index].get(GRB_DoubleAttr_X);
      for (int i = n0; i < n; i++) {
        for (int j = n0; j < n; j++) {
          if (i == j) {
            continue;
          }
          double y_ij = vars[LP_get_variable_index(i, j, n, n0)].get(GRB_DoubleAttr_X);
          solution.insert(std::make_pair(std::make_pair(i, j), y_ij / s));
        }
      }
    }
    return std::make_tuple(objective, solution);
  } else {
    throw std::domain_error("Invalid LP status: " + std::to_string(status));
  }
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
    #pragma omp parallel for
    for (int i = n0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        #pragma omp critical
        {
          GRBLinExpr row =
              vars[LP_get_variable_index(i, j, n, n0)] + vars[LP_get_variable_index(j, i, n, n0)];
          LP->addConstr(row, GRB_LESS_EQUAL, vars[s_index], LP_name("A", { i, j }));
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
              add_transitivity_constraint(LP, vars, n, n0, s_index, i, j, k);
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
        add_transitivity_constraint(LP, vars, n, n0, s_index, i, j, k);
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
    LP->addConstr(row, GRB_EQUAL, 1.0, LP_name("D", {}));
    auto value = LP_solve(LP, vars, n, n0, s_index, get_solution);
    delete LP, delete environment;
    return value;
  } catch (const GRBException &e) {
    throw std::domain_error(
        "LP solver exception code: " + std::to_string(e.getErrorCode())
            + ", message: " + e.getMessage());
  } catch (const std::exception &e) {
    delete LP, delete environment;
    throw e;
  }
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

    // TODO

    auto value = LP_solve(LP, vars, n, n0, s_index, get_solution);
    delete LP, delete environment;
    return value;
  } catch (const GRBException &e) {
    throw std::domain_error(
        "LP solver exception code: " + std::to_string(e.getErrorCode())
            + ", message: " + e.getMessage());
  } catch (const std::exception &e) {
    delete LP, delete environment;
    throw e;
  }
}
