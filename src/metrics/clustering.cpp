#include "clustering.h"

std::pair<int64_t, int64_t> count_triangles(const Graph &G) {
  std::pair<int64_t, int64_t> result{ 0, 0 };

  for (auto &&v : G) {
    result.second += v.size() * (v.size() - 1) / 2;

    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it, 1); jt != v.end(); ++jt) {
        if (G[*it].count(*jt)) {
          result.first++;
        }
      }
    }
  }

  result.first /= 3;

  return result;
}

double clustering_coefficient_two(const Graph &graph) {
  auto triangles = count_triangles(graph);
  return triangles.second > 0 ? static_cast<double>(6 * triangles.first) / triangles.second : 0;
}

std::pair<int64_t, int64_t> count_four_cliques(const Graph &G) {
  std::pair<int64_t, int64_t> result{ 0, 0 };

  for (auto &&v : G) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it, 1); jt != v.end(); ++jt) {
        if (!G[*it].count(*jt)) {
          continue;
        }

        result.second += v.size() - 2;
        for (auto kt = std::next(jt, 1); kt != v.end(); ++kt) {
          if (G[*it].count(*kt) && G[*jt].count(*kt)) {
            result.first++;
          }
        }
      }
    }
  }

  result.first /= 4;

  return result;
}

double clustering_coefficient_three(const Graph &graph) {
  auto four_cliques = count_four_cliques(graph);
  return four_cliques.second > 0
             ? static_cast<double>(3 * (3 + 1) * four_cliques.first) / four_cliques.second
             : 0;
}

std::pair<int64_t, int64_t> count_five_cliques(const Graph &G) {
  std::pair<int64_t, int64_t> result{ 0, 0 };

  for (auto &&v : G) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it, 1); jt != v.end(); ++jt) {
        if (!G[*it].count(*jt)) {
          continue;
        }

        for (auto kt = std::next(jt, 1); kt != v.end(); ++kt) {
          if (!G[*it].count(*kt) || !G[*jt].count(*kt)) {
            continue;
          }
          result.second += v.size() - 3;

          for (auto lt = std::next(kt, 1); lt != v.end(); ++lt) {
            if (G[*it].count(*lt) && G[*jt].count(*lt) && G[*kt].count(*lt)) {
              result.first++;
            }
          }
        }
      }
    }
  }

  result.first /= 5;

  return result;
}

double clustering_coefficient_four(const Graph &graph) {
  auto five_cliques = count_five_cliques(graph);
  return five_cliques.second > 0
             ? static_cast<double>(5 * (5 + 1) * five_cliques.first) / five_cliques.second
             : 0;
}
