#include <set>
#include <vector>

class DAG {
 private:
  std::vector<std::set<int>> H;
  std::vector<int> sources;

 public:
  DAG(const int &size) : H(size), sources(size) { }

  DAG(const DAG &other) : H(other.H.size()), sources(other.sources) {
    for (size_t i = 0; i < H.size(); i++) {
      H[i] = other.H[i];
    }
  }

  void add_edge(const int &u, const int &v) {
    H[u].insert(v), ++sources[v];
  }

  std::set<int> get_neighbors(const int &v) const {
    return H[v];
  }

  std::set<int> get_sources() const {
    std::set<int> out;
    for (size_t i = 0; i < sources.size(); i++) {
      if (is_source(i)) {
        out.insert(i);
      }
    }
    return out;
  }

  bool decrement_source(const int &v) {
    return --sources[v];
  }

  bool is_source(const int &v) const {
    return sources[v] == 0;
  }

  void remove_vertex(const int &v) {  
    for (const auto &u : get_neighbors(v)) {
      decrement_source(u);
    }
  }
};
