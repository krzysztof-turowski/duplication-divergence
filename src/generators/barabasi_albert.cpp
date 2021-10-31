#include "barabasi_albert.h"

BarabasiAlbertParameters::BarabasiAlbertParameters(int const &_m) : m(_m) {
  this->mode = Mode::BA;
}

std::string BarabasiAlbertParameters::to_string() const {
  std::stringstream out;
  out << LONG_NAME.find(this->mode)->second << " ";
  out << "a = " << this->m << " ";
  return out.str();
}

std::string BarabasiAlbertParameters::to_filename() const {
  std::stringstream out;
  out << SHORT_NAME.find(this->mode)->second;
  out << "-" << this->m;
  return out.str();
}

std::string BarabasiAlbertParameters::to_csv() const {
  std::stringstream out;
  out << this->m;
  return out.str();
}

void generate_ba_graph(Graph &G, const int &n, const BarabasiAlbertParameters &params) {
  std::random_device device;
  std::mt19937 generator(device());

  const uint32_t &m_param = params.m;

  const size_t seed_size = G.size();
  if (m_param > seed_size) {
    throw std::invalid_argument("m parameter must not be greater than seed graph size.");
  }
  G.resize(n);

  std::vector<int> degrees(n, 0);
  for (size_t i = 0; i < seed_size; i++) {
    degrees[i] = G[i].size();
  }

  for (int i = seed_size; i < n; i++) {
    std::discrete_distribution<int> distribution(degrees.begin(), degrees.end());
    while (G[i].size() < m_param) {
      const int new_edge = distribution(generator);
      if (G[i].find(new_edge) != G[i].end()) {
        continue;
      }

      degrees[i]++;
      degrees[new_edge]++;
      G[i].insert(new_edge);
      G[new_edge].insert(i);
    }
  }
}
