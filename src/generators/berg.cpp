#include "berg.h"

const int PRECISION_TU = 1;
const int PRECISION_ED = 0;
const int PRECISION_AC = 3;
const int PRECISION_DR = 3;

BergParameters::BergParameters(const double &tu, const double &ed, const double &ac,
    const double &dr, const double &lar, const double &ldr)
    : time_unit(tu), evolution_duration(ed), average_connectivity(ac), duplication_rate(dr),
      link_addition_rate(lar), link_detachment_rate(ldr) {
  this->mode = Mode::BERG;
}

std::string BergParameters::to_string() const {
  std::stringstream out;
  out << this->long_name() << " ";
  out << "tu = " << std::fixed << std::setprecision(PRECISION_TU) << this->time_unit << " ";
  out << "ed = " << std::fixed << std::setprecision(PRECISION_ED) << this->evolution_duration
      << " ";
  out << "ac = " << std::fixed << std::setprecision(PRECISION_AC) << this->average_connectivity
      << " ";
  out << "dr = " << std::fixed << std::setprecision(PRECISION_DR) << this->duplication_rate << " ";
  out << "lar = " << std::fixed << std::setprecision(PRECISION_DR) << this->link_addition_rate
      << " ";
  out << "ldr = " << std::fixed << std::setprecision(PRECISION_DR) << this->link_detachment_rate
      << " ";
  return out.str();
}

std::string BergParameters::to_filename() const {
  std::stringstream out;
  out << this->short_name();
  out << "-" << std::fixed << std::setprecision(PRECISION_TU) << this->time_unit;
  out << "-" << std::fixed << std::setprecision(PRECISION_ED) << this->evolution_duration;
  out << "-" << std::fixed << std::setprecision(PRECISION_AC) << this->average_connectivity;
  out << "-" << std::fixed << std::setprecision(PRECISION_AC) << this->duplication_rate;
  out << "-" << std::fixed << std::setprecision(PRECISION_AC) << this->link_addition_rate;
  out << "-" << std::fixed << std::setprecision(PRECISION_AC) << this->link_detachment_rate;
  return out.str();
}

std::string BergParameters::to_csv() const {
  std::stringstream out;
  out << this->time_unit;
  out << ",";
  out << this->evolution_duration;
  out << ",";
  out << this->average_connectivity;
  out << ",";
  out << this->duplication_rate;
  out << ",";
  out << this->link_addition_rate;
  out << ",";
  out << this->link_detachment_rate;
  return out.str();
}

std::string BergParameters::short_name() const {
  return "BERG";
}

std::string BergParameters::long_name() const {
  return "Berg";
}

void generate_berg_graph(SimpleGraph &G, const BergParameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  std::vector<int> degrees(G.size(), 0);
  int sum_degrees = 0;
  for (size_t i = 0; i < G.size(); i++) {
    degrees[i] = G[i].size();
    sum_degrees += degrees[i];
  }

  for (double t = 0; t <= params.evolution_duration; t += params.time_unit) {
    for (size_t i = 0; i < G.size(); i++) {
      auto &node = G[i];
      if (distribution(generator) <= params.time_unit * params.duplication_rate) {
        G.resize(G.size() + 1);
        degrees.push_back(0);
      }

      if (distribution(generator) <= params.time_unit * params.link_addition_rate) {
        std::discrete_distribution<int> vertex_distribution(degrees.begin(), degrees.end());
        size_t new_node = vertex_distribution(generator);
        for (size_t j = 0; j < G.size() && (new_node == i || node.find(new_node) != node.end());
             j++) {
          new_node = vertex_distribution(generator);
        }

        G[i].insert(new_node);
        G[new_node].insert(i);
        sum_degrees++;
      }

      if (sum_degrees >= params.average_connectivity * degrees.size()
          && distribution(generator) <= params.time_unit * params.link_detachment_rate) {
        std::uniform_int_distribution<int> vertex_distribution(0, G[i].size());
        auto it = G[i].begin();
        std::advance(it, vertex_distribution(generator));
        G[i].erase(*it);
        G[*it].erase(i);
        sum_degrees--;
      }
    }
  }
}
