#include "parameters.h"

const int PRECISION_P = 3, PRECISION_Q = 3, PRECISION_R = 2, WIDTH_R = 6;

const std::map<Mode, std::string> SHORT_NAME = {
  { Mode::PURE_DUPLICATION, "PD" },
  { Mode::PURE_DUPLICATION_CONNECTED, "PDC" },
  { Mode::CHUNG_LU, "CL" },
  { Mode::PASTOR_SATORRAS, "PS" },
  { Mode::STICKY, "STICKY" },
  { Mode::BA, "BA" },
  { Mode::COPY_GRAPH, "COPY" },
  { Mode::TWO_STEP, "2STEP" },
};

const std::map<Mode, std::string> LONG_NAME = {
  { Mode::PURE_DUPLICATION, "Pure duplication" },
  { Mode::PURE_DUPLICATION_CONNECTED, "Pure duplication without isolated vertices" },
  { Mode::CHUNG_LU, "Chung-Lu" },
  { Mode::PASTOR_SATORRAS, "Pastor-Satorras" },
  { Mode::STICKY, "STICKY" },
  { Mode::BA, "Barabasi-Albert" },
  { Mode::COPY_GRAPH, "Copy graph" },
  { Mode::TWO_STEP, "Two step" },
};

const std::map<std::string, Mode> REVERSE_NAME = {
  { "pure_duplication", Mode::PURE_DUPLICATION },
  { "pure_duplication_connected", Mode::PURE_DUPLICATION_CONNECTED },
  { "chung_lu", Mode::CHUNG_LU },
  { "pastor_satorras", Mode::PASTOR_SATORRAS },
  { "sticky", Mode::STICKY },
  { "ba", Mode::BA },
  { "copy_graph", Mode::COPY_GRAPH },
  { "two_step", Mode::TWO_STEP },
};

Parameters::Parameters() : mode(Mode::INVALID), p(nan("")), q(nan("")), r(nan("")) { }

void Parameters::initialize(const std::string &mode_v, char *argv[]) {
  if (mode_v == "pure_duplication") {
    initialize_pure_duplication(std::stod(argv[0]));
  } else if (mode_v == "chung_lu") {
    initialize_chung_lu(std::stod(argv[0]), std::stod(argv[1]));
  } else if (mode_v == "pastor_satorras") {
    initialize_pastor_satorras(std::stod(argv[0]), std::stod(argv[1]));
  } else {
    throw std::invalid_argument("Invalid mode: " + mode_v);
  }
}

void Parameters::initialize_pure_duplication(const double &p_v) {
  this->mode = Mode::PURE_DUPLICATION;
  this->p = p_v;
  this->q = this->r = nan("");
}

void Parameters::initialize_pure_duplication_connected(const double &p_v) {
  this->mode = Mode::PURE_DUPLICATION_CONNECTED;
  this->p = p_v;
  this->q = this->r = nan("");
}

void Parameters::initialize_chung_lu(const double &p_v, const double &q_v) {
  this->mode = Mode::CHUNG_LU;
  this->p = p_v, this->q = q_v;
  this->r = nan("");
}

void Parameters::initialize_pastor_satorras(const double &p_v, const double &r_v) {
  this->mode = Mode::PASTOR_SATORRAS;
  this->p = p_v, this->r = r_v;
  this->q = nan("");
}

std::string Parameters::to_string() const {
  std::stringstream out;
  out << LONG_NAME.find(this->mode)->second << " ";
  if (!std::isnan(this->p)) {
    out << "p = " << std::fixed << std::setprecision(PRECISION_P) << this->p << " ";
  }
  if (!std::isnan(this->q)) {
    out << "q = " << std::fixed << std::setprecision(PRECISION_Q) << this->q << " ";
  }
  if (!std::isnan(this->r)) {
    out << "r = " << std::fixed << std::setw(WIDTH_R) << std::setprecision(PRECISION_R) << this->r
        << " ";
  }
  return out.str();
}

std::string Parameters::to_string(const Parameters &low, const Parameters &high) const {
  std::stringstream out;
  out << LONG_NAME.find(this->mode)->second << " ";
  out << "p_min = " << std::fixed << std::setprecision(PRECISION_P) << low.p << " "
      << "p = " << std::fixed << std::setprecision(PRECISION_P) << this->p << " "
      << "p_max = " << std::fixed << std::setprecision(PRECISION_P) << high.p << " ";
  if (!std::isnan(this->q)) {
    out << "q = " << std::fixed << std::setprecision(PRECISION_Q) << this->q << " ";
  }
  if (!std::isnan(this->r)) {
    out << "r = " << std::fixed << std::setw(WIDTH_R) << std::setprecision(PRECISION_R) << this->r
        << " ";
  }
  return out.str();
}

std::string Parameters::to_filename() const {
  std::stringstream out;
  out << SHORT_NAME.find(this->mode)->second;
  if (!std::isnan(this->p)) {
    out << "-" << std::fixed << std::setprecision(PRECISION_P) << this->p;
  }
  if (!std::isnan(this->q)) {
    out << "-" << std::fixed << std::setprecision(PRECISION_Q) << this->q;
  }
  if (!std::isnan(this->r)) {
    out << "-" << std::fixed << std::setprecision(PRECISION_R) << this->r;
  }
  return out.str();
}

std::string Parameters::to_csv() const {
  std::stringstream out;
  if (!std::isnan(this->p)) {
    out << this->p;
  }
  out << ",";
  if (!std::isnan(this->q)) {
    out << this->q;
  }
  out << ",";
  if (!std::isnan(this->r)) {
    out << this->r;
  }
  return out.str();
}

std::string Parameters::short_name() const {
  return SHORT_NAME.find(this->mode)->second;
}

std::string Parameters::long_name() const {
  return LONG_NAME.find(this->mode)->second;
}
