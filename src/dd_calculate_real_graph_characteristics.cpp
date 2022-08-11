
#include "./dd_stats.h"

int main(int argc, char const *argv[]) {
  if (argc < 2) {
    return 1;
  }
  int mode = FAST | SLOW;
  if (argc > 2) {
    mode = std::stoi(argv[2]);
  }
  const std::string RESULTS_FOLDER = "results/";
  std::string graph_name = argv[1];

  std::ofstream result;
  result.exceptions(std::ifstream::failbit | std::ifstream::badbit);

  SimpleGraph G(read_graph_simple(FILES_FOLDER + graph_name));
  result.open(RESULTS_FOLDER + graph_name);
  count_graph_statistics(G, result, mode);
  result.close();

  return 0;
}
