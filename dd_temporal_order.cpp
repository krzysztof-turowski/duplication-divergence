#include "dd_koala.h"

using namespace std;

typedef Koala::Graph<int, int> Graph;
typedef Koala::Graph<int, int>::PVertex Vertex;

void LP_bound_exact(const int &n, const int &n0, const Parameters &params) {
  // generate g_n0 on n0 vertces and g_n on n vertices using params
  // draw sigma_n uniformly at random and get h_n = sigma_n(g_n)
  // compute P(Pi_n = sigma_n | Pi_n(G_n) = h_n, G_n0 = g_n0) - hint: use Lehmer code to store sigma_n
  // compute p_uv coefficients
  // construct and solve LP
  // export to file and to stdout
}

int main(int argc, char *argv[]) {
  string action(argv[1]), mode(argv[2]);
  int n = stoi(argv[3]), n0 = stoi(argv[4]);
  Parameters params;
  params.initialize(mode, argv + 5);
  if (action == "exact_bound") {
    LP_bound_exact(n, n0, params);
  }
  else {
    assert(0);
  }

  return 0;
}
