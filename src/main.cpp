#include "analysis_driver.hpp"

int main() {
  run_basic_analysis();
  run_cooling_evolution_analysis();
  run_ensemble_analysis();
  run_rilm_analysis();
  run_heated_rilm_analysis();
  run_iilm_analysis();
  run_qmidens_analysis();
  return 0;
}
