/**
 * @file main.cpp
 * @brief Command-line interface for the instanton simulation project.
 *
 * This executable acts as a dispatcher for all high-level numerical
 * experiments defined in analysis_driver.hpp.
 *
 * The program allows selective execution of:
 *
 *   - Full Monte Carlo analysis (path + correlators)
 *   - Cooling evolution studies
 *   - Ensemble averaging
 *   - Random Instanton Liquid Model (RILM)
 *   - Heated RILM
 *   - Interacting Instanton Liquid Model (IILM)
 *   - Non-Gaussian density correction (adiabatic switching)
 *   - Instanton density scan versus eta
 *
 * This design ensures:
 *   - Reproducibility (via fixed RNG seed)
 *   - Modular execution (figure-by-figure workflow)
 *   - No unnecessary recomputation
 */

#include "analysis_driver.hpp"
#include <iostream>
#include <random>
#include <string>
#include <vector>

/**
 * @brief Print usage information for the executable.
 *
 * @param prog Name of the executable (argv[0]).
 *
 * Displays available commands and optional arguments.
 */
static void print_help(const char *prog) {
  std::cout
      << "Usage:\n"
      << "  " << prog << " <command> [options]\n\n"
      << "Commands:\n"
      << "  basic                Run basic MC + cooled config export (Fig. 2, "
         "correlators)\n"
      << "  cooling-evolution    Run density vs cooling sweeps (Fig. 7a)\n"
      << "  ensemble             Run ensemble averages (Fig. 4/6 style)\n"
      << "  rilm                 Run RILM (Fig. 10)\n"
      << "  heated-rilm          Run heated RILM (Fig. 12/13)\n"
      << "  iilm                 Run IILM\n"
      << "  qmidens              Run non-Gaussian density correction\n"
      << "  eta-scan             Run cooling eta scan (writes "
         "data/cooling_eta_scan.csv)\n\n"
      << "Options:\n"
      << "  --seed <int>         RNG seed (default: 12345)\n"
      << "  --etas a,b,c         Only for eta-scan (comma-separated list)\n";
}

/**
 * @brief Parse a comma-separated list of eta values.
 *
 * Example:
 *     "1.0,1.2,1.4"
 *
 * @param s Input string.
 * @return Vector of eta values as doubles.
 */
static std::vector<double> parse_etas(const std::string &s) {
  std::vector<double> etas;
  std::string token;

  for (char ch : s) {
    if (ch == ',') {
      if (!token.empty()) {
        etas.push_back(std::stod(token));
      }
      token.clear();
    } else {
      token.push_back(ch);
    }
  }

  if (!token.empty()) {
    etas.push_back(std::stod(token));
  }

  return etas;
}

/**
 * @brief Entry point of the instanton simulation project.
 *
 * The program expects at least one command argument.
 * Optional flags allow control of:
 *   - RNG seed
 *   - eta scan values
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 *
 * @return 0 on success, non-zero on error.
 */
int main(int argc, char **argv) {

  if (argc < 2) {
    print_help(argv[0]);
    return 1;
  }

  std::string cmd;
  int seed = 12345;
  std::vector<double> etas;

  cmd = argv[1];

  /**
   * Parse optional flags following the command.
   *
   * Supported:
   *   --seed <int>
   *   --etas <comma-separated list>
   */
  for (int i = 2; i < argc; ++i) {
    std::string a = argv[i];

    if (a == "--seed" && i + 1 < argc) {
      seed = std::stoi(argv[++i]);
    } else if (a == "--etas" && i + 1 < argc) {
      etas = parse_etas(argv[++i]);
    } else if (a == "--help" || a == "-h") {
      print_help(argv[0]);
      return 0;
    } else {
      std::cerr << "Unknown option: " << a << "\n";
      print_help(argv[0]);
      return 1;
    }
  }

  /**
   * Initialize random number generator.
   *
   * If a fixed seed is provided, the run becomes fully reproducible.
   */
  std::mt19937 gen(static_cast<std::mt19937::result_type>(seed));

  /**
   * Dispatch command.
   */
  if (cmd == "basic") {
    run_basic_analysis(gen);
  } else if (cmd == "cooling-evolution") {
    run_cooling_evolution_analysis(gen);
  } else if (cmd == "ensemble") {
    run_ensemble_analysis(gen);
  } else if (cmd == "rilm") {
    run_rilm_analysis(gen);
  } else if (cmd == "heated-rilm") {
    run_heated_rilm_analysis(gen);
  } else if (cmd == "iilm") {
    run_iilm_analysis(gen);
  } else if (cmd == "qmidens") {
    run_qmidens_analysis(gen);
  } else if (cmd == "eta-scan") {
    if (etas.empty()) {
      /**
       * Default eta grid if none provided explicitly.
       */
      etas = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
    }
    run_cooling_eta_scan(etas, gen);
  } else {
    std::cerr << "Unknown command: " << cmd << "\n";
    print_help(argv[0]);
    return 1;
  }

  return 0;
}