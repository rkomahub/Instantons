/**
 * @file main.cpp
 * @brief Command-line interface for the instanton simulation project.
 *
 * This executable dispatches high-level numerical analyses defined in
 * analysis_driver.hpp.
 *
 * The program supports:
 *
 *   - Full Euclidean Monte Carlo simulations
 *   - Cooling analyses and instanton counting
 *   - Ensemble averages of correlators
 *   - Random and interacting instanton liquid models
 *   - Adiabatic switching calculations
 *   - Figure-by-figure data generation
 *
 * The interface is designed for:
 *
 *   - Reproducibility through explicit RNG seeds
 *   - Modular execution of individual analyses
 *   - Minimal recomputation during figure production
 */

#include "analysis_driver.hpp"

#include <iostream>
#include <random>
#include <string>
#include <vector>

/**
 * @brief Print command-line usage information.
 *
 * @param prog Program name, typically argv[0].
 *
 * The help message lists available commands, their purpose, and
 * supported optional flags.
 */
static void print_help(const char *prog) {
  std::cout
      << "\n"
      << "Monte Carlo Instantons in the Double Well Potential\n"
      << "--------------------------------------------------\n\n"
      << "Usage:\n"
      << "  " << prog << " <command> [options]\n"
      << "  " << prog << " --help\n\n"
      << "Main commands:\n"
      << "  basic                Basic Monte Carlo path analysis and cooled "
         "export\n"
      << "                       Produces data for typical paths and "
         "correlators\n"
      << "                       Related to Fig. 2 and basic correlator "
         "studies\n\n"
      << "  ensemble             Ensemble-averaged correlators in quantum and "
         "cooled\n"
      << "                       configurations\n"
      << "                       Related to Figs. 4 and 6\n\n"
      << "  cooling-evolution    Instanton density versus number of cooling "
         "sweeps\n"
      << "                       Cooling history on a single Monte Carlo "
         "configuration\n\n"
      << "Figure-oriented commands:\n"
      << "  fig7                 Averaged cooling evolution for eta = 1.4, "
         "1.5, 1.6\n"
      << "                       Related to Fig. 7\n\n"
      << "  fig8                 Instanton density versus eta\n"
      << "                       Combines cooling estimate and non-Gaussian "
         "correction\n"
      << "                       Related to Fig. 8\n\n"
      << "  fig9                 Switching paths for vacuum and one-instanton "
         "sectors\n"
      << "                       Related to Fig. 9\n\n"
      << "  fig11                Gaussian effective potential around one "
         "instanton\n"
      << "                       Related to Fig. 11\n\n"
      << "  fig14                Instanton–anti-instanton interaction and "
         "separation\n"
      << "                       histograms\n"
      << "                       Related to Figs. 14 and 16\n\n"
      << "  fig15                Relaxed streamline family of IA "
         "configurations\n"
      << "                       Related to Fig. 15\n\n"
      << "  fig16                Alias of fig14 for separation-histogram "
         "production\n"
      << "                       Related to Fig. 16\n\n"
      << "  fig17                Monte Carlo history of instanton positions in "
         "IILM\n"
      << "                       Related to Fig. 17\n\n"
      << "Instanton liquid model commands:\n"
      << "  rilm                 Random Instanton Liquid Model correlators\n"
      << "                       Related to Fig. 10\n\n"
      << "  heated-rilm          Heated RILM paths and correlators\n"
      << "                       Related to Figs. 12 and 13\n\n"
      << "  iilm                 Interacting Instanton Liquid Model "
         "analysis\n\n"
      << "Switching / density commands:\n"
      << "  qmidens              Adiabatic switching analysis\n"
      << "                       Used for free-energy and non-Gaussian "
         "corrections\n"
      << "                       Related to Fig. 5 and part of Fig. 8\n\n"
      << "  eta-scan             Cooling-based density scan over a list of eta "
         "values\n"
      << "                       Writes data/cooling_eta_scan.csv\n\n"
      << "Options:\n"
      << "  --seed <int>         Set RNG seed for reproducible runs\n"
      << "                       Default: 12345\n\n"
      << "  --etas a,b,c         Comma-separated eta list\n"
      << "                       Used by: fig8, eta-scan\n"
      << "                       Example: --etas 1.0,1.2,1.4,1.6\n\n"
      << "  -h, --help           Show this help message\n\n"
      << "Examples:\n"
      << "  " << prog << " basic\n"
      << "  " << prog << " fig8 --etas 1.0,1.2,1.4,1.6 --seed 7\n"
      << "  " << prog << " qmidens --seed 12345\n"
      << "  " << prog << " eta-scan --etas 0.8,1.0,1.2,1.4,1.6,1.8,2.0\n\n";
}

/**
 * @brief Parse a comma-separated list of eta values.
 *
 * Example:
 *
 *     "1.0,1.2,1.4"
 *
 * @param s Input string.
 * @return Vector of eta values.
 */
static std::vector<double> parse_etas(const std::string &s) {
  std::vector<double> etas;
  std::string token;

  for (char ch : s) {
    if (ch == ',') {
      if (!token.empty()) {
        etas.push_back(std::stod(token));
        token.clear();
      }
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
 * The program expects one command and optional flags.
 *
 * Supported optional flags:
 *
 *   - --seed <int>
 *   - --etas <comma-separated list>
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return 0 on success, non-zero on error.
 */
int main(int argc, char **argv) {
  if (argc < 2) {
    print_help(argv[0]);
    return 1;
  }

  const std::string first_arg = argv[1];

  if (first_arg == "--help" || first_arg == "-h") {
    print_help(argv[0]);
    return 0;
  }

  std::string cmd = first_arg;
  int seed = 12345;
  std::vector<double> etas;

  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];

    if (arg == "--seed") {
      if (i + 1 >= argc) {
        std::cerr << "Error: missing integer after --seed\n\n";
        print_help(argv[0]);
        return 1;
      }
      seed = std::stoi(argv[++i]);
    } else if (arg == "--etas") {
      if (i + 1 >= argc) {
        std::cerr << "Error: missing eta list after --etas\n\n";
        print_help(argv[0]);
        return 1;
      }
      etas = parse_etas(argv[++i]);
    } else if (arg == "--help" || arg == "-h") {
      print_help(argv[0]);
      return 0;
    } else {
      std::cerr << "Error: unknown option \"" << arg << "\"\n\n";
      print_help(argv[0]);
      return 1;
    }
  }

  std::mt19937 gen(static_cast<std::mt19937::result_type>(seed));

  if (cmd == "basic") {
    run_basic_analysis(gen);
  } else if (cmd == "cooling-evolution") {
    run_cooling_evolution_analysis(gen);
  } else if (cmd == "ensemble") {
    run_ensemble_analysis(gen);
  } else if (cmd == "fig7") {
    run_fig7_analysis(gen);
  } else if (cmd == "fig8") {
    if (etas.empty()) {
      etas = {1.0, 1.2, 1.4, 1.6, 1.8};
    }
    run_fig8_analysis(etas, gen);
  } else if (cmd == "fig9") {
    run_fig9_analysis(gen);
  } else if (cmd == "fig11") {
    run_fig11_analysis(1.4);
  } else if (cmd == "rilm") {
    run_rilm_analysis(gen);
  } else if (cmd == "heated-rilm") {
    run_heated_rilm_analysis(gen);
  } else if (cmd == "fig14") {
    run_fig14_analysis(gen);
  } else if (cmd == "fig15") {
    run_fig15_analysis(gen);
  } else if (cmd == "fig16") {
    run_fig16_analysis(gen);
  } else if (cmd == "iilm") {
    run_iilm_analysis(gen);
  } else if (cmd == "fig17") {
    run_fig17_analysis(gen);
  } else if (cmd == "qmidens") {
    run_qmidens_analysis(gen);
  } else if (cmd == "eta-scan") {
    if (etas.empty()) {
      etas = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
    }
    run_cooling_eta_scan(etas, gen);
  } else {
    std::cerr << "Error: unknown command \"" << cmd << "\"\n\n";
    print_help(argv[0]);
    return 1;
  }

  return 0;
}