#include "analysis/analysis_driver.hpp"

#include "analysis/fig10_rilm.hpp"
#include "analysis/fig12_13_heating.hpp"
#include "analysis/fig14_16.hpp"
#include "analysis/fig15.hpp"
#include "analysis/fig17_iilm.hpp"
#include "analysis/fig7.hpp"
#include "analysis/fig8.hpp"
#include "analysis/fig9.hpp"

#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @file main.cpp
 * @brief Command-line entry point for the instanton simulation executable.
 *
 * Parses command-line options, initializes the random number generator,
 * and dispatches the selected numerical analysis.
 */

/**
 * @brief Parsed command-line options.
 */
struct CliOptions {
  std::string command;
  int seed = 12345;
  std::vector<double> etas;
};

/**
 * @brief Print available commands and command-line options.
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
 * @param s String of the form "1.0,1.2,1.4".
 * @return Parsed eta values.
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
 * @brief Parse command-line arguments.
 *
 * Recognized options:
 *   - --seed <int>
 *   - --etas a,b,c
 *
 * @throws std::runtime_error on invalid options.
 */
static CliOptions parse_cli_options(int argc, char **argv) {
  CliOptions opt;

  if (argc < 2) {
    return opt;
  }

  opt.command = argv[1];

  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];

    if (arg == "--seed") {
      if (i + 1 >= argc) {
        throw std::runtime_error("missing integer after --seed");
      }
      opt.seed = std::stoi(argv[++i]);
    } else if (arg == "--etas") {
      if (i + 1 >= argc) {
        throw std::runtime_error("missing eta list after --etas");
      }
      opt.etas = parse_etas(argv[++i]);
    } else if (arg == "--help" || arg == "-h") {
      opt.command = "--help";
    } else {
      throw std::runtime_error("unknown option: " + arg);
    }
  }

  return opt;
}

/**
 * @brief Program entry point.
 *
 * @return Zero on success, non-zero on invalid input.
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

  CliOptions opt;

  try {
    opt = parse_cli_options(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n\n";
    print_help(argv[0]);
    return 1;
  }

  std::string cmd = opt.command;
  int seed = opt.seed;
  std::vector<double> etas = opt.etas;

  // ---------------- initialize RNG ----------------

  std::mt19937 gen(static_cast<std::mt19937::result_type>(seed));

  // ---------------- command dispatch ----------------

  std::unordered_map<std::string, std::function<void()>> commands;

  commands["basic"] = [&]() { run_basic_analysis(gen); };
  commands["cooling-evolution"] = [&]() {
    run_cooling_evolution_analysis(gen);
  };
  commands["ensemble"] = [&]() { run_ensemble_analysis(gen); };
  commands["fig7"] = [&]() { run_fig7_analysis(gen); };

  commands["fig8"] = [&]() {
    if (etas.empty()) {
      etas = {1.0, 1.2, 1.4, 1.6, 1.8};
    }
    run_fig8_analysis(etas, gen);
  };

  commands["fig9"] = [&]() { run_fig9_analysis(gen); };
  commands["fig11"] = [&]() { run_fig11_analysis(1.4); };
  commands["rilm"] = [&]() { run_rilm_analysis(gen); };
  commands["heated-rilm"] = [&]() { run_heated_rilm_analysis(gen); };
  commands["fig14"] = [&]() { run_fig14_analysis(gen); };
  commands["fig15"] = [&]() { run_fig15_analysis(gen); };
  commands["fig16"] = [&]() { run_fig16_analysis(gen); };
  commands["iilm"] = [&]() { run_iilm_analysis(gen); };
  commands["fig17"] = [&]() { run_fig17_analysis(gen); };
  commands["qmidens"] = [&]() { run_qmidens_analysis(gen); };

  commands["eta-scan"] = [&]() {
    if (etas.empty()) {
      etas = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
    }
    run_cooling_eta_scan(etas, gen);
  };

  // ----------------- execute command ----------------

  auto it = commands.find(cmd);

  if (it != commands.end()) {
    it->second();
  } else {
    std::cerr << "Error: unknown command \"" << cmd << "\"\n\n";
    print_help(argv[0]);
    return 1;
  }

  return 0;
}