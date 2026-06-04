#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct ParameterConfig {
  std::string name;
  double generator = 1.0;
};

std::string ProcessName(const std::string& parameter) {
  const auto pos = parameter.find('_');
  return pos == std::string::npos ? parameter : parameter.substr(0, pos);
}

std::string CsvEscape(const std::string& value) {
  if (value.find_first_of(",\"\n") == std::string::npos) {
    return value;
  }
  std::string escaped = "\"";
  for (char c : value) {
    if (c == '"') escaped += "\"\"";
    else escaped += c;
  }
  escaped += "\"";
  return escaped;
}

std::vector<ParameterConfig> ReadPDSPFitModel(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("Could not open xsec config: " + path);
  }

  std::regex parameter_re("^\\s*ParameterName:\\s*(\\S+)\\s*$");
  std::regex generator_re("^\\s*Generator:\\s*([-+0-9.eE]+)\\s*$");
  std::smatch match;
  std::vector<ParameterConfig> parameters;
  int current = -1;

  std::string line;
  while (std::getline(input, line)) {
    if (std::regex_match(line, match, parameter_re)) {
      parameters.push_back({match[1].str(), 1.0});
      current = static_cast<int>(parameters.size()) - 1;
    } else if (current >= 0 && std::regex_match(line, match, generator_re)) {
      parameters[current].generator = std::stod(match[1].str());
    }
  }

  return parameters;
}

std::set<std::string> SplitCsv(const std::string& raw) {
  std::set<std::string> values;
  std::stringstream stream(raw);
  std::string item;
  while (std::getline(stream, item, ',')) {
    const auto first = item.find_first_not_of(" \t");
    const auto last = item.find_last_not_of(" \t");
    if (first != std::string::npos) {
      values.insert(item.substr(first, last - first + 1));
    }
  }
  return values;
}

bool ShouldSummarise(const ParameterConfig& parameter, const std::set<std::string>& targets) {
  if (targets.empty()) return std::abs(parameter.generator - 1.0) > 1e-12;
  return targets.count(parameter.name) > 0 || targets.count(ProcessName(parameter.name)) > 0;
}

} // namespace

void summarise_pdsp_fit_recovery(const char* fit_file,
                                 const char* xsec_config,
                                 const char* output_csv,
                                 const char* label = "",
                                 const char* targets_csv = "",
                                 int burn_in = 10000,
                                 double tolerance = 0.10) {
  const std::vector<ParameterConfig> parameters = ReadPDSPFitModel(xsec_config);
  const std::set<std::string> targets = SplitCsv(targets_csv);

  TFile file(fit_file, "READ");
  if (file.IsZombie()) {
    throw std::runtime_error(std::string("Could not open fit file: ") + fit_file);
  }

  auto* tree = dynamic_cast<TTree*>(file.Get("posteriors"));
  if (tree == nullptr) {
    throw std::runtime_error(std::string("Could not find posteriors tree in: ") + fit_file);
  }

  const Long64_t entries = tree->GetEntries();
  const Long64_t first_entry = std::min<Long64_t>(std::max(0, burn_in), entries);
  const Long64_t used_entries = entries - first_entry;
  if (used_entries <= 0) {
    throw std::runtime_error("Burn-in removes all posterior entries");
  }

  std::ofstream output(output_csv);
  if (!output) {
    throw std::runtime_error(std::string("Could not write summary CSV: ") + output_csv);
  }

  output << "label,parameter,process,injected,posterior_mean,posterior_rms,"
         << "abs_error,rel_error,pull,entries_used,tolerance,pass\n";

  for (std::size_t i = 0; i < parameters.size(); ++i) {
    const auto& parameter = parameters[i];
    if (!ShouldSummarise(parameter, targets)) continue;

    const std::string branch_name = "param_" + std::to_string(i);
    if (tree->GetBranch(branch_name.c_str()) == nullptr) {
      std::cerr << "Missing branch " << branch_name << " for " << parameter.name << std::endl;
      continue;
    }

    double value = 0.0;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus(branch_name.c_str(), 1);
    tree->SetBranchAddress(branch_name.c_str(), &value);

    double sum = 0.0;
    double sum2 = 0.0;
    for (Long64_t entry = first_entry; entry < entries; ++entry) {
      tree->GetEntry(entry);
      sum += value;
      sum2 += value * value;
    }

    const double mean = sum / used_entries;
    const double variance = std::max(0.0, sum2 / used_entries - mean * mean);
    const double rms = std::sqrt(variance);
    const double abs_error = mean - parameter.generator;
    const double rel_error = parameter.generator == 0.0 ? 0.0 : abs_error / parameter.generator;
    const double pull = rms == 0.0 ? 0.0 : abs_error / rms;
    const bool pass = std::abs(rel_error) <= tolerance;

    output << CsvEscape(label) << ','
           << CsvEscape(parameter.name) << ','
           << CsvEscape(ProcessName(parameter.name)) << ','
           << std::setprecision(12) << parameter.generator << ','
           << mean << ','
           << rms << ','
           << abs_error << ','
           << rel_error << ','
           << pull << ','
           << used_entries << ','
           << tolerance << ','
           << (pass ? "true" : "false") << '\n';
  }
}
