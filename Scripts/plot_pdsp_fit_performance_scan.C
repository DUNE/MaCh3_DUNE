#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Row {
  std::string scan_target;
  std::string parameter;
  std::string process;
  double injected_value = 0.0;
  double injected = 0.0;
  double rel_error = 0.0;
  double posterior_rms = 0.0;
  bool pass = false;
};

std::vector<std::string> ParseCsvLine(const std::string& line) {
  std::vector<std::string> fields;
  std::string current;
  bool quoted = false;
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (quoted) {
      if (c == '"' && i + 1 < line.size() && line[i + 1] == '"') {
        current += '"';
        ++i;
      } else if (c == '"') {
        quoted = false;
      } else {
        current += c;
      }
    } else if (c == '"') {
      quoted = true;
    } else if (c == ',') {
      fields.push_back(current);
      current.clear();
    } else {
      current += c;
    }
  }
  fields.push_back(current);
  return fields;
}

std::map<std::string, int> HeaderIndex(const std::vector<std::string>& header) {
  std::map<std::string, int> index;
  for (std::size_t i = 0; i < header.size(); ++i) {
    index[header[i]] = static_cast<int>(i);
  }
  return index;
}

std::string Get(const std::vector<std::string>& fields,
                const std::map<std::string, int>& index,
                const std::string& name) {
  const auto iter = index.find(name);
  if (iter == index.end() || iter->second >= static_cast<int>(fields.size())) {
    return "";
  }
  return fields[iter->second];
}

std::vector<Row> ReadRows(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("Could not open performance summary CSV: " + path);
  }

  std::string line;
  if (!std::getline(input, line)) {
    throw std::runtime_error("Empty performance summary CSV: " + path);
  }
  const auto index = HeaderIndex(ParseCsvLine(line));

  std::vector<Row> rows;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    const auto fields = ParseCsvLine(line);
    Row row;
    row.scan_target = Get(fields, index, "scan_target");
    row.parameter = Get(fields, index, "parameter");
    row.process = Get(fields, index, "process");
    row.injected_value = std::stod(Get(fields, index, "injected_value"));
    row.injected = std::stod(Get(fields, index, "injected"));
    row.rel_error = std::stod(Get(fields, index, "rel_error"));
    row.posterior_rms = std::stod(Get(fields, index, "posterior_rms"));
    row.pass = Get(fields, index, "pass") == "true";
    rows.push_back(row);
  }
  return rows;
}

std::string Sanitise(const std::string& raw) {
  std::string output;
  for (char c : raw) {
    if (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-') {
      output += c;
    } else {
      output += '_';
    }
  }
  return output.empty() ? "all" : output;
}

std::vector<std::string> UniqueTargets(const std::vector<Row>& rows, const std::string& requested) {
  std::set<std::string> targets;
  for (const auto& row : rows) {
    if (!requested.empty() && row.scan_target != requested) continue;
    targets.insert(row.scan_target);
  }
  return {targets.begin(), targets.end()};
}

int Colour(int index) {
  static const std::vector<int> colours = {
    kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7,
    kCyan + 2, kViolet + 1, kTeal + 3, kPink + 7, kAzure + 7
  };
  return colours[index % colours.size()];
}

void PlotTarget(const std::vector<Row>& rows,
                const std::string& target,
                const std::string& output_dir,
                double tolerance,
                TFile& root_file) {
  std::map<std::string, std::vector<Row>> by_parameter;
  double xmin = 1e9;
  double xmax = -1e9;
  double ymax = tolerance * 1.25;

  for (const auto& row : rows) {
    if (row.scan_target != target) continue;
    by_parameter[row.parameter].push_back(row);
    xmin = std::min(xmin, row.injected_value);
    xmax = std::max(xmax, row.injected_value);
    const double y = std::abs(row.rel_error);
    const double yerr = row.injected == 0.0 ? 0.0 : row.posterior_rms / row.injected;
    ymax = std::max(ymax, (y + yerr) * 1.25);
  }

  if (by_parameter.empty()) return;
  if (xmin == xmax) {
    xmin -= 0.05;
    xmax += 0.05;
  }

  root_file.cd();
  auto* canvas = new TCanvas(("c_" + Sanitise(target)).c_str(), target.c_str(), 1100, 800);
  canvas->SetGrid();
  auto* multigraph = new TMultiGraph();
  auto* legend = new TLegend(0.62, 0.58, 0.90, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  int index = 0;
  for (auto& [parameter, points] : by_parameter) {
    std::sort(points.begin(), points.end(), [](const Row& a, const Row& b) {
      return a.injected_value < b.injected_value;
    });

    auto* graph = new TGraphErrors(static_cast<int>(points.size()));
    graph->SetName(("g_" + Sanitise(target) + "_" + Sanitise(parameter)).c_str());
    graph->SetTitle(parameter.c_str());
    graph->SetLineColor(Colour(index));
    graph->SetMarkerColor(Colour(index));
    graph->SetMarkerStyle(20 + (index % 10));
    graph->SetLineWidth(2);

    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
      const double y = std::abs(points[i].rel_error);
      const double yerr = points[i].injected == 0.0 ? 0.0 : points[i].posterior_rms / points[i].injected;
      graph->SetPoint(i, points[i].injected_value, y);
      graph->SetPointError(i, 0.0, yerr);
    }
    multigraph->Add(graph, "LP");
    legend->AddEntry(graph, parameter.c_str(), "lp");
    ++index;
  }

  multigraph->SetTitle((target + " fit recovery degradation;Injected Generator normalisation;|posterior mean - injected| / injected").c_str());
  multigraph->Draw("A");
  multigraph->GetXaxis()->SetLimits(xmin - 0.02, xmax + 0.02);
  multigraph->SetMinimum(0.0);
  multigraph->SetMaximum(ymax);
  multigraph->GetXaxis()->SetTitleOffset(1.2);
  multigraph->GetYaxis()->SetTitleOffset(1.35);

  auto* line = new TLine(xmin - 0.02, tolerance, xmax + 0.02, tolerance);
  line->SetLineColor(kRed + 1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw("same");
  legend->AddEntry(line, "degradation threshold", "l");
  legend->AddEntry((TObject*)nullptr, "error bars: posterior RMS / injected", "");
  legend->Draw();

  const std::string base = output_dir + "/fit_performance_degradation_" + Sanitise(target);
  canvas->Write();
  multigraph->Write(("mg_" + Sanitise(target)).c_str());
  canvas->SaveAs((base + ".pdf").c_str());
  canvas->SaveAs((base + ".png").c_str());
}

} // namespace

void plot_pdsp_fit_performance_scan(const char* summary_csv = "PDSPFitPerformanceScan/performance_summary.csv",
                                    const char* output_dir = "PDSPFitPerformanceScan/plots",
                                    double tolerance = 0.10,
                                    const char* target = "") {
  gStyle->SetOptStat(0);
  std::filesystem::create_directories(output_dir);

  const auto rows = ReadRows(summary_csv);
  const auto targets = UniqueTargets(rows, target);
  if (targets.empty()) {
    throw std::runtime_error(std::string("No rows found for requested target: ") + target);
  }

  TFile root_file((std::string(output_dir) + "/fit_performance_degradation.root").c_str(), "RECREATE");
  for (const auto& scan_target : targets) {
    PlotTarget(rows, scan_target, output_dir, tolerance, root_file);
  }
  root_file.Close();
}
