#pragma once

#include "Eigen/Core"

class SampleSpleigen : public SplineBase {

  struct Spleigen {
    std::vector<double> knot_points;
    Eigen::ArrayXXd coeffs;

    // from
    // https://people.clas.ufl.edu/kees/files/CubicSplines.pdf
    Spleigen(Eigen::ArrayXd const &x, Eigen::ArrayXXd const &y)
        : knot_points(x.data(), x.data() + x.size()) {
      if (y.cols() != x.size()) {
        throw std::runtime_error(
            "Spleigen instantiated with " + std::to_string(y.size()) +
            " knots, but "
            "attemping to build spline with yvals matrix " +
            std::to_string(x.cols()) + " columns wide.");
      }

      int n = x.size();

      coeffs = Eigen::ArrayXXd::Zero(y.rows(), 4 * (n - 1));

      Eigen::RowArrayXd h = (x.tail(n - 1) - x.head(n - 1));

      Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);

      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          if (i == j) {

            if ((i == 0) || ((i + 1) == n)) {
              A(i, j) = 1;
            } else {
              A(i, j) = 2.0 * (h(i - 1) + h(i));
            }

          } else if (((i + 1) == j) && (i != 0)) { // upper diag
            A(i, j) = h(i);
          } else if (((i - 1) == j) && ((i + 1) != n)) { // lower diag
            A(i, j) = h(i - 1);
          }
        }
      }

      Eigen::RowArrayXd b, c, d;

      Eigen::ColArrayXd alpha = Eigen::RowArrayXd::Zero(n);

      for (int row = 0; row < y.rows(); ++row) {

        // solve Ac = alpha for c, where alpha
        alpha.segment(1, n - 2) =
            (3.0 / h.tail(n - 2)) *
                (y.row(row).tail(n - 2) - y.row(row).segment(1, n - 2)) -
            (3.0 / h.head(n - 2)) *
                (y.row(row).segment(1, n - 2) - y.row(row).head(n - 2));

        // solve for c
        c = A.colPivHouseholderQr().solve(alpha.matrix()).array();

        d = (c.tail(n - 1) - c.head(n - 1)) / (3.0 * h);
        b = ((y.row(row).tail(n - 1) - y.row(row).head(n - 1)) / h) -
            ((h / 3.0) * (2.0 * c.head(n - 1) + c.tail(n - 1)));

        coeffs.row(row).segment(0, n - 1) = y.row(row).head(n - 1);
        coeffs.row(row).segment(n - 1, n - 1) = b;
        coeffs.row(row).segment(2 * (n - 1), n - 1) = c.head(n - 1).cast<P>();
        coeffs.row(row).segment(3 * (n - 1), n - 1) = d;
      }
    }
  };

  std::string SampleName;
  std::vector<Spleigen> param_blocks;

  Eigen::ArrayXd result;

  Eigen::ArrayXd const &Evaluate(Eigen::ArrayXd const &param_vals) {
    result.setOnes();
    for (size_t ip = 0; ip < param_blocks.size(); ++ip) {
      size_t i = std::distance(
          param_blocks[ip].knot_points.begin(),
          std::upper_bound(param_vals[ip], param_blocks[ip].knot_points.begin(),
                           param_blocks[ip].knot_points.end()));

      double past_knot = val - knot_points[i];

      result *= coeffs.col(i) + past_knot * coeffs.col((num_knots - 1) + i) +
                std::pow(past_knot, 2) * coeffs.col(2 * (num_knots - 1) + i) +
                std::pow(past_knot, 3) * coeffs.col(3 * (num_knots - 1) + i);
    }
    return result;
  }
};

SampleSpleigen ConstructSpleigen(std::string const &file_name,
                                 std::string const &sample_name) {
  std::unique_ptr<TFile> fin(TFile::Open(file_name.c_str(), "OPEN"));

  if (!fin) {
    throw std::runtime_error(std::string("Failed to open root file: ") +
                             file_name);
  }

  std::vector<std::string> *sample_list = nullptr;
  fin->GetObject("sample_list", sample_list);

  for (auto const &sample : *sample_list) {

    if (sample_name != sample) {
      continue;
    }

    SampleSpleigen spl{sample_name};

    std::vector<std::vector<double>> *bin_edges = nullptr;
    fin->GetObject(sample + "/bin_edges", bin_edges);

    std::vector<std::string> *parameter_list = nullptr;
    fin->GetObject(sample + "/parameter_list", parameter_list);

    for (auto const &param : *parameter_list) {

      std::vector<double> *knot_list = nullptr;
      fin->GetObject(sample + "/" + param + "/knot_list", knot_list);

      Eigen::ArrayXXd y =
          Eigen::ArrayXXd::Zero(bin_edges->size(), knot_list->size());

      for (size_t kti = 0; kti < knot_list->size(); ++kti) {
        auto knot_prediction = fin->Get<TH1>(
            sample + "/" + param + "/knot_prediction_" + std::to_string(kti));

        y.col(kti) = Eigen::Map<ArrayXd>(knot_prediction->GetArray() + 1,
                                         knot_prediction->GetNbins());
      }

      spl.param_blocks.emplace_back(
          Eigen::Map<ArrayXd>(knot_list->data(), knot_list->size()), y);
    }

    return spl;
  }

  throw std::runtime_error(
      std::string("Failed to find folder for requested sample: ") +
      sample_name + " in file: " + file_name);
}
