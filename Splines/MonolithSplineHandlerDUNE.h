#pragma once

#include "Splines/SplineMonolith.h"

struct SplineHeader {
    std::string name;
    std::vector<double> knots;
    bool isCorrection;
    double cv = 0.0;
};

class MonolithSplineHandlerDUNE : virtual public SMonolith {
 public:
  MonolithSplineHandlerDUNE(const std::vector<SplineParameter>& splinePars, 
    const std::vector<uint>& eventIndices, const std::string& spline_filename);
  virtual ~MonolithSplineHandlerDUNE();

 private:
  void InitFromFile(std::string &spline_filename);
  MonolithSplineHandlerDUNE(std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> initParams);
  static std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> GetInitParamsFromConfig(
    const std::vector<SplineParameter>& splinePars, const std::vector<uint>& eventIndices, const std::string& spline_filename
  );
  static std::vector<struct SplineHeader> GetSplineParametersFromFile(TFile *f, const std::string& treeName);

};

//Forced to write a new version as the M3 core one delected the splines it is passed................
class TSpline3_redDUNE : public TSpline3_red {
public:
    TSpline3_redDUNE(const TSpline3* spline) : TSpline3_red() {
        nPoints = spline->GetNp();
        Par = new M3::float_t*[nPoints];
        XPos = new M3::float_t[nPoints];
        YResp = new M3::float_t[nPoints];

        for (int i = 0; i < nPoints; ++i) {
            // 3 is the size of the TSpline3 coefficients
            Par[i] = new M3::float_t[3];
            double x = -999.99, y = -999.99, b = -999.99, c = -999.99, d = -999.99;
            spline->GetCoeff(i, x, y, b, c, d);
            XPos[i]   = M3::float_t(x);
            YResp[i]  = M3::float_t(y);
            Par[i][0] = M3::float_t(b);
            Par[i][1] = M3::float_t(c);
            Par[i][2] = M3::float_t(d);
      }
    }

    //Adding a copy constructor because we cannot reuse the same object as M3 core once again loves to delete objects it's given...
    TSpline3_redDUNE(const TSpline3_redDUNE& other) : TSpline3_red() {
        nPoints = other.nPoints;
        Par = new M3::float_t*[nPoints];
        XPos = new M3::float_t[nPoints];
        YResp = new M3::float_t[nPoints];

        for (int i = 0; i < nPoints; ++i) {
            Par[i] = new M3::float_t[3];
            XPos[i] = other.XPos[i];
            YResp[i] = other.YResp[i];
            Par[i][0] = other.Par[i][0];
            Par[i][1] = other.Par[i][1];
            Par[i][2] = other.Par[i][2];
        }
    }
};

class ReusableSpline : public TSpline3 {
public:
    // Inherit standard constructors
    ReusableSpline(const std::vector<double>& knots) : TSpline3() {
        fNp = knots.size();
        fPoly = new TSplinePoly3[fNp];
        for (int i = 0; i < fNp; ++i) {
            fPoly[i].X() = knots[i];
        }
    }

    // Expose the protected BuildCoeff method
    void Refresh() {
        this->BuildCoeff();
    }

    void SetVariationWeights(const std::vector<double>& weights) {
        if(weights.size() != fNp) {
            MACH3LOG_ERROR("Size of weights does not match number of nodes in spline.");
            throw MaCh3Exception(__FILE__, __LINE__);
        }
        for (int i = 0; i < fNp; ++i) {
            fPoly[i].Y() = weights[i];
        }
    }

    void PrintPoly() const {
        for (int i = 0; i < fNp; ++i) {
            MACH3LOG_INFO("Knot {}: X = {}, Y = {}", i, fPoly[i].X(), fPoly[i].Y());
        }
    }

};