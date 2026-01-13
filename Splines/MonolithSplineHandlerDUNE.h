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

    bool IsFlat() const {
        for (int i = 0; i < nPoints; ++i) {
            if (std::abs(YResp[i] - 1.0) > 1e-9) {
                return false;
            }
        }
        return true;
    }
};

class ReusableSpline : public TSpline3 {
public:
    // Inherit standard constructors
    ReusableSpline(const std::vector<double>& knots) : TSpline3(){
        fValBeg = 0;
        fValEnd = 0;
        fBegCond = 0;
        fEndCond = 0;
        fNp = knots.size();
        fPoly = new TSplinePoly3[fNp];
        for (int i = 0; i < fNp; ++i) {
            fPoly[i].X() = knots[i];
            fPoly[i].Y() = 1.0;  // Initialize Y to 1.0 (neutral weight)
        }
        fXmin = fPoly[0].X();
        fXmax = fPoly[fNp-1].X();
        fName = "ReusableSpline";
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
            double x = 0, y = 0, b = 0, c = 0, d = 0;
            this->GetCoeff(i, x, y, b, c, d);
            MACH3LOG_INFO("Knot {}: X = {:.2f}, Y = {:.2f}, b = {:.2f}, c = {:.2f}, d = {:.2f}", i, x, y, b, c, d);
        }
    }

};