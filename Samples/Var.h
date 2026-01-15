#pragma once

#include <functional>
#include <memory>

// Forward declarations to avoid including full CAF headers
namespace caf {
    template <typename T> class Proxy;
    class StandardRecord;
    using StandardRecordProxy = Proxy<StandardRecord>;
}

namespace dune {

// Forward declaration
class Cut;

/// Variable that extracts a double value from a StandardRecord
class Var {
public:
    /// Function type for variable extraction
    typedef std::function<double(const caf::StandardRecordProxy*)> VarFunc_t;

    /// Default constructor (creates invalid variable)
    Var() : fFunc(nullptr) {}

    /// Construct from lambda or function
    explicit Var(const VarFunc_t& func) : fFunc(func) {}

    /// Evaluate the variable for a given record
    double operator()(const caf::StandardRecordProxy* sr) const {
        return fFunc(sr);
    }

    // Comparison operators that return Cuts
    Cut operator<(double val) const;
    Cut operator>(double val) const;
    Cut operator<=(double val) const;
    Cut operator>=(double val) const;
    Cut operator==(double val) const;
    Cut operator!=(double val) const;

    // Arithmetic operators (Var op Var)
    Var operator+(const Var& other) const;
    Var operator-(const Var& other) const;
    Var operator*(const Var& other) const;
    Var operator/(const Var& other) const;
    
    // Arithmetic operators with scalars (Var op double)
    Var operator+(double val) const;
    Var operator-(double val) const;
    Var operator*(double val) const;
    Var operator/(double val) const;

protected:
    VarFunc_t fFunc;
};

} // namespace dune
