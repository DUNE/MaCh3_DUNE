#pragma once

#include <functional>
#include <memory>

// Forward declarations to avoid including full CAF headers
namespace caf {
    template <typename T> class Proxy;
    class StandardRecord;
    using StandardRecordProxy = Proxy<StandardRecord>;
}

#include "Var.h"

namespace dune {

/// Boolean selection criteria that can be applied to StandardRecords
class Cut {
public:
    /// Function type for cut evaluation
    typedef std::function<bool(const caf::StandardRecordProxy*)> CutFunc_t;

    /// Default constructor (creates always-true cut)
    Cut() : fFunc([](const caf::StandardRecordProxy*) { return true; }) {}

    /// Construct from lambda or function
    explicit Cut(const CutFunc_t& func) : fFunc(func) {}

    /// Evaluate the cut for a given record
    bool operator()(const caf::StandardRecordProxy* sr) const {
        return fFunc(sr);
    }

    // Boolean combination operators
    Cut operator&&(const Cut& other) const;
    Cut operator||(const Cut& other) const;
    Cut operator!() const;

protected:
    CutFunc_t fFunc;

    friend class Var;
};

// Implementation of boolean operators
inline Cut Cut::operator&&(const Cut& other) const {
    return Cut([*this, other](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) && other(sr);
    });
}

inline Cut Cut::operator||(const Cut& other) const {
    return Cut([*this, other](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) || other(sr);
    });
}

inline Cut Cut::operator!() const {
    return Cut([*this](const caf::StandardRecordProxy* sr) {
        return !(*this)(sr);
    });
}

// Implementation of Var comparison operators
inline Cut Var::operator<(double val) const {
    return Cut([*this, val](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) < val;
    });
}

inline Cut Var::operator>(double val) const {
    return Cut([*this, val](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) > val;
    });
}

inline Cut Var::operator<=(double val) const {
    return Cut([*this, val](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) <= val;
    });
}

inline Cut Var::operator>=(double val) const {
    return Cut([*this, val](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) >= val;
    });
}

inline Cut Var::operator==(double val) const {
    return Cut([*this, val](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) == val;
    });
}

inline Cut Var::operator!=(double val) const {
    return Cut([*this, val](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) != val;
    });
}

// Implementation of Var arithmetic operators
inline Var Var::operator+(const Var& other) const {
    return Var([*this, other](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) + other(sr);
    });
}

inline Var Var::operator-(const Var& other) const {
    return Var([*this, other](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) - other(sr);
    });
}

inline Var Var::operator*(const Var& other) const {
    return Var([*this, other](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) * other(sr);
    });
}

inline Var Var::operator/(const Var& other) const {
    return Var([*this, other](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) / other(sr);
    });
}

// Implementation of Var arithmetic operators with scalars
inline Var Var::operator+(double val) const {
    return Var([*this, val](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) + val;
    });
}

inline Var Var::operator-(double val) const {
    return Var([*this, val](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) - val;
    });
}

inline Var Var::operator*(double val) const {
    return Var([*this, val](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) * val;
    });
}

inline Var Var::operator/(double val) const {
    return Var([*this, val](const caf::StandardRecordProxy* sr) {
        return (*this)(sr) / val;
    });
}

// Allow scalar arithmetic with constant on the left (double op Var)
inline Var operator+(double val, const Var& var) { return var + val; }
inline Var operator-(double val, const Var& var) {
    return Var([val, var](const caf::StandardRecordProxy* sr) {
        return val - var(sr);
    });
}
inline Var operator*(double val, const Var& var) { return var * val; }
inline Var operator/(double val, const Var& var) {
    return Var([val, var](const caf::StandardRecordProxy* sr) {
        return val / var(sr);
    });
}

// Allow comparisons with constant on the left
inline Cut operator>(double val, const Var& var) { return var < val; }
inline Cut operator<(double val, const Var& var) { return var > val; }
inline Cut operator>=(double val, const Var& var) { return var <= val; }
inline Cut operator<=(double val, const Var& var) { return var >= val; }
inline Cut operator==(double val, const Var& var) { return var == val; }
inline Cut operator!=(double val, const Var& var) { return var != val; }

} // namespace dune
