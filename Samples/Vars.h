#pragma once

#include "Var.h"
#include "Cut.h"
#include "CAFIncludes.h"

#include <map>
#include <string>
#include <cmath>
#include <TVector3.h>
#include "Manager/MaCh3Logger.h"
#include "Manager/MaCh3Exception.h"

namespace dune {

// ============================================================================
// Neutrino Variables
// ============================================================================

/// True neutrino energy
const Var kTrueNuE([](const caf::StandardRecordProxy* sr) {
    return static_cast<double>(sr->mc.nu[0].E);
});

/// True neutrino PDG code  
const Var kTrueNuPDG([](const caf::StandardRecordProxy* sr) {
    return static_cast<double>(sr->mc.nu[0].pdg);
});

/// Is CC interaction
const Var kIsCC([](const caf::StandardRecordProxy* sr) {
    return static_cast<double>(sr->mc.nu[0].iscc);
});

/// Interaction mode
const Var kMode([](const caf::StandardRecordProxy* sr) {
    return static_cast<double>(sr->mc.nu[0].mode);
});

/// True neutrino cosine zenith
const Var kTrueCosZ([](const caf::StandardRecordProxy* sr) {
    TVector3 TrueNuMomentumVector(sr->mc.nu[0].momentum.x,
                                   sr->mc.nu[0].momentum.y,
                                   sr->mc.nu[0].momentum.z);
    TrueNuMomentumVector = TrueNuMomentumVector.Unit();
    return TrueNuMomentumVector.Z();
});

// ============================================================================
// Reconstruction Variables
// ============================================================================

/// Number of Pandora slices
const Var kNPandoraSlices([](const caf::StandardRecordProxy* sr) {
    return static_cast<double>(sr->common.ixn.npandora);
});

/// Reconstructed energy (e-like)
const Var kRecoEnu_ELike([](const caf::StandardRecordProxy* sr) {
    if (sr->common.ixn.npandora != 1) return -999.0;
    return static_cast<double>(sr->common.ixn.pandora[0].Enu.e_calo);
});

/// Reconstructed energy (mu-like)
const Var kRecoEnu_MuLike([](const caf::StandardRecordProxy* sr) {
    if (sr->common.ixn.npandora != 1) return -999.0;
    return static_cast<double>(sr->common.ixn.pandora[0].Enu.lep_calo);
});

/// Reconstructed energy (nc-like)
const Var kRecoEnu_NCLike([](const caf::StandardRecordProxy* sr) {
    if (sr->common.ixn.npandora != 1) return -999.0;
    return static_cast<double>(sr->common.ixn.pandora[0].Enu.calo);
});

/// Reconstructed cosine Z (e-like)
const Var kRecoCosZ_ELike([](const caf::StandardRecordProxy* sr) {
    if (sr->common.ixn.npandora != 1) return -999.0;
    TVector3 recoDir(sr->common.ixn.pandora[0].dir.heshw.x,
                     sr->common.ixn.pandora[0].dir.heshw.y,
                     sr->common.ixn.pandora[0].dir.heshw.z);
    return recoDir.Unit().Z();
});

/// Reconstructed cosine Z (mu-like)
const Var kRecoCosZ_MuLike([](const caf::StandardRecordProxy* sr) {
    if (sr->common.ixn.npandora != 1) return -999.0;
    TVector3 recoDir(sr->common.ixn.pandora[0].dir.lngtrk.x,
                     sr->common.ixn.pandora[0].dir.lngtrk.y,
                     sr->common.ixn.pandora[0].dir.lngtrk.z);
    return recoDir.Unit().Z();
});

// ============================================================================
// Variable Registry for YAML Configuration
// ============================================================================

/// Map of variable names to Var objects for use in YAML configs
inline const std::map<std::string, Var>& GetVarRegistry() {
    static const std::map<std::string, Var> registry = {
        {"TrueNuE", kTrueNuE},
        {"TrueNuPDG", kTrueNuPDG},
        {"IsCC", kIsCC},
        {"Mode", kMode},
        {"TrueCosZ", kTrueCosZ},
        {"NPandoraSlices", kNPandoraSlices},
        {"RecoEnu_ELike", kRecoEnu_ELike},
        {"RecoEnu_MuLike", kRecoEnu_MuLike},
        {"RecoCosZ_ELike", kRecoCosZ_ELike},
        {"RecoCosZ_MuLike", kRecoCosZ_MuLike}
    };
    return registry;
}

/// Get a variable by name from the registry
inline Var GetVarByName(const std::string& name) {
    const auto& registry = GetVarRegistry();
    auto it = registry.find(name);
    if (it == registry.end()) {
        MACH3LOG_ERROR("Variable '{}' not found in registry", name);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    return it->second;
}

// ============================================================================
// Variable Expression Parser for YAML
// ============================================================================

/// Helper to check if a string is a number
inline bool IsNumber(const std::string& str) {
    if (str.empty()) return false;
    size_t start = 0;
    if (str[0] == '-' || str[0] == '+') start = 1;
    if (start >= str.length()) return false;
    
    bool hasDot = false;
    bool hasE = false;
    for (size_t i = start; i < str.length(); ++i) {
        if (str[i] == '.') {
            if (hasDot || hasE) return false;
            hasDot = true;
        } else if (str[i] == 'e' || str[i] == 'E') {
            if (hasE || i == start) return false;
            hasE = true;
            if (i + 1 < str.length() && (str[i+1] == '+' || str[i+1] == '-')) ++i;
        } else if (!std::isdigit(str[i])) {
            return false;
        }
    }
    return true;
}

/// Helper to trim whitespace
inline std::string TrimVar(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

/// Find matching closing parenthesis
inline size_t FindMatchingParenVar(const std::string& expr, size_t openPos) {
    int depth = 1;
    for (size_t i = openPos + 1; i < expr.length(); ++i) {
        if (expr[i] == '(') depth++;
        else if (expr[i] == ')') {
            depth--;
            if (depth == 0) return i;
        }
    }
    return std::string::npos;
}

/// Parse a variable expression from a string
/// Supports: "VarName", "123.4", "Var1 + Var2", "Var1 * 2.0", "Var1 / Var2", etc.
inline Var ParseVar(const std::string& expression) {
    std::string expr = TrimVar(expression);
    
    // Handle empty expression
    if (expr.empty()) {
        MACH3LOG_ERROR("Empty variable expression");
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    // Handle parentheses
    if (expr[0] == '(') {
        size_t closePos = FindMatchingParenVar(expr, 0);
        if (closePos == std::string::npos) {
            MACH3LOG_ERROR("Mismatched parentheses in variable expression: {}", expression);
            throw MaCh3Exception(__FILE__, __LINE__);
        }
        
        // If the entire expression is in parentheses, remove them
        if (closePos == expr.length() - 1) {
            return ParseVar(expr.substr(1, closePos - 1));
        }
    }
    
    // Find operators (scan from left to right, respecting parentheses)
    // Order: + and - (lowest precedence), then * and / (higher precedence)
    int parenDepth = 0;
    
    // First pass: look for + and -
    for (size_t i = 1; i < expr.length(); ++i) { // Start at 1 to skip unary minus
        if (expr[i] == '(') parenDepth++;
        else if (expr[i] == ')') parenDepth--;
        else if (parenDepth == 0) {
            if (expr[i] == '+') {
                std::string left = TrimVar(expr.substr(0, i));
                std::string right = TrimVar(expr.substr(i + 1));
                return ParseVar(left) + ParseVar(right);
            } else if (expr[i] == '-' && i > 0 && expr[i-1] != 'e' && expr[i-1] != 'E') {
                // Make sure it's not part of scientific notation
                std::string left = TrimVar(expr.substr(0, i));
                std::string right = TrimVar(expr.substr(i + 1));
                return ParseVar(left) - ParseVar(right);
            }
        }
    }
    
    // Second pass: look for * and /
    parenDepth = 0;
    for (size_t i = 1; i < expr.length(); ++i) {
        if (expr[i] == '(') parenDepth++;
        else if (expr[i] == ')') parenDepth--;
        else if (parenDepth == 0) {
            if (expr[i] == '*') {
                std::string left = TrimVar(expr.substr(0, i));
                std::string right = TrimVar(expr.substr(i + 1));
                return ParseVar(left) * ParseVar(right);
            } else if (expr[i] == '/') {
                std::string left = TrimVar(expr.substr(0, i));
                std::string right = TrimVar(expr.substr(i + 1));
                return ParseVar(left) / ParseVar(right);
            }
        }
    }
    
    // If we get here, it should be either a constant or a variable name
    if (IsNumber(expr)) {
        double value = std::stod(expr);
        return Var([value](const caf::StandardRecordProxy*) { return value; });
    }
    
    // Otherwise, it's a variable name
    return GetVarByName(expr);
}

} // namespace dune
