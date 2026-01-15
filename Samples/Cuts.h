#pragma once

#include "Cut.h"
#include "Vars.h"
#include <map>
#include <string>
#include <cmath>

namespace dune {

// ============================================================================
// Common Cuts
// ============================================================================

/// Accept all events
const Cut kNoCut([](const caf::StandardRecordProxy*) {
    return true;
});

/// Require exactly one Pandora slice
const Cut kOnePandoraSlice = (kNPandoraSlices == 1.0);

/// Require no NaN in reconstructed energy (e-like)
const Cut kValidRecoEnu_ELike([](const caf::StandardRecordProxy* sr) {
    if (sr->common.ixn.npandora != 1) return false;
    return !std::isnan(sr->common.ixn.pandora[0].Enu.e_calo);
});

/// Require no NaN in reconstructed energy (mu-like)
const Cut kValidRecoEnu_MuLike([](const caf::StandardRecordProxy* sr) {
    if (sr->common.ixn.npandora != 1) return false;
    return !std::isnan(sr->common.ixn.pandora[0].Enu.lep_calo);
});

/// Require no NaN in reconstructed cosZ (e-like)
const Cut kValidRecoCosZ_ELike([](const caf::StandardRecordProxy* sr) {
    if (sr->common.ixn.npandora != 1) return false;
    TVector3 recoDir(sr->common.ixn.pandora[0].dir.heshw.x,
                     sr->common.ixn.pandora[0].dir.heshw.y,
                     sr->common.ixn.pandora[0].dir.heshw.z);
    double cosZ = -recoDir.Unit().Y();
    return !std::isnan(cosZ);
});

/// Require no NaN in reconstructed cosZ (mu-like)
const Cut kValidRecoCosZ_MuLike([](const caf::StandardRecordProxy* sr) {
    if (sr->common.ixn.pandora.size() != 1) return false;
    TVector3 recoDir(sr->common.ixn.pandora[0].dir.lngtrk.x,
                     sr->common.ixn.pandora[0].dir.lngtrk.y,
                     sr->common.ixn.pandora[0].dir.lngtrk.z);
    double cosZ = -recoDir.Unit().Y();
    return !std::isnan(cosZ);
});

const Cut kIsMuLike([](const caf::StandardRecordProxy* sr) {
    if (sr->common.ixn.pandora.size() != 1) return false;
    return (
        (sr->common.ixn.pandora[0].nuhyp.cvn.numu > sr->common.ixn.pandora[0].nuhyp.cvn.nue) &&
        (sr->common.ixn.pandora[0].nuhyp.cvn.numu > sr->common.ixn.pandora[0].nuhyp.cvn.nc)
    );
});

const Cut kIsELike([](const caf::StandardRecordProxy* sr) {
    if (sr->common.ixn.pandora.size() != 1) return false;
    return (
        (sr->common.ixn.pandora[0].nuhyp.cvn.nue > sr->common.ixn.pandora[0].nuhyp.cvn.numu) &&
        (sr->common.ixn.pandora[0].nuhyp.cvn.nue > sr->common.ixn.pandora[0].nuhyp.cvn.nc)
    );
});

/// Standard quality cuts for e-like selection
const Cut kStandardQuality_ELike = kOnePandoraSlice && kValidRecoEnu_ELike && kValidRecoCosZ_ELike;

/// Standard quality cuts for mu-like selection
const Cut kStandardQuality_MuLike = kOnePandoraSlice && kValidRecoEnu_MuLike && kValidRecoCosZ_MuLike;

// ============================================================================
// Cut Registry for YAML Configuration
// ============================================================================

/// Map of cut names to Cut objects for use in YAML configs
inline const std::map<std::string, Cut>& GetCutRegistry() {
    static const std::map<std::string, Cut> registry = {
        {"NoCut", kNoCut},
        {"OnePandoraSlice", kOnePandoraSlice},
        {"ValidRecoEnu_ELike", kValidRecoEnu_ELike},
        {"ValidRecoEnu_MuLike", kValidRecoEnu_MuLike},
        {"ValidRecoCosZ_ELike", kValidRecoCosZ_ELike},
        {"ValidRecoCosZ_MuLike", kValidRecoCosZ_MuLike},
        {"StandardQuality_ELike", kStandardQuality_ELike},
        {"StandardQuality_MuLike", kStandardQuality_MuLike},
        {"IsELike", kIsELike},
        {"IsMuLike", kIsMuLike}
    };
    return registry;
}

/// Get a cut by name from the registry
inline Cut GetCutByName(const std::string& name) {
    const auto& registry = GetCutRegistry();
    auto it = registry.find(name);
    if (it == registry.end()) {
        MACH3LOG_ERROR("Cut '{}' not found in registry", name);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    return it->second;
}

// ============================================================================
// Cut Expression Parser for YAML
// ============================================================================

/// Helper to trim whitespace from both ends of a string
inline std::string Trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

/// Find matching closing parenthesis
inline size_t FindMatchingParen(const std::string& expr, size_t openPos) {
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

/// Parse a cut expression from a string
/// Supports: "CutName", "Cut1 && Cut2", "Cut1 || Cut2", "!CutName", "(Cut1 || Cut2) && Cut3"
inline Cut ParseCut(const std::string& expression) {
    std::string expr = Trim(expression);
    
    // Handle empty expression
    if (expr.empty()) {
        MACH3LOG_ERROR("Empty cut expression");
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    // Handle NOT operator (!)
    if (expr[0] == '!') {
        std::string rest = Trim(expr.substr(1));
        return !ParseCut(rest);
    }
    
    // Handle parentheses
    if (expr[0] == '(') {
        size_t closePos = FindMatchingParen(expr, 0);
        if (closePos == std::string::npos) {
            MACH3LOG_ERROR("Mismatched parentheses in cut expression: {}", expression);
            throw MaCh3Exception(__FILE__, __LINE__);
        }
        
        // If the entire expression is in parentheses, remove them
        if (closePos == expr.length() - 1) {
            return ParseCut(expr.substr(1, closePos - 1));
        }
    }
    
    // Find OR operator (||) - lowest precedence
    int parenDepth = 0;
    for (size_t i = 0; i < expr.length() - 1; ++i) {
        if (expr[i] == '(') parenDepth++;
        else if (expr[i] == ')') parenDepth--;
        else if (parenDepth == 0 && expr[i] == '|' && expr[i + 1] == '|') {
            std::string left = Trim(expr.substr(0, i));
            std::string right = Trim(expr.substr(i + 2));
            return ParseCut(left) || ParseCut(right);
        }
    }
    
    // Find AND operator (&&) - higher precedence than OR
    parenDepth = 0;
    for (size_t i = 0; i < expr.length() - 1; ++i) {
        if (expr[i] == '(') parenDepth++;
        else if (expr[i] == ')') parenDepth--;
        else if (parenDepth == 0 && expr[i] == '&' && expr[i + 1] == '&') {
            std::string left = Trim(expr.substr(0, i));
            std::string right = Trim(expr.substr(i + 2));
            return ParseCut(left) && ParseCut(right);
        }
    }
    
    // If we get here, it should be a simple cut name
    return GetCutByName(expr);
}

} // namespace dune
