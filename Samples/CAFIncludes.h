#pragma once

// This header provides a single point of inclusion for CAF (Common Analysis Format) headers
// to avoid multiple definition errors from non-inline functions in BasicTypesProxy.txx

#ifndef DUNE_CAF_INCLUDES_GUARD
#define DUNE_CAF_INCLUDES_GUARD

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#pragma GCC diagnostic pop

#endif // DUNE_CAF_INCLUDES_GUARD
