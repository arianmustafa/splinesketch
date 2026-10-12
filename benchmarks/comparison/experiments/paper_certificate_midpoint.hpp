#pragma once

#include <splinesketch/certificate_midpoint.hpp>

// Preserve the experiment entry points while sharing the public implementation.
namespace splinesketch::experimental {
using ::splinesketch::CertificateMidpointEstimate;
using ::splinesketch::certificate_midpoint;
using ::splinesketch::certificate_midpoint_rank;
using ::splinesketch::certificate_midpoint_uniform_error;
using ::splinesketch::certificate_midpoint_max_error;
using ::splinesketch::certificate_midpoint_select;
namespace midpoint_detail = ::splinesketch::detail::midpoint;
} // namespace splinesketch::experimental
