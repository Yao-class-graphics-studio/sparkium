#pragma once
#include "cstdint"
namespace sparks {
struct ObjectInfo {
  uint32_t vertex_offset;
  uint32_t index_offset;
  uint32_t num_faces;
};
}  // namespace sparks
