#pragma once
#include "random"
#include "sparks/assets/scene.h"
#include "sparks/renderer/renderer_settings.h"

namespace sparks {
class PathTracer {
 public:
  PathTracer(const RendererSettings *render_settings, const Scene *scene);
  [[nodiscard]] glm::vec3 SampleRay(glm::vec3 origin,
                                    glm::vec3 direction,
                                    int x,
                                    int y,
                                    int sample) const;
  glm::vec3 SampleRay_test(glm::vec3 origin,
                      glm::vec3 direction,
                      int x,
                      int y,
                      int sample,std::mt19937 rd) const;

 private:
  const RendererSettings *render_settings_{};
  const Scene *scene_{};

};
}  // namespace sparks
