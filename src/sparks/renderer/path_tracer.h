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
                                    int sample,
                                    bool moreBounces = true,
                                    float weight = 1.0f);
  void SampleFromLight(glm::vec3 &res, glm::vec3 &norm, float &pdf);
  void SampleFromCosine(glm::vec3 &res, float &pdf, glm::vec3 localNorm);
  void MultipleImportanceSample(glm::vec3 &res, float &pdf);
 private:
  const RendererSettings *render_settings_{};
  const Scene *scene_{};
  std::mt19937 rd;
  std::uniform_real_distribution<float> uniform;
  const float continueProb = 0.8f;
};
}  // namespace sparks
