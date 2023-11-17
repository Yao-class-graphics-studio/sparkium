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
                                    float weight = 1.0f,
                                    int bounces = 0);
  void SampleFromLight(glm::vec3 &res, glm::vec3 &norm, float &area, int except);
  void SampleFromCosine(glm::vec3 &res, float &pdf, glm::vec3 localNorm);
  void Sample(glm::vec3 &res, float &pdf, glm::vec3 in, glm::vec3 localNorm, Material &material);
  float getPdfByMaterial(glm::vec3 sample, glm::vec3 in, glm::vec3 norm, Material &material);
  float getPdfByLight(glm::vec3 pos, glm::vec3 sample, float area);
  glm::vec3 getBRDF(Material &material, HitRecord &hit, glm::vec3 norm, glm::vec3 in, glm::vec3 out);
 private:
  const RendererSettings *render_settings_{};
  const Scene *scene_{};
  std::mt19937 rd;
  std::uniform_real_distribution<float> uniform;
  const float continueProb = 0.8f;
};
}  // namespace sparks
