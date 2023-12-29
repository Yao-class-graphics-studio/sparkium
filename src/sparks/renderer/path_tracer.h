#pragma once
#include "random"
#include "sparks/assets/scene.h"
#include "sparks/assets/bssrdf.h"
#include "sparks/renderer/renderer_settings.h"
#include <stack>

namespace sparks {
class PathTracer {
 public:
  PathTracer(const RendererSettings *render_settings, const Scene *scene);
  [[nodiscard]] glm::vec3 SampleRay(glm::vec3 origin,
                                    glm::vec3 direction,
                                    int x,
                                    int y,
                                    int sample,
                                    int bounces = 0,
                                    float currentRatio = 1.0f);
  void SampleFromLight(glm::vec3 &res, glm::vec3 &norm, float &area, int except);
  float getPdfByLight(glm::vec3 pos, glm::vec3 sample, float area);
  int shadowRay(glm::vec3 pos, glm::vec3 &dir, glm::vec3 sample, glm::vec3 &throughput);
  glm::vec3 directIllumination(glm::vec3 pos, float &lightPdf, glm::vec3 &dir, float &lightArea, int except);
  void initialMedium(glm::vec3 pos);
  void updateStack(std::stack<Medium*> &st, HitRecord hit, glm::vec3 sample);
 private:
  const RendererSettings *render_settings_{};
  const Scene *scene_{};
  std::mt19937 rd;
  std::uniform_real_distribution<float> uniform;
  const float continueProb = 0.85f;
  std::stack<Medium*> mediumStack;
};
}  // namespace sparks
