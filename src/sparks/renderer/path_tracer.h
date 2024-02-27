#pragma once
#include "random"
#include "sparks/assets/scene.h"
#include "sparks/assets/bssrdf.h"
#include "sparks/renderer/renderer_settings.h"

namespace sparks {

// Jingyi Lyu
// Stack for nested medium
class MyStack {
  int cnt;
  Medium *data[10];
public:
  bool flag;
  MyStack() {
    cnt = 0;
    for(int i = 0; i <= 9; i++)
      data[i] = nullptr;
    flag = false;
  }
  
  ~MyStack() {
    for(int i = 0; i <= 9; i++)
      data[i] = nullptr;
  }
  
  MyStack(const MyStack &t) {
    cnt = t.cnt;
    flag = true;
    for(int i = 0; i <= 9; i++)
      data[i] = t.data[i];
  }

  void push(Medium* newMedium) {
    data[++cnt] = newMedium;
    // std::cerr << cnt << std::endl;
  }

  void pop() {
    data[cnt--] = nullptr;
    // std::cerr << cnt << std::endl;
  }

  Medium* top() {
    if(cnt < 0)
      return nullptr;
    return data[cnt];
  }

  bool empty() {
    return cnt <= 0;
  }

  int size() {
    return cnt;
  }

  void clear() {
    cnt = 0;
    for(int i = 0; i <= 9; i++)
      data[i] = nullptr;
  }
};

class PathTracer {
 public:
  PathTracer(const RendererSettings *render_settings, const Scene *scene);
  [[nodiscard]] glm::vec3 SampleRay(glm::vec3 origin,
                                    glm::vec3 direction,
                                    int x,
                                    int y,
                                    float sampleTime,
                                    int sample,
                                    int bounces = 0,
                                    bool initialized = false,
                                    float currentRatio = 1.0f);
  void SampleFromLight(glm::vec3 &res, glm::vec3 &norm, float &area, int except);
  float getPdfByLight(glm::vec3 pos, glm::vec3 sample, float area, float sampleTime);
  int shadowRay(glm::vec3 pos, glm::vec3 &dir, glm::vec3 sample, glm::vec3 &throughput, float sampleTime);
  glm::vec3 directIllumination(glm::vec3 pos, float &lightPdf, glm::vec3 &dir, float &lightArea, int except, float sampleTime);
  void initialMedium(glm::vec3 pos, float sampleTime);
  void updateStack(MyStack &st, HitRecord hit, bool penetrate);
 private:
  const RendererSettings *render_settings_{};
  const Scene *scene_{};
  std::mt19937 rd;
  std::uniform_real_distribution<float> uniform;
  const float continueProb = 0.85f;
  MyStack mediumStack;
};

}  // namespace sparks
