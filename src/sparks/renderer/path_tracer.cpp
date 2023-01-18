#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"

namespace sparks {
PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene) {
  render_settings_ = render_settings;
  scene_ = scene;
}

glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample) const {
  std::mt19937 rd(sample ^ x ^ y);
  return SampleRay_test(origin, direction, x, y, sample,rd);
  glm::vec3 throughput{1.0f};
  glm::vec3 radiance{0.0f};
  HitRecord hit_record;
  const int max_bounce = render_settings_->num_bounces;
  for (int i = 0; i < max_bounce; i++) {
    auto t = scene_->TraceRay(origin, direction, 0.5, 1e4f, &hit_record);
    if (t > 0.0f) {
      auto &material =
          scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
      if (material.material_type == MATERIAL_TYPE_EMISSION) {
        radiance += throughput * material.emission * material.emission_strength;
        break;
      } else {
        throughput *=
            material.albedo_color *
            glm::vec3{scene_->GetTextures()[material.albedo_texture_id].Sample(
                hit_record.tex_coord)};
        origin = hit_record.position;
        direction = scene_->GetEnvmapLightDirection();
        radiance += throughput * scene_->GetEnvmapMinorColor();
        throughput *=
            std::max(glm::dot(direction, hit_record.normal), 0.0f) * 2.0f;
        if (scene_->TraceRay(origin, direction, 1e-3f, 1e4f, nullptr) < 0.0f) {
          radiance += throughput * scene_->GetEnvmapMajorColor();
        }
        break;
      }
    } else {
      radiance += throughput * glm::vec3{scene_->SampleEnvmap(direction)};
      break;
    }
  }
  return radiance;
}
glm::vec3 PathTracer::SampleRay_test(glm::vec3 origin,
                                     glm::vec3 direction,
                                     int x,
                                     int y,
                                     int sample,std::mt19937 rd) const {
  
  glm::vec3 throughput{1.0f};
  glm::vec3 radiance{0.0f};
  bool specularBounce = false;
  int bounces;
  float etaScale = 1;
  std::random_device raad;
  rd=std::mt19937((sample * 182741431) ^ (x * 1239239423) ^ (y * 129423741));
  //std::mt19937 rd(sample*sample*sample^x);
  float time = std::uniform_real_distribution<float>(0.0f, scene_->GetCamera().shutter_speed)(rd);
  for (bounces = 0;; ++bounces) {
    HitRecord hit_record;
    auto t = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit_record,time);
    bool intersected = (t > 0.0f);
    //return glm::vec3(scene_->GetEntity(hit_record.hit_entity_id).GetMaterial().IsEmission());
    if (bounces == 0 || specularBounce) {
      if (intersected) {
        auto &material =
            scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
        radiance += throughput * material.emission *
                        material.emission_strength;
      } else {
        glm::vec4 env_sample = scene_->SampleEnvmap(direction);
        env_sample[3] = 0;
        radiance += throughput * glm::vec3{env_sample};
      }
    }
    if (!intersected || bounces >= render_settings_->num_bounces)
      break;
    // sample a light from certain distribution
    auto &material = scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
    if (material.material_type != MATERIAL_TYPE_SPECULAR) {
      glm::vec3 light = scene_->SampleLight(direction,hit_record,rd,time);
      glm::vec3 Ld = throughput * light;
      radiance += Ld;
    } else {
      bounces -= 1;
      specularBounce = 1;
    }
    
    glm::vec3 wo = -direction, wi;
    float pdf;
    glm::vec3 f = material.Sample_f(hit_record,rd,wo,&wi, &pdf);
    glm::vec3 printvec = hit_record.tangent;
    if (f == glm::vec3(0.0f) || pdf == 0.f)
      break;
    throughput *= f * abs(dot(wi, hit_record.normal)) / pdf;
    origin = hit_record.position;
    direction = wi;
    if (bounces > 3) {
      float q = 0.5f;
      if (std::uniform_real_distribution<float>(0.0f, 1.0f)(rd) < q)
        break;
      throughput /= 1 - q;
    }

  }
  //return glm::vec3(0.0f);
  radiance = glm::clamp(radiance, glm::vec3(0.0f),
                  glm::vec3(scene_->GetCamera().GetClamp()));
  return radiance;
}
}  // namespace sparks
