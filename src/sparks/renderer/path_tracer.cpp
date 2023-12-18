#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"
#include <glm/glm.hpp>
#include <glm/gtc/random.hpp>

namespace sparks {
PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene) {
  render_settings_ = render_settings;
  scene_ = scene;
}
bool PathTracer::RR(double p) const{
  static std::default_random_engine e(time(NULL));
  static std::uniform_real_distribution<double> u1(0, 1);
  double rnd = u1(e);
  if (rnd < p)
    return true;
  else
    return false;
}
glm::vec3 PathTracer::Shade(HitRecord intersection, glm::vec3 wo, int depth) const {
  const float pi = 3.1415926;
  glm::vec3 L_dir(0, 0, 0), L_ind(0, 0, 0);
  auto &material = scene_->GetEntity(intersection.hit_entity_id).GetMaterial();
  if (material.material_type == MATERIAL_TYPE_EMISSION) {
    return material.emission * material.emission_strength;
  }
  static std::default_random_engine e(time(NULL));
  int num_ob = scene_->GetEntityCount();
  for (int i = 0; i < num_ob; i++) {
    auto &object = scene_->GetEntity(i);
    auto &material_i = scene_->GetEntity(i).GetMaterial();
    if (material_i.material_type == MATERIAL_TYPE_EMISSION) {
      double total_area = 0;
      glm::vec3 xl, vn;
      const Mesh *mesh_ptr = reinterpret_cast<const Mesh *>(object.GetModel());
      auto &vertices = mesh_ptr->Get_vertices();
      auto &indices = mesh_ptr->Get_indices();
      int n_triangle_face = indices.size() / 3;
      double *triangle_area = new double[n_triangle_face];
      for (int j = 0; j < n_triangle_face; j++) {
        total_area +=
            reinterpret_cast<const Mesh *>(object.GetModel())->area(j);
        triangle_area[j] = total_area;
      }
      static std::uniform_real_distribution<double> u1(0, total_area);
      double rnd = u1(e);
      for (int j = 0; j < n_triangle_face; j++) {
        if (rnd < triangle_area[j] && (j == 0 || rnd >= triangle_area[j - 1])) {
          auto v1 = vertices[indices[3 * j]], v2 = vertices[indices[3 * j + 1]],
               v3 = vertices[indices[3 * j + 2]];

          static std::uniform_real_distribution<double> u2(0, 1);
          double rnd1 = u2(e), rnd2 = u2(e), rnd3 = u2(e);
          float p1 = rnd1 / (rnd1 + rnd2 + rnd3),
                p2 = rnd2 / (rnd1 + rnd2 + rnd3),
                p3 = rnd3 / (rnd1 + rnd2 + rnd3);
          xl = v1.position * p1 + v2.position * p2 + v3.position * p3;
          vn = v1.normal * p1 + v2.normal * p2 + v3.normal * p3;
          break;
        }
      }
      delete[] triangle_area;
      glm::vec3 wi = glm::normalize(xl - intersection.position);
      float visible = 1;
      HitRecord hit_light;
      auto vt =
          scene_->TraceRay(intersection.position, wi, 1e-3f, 1e4f, &hit_light);
      if (hit_light.hit_entity_id != i)
        visible = 0;
      intersection.geometry_normal = glm::normalize(intersection.geometry_normal);
      auto dots = glm::dot(intersection.geometry_normal, wi);
      if (dots > 0) {
        float pdf_light = float(1) / total_area;
        float cos_theta = abs(glm::dot(wi, glm::normalize(vn)));
        float cos_theta_hat = abs(glm::dot(wi, intersection.geometry_normal));
        float dist = glm::distance(xl, intersection.position);
        float dist2 = dist * dist;
        glm::vec3 intensity = material.emission * material.emission_strength *
                              cos_theta * cos_theta_hat / dist2 / pdf_light *
                              visible;
        L_dir += intensity * material.albedo_color * dots / pi;
      }

    }
  }
  auto light_direction = scene_->GetEnvmapLightDirection();
  auto t_l = scene_->TraceRay(intersection.position, light_direction, 1e-3f,
                              1e4f, nullptr);
  if (t_l < 0.0f) {
    L_dir +=
        scene_->GetEnvmapMajorColor() *
        std::max(glm::dot(light_direction, intersection.geometry_normal),
                 0.0f) * 2.0f;
  }
  float p_rr = 0.6;
  if (RR(p_rr)) {
    glm::vec3 next_ray = glm::ballRand(1.0f);
    if (glm::dot(next_ray, intersection.geometry_normal) < 0)
      next_ray = -next_ray;
    HitRecord next_intersection;
    auto t_n = scene_->TraceRay(intersection.position, next_ray, 1e-3f, 1e4f,
                              &next_intersection);
    auto throughput = material.albedo_color * glm::vec3{scene_->GetTextures()[material.albedo_texture_id].Sample(intersection.tex_coord)};
    if (t_n > 0 && scene_->GetEntity(next_intersection.hit_entity_id)
                           .GetMaterial()
                           .material_type != MATERIAL_TYPE_EMISSION) {
      glm::vec3 intensity = Shade(next_intersection, -next_ray, depth+1) / p_rr;
      L_ind += throughput * intensity;
    } else {
      L_ind += throughput * glm::vec3{scene_->SampleEnvmap(next_ray)};
    }
  }
  return L_dir + L_ind;
 
  
}
glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample) const {
  const float num_sample_per_pixel = 16;
  glm::vec3 current_radiance(0, 0, 0);
  HitRecord intersection;
  auto t = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &intersection);
  if (t > 0) {
    for (int k = 0; k < num_sample_per_pixel; k++) {
      glm::vec3 radiance = sparks::PathTracer::Shade(intersection, -direction, 0);
      current_radiance += radiance / num_sample_per_pixel;
    }
  } else {
    current_radiance = glm::vec3{scene_->SampleEnvmap(direction)};
  }
  return current_radiance;
}
}  // namespace sparks
