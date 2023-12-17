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
glm::vec3 PathTracer::Shade(glm::vec3 origin, glm::vec3 direction) const {
  HitRecord hit_record;
  auto t = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit_record);
  if (t > 0.0f) {
    int num_ob = scene_->GetEntityCount();
    auto &material = scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
    if (material.material_type == MATERIAL_TYPE_EMISSION) {
      return material.emission * material.emission_strength;
    }
    glm::vec3 L_dir(0, 0, 0), L_ind(0, 0, 0);
    origin = hit_record.position;
    static std::default_random_engine e(time(NULL));
    for (int i = 0; i < num_ob; i++) {
      auto &object = scene_->GetEntity(i);
      auto &material =
          scene_->GetEntity(i).GetMaterial();
      if (material.material_type == MATERIAL_TYPE_EMISSION) {
        double total_area = 0;
        glm::vec3 xl, vn;
        const Mesh * mesh_ptr = reinterpret_cast<const Mesh *>(object.GetModel());
        auto &vertices = mesh_ptr->Get_vertices();
        auto &indices = mesh_ptr->Get_indices();
        int n_triangle_face = indices.size() / 3;
        double *triangle_area = new double[n_triangle_face];
        for (int j = 0; j < n_triangle_face; j++) {
          total_area += reinterpret_cast<const Mesh *>(object.GetModel())->area(j);
          triangle_area[j] = total_area;
        }
        static std::uniform_real_distribution<double> u1(0, total_area);
        double rnd = u1(e);
        for (int j = 0; j < n_triangle_face; j++) {
          if (rnd < triangle_area[j]) {
            // sample a vertex on triangle mesh j
            auto v1 = vertices[indices[3 * j]],
                 v2 = vertices[indices[3 * j + 1]],
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
        glm::vec3 new_direction = xl - hit_record.position;
        new_direction = glm::normalize(new_direction);
        float visible = 1;
        HitRecord intersection;
        auto vt = scene_->TraceRay(xl, vn, 1e-3f, 1e4f, &intersection);
        if (intersection.hit_entity_id != i)
          visible = 0;
        hit_record.normal = glm::normalize(hit_record.normal);
        auto dots = glm::dot(new_direction , hit_record.normal);
        if (dots > 0) {
          float pdf_light = double(1) / total_area;
          float cos_theta = abs(glm::dot(glm::normalize(new_direction) , glm::normalize(vn)));
          float cos_theta_hat = abs(glm::dot(glm::normalize(new_direction), hit_record.normal));
          float dist = glm::distance(xl, hit_record.position);
          float dist2 = dist * dist;
          glm::vec3 intensity = material.emission * material.emission_strength * cos_theta * cos_theta_hat /
                             dist2 / pdf_light * visible;
          L_dir += intensity * material.albedo_color;       
        }

      }
      
    }
    float p_rr = 0.5;
    if (RR(p_rr)) {
      glm::vec3 next_ray = glm::ballRand(1.0f);
      if (glm::dot(next_ray, hit_record.normal) < 0)
        next_ray = -next_ray;
      HitRecord intersection;
      auto t = scene_->TraceRay(hit_record.position, next_ray, 1e-3f, 1e4f, &intersection);
      if (t > 0) {
        glm::vec3 intensity = Shade(hit_record.position, next_ray) / p_rr;
        auto &material =
            scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
        L_ind += intensity * material.albedo_color;
      }
    }
    return L_dir + L_ind;

  } else {
    return glm::vec3{scene_->SampleEnvmap(direction)};
  }
 
  
}
glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample) const {
  const int num_sample_per_pixel = 1;
  glm::vec3 current_radiance(0, 0, 0);
  for (int k = 0; k < num_sample_per_pixel; k++) {
    glm::vec3 radiance = sparks::PathTracer::Shade(origin, direction);
    current_radiance.x += radiance.x / num_sample_per_pixel;
    current_radiance.y += radiance.y / num_sample_per_pixel;
    current_radiance.z += radiance.z / num_sample_per_pixel;
  }

  return current_radiance;
}
}  // namespace sparks
