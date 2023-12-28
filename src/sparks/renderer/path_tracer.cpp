#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <memory>

#include <cassert>

namespace sparks {
PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene) {
  render_settings_ = render_settings;
  scene_ = scene;
  rd = std::mt19937(std::random_device()());
  uniform = std::uniform_real_distribution<float>(0, 1);
}


// glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
//                                 glm::vec3 direction,
//                                 int x,
//                                 int y,
//                                 int sample,
//                                 bool moreBounces,
//                                 float weight) {
//   glm::vec3 throughput{1.0f};
//   glm::vec3 radiance{0.0f};
//   HitRecord hit_record;
//   const int max_bounce = render_settings_->num_bounces;
//   std::mt19937 rd(sample ^ x ^ y);
//   for (int i = 0; i < max_bounce; i++) {
//     auto t = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit_record);
//     if (t > 0.0f) {
//       auto &material =
//           scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
//       if (material.material_type == MATERIAL_TYPE_EMISSION) {
//         radiance += throughput * material.emission * material.emission_strength;
//         break;
//       } else {
//         throughput *=
//             material.albedo_color *
//             glm::vec3{scene_->GetTextures()[material.albedo_texture_id].Sample(
//                 hit_record.tex_coord)};
//         origin = hit_record.position;
//         direction = scene_->GetEnvmapLightDirection();
//         radiance += throughput * scene_->GetEnvmapMinorColor();
//         throughput *=
//             std::max(glm::dot(direction, hit_record.normal), 0.0f) * 2.0f;
//         if (scene_->TraceRay(origin, direction, 1e-3f, 1e4f, nullptr) < 0.0f) {
//           radiance += throughput * scene_->GetEnvmapMajorColor();
//         }
//         break;
//       }
//     } else {
//       radiance += throughput * glm::vec3{scene_->SampleEnvmap(direction)};
//       break;
//     }
//   }
//   return radiance;
// }


// sample from light source, returns sample result(position on light source), pdf, normal vector and light
void PathTracer::SampleFromLight(glm::vec3 &res, glm::vec3 &norm, float &area, int except) {
  const std::vector<Entity> &lightSources = scene_->GetEntities(); // all of them shall be meshes
  std::vector<float> prefixSum;
  std::vector<std::tuple<glm::vec3, glm::vec3, glm::vec3>> triangles; 
  std::vector<Material> materials;
  float totalArea = 0;
  // calculate total area
  for (int i = 0; i < lightSources.size(); i++) {
    const Entity &object = lightSources[i];
    if (object.GetMaterial().material_type != MATERIAL_TYPE_EMISSION || i == except)
      continue;
    const std::vector<Vertex> &vertices = object.GetModel()->GetVertices();
    const std::vector<uint32_t> &indices = object.GetModel()->GetIndices();
    glm::mat4 transform = object.GetTransformMatrix();
    // std::cerr << "yeah" << std::endl;
    for(int i = 0; i < indices.size(); i += 3) {
      int j = i + 1, k = i + 2;
      glm::vec4 v0(vertices[indices[i]].position, 1.0f);
      glm::vec4 v1(vertices[indices[j]].position, 1.0f);
      glm::vec4 v2(vertices[indices[k]].position, 1.0f);
      glm::vec3 u0 = glm::vec3(transform * v0);
      glm::vec3 u1 = glm::vec3(transform * v1);
      glm::vec3 u2 = glm::vec3(transform * v2);
      float nowArea = glm::length(glm::cross(u2 - u0, u1 - u0)) / 2;
      totalArea += nowArea;
      prefixSum.emplace_back(totalArea);
      triangles.emplace_back(std::make_tuple(u0, u1, u2));
      materials.emplace_back(object.GetMaterial());
    }
  }
  if(totalArea == 0) {
    area = 0;
    return;
  }
  // uniformly sample
  float sample = uniform(rd) * totalArea;
  int totalTriangles = triangles.size();
  int triangleFound = std::lower_bound(prefixSum.begin(), prefixSum.end(), sample) - prefixSum.begin();
  float t = uniform(rd);
  float x = 1 - sqrt(1 - t); 
  float y = uniform(rd) * (1 - x);
  glm::vec3 v[3];
  std::tie(v[0], v[1], v[2]) = triangles[triangleFound];
  res = x * v[0] + y * v[1] + (1 - x - y) * v[2];
  norm = glm::normalize(glm::cross(v[2] - v[0], v[1] - v[0]));
  area = totalArea;
}

float PathTracer::getPdfByLight(glm::vec3 pos, glm::vec3 sample, float area) {
  HitRecord tmp;
  float t = scene_->TraceRay(pos, sample, 1e-3f, 1e4f, &tmp);
  if(t > 0.0f && scene_->GetEntity(tmp.hit_entity_id).GetMaterial().material_type == MATERIAL_TYPE_EMISSION)
    return 1 / area * t * t / std::fabs(glm::dot(tmp.normal, sample));
  else  
    return 0.0f;
}

glm::vec3 PathTracer::directIllumination(glm::vec3 pos, float &pdf, glm::vec3 &dir, float &lightArea, int except, Medium *currentMedium) {
  glm::vec3 directSample, directNorm, directDir;
  glm::vec3 res{0.0f};
  lightArea = 0.0f;
  SampleFromLight(directSample, directNorm, lightArea, except);
  if(lightArea != 0) {
    directDir = glm::normalize(directSample - pos);
    HitRecord tmpHit;
    glm::vec3 currentPos = pos, throughput = glm::vec3{1.0f};
    Medium *tmpMedium = currentMedium; 
    while(true) { // repeatly accumulate all attenuation into the throughput
      float check = scene_->TraceRay(currentPos, directDir, 1e-3f, 1e4f, &tmpHit);
      if(check < 0) {
        pdf = 0.0f;
        return glm::vec3{0.0f};
      }
      Material isectMaterial = scene_->GetEntity(tmpHit.hit_entity_id).GetMaterial();
      if(tmpMedium != nullptr) {
          throughput *= tmpMedium->Tr(currentPos, directDir, check, rd, uniform); // direct light get attenuated while passing through scattering medium
          // if(pos[1] == 0.0f && glm::floor(pos[0]) == 200 && glm::floor(pos[2]) == 200)
          //   std::cerr << throughput[0] << std::endl;
      }
      if(glm::length(currentPos + check * directDir - directSample) < 0.01f) {
        pdf = 1 / lightArea * glm::dot(directSample - pos, directSample - pos) / std::fabs(glm::dot(directNorm, -directDir));
        res = isectMaterial.emission * isectMaterial.emission_strength / pdf * throughput;
        break;
      } else if(!isectMaterial.false_surface) { 
        // false surface denotes those boundaries of participating medium
        // they determine region of medium, do not effect light ray directly
        dir = directDir;
        pdf = 0.0f;
        return glm::vec3{0.0f};
      }
      currentPos = tmpHit.position + 3e-5f * directDir;
      if(tmpHit.front_face)
        tmpMedium = isectMaterial.medium;
      else
        tmpMedium = nullptr;
    }
  } else {
    pdf = 0.0f;
    directDir = glm::vec3{0.0f};
  }
  dir = directDir;
  return res;
}

static float PowerHeuristic(float pdf1, float pdf2) {
  if (pdf1 == 0.0f)
    return 0.0f;
  return pdf1 / (pdf1 + pdf2);
  return pdf1 * pdf1 / (pdf1 * pdf1 + pdf2 * pdf2);
}

// path trace main algorithm
glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample,
                                int bounces,
                                Medium *currentMedium, 
                                float currentRatio) {           
  const float pdfClamp = 0.05f;
  glm::vec3 emission{0.0f}, direct{0.0f}, env{0.0f}, incident{0.0f};
  HitRecord hit;
  glm::vec3 volSample, volWeight = glm::vec3{1.0f};
  float volPdf, intersection;
  bool volFlag = false;
  glm::vec3 directSample, directDir = glm::vec3(0, 0, 1);
  float directPdf = 0, lightArea=0;
  glm::vec3 directNorm;
  float actualContinueProb = bounces > 2 ? continueProb : 1;
  float russianRoulette = uniform(rd);
  const std::function<float(const glm::vec3&, const glm::vec3&, float, float, HitRecord*)> TraceRayMethod = 
      std::bind(&Scene::TraceRay, scene_, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
  intersection = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit);
  if(currentMedium != nullptr) {
    volWeight = currentMedium->Sample(origin, direction, volSample, volPdf, volFlag, intersection, rd, uniform);
  }
  // If ray is not in vaccum, do medium sampling
  if(volFlag) {
    // directIllumination evaluates direct light strength ingoing
    // as for normal reflective material, only bsdf and cos(theta_i) need to be further multiplied
    direct = directIllumination(volSample, directPdf, directDir, lightArea, -1, currentMedium); 
    float p = currentMedium->p(-direction, directDir);
    direct = glm::vec3(p) * direct * PowerHeuristic(directPdf, p) * volWeight;
    if(russianRoulette >= actualContinueProb || bounces >= render_settings_->num_bounces)
      return direct;
    glm::vec3 incidentSample;
    volWeight *= currentMedium->Sample_p(-direction, incidentSample, rd, uniform) / actualContinueProb;
    return direct + volWeight * SampleRay(volSample, incidentSample, x, y, sample, bounces + 1, currentMedium) + currentMedium->getEmission();
  } else { // otherwise, do normal sampling
    auto FirstBounceClamp = [bounces](const glm::vec3 &L) -> glm::vec3 {
      return bounces == 0 ? glm::clamp(L, 0.0f, 1.0f) : L;
    };
    if(intersection <= 0.0f)
      return FirstBounceClamp(glm::vec3(scene_->SampleEnvmap(direction)));
    Material material = scene_->GetEntity(hit.hit_entity_id).GetMaterial();
    std::unique_ptr<BSDF> bsdf(material.ComputeBSDF(hit, scene_)); 
    std::unique_ptr<BSSRDF> bssrdf(material.ComputeBSSRDF(hit, direction, scene_));
    //sample incident
    glm::vec3 incidentSample;
    float incidentPdf;
    BxDFType sampledType;
    glm::vec3 brdf;
    if (russianRoulette > actualContinueProb || bounces >= render_settings_->num_bounces) {
      incidentPdf = 0;
      sampledType = BxDFType(0);
      brdf = glm::vec3{0.0f};
      incidentSample = glm::vec3{0.0f};
      incident = glm::vec3{0.0f};
    } else {
      // Coming into an entiy with BSSRDF
      if(bssrdf != nullptr) {
        std::cerr << "???" << std::endl;
        glm::vec3 ssSample, ssDir, ssNormal;
        float rawPdf = 0.0f, fullPdf = 0.0f;
        glm::vec3 ss = bssrdf->Sample_S(TraceRayMethod, ssSample, ssDir, ssNormal, rawPdf, fullPdf,
                                        hit.hit_entity_id, rd, uniform);
        if(ss != glm::vec3{0.0f} && rawPdf != 0.0f && fullPdf != 0.0f) {
          // direct illumination on the point ray leaving
          glm::vec3 ssDirect{0.0f}, ssDirectDir{0.0f};
          float ssDirectPdf;
          ssDirect = directIllumination(ssSample, ssDirectPdf, ssDirectDir, lightArea, -1);
          float cosTheta = glm::dot(ssDirectDir, ssNormal);
          float surfacePdf = glm::sqrt(1 - cosTheta * cosTheta) * cosTheta * INV_PI;
          if(cosTheta == 0.0f)
            std::cerr << glm::length(ssDirectDir) << std::endl;
          if(cosTheta < 1.0f && cosTheta > 0.0f) {
            incident += ss * ssDirect * cosTheta * INV_PI / (surfacePdf * rawPdf); // we assume the brdf at the point is lambertian
          }
          // incident illumination on the point ray leaving
          incident += ss * SampleRay(ssSample + 3e-5f * ssDir, ssDir, x, y, sample, bounces + 1, nullptr) *
                      glm::dot(ssDir, ssNormal) * INV_PI / fullPdf;
          incident /= actualContinueProb;
          // reflection on surface of bssrdf material is not implemented yet
          // it should be identical to other material (in fact lambertian material, as assumed)
          // but unknown bug exists so I cannot reuse the existing code
          // TODO
          incident = glm::clamp(incident, 0.0f, 1.0f); 
          // it shouldn't appear since we don't have clamp in procedure. but removing it produces BLOODY scene in translucent example.
          // BSSRDF is really confusing for debugging and parameter setting. I will try to fix it (or at least understand it) later. 
        }
      } else {
        // incident illumination
        brdf = bsdf->Sample_f(-direction, incidentSample,
                              glm::vec2{uniform(rd), uniform(rd)}, incidentPdf,
                              BxDFType(BSDF_ALL), sampledType);
        assert(!std::isnan(incidentPdf));
        assert(!std::isinf(incidentPdf));
        glm::vec3 pos = hit.position + 3e-5f * incidentSample;
        glm::vec3 norm = material.GetShaderNormal(hit, scene_);
        if (incidentPdf > 0.0f) {
          float incidentPdfX = getPdfByLight(hit.position, incidentSample, lightArea);
          assert(!std::isinf(incidentPdf)); 
          float incidentRatio = incidentPdf == 0 ? 0 : PowerHeuristic(incidentPdf, incidentPdfX);
          bool penetration = glm::dot(hit.geometry_normal, incidentSample) < 0.0f;
          Medium *newMedium = penetration ? (hit.front_face ? material.medium : nullptr) : currentMedium; 
          incidentPdf = glm::clamp(incidentPdf, pdfClamp, 1e10f);
          incident = SampleRay(hit.position + 3e-5f * incidentSample,
                              incidentSample, x, y, hit.position[1] == 330.0f ? -1 : sample, bounces + 1, newMedium, incidentRatio) *
                      brdf *
                      std::fabs(glm::dot(norm, incidentSample)) /
                      (actualContinueProb * incidentPdf);
        }
        else
          incident = glm::vec3{0.0f};
        if (sampledType & BSDF_SPECULAR) {
          return FirstBounceClamp(incident);
        }
        // emission
        emission = material.material_type == MATERIAL_TYPE_EMISSION
                          ? material.emission * material.emission_strength * currentRatio
                          : glm::vec3{0.0f};
        // std::cerr << currentRatio << std::endl;


        // direct illumination (from environment map)
      
          // minor color used for diffuse reflection only
        env += bsdf->f(-direction, norm, BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)) * scene_->GetEnvmapMinorColor();
        glm::vec3 envDirection = scene_->GetEnvmapLightDirection();
        if (scene_->TraceRay(pos, envDirection, 1e-3f, 1e4f, nullptr) < 0.0f) {
        env += bsdf->f(-direction, envDirection,
                        BxDFType(BSDF_ALL)) *
                scene_->GetEnvmapMajorColor() *
                std::fabs(glm::dot(envDirection, norm)) * 2.0f;
        }

        // direct illumination (from emitting entities)
        direct = directIllumination(hit.position, directPdf, directDir, lightArea, hit.hit_entity_id);
        if(direct != glm::vec3{0.0f}) {
          direct *= bsdf->f(-direction, directDir, BxDFType(BSDF_ALL)) * std::fabs(glm::dot(norm, directDir));
          float directPdfX = bsdf->Pdf(-direction, directDir);
          float directRatio = directPdf == 0 ? 0 : PowerHeuristic(directPdf, directPdfX);
          direct *= directRatio;
          assert(!std::isnan(direct[0]));
        } else {
          directPdf = 0;
        }
        // incident illumination
        
        // float incidentPdfX = getPdfByLight(pos, incidentSample, lightArea);
        // assert(!std::isinf(incidentPdf));
        // float incidentRatio = incidentPdf == 0 ? 0 : PowerHeuristic(incidentPdf, incidentPdfX);
        // //if (std::isnan(incidentRatio)) {
        // //  std::cerr << incidentPdf << " " << incidentPdfX << std::endl;
        // //}
        // assert(!std::isnan(incidentRatio));
        // incident *= incidentRatio;
      }
    }
    if (sampledType & BSDF_SPECULAR) {
      //printf("incident %f %f %f\n", incident[0], incident[1], incident[2]);
      return FirstBounceClamp(incident);
    }
    glm::vec3 L = (emission + direct + env + incident) * volWeight;
      //printf("%f %f %f\n", incident[0], incident[1], incident[2]);
    assert(!std::isnan(L[0]) && !std::isnan(L[1]) && !std::isnan(L[2]));
    //if (bounces == 0)
    //  L = glm::clamp(L, 0.0f, 1.0f);
    return FirstBounceClamp(L);
  }

}
}  // namespace sparks
