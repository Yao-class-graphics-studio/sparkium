#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <memory>
#include <stack>

#include <cassert>
#include <glm/gtc/type_ptr.hpp>

namespace sparks {

glm::vec3 shift(glm::vec3 pos, glm::vec3 norm, glm::vec3 dir) {
  float cos0 = glm::dot(norm, dir);
  if(std::fabs(cos0) < 1e-6f) 
    return pos + 3e-4f * norm * (cos0 > 0 ? 1.0f : -1.0f);
  else
    return pos + 3e-4f / std::fabs(cos0) * dir;
}

PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene) {
  render_settings_ = render_settings;
  scene_ = scene;
  rd = std::mt19937(std::random_device()());
  uniform = std::uniform_real_distribution<float>(0, 1);
}

// Jingyi Lyu
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

// Jingyi Lyu
float PathTracer::getPdfByLight(glm::vec3 pos, glm::vec3 sample, float area, float sampleTime) {
  HitRecord tmp;
  float t = scene_->TraceRay(pos, sample, 1e-3f, 1e4f, &tmp, sampleTime);
  if(t > 0.0f && scene_->GetEntity(tmp.hit_entity_id).GetMaterial().material_type == MATERIAL_TYPE_EMISSION)
    return 1 / area * t * t / std::fabs(glm::dot(tmp.normal, sample));
  else  
    return 0.0f;
}

// Jingyi Lyu
int PathTracer::shadowRay(glm::vec3 pos, glm::vec3 &dir, glm::vec3 sample, glm::vec3 &throughput, float sampleTime) {
  MyStack tmpStack = mediumStack;
  HitRecord tmpHit;
  glm::vec3 currentPos = pos;
  throughput = glm::vec3{1.0f};
  glm::vec3 norm;
  while(true) {
    Medium *tmpMedium = tmpStack.top();
    float check = scene_->TraceRay(currentPos, dir, 1e-3f, 1e4f, &tmpHit, sampleTime);
    if(check < 0) {
      if(tmpMedium != nullptr) 
        throughput *= tmpMedium->Tr(currentPos, dir, 1e10f, rd, uniform);
      return 1;
    }
    const Material &material = scene_->GetEntity(tmpHit.hit_entity_id).GetMaterial();
    norm = material.GetShaderNormal(tmpHit, scene_);
    if(tmpMedium != nullptr)
      throughput *= tmpMedium->Tr(currentPos, dir, check, rd, uniform);
    if(glm::length(currentPos + check * dir - sample) < 0.01f) {
      throughput *= material.emission * material.emission_strength;
      if(std::fabs(glm::dot(dir, tmpHit.geometry_normal)) < material.cone)
        throughput *= 0.0f;
      return 0;
    } else if (!material.false_surface) 
      return 2;
    currentPos = shift(tmpHit.position, norm, dir);
    updateStack(tmpStack, tmpHit, true);
  }
}

// Jingyi Lyu
glm::vec3 PathTracer::directIllumination(glm::vec3 pos, float &pdf, glm::vec3 &dir, float &lightArea, int except, float sampleTime) {
  // return glm::vec3{1.0f};
  glm::vec3 directSample, directNorm, directDir;
  glm::vec3 res{0.0f};
  lightArea = 0.0f;
  SampleFromLight(directSample, directNorm, lightArea, except);
  if(lightArea != 0) {
    directDir = glm::normalize(directSample - pos);
    glm::vec3 throughput;
    int flag = shadowRay(pos, directDir, directSample, throughput, sampleTime);
    if(flag != 0) {
      pdf = 0.0f;
    } else {
      float cosCone = std::fabs(glm::dot(directNorm, -directDir));
      pdf = 1 / lightArea * glm::dot(directSample - pos, directSample - pos) / cosCone;
      res = throughput / pdf;
    }
  } else {
    pdf = 0.0f;
    directDir = glm::vec3{0.0f};
  }
  dir = directDir;
  return res;
}

// Jingyi Lyu
// To process multiple nested medium, use a stack
void PathTracer::updateStack(MyStack &st, HitRecord hit, bool penetrate) {
  const Material material = scene_->GetEntity(hit.hit_entity_id).GetMaterial();
  // std::cerr << hit.position[1] << std::endl;
  if(!material.false_surface || !penetrate)
    return;
  if(hit.front_face)
    st.push(material.medium);
  else
    st.pop();
}

// Jingyi Lyu
void PathTracer::initialMedium(glm::vec3 pos, float sampleTime) {
  glm::vec3 currentPos(0.0f, 1e3f, 0.0f);
  glm::vec3 dir = glm::normalize(pos - currentPos);
  mediumStack.clear();
  mediumStack.push(nullptr);
  int cnt = 0;
  while(true) {
    HitRecord hit;
    int t = scene_->TraceRay(currentPos, dir, 1e-3f, 1e4f, &hit, sampleTime);
    int dis = glm::length(pos - currentPos);
    if(t < 0 || t > dis)
      break;
    const Material &material = scene_->GetEntity(hit.hit_entity_id).GetMaterial();
    glm::vec3 norm = material.GetShaderNormal(hit, scene_);
    if(!material.false_surface && hit.front_face)
      cnt++;
    if(cnt == 0)
      updateStack(mediumStack, hit, true);
    currentPos = hit.position + 1e-3f / std::fabs(glm::dot(norm, dir)) * dir;
    dir = glm::normalize(pos - currentPos);
    if(!material.false_surface && !hit.front_face)
      cnt--;
  }
  mediumStack.flag = true;
}

// Jingyi Lyu + Shengquan Du
static float PowerHeuristic(float pdf1, float pdf2) {
  if (pdf1 == 0.0f)
    return 0.0f;
  return pdf1 / (pdf1 + pdf2);
  return pdf1 * pdf1 / (pdf1 * pdf1 + pdf2 * pdf2);
}

// Jingyi Lyu implemented the basic algorithm first
// Then slightly reconstructed by Shengquan Du
// Both have added new features after reconstruction
// path trace main algorithm
glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                float sampleTime,
                                int sample,
                                int bounces,
                                bool initialized,
                                float currentRatio) {    
  // Initialization: Jingyi Lyu + Shengquan Du       
  if(!initialized)
    initialMedium(origin, sampleTime);
  Medium *currentMedium = mediumStack.top();
  const float pdfClamp = 0.1f;
  glm::vec3 emission{0.0f}, direct{0.0f}, env{0.0f}, incident{0.0f};
  glm::vec3 envDirection = scene_->GetSceneLightDirection();
  glm::vec3 envLight = scene_->GetSceneLight(), envThroughput{1.0f};
  HitRecord hit;
  glm::vec3 volSample, volWeight = glm::vec3{1.0f};
  float volPdf, intersection;
  bool volFlag = false;
  glm::vec3 directSample, directDir = glm::vec3(0, 0, 1);
  float directPdf = 0, lightArea = 0;
  glm::vec3 directNorm;
  float actualContinueProb = bounces > 2 ? continueProb : 1;
  float russianRoulette = uniform(rd);
  const std::function<float(const glm::vec3&, const glm::vec3&, float, float, HitRecord*)> TraceRayMethod = 
      std::bind(&Scene::TraceRay, scene_, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, sampleTime);
  auto FirstBounceClamp = [bounces](const glm::vec3 &L) -> glm::vec3 {
    //return L;
    if (bounces == 0)
      return glm::clamp(L, 0.0f, 10.0f);
    return L;
  };
  intersection = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit, sampleTime);
  if(currentMedium != nullptr) {
    volWeight = currentMedium->Sample(origin, direction, volSample, volPdf, volFlag, intersection, rd, uniform);
  }
  // Medium sampling: Jingyi Lyu
  // If ray is not in vaccum, do medium sampling
  if(volFlag) {
    // directIllumination evaluates direct light strength ingoing
    // as for normal reflective material, only bsdf and cos(theta_i) need to be further multiplied
    direct = directIllumination(volSample, directPdf, directDir, lightArea, -1, sampleTime) * volWeight; 
    float p = currentMedium->p(-direction, directDir);
    direct = glm::vec3(p) * direct * PowerHeuristic(directPdf, p);
    if(scene_->use_scene_light) {
      int flag = shadowRay(volSample, envDirection, glm::vec3{0.0f}, envThroughput, sampleTime);
      // if(bounces == 0)
      //   std::cerr << envThroughput[0] * volWeight[0] << std::endl;
      if (flag == 1) {
        env += currentMedium->p(-direction, envDirection) * scene_->GetSceneLight() * envThroughput * volWeight;
      }
    }
    if(russianRoulette >= actualContinueProb || bounces >= render_settings_->num_bounces)
      return (direct + env) * volWeight;
    glm::vec3 incidentSample;
    volWeight *= currentMedium->Sample_p(-direction, incidentSample, rd, uniform) / actualContinueProb;
    return direct + env + SampleRay(volSample, incidentSample, x, y, sampleTime, sample, bounces + 1, true) + currentMedium->getEmission() * volWeight;
  } else { // otherwise, do normal sampling
    if (intersection <= 0.0f)
      return scene_->SampleEnvmap(direction);
    Material material = scene_->GetEntity(hit.hit_entity_id).GetMaterial();
    std::unique_ptr<BSDF> bsdf(material.ComputeBSDF(hit, scene_)); 
    std::unique_ptr<BSSRDF> bssrdf(material.ComputeBSSRDF(hit, direction, scene_));
    // sample incident
    glm::vec3 incidentSample;
    float incidentPdf;
    BxDFType sampledType = BxDFType(0);
    glm::vec3 brdf;
    if (russianRoulette > actualContinueProb || bounces >= render_settings_->num_bounces) {
      incidentPdf = 0;
      sampledType = BxDFType(0);
      brdf = glm::vec3{0.0f};
      incidentSample = glm::vec3{0.0f};
      incident = glm::vec3{0.0f};
    } else {
      // Coming into an entiy with BSSRDF. Jingyi Lyu
      if(bssrdf != nullptr) {
        glm::vec3 ssSample, ssDir, ssNormal;
        float rawPdf = 0.0f, fullPdf = 0.0f;
        glm::vec3 ss = bssrdf->Sample_S(TraceRayMethod, ssSample, ssDir, ssNormal, rawPdf, fullPdf,
                                        hit.hit_entity_id, rd, uniform);
        if(glm::length(ss) > 1e-6f && rawPdf > 1e-6f && fullPdf > 1e-6f) {
          assert(!std::isnan(ss[0]));
          assert(!std::isnan(ss[1]));
          assert(!std::isnan(ss[2]));
          // direct illumination on the point ray leaving
          glm::vec3 ssDirect{0.0f}, ssDirectDir{0.0f};
          float ssDirectPdf;
          ssDirect = directIllumination(ssSample, ssDirectPdf, ssDirectDir, lightArea, -1, sampleTime);
          float cosTheta = glm::dot(ssDirectDir, ssNormal);
          float surfacePdf = glm::sqrt(1 - cosTheta * cosTheta) * cosTheta * INV_PI;
          if(cosTheta < 1.0f && cosTheta > 0.0f) {
            incident += ss * ssDirect * cosTheta * INV_PI / (surfacePdf * rawPdf); 

            assert(!std::isnan(incident[0]));
            // we assume the brdf at the point is lambertian
          }
          // incident illumination on the point ray leaving
          incident += ss * SampleRay(ssSample + 1e-3f * ssDir, ssDir, x, y, sampleTime, sample, bounces + 1, true) *
                      glm::dot(ssDir, ssNormal) * INV_PI / fullPdf;

          assert(!std::isnan(incident[0]));
          incident /= actualContinueProb;
          incident = glm::clamp(incident, 0.0f, 1.0f); 
        }
      } else { // Normal surface. Jingyi Lyu + Shengquan Du
        brdf = bsdf->Sample_f(-direction, incidentSample,
                              glm::vec2{uniform(rd), uniform(rd)}, incidentPdf,
                              BxDFType(BSDF_ALL), sampledType);
        
        glm::vec3 norm = material.GetShaderNormal(hit, scene_);
        glm::vec3 pos = hit.position + 3e-5f * norm;
        // direct illumination (from emitting entities). Jingyi Lyu + Shengquan Du
        direct = directIllumination(hit.position + 3e-5f * norm, directPdf, directDir, lightArea, hit.hit_entity_id, sampleTime);
        if(direct != glm::vec3{0.0f}) {
          direct *= bsdf->f(-direction, directDir, BxDFType(BSDF_ALL)) * std::fabs(glm::dot(norm, directDir));
          float directPdfX = bsdf->Pdf(-direction, directDir);
          float directRatio = directPdf == 0 ? 0 : PowerHeuristic(directPdf, directPdfX);
          direct *= directRatio;
          assert(!std::isnan(direct[0]));
        } else {
          directPdf = 0;
        }

        // incident illumination. Jingyi Lyu + Shengquan Du
        assert(!std::isnan(incidentPdf));
        assert(!std::isinf(incidentPdf));
        bool penetrate = glm::dot(norm, incidentSample) * glm::dot(norm, -direction) < 0.0f;
        if (incidentPdf > 0.0f) {
          float incidentPdfX = getPdfByLight(hit.position, incidentSample, lightArea, sampleTime);
          assert(!std::isinf(incidentPdf)); 
          float incidentRatio = incidentPdf == 0 ? 0 : PowerHeuristic(incidentPdf, incidentPdfX);
          updateStack(mediumStack, hit, penetrate);
          incidentPdf = glm::clamp(incidentPdf, pdfClamp, 1e10f);
          incident = SampleRay(shift(hit.position, norm, incidentSample),
                              incidentSample, x, y, sampleTime, sample, bounces + (material.false_surface ? 0 : 1), true, incidentRatio);
          incident = incident * brdf *
                      std::fabs(glm::dot(norm, incidentSample)) /
                      (actualContinueProb * incidentPdf);
          
        }
        else
          incident = glm::vec3{0.0f};
        if ((sampledType & BSDF_SPECULAR) || material.false_surface) {
          if (sampledType & BSDF_SPECULAR)
          return FirstBounceClamp(incident);
        }
        // emission. Jingyi Lyu
        emission = material.material_type == MATERIAL_TYPE_EMISSION
                          ? material.emission * material.emission_strength * currentRatio
                          : glm::vec3{0.0f};
        // direct illumination (from environment). Jingyi Lyu
        if(scene_->use_scene_light) {
          int flag = shadowRay(pos, envDirection, glm::vec3{0.0f}, envThroughput, sampleTime);
          if (flag == 1) {
            env += bsdf->f(-direction, envDirection, BxDFType(BSDF_ALL)) * scene_->GetSceneLight() *
                   std::fabs(glm::dot(envDirection, norm)) * envThroughput;
          }
        }
        
      }
    }
    if (sampledType & BSDF_SPECULAR) {
      return FirstBounceClamp(incident);
    }
    glm::vec3 L = (emission + direct + env + incident) * volWeight;
    assert(!std::isnan(L[0]) && !std::isnan(L[1]) && !std::isnan(L[2]));
    if(material.material_type == MATERIAL_TYPE_PRINCIPLED)
      L = clamp(L, 0.0f, 1.5f);
    return FirstBounceClamp(L);
  }

}
}  // namespace sparks
