#include "sparks/assets/entity.h"
#include "sparks/assets/model.h"
#include "sparks/assets/scene.h"
namespace sparks {

const Model *Entity::GetModel() const {
  return model_.get();
}

glm::mat4 &Entity::GetTransformMatrix() {
  return transform_;
}

const glm::mat4 &Entity::GetTransformMatrix() const {
  return transform_;
}

Material &Entity::GetMaterial() {
  return material_;
}

const Material &Entity::GetMaterial() const {
  return material_;
}

const std::string &Entity::GetName() const {
  return name_;
}
glm::vec3 Entity::Sample_Li(HitRecord &hit_record,
                    std::mt19937 &rd,
                    glm::vec3 *wi,
                    float *pdf,float time) const{
  //printf("hello");
  glm::vec3 normal(0,0,1);
  glm::mat4 transform = transform_;
  if (time > 0 && duration != 0.0f)
    transform = transform*matrix_interpolation(anime_transform_, time/duration);
  glm::vec3 p_sample = GetModel()->Sample(rd,transform_,pdf,&normal);
  
   
  *wi = glm::normalize(p_sample - hit_record.position);
  *pdf=(*pdf)*glm::dot(p_sample - hit_record.position, p_sample - hit_record.position) /
           abs(glm::dot(normal, -*wi));
  if (dot(normal, -*wi) == 0.f)
    *pdf = 0;
  return dot(normal, -*wi) > 0.f
             ? GetMaterial().emission * GetMaterial().emission_strength
             : glm::vec3(0.f);


}
}  // namespace sparks
