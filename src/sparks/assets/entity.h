#pragma once
#include "memory"
#include "sparks/assets/material.h"
#include "sparks/assets/mesh.h"
#include "sparks/assets/model.h"
#include <glm/gtx/matrix_decompose.hpp>
#include<glm/gtc/quaternion.hpp>
#include "glm/gtc/matrix_transform.hpp"

namespace sparks {
class Entity {
 public:
  template <class ModelType>
  Entity(const ModelType &model,
         const Material &material,
         const glm::mat4 &transform = glm::mat4{1.0f}) {
    model_ = std::make_unique<ModelType>(model);
    material_ = material;
    transform_ = transform;
    name_ = model_->GetDefaultEntityName();
  }

  template <class ModelType>
  Entity(const ModelType &model,
         const Material &material,
         const glm::mat4 &transform,
         const std::string &name) {
    model_ = std::make_unique<ModelType>(model);
    material_ = material;
    transform_ = transform;
    name_ = name;
  }
  [[nodiscard]] const Model *GetModel() const;
  [[nodiscard]] glm::mat4 &GetTransformMatrix();
  [[nodiscard]] const glm::mat4 &GetTransformMatrix() const;
  [[nodiscard]] Material &GetMaterial();
  [[nodiscard]] const Material &GetMaterial() const;
  [[nodiscard]] const std::string &GetName() const;
  [[nodiscard]] glm::vec3 Sample_Li(HitRecord &hit_record,
                      std::mt19937 &rd,
                      glm::vec3 *wi,
                      float *pdf,float time=-1.0f)const;
//  glm::mat4 trans = glm::translate(glm::mat4{1.0f}, glm::vec3(0,0,0));
//  glm::mat4 rot = glm::mat4_cast(glm::rotate(glm::quat(1,0,0,0),3.14f,glm::vec3(1,0,0)));
//  glm::mat4 scal = glm::scale(glm::mat4{1.0f}, glm::vec3(1,1,1));
  glm::mat4 anime_transform_{1.0f};
float duration{0.0f};

 private:
  std::unique_ptr<Model> model_;
  Material material_{};
  glm::mat4 transform_{1.0f};
  std::string name_;
  
};
}  // namespace sparks
