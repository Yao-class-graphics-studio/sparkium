﻿#include "sparks/assets/scene.h"

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "imgui.h"
#include "sparks/assets/accelerated_mesh.h"
#include "sparks/util/util.h"

/*
#include "assimp/Importer.hpp"     // C++ importer interface
#include "assimp/scene.h"           // Output data structure
#include "assimp/postprocess.h"
*/
#include "glm/gtc/constants.hpp"
namespace sparks {
Scene::Scene() {
  AddTexture(Texture(1, 1, glm::vec4{1.0f}, SAMPLE_TYPE_LINEAR), "Pure White");
  envmap_id_ = AddTexture(Texture(1, 1, glm::vec4{0.0f}, SAMPLE_TYPE_LINEAR),
                          "Pure Black");
  Texture envmap;
  Texture::Load(u8"../../textures/envmap_clouds_4k.hdr", envmap);
  envmap.SetSampleType(SAMPLE_TYPE_LINEAR);
  AddTexture(envmap, "Clouds");
}

int Scene::AddTexture(const Texture &texture, const std::string &name) {
  textures_.push_back(texture);
  texture_names_.push_back(name);
  return int(textures_.size() - 1);
}

const std::vector<Texture> &Scene::GetTextures() const {
  return textures_;
}

int Scene::GetTextureCount() const {
  return int(textures_.size());
}

void Scene::Clear() {
  textures_.clear();
  entities_.clear();
  camera_ = Camera{};
}

Entity &Scene::GetEntity(int entity_index) {
  return entities_[entity_index];
}

const Entity &Scene::GetEntity(int entity_index) const {
  return entities_[entity_index];
}

std::vector<Entity> &Scene::GetEntities() {
  return entities_;
}

const std::vector<Entity> &Scene::GetEntities() const {
  return entities_;
}

int Scene::GetEntityCount() const {
  return int(entities_.size());
}

Camera &Scene::GetCamera() {
  return camera_;
}

const Camera &Scene::GetCamera() const {
  return camera_;
}

void Scene::SetCamera(const Camera &camera) {
  camera_ = camera;
}

glm::mat4 Scene::GetCameraToWorld() const {
  return glm::translate(glm::mat4{1.0f}, camera_position_) *
         ComposeRotation(camera_pitch_yaw_roll_) *
         glm::scale(glm::mat4{1.0f}, glm::vec3{1.0f, 1.0f, 1.0f});
}

void Scene::SetCameraToWorld(const glm::mat4 &camera_to_world) {
  camera_pitch_yaw_roll_ = DecomposeRotation(
      camera_to_world *
      glm::scale(glm::mat4{1.0f}, glm::vec3{1.0f, 1.0f, 1.0f}));
  camera_position_ = camera_to_world[3];
}

int &Scene::GetEnvmapId() {
  return envmap_id_;
}

const int &Scene::GetEnvmapId() const {
  return envmap_id_;
}

float &Scene::GetEnvmapOffset() {
  return envmap_offset_;
}

const float &Scene::GetEnvmapOffset() const {
  return envmap_offset_;
}

glm::vec3 &Scene::GetCameraPosition() {
  return camera_position_;
}
const glm::vec3 &Scene::GetCameraPosition() const {
  return camera_position_;
}
float &Scene::GetCameraSpeed() {
  return camera_speed_;
}
const float &Scene::GetCameraSpeed() const {
  return camera_speed_;
}
glm::vec3 &Scene::GetCameraPitchYawRoll() {
  return camera_pitch_yaw_roll_;
}
const glm::vec3 &Scene::GetCameraPitchYawRoll() const {
  return camera_pitch_yaw_roll_;
}

void Scene::UpdateEnvmapConfiguration() {
  const auto &scene = *this;
  auto envmap_id = scene.GetEnvmapId();
  auto &envmap_texture = scene.GetTextures()[envmap_id];
  auto buffer = envmap_texture.GetBuffer();
  envmap_tot_weight_ = 0;
  envmap_minor_color_ = glm::vec3{0.0f};
  envmap_major_color_ = glm::vec3{0.0f};
  envmap_cdf_.resize(envmap_texture.GetWidth() * envmap_texture.GetHeight());

  std::vector<float> sample_scale_(envmap_texture.GetHeight() + 1);
  auto inv_width = 1.0f / float(envmap_texture.GetWidth());
  auto inv_height = 1.0f / float(envmap_texture.GetHeight());
  for (int i = 0; i <= envmap_texture.GetHeight(); i++) {
    float x = float(i) * glm::pi<float>() * inv_height;
    sample_scale_[i] = (-std::cos(x) + 1.0f) * 0.5f;
  }

  auto width_height = envmap_texture.GetWidth() * envmap_texture.GetHeight();
  float total_weight = 0.0f;
  float major_strength = -1.0f;
  for (int y = 0; y < envmap_texture.GetHeight(); y++) {
    auto scale = sample_scale_[y + 1] - sample_scale_[y];

    auto theta = (float(y) + 0.5f) * inv_height * glm::pi<float>();
    auto sin_theta = std::sin(theta);
    auto cos_theta = std::cos(theta);

    for (int x = 0; x < envmap_texture.GetWidth(); x++) {
      auto phi = (float(x) + 0.5f) * inv_width * glm::pi<float>() * 2.0f;
      auto sin_phi = std::sin(phi);
      auto cos_phi = std::cos(phi);

      auto i = y * envmap_texture.GetWidth() + x;
      auto color = glm::vec3{buffer[i]};
      auto minor_color = glm::clamp(color, 0.0f, 1.0f);
      auto major_color = color - minor_color;
      envmap_major_color_ += major_color * (scale * inv_width);
      envmap_minor_color_ += minor_color * (scale * inv_width);
      envmap_tot_weight_ += length(color) * scale * 2 * glm::pi<float>()*inv_width*radius;
      color *= scale;

      auto strength = std::max(color.x, std::max(color.y, color.z));
      if (strength > major_strength) {
        envmap_light_direction_ = {sin_theta * sin_phi, cos_theta,
                                   -sin_theta * cos_phi};
        major_strength = strength;
      }

      total_weight += strength * scale;
      envmap_cdf_[i] = total_weight;
    }
  }
  if (light_distribution.size() != 0) {
    light_distribution[0] = envmap_tot_weight_;
    float tot = 0;
    for (auto dist : light_distribution)
      tot += dist;
    if (tot != 0)
      for (int i = 0; i < light_distribution.size(); ++i)
        light_distribution[i] /= tot;

  }
    
  auto inv_total_weight = 1.0f / total_weight;
  for (auto &v : envmap_cdf_) {
    v *= inv_total_weight;
  }
}
glm::vec3 Scene::GetEnvmapLightDirection() const {
  float sin_offset = std::sin(envmap_offset_);
  float cos_offset = std::cos(envmap_offset_);
  return {cos_offset * envmap_light_direction_.x +
              sin_offset * envmap_light_direction_.z,
          envmap_light_direction_.y,
          -sin_offset * envmap_light_direction_.x +
              cos_offset * envmap_light_direction_.z};
}
const glm::vec3 &Scene::GetEnvmapMinorColor() const {
  return envmap_minor_color_;
}
const glm::vec3 &Scene::GetEnvmapMajorColor() const {
  return envmap_major_color_;
}
const std::vector<float> &Scene::GetEnvmapCdf() const {
  return envmap_cdf_;
}

float Scene::TraceRay(const glm::vec3 &origin,
                      const glm::vec3 &direction,
                      float t_min,
                      float t_max,
                      HitRecord *hit_record,float time) const {// hit in world corrdinate
  float result = -1.0f;
  HitRecord local_hit_record;
  float local_result;
  for (int entity_id = 0; entity_id < entities_.size(); entity_id++) {
    auto &entity = entities_[entity_id];
    auto transform = entity.GetTransformMatrix();
    if (time > 0 && entity.duration!=0.f)
      transform =transform*matrix_interpolation(entity.anime_transform_, time / entity.duration) ;
    auto inv_transform = glm::inverse(transform);
    auto transformed_direction =
        glm::vec3{inv_transform * glm::vec4{direction, 0.0f}};
    auto transformed_direction_length = glm::length(transformed_direction);
    if (transformed_direction_length < 1e-6) {
      continue;
    }
    local_result = entity.GetModel()->TraceRay(
        inv_transform * glm::vec4{origin, 1.0f},
        transformed_direction / transformed_direction_length, t_min,
        hit_record ? &local_hit_record : nullptr);
    local_result /= transformed_direction_length;
    if (local_result > t_min && local_result < t_max &&
        (result < 0.0f || local_result < result)) {
      result = local_result;
      if (hit_record) {
        local_hit_record.position =
            transform * glm::vec4{local_hit_record.position, 1.0f};
        local_hit_record.normal = glm::transpose(inv_transform) *
                                  glm::vec4{local_hit_record.normal, 0.0f};
        local_hit_record.tangent =
            transform * glm::vec4{local_hit_record.tangent, 0.0f};
        local_hit_record.geometry_normal =
            glm::transpose(inv_transform) *
            glm::vec4{local_hit_record.geometry_normal, 0.0f};
        *hit_record = local_hit_record;
        hit_record->hit_entity_id = entity_id;
      }
    }
  }
  if (hit_record) {
    hit_record->geometry_normal = glm::normalize(hit_record->geometry_normal);
    hit_record->normal = glm::normalize(hit_record->normal);
    hit_record->tangent = glm::normalize(hit_record->tangent);
  }
  return result;
}

glm::vec4 Scene::SampleEnvmap(const glm::vec3 &direction) const {
  float x = envmap_offset_;
  float y = acos(direction.y) * INV_PI;
  if (glm::length(glm::vec2{direction.x, direction.y}) > 1e-4) {
    x += glm::atan(direction.x, -direction.z);
  }
  x *= INV_PI * 0.5;
  return textures_[envmap_id_].Sample(glm::vec2{x, y});
}

const Texture &Scene::GetTexture(int texture_id) const {
  return textures_[texture_id];
}

std::vector<const char *> Scene::GetTextureNameList() const {
  std::vector<const char *> result;
  result.reserve(texture_names_.size());
  for (const auto &texture_name : texture_names_) {
    result.push_back(texture_name.data());
  }
  return result;
}

std::vector<const char *> Scene::GetEntityNameList() const {
  std::vector<const char *> result;
  result.reserve(entities_.size());
  for (const auto &entity : entities_) {
    result.push_back(entity.GetName().data());
  }
  return result;
}

bool Scene::TextureCombo(const char *label, int *current_item) const {
  return ImGui::Combo(label, current_item, GetTextureNameList().data(),
                      textures_.size());
}

bool Scene::EntityCombo(const char *label, int *current_item) const {
  return ImGui::Combo(label, current_item, GetEntityNameList().data(),
                      entities_.size());
}

int Scene::LoadTexture(const std::string &file_path) {
  Texture texture;
  if (Texture::Load(file_path, texture)) {
    return AddTexture(texture, PathToFilename(file_path));
  } else {
    LAND_WARN("[Sparks] Load Texture \"{}\" failed.", file_path);
    return 0;
  }
}

int Scene::LoadObjMesh(const std::string &file_path) {
  AcceleratedMesh mesh;
  if (Mesh::LoadObjFile(file_path, mesh)) {
    mesh.BuildAccelerationStructure();
    return AddEntity(mesh, Material{}, glm::mat4{1.0f},
                     PathToFilename(file_path));
  } else {
    return -1;
  }
}

Scene::Scene(const std::string &filename) : Scene() {
  auto doc = std::make_unique<tinyxml2::XMLDocument>();
  doc->LoadFile(filename.c_str());
  tinyxml2::XMLElement *rootElement = doc->RootElement();

  glm::mat4 camera_to_world = glm::inverse(
      glm::lookAt(glm::vec3{2.0f, 1.0f, 3.0f}, glm::vec3{0.0f, 0.0f, 0.0f},
                  glm::vec3{0.0f, 1.0f, 0.0f}));

  for (tinyxml2::XMLElement *child_element = rootElement->FirstChildElement();
       child_element; child_element = child_element->NextSiblingElement()) {
    std::string element_type{child_element->Value()};
    if (element_type == "envmap") {
      std::string envmap_type = child_element->FindAttribute("type")->Value();
      if (envmap_type == "file") {
        std::string envmap_filename =
            child_element->FindAttribute("value")->Value();
        Texture envmap;
        Texture::Load(envmap_filename, envmap);
        envmap_id_ = AddTexture(envmap, PathToFilename(envmap_filename));
      } else if (envmap_type == "color") {
        glm::vec3 color =
            StringToVec3(child_element->FindAttribute("value")->Value());
        Texture envmap(1, 1, glm::vec4{color, 1.0f});
        envmap_id_ = AddTexture(envmap, "Environment Map");
      }
    } else if (element_type == "camera") {
      camera_to_world =
          XmlTransformMatrix(child_element->FirstChildElement("transform"));
      float fov = 60.0f;
      float aperture = 0.0f;
      float focal_length = 3.0f;
      auto grandchild_element = child_element->FirstChildElement("fov");
      if (grandchild_element) {
        fov = std::stof(grandchild_element->FindAttribute("value")->Value());
      }
      grandchild_element = child_element->FirstChildElement("speed");
      if (grandchild_element) {
        camera_speed_ =
            std::stof(grandchild_element->FindAttribute("value")->Value());
      }
      grandchild_element = child_element->FirstChildElement("aperture");
      if (grandchild_element) {
        aperture =
            std::stof(grandchild_element->FindAttribute("value")->Value());
      }
      grandchild_element = child_element->FirstChildElement("focal_length");
      if (grandchild_element) {
        focal_length =
            std::stof(grandchild_element->FindAttribute("value")->Value());
      }
      camera_ = Camera(fov, aperture, focal_length);
    } else if (element_type == "model") {
      Mesh mesh = Mesh(child_element);
      Material material{};

      auto grandchild_element = child_element->FirstChildElement("material");
      if (grandchild_element) {
        material = Material(this, grandchild_element);
      }

      glm::mat4 transformation = XmlComposeTransformMatrix(child_element);

      auto name_attribute = child_element->FindAttribute("name");
      if (name_attribute) {
        AddEntity(AcceleratedMesh(mesh), material, transformation,
                  std::string(name_attribute->Value()));
      } else {
        AddEntity(AcceleratedMesh(mesh), material, transformation);
      }
    } else {
      LAND_ERROR("Unknown Element Type: {}", child_element->Value());
    }
  }

  SetCameraToWorld(camera_to_world);
  UpdateEnvmapConfiguration();

  float total_strength = 0;
  light_id.push_back(-1);
  light_distribution.push_back(envmap_tot_weight_);
  total_strength += envmap_tot_weight_;
  for (int i = 0; i < entities_.size(); ++i) {
    if (entities_[i].GetMaterial().IsEmission()) {
      light_id.push_back(i);
      float light_size = entities_[i].GetMaterial().emission_strength *
                         glm::pi<float>() *
                         entities_[i]
                             .GetModel()
                             ->GetAABB(entities_[i].GetTransformMatrix())
                             .SurfaceArea();
      light_distribution.push_back(light_size);
      total_strength += light_size;
    }
    if (entities_[i].GetName() == std::string("Lucy")) {
      printf("Success!");
      glm::mat4 trans = glm::translate(glm::mat4{1.0f}, glm::vec3(0,-0.2, 0));
      glm::mat4 rot =
      glm::mat4_cast(glm::rotate(glm::quat(1,0,0,0),3.14f*0.8f,glm::vec3(1,1,-0.5)));
      glm::mat4 scal = glm::scale(glm::mat4{1.0f}, glm::vec3(1,1,1));
      entities_[i].anime_transform_ = trans * rot * scal;
      entities_[i].duration = 0.3f;
    }
  }

  if (total_strength != 0) {
    for (int i = 0; i < light_distribution.size(); ++i) {
      light_distribution[i] /= total_strength;
      //sprintf("" light_distribution[i]);
    }
  }
}
/*
Entity Scene::processMesh(aiMesh *mesh, const aiScene *scene,glm::mat4 transform) {
  std::vector<Vertex> vertices;
  std::vector<unsigned int> indices;
  for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
    Vertex vertex;

    glm::vec3 vector;
    vector.x = mesh->mVertices[i].x;
    vector.y = mesh->mVertices[i].y;
    vector.z = mesh->mVertices[i].z;
    vertex.position = vector;
    vector.x = mesh->mNormals[i].x;
    vector.y = mesh->mNormals[i].y;
    vector.z = mesh->mNormals[i].z;
    vertex.normal = vector;
    if (mesh->mTextureCoords[0])
    {
      glm::vec2 vec;
      //only use tex coords 0
      vec.x = mesh->mTextureCoords[0][i].x;
      vec.y = mesh->mTextureCoords[0][i].y;
      vertex.tex_coord = vec;
    } else
      vertex.tex_coord = glm::vec2(0.0f, 0.0f);
    vector.x = mesh->mTangents[i].x;
    vector.y = mesh->mTangents[i].y;
    vector.z = mesh->mTangents[i].z;
    vertex.tangent = vector;
    vertices.push_back(vertex);
  }
  for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
    aiFace face = mesh->mFaces[i];
    for (unsigned int j = 0; j < face.mNumIndices; j++)
      indices.push_back(face.mIndices[j]);
  }
  aiMaterial *material = scene->mMaterials[mesh->mMaterialIndex];
  Material new_material;
  new_material.
  std::vector<Texture> diffuseMaps =
      loadMaterialTextures(material, aiTextureType_DIFFUSE, "texture_diffuse");
  textures_.insert(textures.end(), diffuseMaps.begin(), diffuseMaps.end());
  std::vector<Texture> specularMaps = loadMaterialTextures(
      material, aiTextureType_SPECULAR, "texture_specular");
  textures_.insert(textures.end(), specularMaps.begin(), specularMaps.end());
  std::vector<Texture> normalMaps =
      loadMaterialTextures(material, aiTextureType_HEIGHT, "texture_normal");
  textures_.insert(textures.end(), normalMaps.begin(), normalMaps.end());
  std::vector<Texture> heightMaps =
      loadMaterialTextures(material, aiTextureType_AMBIENT, "texture_height");
  textures_.insert(textures.end(), heightMaps.begin(), heightMaps.end());

  // 返回从提取的网格数据创建的网格对象
  return Entity(Mesh(vertices, indices),new_material,transform);
}
std::vector<int> Scene::loadMaterialTextures(aiMaterial *mat,
                                     aiTextureType type,
                                     std::string typeName) {
  std::vector<int> tex_id;
  for (unsigned int i = 0; i < mat->GetTextureCount(type); i++) {
    aiString str;
    mat->GetTexture(type, i, &str);
    // 检查之前是否加载了纹理，如果是，则继续下一次迭代：跳过加载新纹理
    bool skip = false;
    for (unsigned int j = 0; j < textures_.size(); j++) {
      if (std::strcmp(texture_names_[i].c_str(), str.C_Str()) == 0) {
        tex_id.push_back(j);
        skip = true;
        break;  // 已加载具有相同文件路径的纹理，继续下一个（优化）。
      }
    }
    if (!skip) {  // 如果尚未加载纹理，请加载它
      Texture texture;
      Texture::Load(texture)
      texture.id = TextureFromFile(str.C_Str(), this->directory);
      texture.type = typeName;
      texture.path = str.C_Str();
      textures.push_back(texture);
      textures_loaded.push_back(
          texture);  //将其存储为整个模型加载的纹理，以确保我们不会加载重复纹理。
      AddTexture
    }
  }
  return textures;
}
unsigned int TextureFromFile(const char *path,
                             const string &directory,
                             bool gamma) {
  string filename = string(path);
  filename = directory + '/' + filename;

  unsigned int textureID;
  glGenTextures(1, &textureID);

  int width, height, nrComponents;
  unsigned char *data =
      stbi_load(filename.c_str(), &width, &height, &nrComponents, 0);
  if (data) {
    GLenum format;
    if (nrComponents == 1)
      format = GL_RED;
    else if (nrComponents == 3)
      format = GL_RGB;
    else if (nrComponents == 4)
      format = GL_RGBA;

    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format,
                 GL_UNSIGNED_BYTE, data);
    glGenerateMipmap(GL_TEXTURE_2D);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                    GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    stbi_image_free(data);
  } else {
    std::cout << "纹理无法从此路径加载: " << path << std::endl;
    stbi_image_free(data);
  }
  return textureID;
}
void Scene::FindMesh(aiScene *scene, aiNode *node,glm::mat4 transform) {
  std::cerr << node->mName.C_Str() << std::endl;
  for (unsigned int i = 0; i < node->mNumMeshes; i++) {
    aiMesh *mesh = scene->mMeshes[node->mMeshes[i]];
    entities_.emplace_back(processMesh(mesh, scene,transform));
  }
  for (int i = 0; i < node->mNumChildren; ++i) {
    FindMesh(scene, node->mChildren[i]);
  }
}
bool Scene::LoadAssimp(const std::string &file_path,glm::mat4 transform){
  // Create an instance of the Importer class
  Assimp::Importer importer;

  // And have it read the given file with some example postprocessing
  // Usually - if speed is not the most important aspect for you - you'll
  // probably to request more postprocessing than we do in this example.
  aiScene *scene = importer.ReadFile(
      file_path, aiProcess_CalcTangentSpace | aiProcess_Triangulate |
                 aiProcess_JoinIdenticalVertices | aiProcess_SortByPType);

  // If the import failed, report it
  if (nullptr != scene) {
    LAND_WARN(importer.GetErrorString());
    return false;
  }
  auto root = scene->mRootNode;
  FindMesh(scene, root,transform);
  
  // Now we can access the file's contents.
  return true;


}*/
inline float PowerHeuristic(int nf, float fPdf, int ng, float gPdf) {
  float f = nf * fPdf, g = ng * gPdf;
  return (f * f) / (f * f + g * g);
}
glm::vec3 Scene::SampleLight(glm::vec3 direction,
                             HitRecord &hit_record,
                             std::mt19937 &rd,float time) const {
  float test = std::uniform_real_distribution<float>(0.0f, 1.0f)(rd);
  int index = 0;
  for (index = 0; index < light_distribution.size(); ++index) {
    test -= light_distribution[index];
    if (test <= 0) {
      break;
    }
  }
  if (index == light_distribution.size())  // boundary error
    return glm::vec3(0.0f);
  //printf("%d", index);
  float discrete_pdf = light_distribution[index];
  index = light_id[index];
  glm::vec3 Ld(0.f);
  glm::vec3 wi;
  float lightpdf = 0, scatterpdf = 0;
  bool visible = false;
  glm::vec3 Li;
  auto &material = GetEntity(hit_record.hit_entity_id).GetMaterial();
  //index = 0;
  if (index == -1) {
    Li = SampleEnvmap_Li(rd, &wi, &lightpdf);
  } else {
    Li = entities_[index].Sample_Li(hit_record, rd, &wi, &lightpdf);
  }
  HitRecord vis_record;
  glm::vec3 tp_direction = wi;
  if (lightpdf > 0 && Li != glm::vec3(0.f)) {
    // Compute BSDF or phase function's value for light sample
    glm::vec3 f;
    f = material.f(hit_record, -direction, wi) *
        abs(dot(wi, hit_record.normal));
    //f = material.f(hit_record, -direction, wi);
    scatterpdf = material.pdf(hit_record, -direction, wi);
    //return glm::vec3(scatterpdf);
    //return material.f(hit_record, -direction, wi);
    if (f != glm::vec3(0.f)) {
      // Compute effect of visibility for light source sample
      vis_record.hit_entity_id = -1;
      float t = TraceRay(hit_record.position, tp_direction, 1e-3f, 1e4f, &vis_record,time);
      //print(visible)
      //if (index == -1)
        //return 3.0f*f*Li/lightpdf;
      if (vis_record.hit_entity_id == index) {
        visible = true;
      }
      if (!visible) {
        Li = glm::vec3(0.f);
      }
      //return glm::vec3(3.0)*f;
      // Add light's contribution to reflected radiance
      if (Li != glm::vec3(0.0f)) {
        if (false)  // Is point light
          Ld += f * Li / lightpdf;
        else {
          float weight = PowerHeuristic(1, lightpdf, 1, scatterpdf);
          //weight = 0;//debug------------------------------
          Ld += f * Li * weight / lightpdf;
        }
      }
    }
  }
  //return Ld;
  // Sample BSDF with multiple importance sampling
  lightpdf=0;
  if (true) {  //! Is Delta Light
    glm::vec3 f;
    bool sampledSpecular = false;
    // Sample scattered direction for surface interactions
    f = material.Sample_f(hit_record, rd, -direction, &wi, &scatterpdf);
    f *= abs(dot(wi, hit_record.normal));
    //return f;
    sampledSpecular = material.IsSpecular();
    float area, normal, square;
    if (f != glm::vec3(0.0f) && scatterpdf > 0) {
      // Account for light contributions along sampled direction _wi_
      float weight = 1;
      // Find intersection and compute transmittance
      HitRecord light_record;
      glm::vec3 new_origin = hit_record.position;
      glm::vec3 new_direction = wi;
      float t = TraceRay(new_origin, new_direction, 1e-3f, 1e4f, &light_record,time);
      //return glm::vec3(0, 0, float(light_record.hit_entity_id == 2));
      // Add light contribution from material sampling
      glm::vec3 Li(0.f);
      if (t > 0.0f) {
        auto &new_material = GetEntity(light_record.hit_entity_id).GetMaterial();
        //return new_material.albedo_color;
        Li = (light_record.hit_entity_id == index && light_record.front_face)
            ? new_material.emission * new_material.emission_strength
                 : glm::vec3(0.f);
      } else {
        //return glm::vec3(0.0, 1.0, 1.0);
        if (index==-1)
          Li = SampleEnvmap(wi);
      }
      if (Li != glm::vec3(0.0f)) {
        if (!sampledSpecular) {
          if (index == -1) {
            lightpdf = EnvMapPdfLi(hit_record, wi);
          } else {
            glm::mat4 transform =
                (glm::mat4)entities_[index].GetTransformMatrix();
            if (time > 0 && entities_[index].duration != 0.0f)
              transform = transform*matrix_interpolation(
                  entities_[index].anime_transform_, time / entities_[index].duration);
            lightpdf = dot(light_record.position - hit_record.position,
                           light_record.position - hit_record.position) /
                       abs(dot(light_record.normal, -wi));
            int face_id = light_record.index_id;
            auto &indices = entities_[index].GetModel()->GetIndices();
            auto &vertices = entities_[index].GetModel()->GetVertices();
            int num_faces = indices.size() / 3;
            glm::vec3 p0 = glm::vec3{
                transform *
                glm::vec4{vertices[indices[face_id * 3]].position, 1.0f}};
            glm::vec3 p1 = glm::vec3{
                transform *
                glm::vec4{vertices[indices[face_id * 3 + 1]].position, 1.0f}};
            glm::vec3 p2 = glm::vec3{
                transform *
                glm::vec4{vertices[indices[face_id * 3 + 2]].position, 1.0f}};
            lightpdf /=
                0.5 * glm::length(glm::cross(p1 - p0, p2 - p0)) * num_faces;
          }
          // printf("prev:%f now:%f\n", tp, lightpdf);
          if (lightpdf == 0)
            return Ld / discrete_pdf;
          weight = PowerHeuristic(1, scatterpdf, 1, lightpdf);
        }
        //return normalize(hit_record.position);
        //weight = 1;  // debug-------------------------------------------------
        Ld += f * Li * weight / scatterpdf;
      }
    }
  }
  return Ld / discrete_pdf;
}
glm::vec3 Scene::SampleEnvmap_Li(std::mt19937 &rd,
                          glm::vec3 *wi,
                          float *lightpdf) const {
  float z = std::uniform_real_distribution<float>(-1.0f, 1.0f)(rd);
  float r = std::sqrt(std::max((float)0, (float)1. - z * z));
  float phi =
      2 * glm::pi<float>() * std::uniform_real_distribution<float>(0.0f, 1.0f)(rd);
  *lightpdf = glm::one_over_two_pi<float>()/2;
  *wi = glm::vec3(r * std::cos(phi), r * std::sin(phi), z);
  return SampleEnvmap(*wi);

}
float Scene::EnvMapPdfLi(HitRecord &hit_record, glm::vec3 wi) const {
  return glm::one_over_two_pi<float>()/2;
}
}  // namespace sparks
