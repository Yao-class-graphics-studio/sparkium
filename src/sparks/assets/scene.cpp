#include "sparks/assets/scene.h"

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "imgui.h"
#include "sparks/assets/accelerated_mesh.h"
#include "sparks/util/util.h"

// #define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

namespace sparks {

Scene::Scene() {
  AddTexture(Texture(1, 1, glm::vec4{1.0f}, SAMPLE_TYPE_LINEAR), "Pure White");
  AddTexture(Texture(1, 1, glm::vec4{0.0f}, SAMPLE_TYPE_LINEAR), "Pure Black");
  Texture envmap;
  Texture::Load(u8"../../textures/envmap_clouds_4k.hdr", envmap);
  envmap.SetSampleType(SAMPLE_TYPE_LINEAR);
  envmap_id_ = AddTexture(envmap, "Clouds");
  // AddEntity(
  //     AcceleratedMesh({{{-1.0f, 0.0f, 1.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f}},
  //                      {{-1.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 1.0f}},
  //                      {{1.0f, 0.0f, 1.0f}, {0.0f, 1.0f, 0.0f}, {1.0f, 0.0f}},
  //                      {{1.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f}, {1.0f, 1.0f}}},
  //                     {0, 1, 2, 2, 1, 3}),
  //     Material{}, glm::mat4{1.0f});
  SetCameraToWorld(glm::inverse(glm::lookAt(glm::vec3{0.278f, 0.273f, -0.800f},
                                            glm::vec3{0.278f, 0.273f, 0.2f},
                                            glm::vec3{0.0f, 1.0f, 0.0f})));

  // Texture texture;
  // Texture::Load("../../textures/earth.jpg", texture);
  // AddEntity(AcceleratedMesh(Mesh::Sphere(glm::vec3{0.0f, 0.0f, 0.0f}, 0.5f)),
  //           Material{glm::vec3{1.0f}, 0, glm::vec3{1.0f}, AddTexture(texture, "Earth")},
  //           glm::translate(glm::mat4{1.0f}, glm::vec3{0.0f, 0.5f, 0.0f}));
  // AddEntity(AcceleratedMesh(Mesh::Sphere(glm::vec3{0.0f, 0.0f, 0.0f}, 0.5f)),
  //           Material{glm::vec3{1.0f}, 0, glm::vec3{1.0f, 0.0f, 0.0f}},
  //           glm::translate(glm::mat4{1.0f}, glm::vec3{0.0f, 0.5f, 0.0f}));
  LoadObjFile(
      "D:/Seafile/Yao Class/Semester 5/Graphics/models/nashida/nashida.obj", glm::mat4{0.05f});
  // LoadObjFile(
  //     "D:/Seafile/Yao Class/Semester "
  //     "5/Graphics/models/cornellbox/cornell_box.obj",
  //     glm::translate(glm::mat4{0.001f}, glm::vec3{0.0f, 0.0f, 0.0f}));
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
                      HitRecord *hit_record) const {
  float result = -1.0f;
  HitRecord local_hit_record;
  float local_result;
  for (int entity_id = 0; entity_id < entities_.size(); entity_id++) {
    auto &entity = entities_[entity_id];
    auto &transform = entity.GetTransformMatrix();
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

int Scene::LoadObjFile(const std::string &file_path, const glm::mat4 &transform = glm::mat4{1.0f}) {
  auto cal_tangent = [](Vertex &v0, Vertex &v1, Vertex &v2) {
    glm::vec3 E0 = v1.position - v0.position;
    glm::vec3 E1 = v2.position - v1.position;
    glm::vec3 E2 = v0.position - v2.position;
    glm::vec2 D0 = v1.tex_coord - v0.tex_coord;
    glm::vec2 D1 = v2.tex_coord - v1.tex_coord;
    glm::vec2 D2 = v0.tex_coord - v2.tex_coord;
    glm::vec3 T0 = (D0.y * E2 - D2.y * E0) / (D0.y * D2.x - D2.y * D0.x);
    glm::vec3 T1 = (D1.y * E0 - D0.y * E1) / (D1.y * D0.x - D0.y * D1.x);
    glm::vec3 T2 = (D2.y * E1 - D1.y * E2) / (D2.y * D1.x - D1.y * D2.x);
    v0.tangent = glm::normalize(T0 - glm::dot(T0, v0.normal) * v0.normal);
    v1.tangent = glm::normalize(T1 - glm::dot(T1, v1.normal) * v1.normal);
    v2.tangent = glm::normalize(T2 - glm::dot(T2, v2.normal) * v2.normal);
  };
  auto default_tangent = [](Vertex &v0, Vertex &v1, Vertex &v2) {
    glm::vec3 T0 = v1.position - v0.position;
    v0.tangent = glm::normalize(T0 - glm::dot(T0, v0.normal) * v0.normal);
    glm::vec3 T1 = v2.position - v1.position;
    v1.tangent = glm::normalize(T1 - glm::dot(T1, v1.normal) * v1.normal);
    glm::vec3 T2 = v0.position - v2.position;
    v2.tangent = glm::normalize(T2 - glm::dot(T2, v2.normal) * v2.normal);
  };

  int entity_cnt = 0;
  std::string dir_path;
  size_t pos = file_path.find_last_of("///");
  if (pos != std::string::npos) {
    dir_path = file_path.substr(0, pos);
  }

  tinyobj::ObjReaderConfig reader_config;
  // reader_config.mtl_search_path = "./";  // Path to material files

  tinyobj::ObjReader reader;

  if (!reader.ParseFromFile(file_path, reader_config)) {
    if (!reader.Error().empty()) {
      LAND_WARN("[Load obj, ERROR]: {}", reader.Error());
    }
    return false;
  }

  if (!reader.Warning().empty()) {
    LAND_WARN("{}", reader.Warning());
  }

  auto &attrib = reader.GetAttrib();
  auto &shapes = reader.GetShapes();
  auto &materials = reader.GetMaterials();

  // std::cout << "shapes: " << shapes.size() << std::endl;

  // Loop over shapes
  std::vector<Vertex> vertices;
  std::vector<uint32_t> indices;
  for (size_t s = 0; s < shapes.size(); s++) {
    // std::cout << "shape: " << s << std::endl;
    // Loop over faces(polygon)
    int last_material_id = -1;
    size_t index_offset = 0;
    for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
      int material_id = shapes[s].mesh.material_ids[f];
      // std::cout << "face: " << f <<  " material_id: " << material_id << std::endl;
      if (last_material_id == -1) {
        last_material_id = material_id;
      } else if (material_id != last_material_id) {
        // std::cout << "last:" << last_material_id << " cur: " << material_id << std::endl;

        Mesh mesh = Mesh(vertices, indices);
        mesh.MergeVertices();

        vertices.clear();
        indices.clear();

        Material material{};

        tinyobj::material_t mtr = materials[last_material_id];
        material.ambient.r = mtr.ambient[0];
        material.ambient.g = mtr.ambient[1];
        material.ambient.b = mtr.ambient[2];
        Texture ambient_texture;
        if (Texture::Load(dir_path + "/" + mtr.ambient_texname, ambient_texture)) {
          material.ambient_texture_id = AddTexture(ambient_texture, mtr.ambient_texname);
          material.ambient = glm::vec3{1.0f};
        }

        material.diffuse.r = mtr.diffuse[0];
        material.diffuse.g = mtr.diffuse[1];
        material.diffuse.b = mtr.diffuse[2];
        Texture diffuse_texture;
        if (Texture::Load(dir_path + "/" + mtr.diffuse_texname, diffuse_texture)) {
          material.diffuse_texture_id = AddTexture(diffuse_texture, mtr.diffuse_texname);
          material.diffuse = glm::vec3{1.0f};
        }

        material.specular.r = mtr.specular[0];
        material.specular.g = mtr.specular[1];
        material.specular.b = mtr.specular[2];
        Texture specular_texture;
        if (Texture::Load(dir_path + "/" + mtr.specular_texname, specular_texture)) {
          material.specular_texture_id = AddTexture(specular_texture, mtr.specular_texname);
          material.specular = glm::vec3{1.0f};
        }

        material.transmittance.r = mtr.transmittance[0];
        material.transmittance.g = mtr.transmittance[1];
        material.transmittance.b = mtr.transmittance[2];

        material.emission.r = mtr.emission[0];
        material.emission.g = mtr.emission[1];
        material.emission.b = mtr.emission[2];

        int cnt = 0;
        if (mtr.emission[0] > 0 || mtr.emission[1] > 0 || mtr.emission[2] > 0) {
          material.material_type = MATERIAL_TYPE_EMISSION;
          cnt++;
        }
        if (mtr.transmittance[0] > 0 || mtr.transmittance[1] > 0 || mtr.transmittance[2] > 0) {
          material.material_type = MATERIAL_TYPE_TRANSMISSIVE;
          cnt++;
        }
        // std::cout << "face " << f << " specular: " << material.specular_texture_id << std::endl;
        if (material.specular_texture_id > 0 || mtr.specular[0] > 0 || mtr.specular[1] > 0 || mtr.specular[2] > 0) {
          material.material_type = MATERIAL_TYPE_SPECULAR;
          cnt++;
        }
        // std::cout << "face " << f << " diffuse: " << material.diffuse_texture_id << std::endl;
        if (material.diffuse_texture_id > 0 || mtr.diffuse[0] > 0 || mtr.diffuse[1] > 0 || mtr.diffuse[2] > 0) {
          material.material_type = MATERIAL_TYPE_LAMBERTIAN;
          cnt++;
        }
        if (cnt >= 2) {
          material.material_type = MATERIAL_TYPE_PRINCIPLED;
        }
        // std::cout << "material_type is: " << material.material_type << std::endl;
        AddEntity(AcceleratedMesh(mesh), material, transform);
        entity_cnt++;
        last_material_id = material_id;
      }
      size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);

      // Loop over vertices in the face.
      std::vector<Vertex> face_vertices;
      for (size_t v = 0; v < fv; v++) {
        Vertex vertex{};
        // access to vertex
        tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
        tinyobj::real_t vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
        tinyobj::real_t vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
        tinyobj::real_t vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];
        vertex.position = {vx, vy, vz};
        // Check if `normal_index` is zero or positive. negative = no normal
        // data
        if (idx.normal_index >= 0) {
          tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
          tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
          tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];
          vertex.normal = {nx, ny, nz};
        } else {
          vertex.normal = {0.0f, 0.0f, 0.0f};
        }

        // Check if `texcoord_index` is zero or positive. negative = no texcoord
        // data
        if (idx.texcoord_index >= 0) {
          tinyobj::real_t tx =
              attrib.texcoords[2 * size_t(idx.texcoord_index) + 0];
          tinyobj::real_t ty =
              attrib.texcoords[2 * size_t(idx.texcoord_index) + 1];
          vertex.tex_coord = {tx, ty};
        }
        else {
          vertex.tex_coord = {-1.0f, -1.0f};
        }
        face_vertices.push_back(vertex);
      }

      for (int i = 1; i < face_vertices.size() - 1; i++) {
        Vertex v0 = face_vertices[0];
        Vertex v1 = face_vertices[i];
        Vertex v2 = face_vertices[i + 1];
        auto geometry_normal = glm::normalize(
            glm::cross(v2.position - v0.position, v1.position - v0.position));
        if (v0.normal == glm::vec3{0.0f, 0.0f, 0.0f}) {
          v0.normal = geometry_normal;
        } else if (glm::dot(geometry_normal, v0.normal) < 0.0f) {
          v0.normal = -v0.normal;
        }
        if (v1.normal == glm::vec3{0.0f, 0.0f, 0.0f}) {
          v1.normal = geometry_normal;
        } else if (glm::dot(geometry_normal, v1.normal) < 0.0f) {
          v1.normal = -v1.normal;
        }
        if (v2.normal == glm::vec3{0.0f, 0.0f, 0.0f}) {
          v2.normal = geometry_normal;
        } else if (glm::dot(geometry_normal, v2.normal) < 0.0f) {
          v2.normal = -v2.normal;
        }
        if (v0.tex_coord == glm::vec2{-1.0f, -1.0f}) {
          v0.tex_coord = glm::vec2{0.0f, 0.0f};
          default_tangent(v0, v1, v2);
        }
        else {
          cal_tangent(v0, v1, v2);
        }
        indices.push_back(vertices.size());
        indices.push_back(vertices.size() + 1);
        indices.push_back(vertices.size() + 2);
        vertices.push_back(v0);
        vertices.push_back(v1);
        vertices.push_back(v2);
      }

      index_offset += fv;
    }
    Mesh mesh = Mesh(vertices, indices);
    mesh.MergeVertices();

    vertices.clear();
    indices.clear();

    Material material{};

    tinyobj::material_t mtr = materials[last_material_id];
    material.ambient.r = mtr.ambient[0];
    material.ambient.g = mtr.ambient[1];
    material.ambient.b = mtr.ambient[2];
    Texture ambient_texture;
    if (Texture::Load(dir_path + "/" + mtr.ambient_texname, ambient_texture)) {
      material.ambient_texture_id = AddTexture(ambient_texture, mtr.ambient_texname);
      material.ambient = glm::vec3(1.0f);
    }

    material.diffuse.r = mtr.diffuse[0];
    material.diffuse.g = mtr.diffuse[1];
    material.diffuse.b = mtr.diffuse[2];
    Texture diffuse_texture;
    if (Texture::Load(dir_path + "/" + mtr.diffuse_texname, diffuse_texture)) {
      material.diffuse_texture_id = AddTexture(diffuse_texture, mtr.diffuse_texname);
      material.diffuse = glm::vec3(1.0f);
    }

    material.specular.r = mtr.specular[0];
    material.specular.g = mtr.specular[1];
    material.specular.b = mtr.specular[2];
    Texture specular_texture;
    if (Texture::Load(dir_path + "/" + mtr.specular_texname, specular_texture)) {
      material.specular_texture_id = AddTexture(specular_texture, mtr.specular_texname);
      material.specular = glm::vec3(1.0f);
    }

    material.transmittance.r = mtr.transmittance[0];
    material.transmittance.g = mtr.transmittance[1];
    material.transmittance.b = mtr.transmittance[2];

    material.emission.r = mtr.emission[0];
    material.emission.g = mtr.emission[1];
    material.emission.b = mtr.emission[2];

    int cnt = 0;
    if (mtr.emission[0] > 0 || mtr.emission[1] > 0 || mtr.emission[2] > 0) {
      material.material_type = MATERIAL_TYPE_EMISSION;
      cnt++;
      // std::cout << "emission" << std::endl;
    }
    if (mtr.transmittance[0] > 0 || mtr.transmittance[1] > 0 || mtr.transmittance[2] > 0) {
      material.material_type = MATERIAL_TYPE_TRANSMISSIVE;
      cnt++;
      // std::cout << "transmittance" << std::endl;
    }
    if (material.specular_texture_id > 0 || mtr.specular[0] > 0 || mtr.specular[1] > 0 || mtr.specular[2] > 0) {
      material.material_type = MATERIAL_TYPE_SPECULAR;
      cnt++;
      // std::cout << "specular " << material.specular_texture_id << std::endl;
    }
    if (material.diffuse_texture_id > 0 || mtr.diffuse[0] > 0 || mtr.diffuse[1] > 0 || mtr.diffuse[2] > 0) {
      material.material_type = MATERIAL_TYPE_LAMBERTIAN;
      cnt++;
      // std::cout << "diffuse " << material.diffuse_texture_id << std::endl;
    }
    if (cnt >= 2) {
      material.material_type = MATERIAL_TYPE_PRINCIPLED;
    }
    // std::cout << "material_type is: " << material.material_type << std::endl;
    AddEntity(AcceleratedMesh(mesh), material, transform);
    entity_cnt++;
    // std::cout << entity_cnt << std::endl;
    // Mesh mesh = Mesh(vertices, indices);
    // mesh.MergeVertices();
    // Material material{};
    // int material_id = -1;
    // if (shapes[s].mesh.material_ids.size() > 0)
    // {
    //   material_id = shapes[s].mesh.material_ids[0];
    // }
    // if (material_id >= 0) {
    //   tinyobj::material_t mtr = materials[material_id];
    //   material.ambient.r = mtr.ambient[0];
    //   material.ambient.g = mtr.ambient[1];
    //   material.ambient.b = mtr.ambient[2];
    //   Texture ambient_texture;
    //   if (Texture::Load(dir_path + "/" + mtr.ambient_texname, ambient_texture));
    //     material.ambient_texture_id = AddTexture(ambient_texture, mtr.ambient_texname);

    //   material.diffuse.r = mtr.diffuse[0];
    //   material.diffuse.g = mtr.diffuse[1];
    //   material.diffuse.b = mtr.diffuse[2];
    //   Texture diffuse_texture;
    //   if (Texture::Load(dir_path + "/" + mtr.diffuse_texname, diffuse_texture));
    //     material.diffuse_texture_id = AddTexture(diffuse_texture, mtr.diffuse_texname);

    //   material.specular.r = mtr.specular[0];
    //   material.specular.g = mtr.specular[1];
    //   material.specular.b = mtr.specular[2];
    //   Texture specular_texture;
    //   if (Texture::Load(dir_path + "/" + mtr.specular_texname, specular_texture));
    //     material.specular_texture_id = AddTexture(specular_texture, mtr.specular_texname);

    //   material.transmittance.r = mtr.transmittance[0];
    //   material.transmittance.g = mtr.transmittance[1];
    //   material.transmittance.b = mtr.transmittance[2];

    //   material.emission.r = mtr.emission[0];
    //   material.emission.g = mtr.emission[1];
    //   material.emission.b = mtr.emission[2];

    //   int cnt = 0;
    //   if (mtr.emission[0] > 0 || mtr.emission[1] > 0 || mtr.emission[2] > 0) {
    //     material.material_type = MATERIAL_TYPE_EMISSION;
    //     cnt++;
    //   }
    //   if (mtr.transmittance[0] > 0 || mtr.transmittance[1] > 0 || mtr.transmittance[2] > 0) {
    //     material.material_type = MATERIAL_TYPE_TRANSMISSIVE;
    //     cnt++;
    //   }
    //   if (material.specular_texture_id > 0 || mtr.specular[0] > 0 || mtr.specular[1] > 0 || mtr.specular[2] > 0) {
    //     material.material_type = MATERIAL_TYPE_SPECULAR;
    //     cnt++;
    //   }
    //   if (material.diffuse_texture_id > 0 || mtr.diffuse[0] > 0 || mtr.diffuse[1] > 0 || mtr.diffuse[2] > 0) {
    //     material.material_type = MATERIAL_TYPE_LAMBERTIAN;
    //     cnt++;
    //   }
    //   if (cnt > 2) {
    //     material.material_type = MATERIAL_TYPE_PRINCIPLED;
    //   }
    // }
    // std::cout << "material_type is: " << material.material_type << std::endl;
    // AddEntity(AcceleratedMesh(mesh), material, glm::mat4{0.05f});
  }
  return entity_cnt;
}

}  // namespace sparks
