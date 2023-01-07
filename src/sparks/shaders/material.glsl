
struct Material {
  // vec3 ambient;
  // int ambient_texture_id;
  vec3 diffuse;
  int diffuse_texture_id;

  vec3 specular;
  int specular_texture_id;

  vec3 opacity;
  int opacity_texture_id;

  vec3 emission;
  int emission_texture_id;

  vec3 transmittance;
  float ior;

  float roughness;
  int roughness_texture_id;
  float metallic;
  int metallic_texture_id;

  float sheen;
  int sheen_texture_id;
  float clearcoat_thickness;
  float clearcoat_roughness;

  float anisotropy;
  float anisotropy_rotation;
  int normal_texture_id;
  uint material_type;
};

#define MATERIAL_TYPE_LAMBERTIAN 0
#define MATERIAL_TYPE_SPECULAR 1
#define MATERIAL_TYPE_TRANSMISSIVE 2
#define MATERIAL_TYPE_PRINCIPLED 3
#define MATERIAL_TYPE_EMISSION 4
