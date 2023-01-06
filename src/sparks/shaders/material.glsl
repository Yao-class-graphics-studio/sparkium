
struct Material {
  vec3 ambient;
  int ambient_texture_id;
  vec3 diffuse;
  int diffuse_texture_id;
  vec3 specular;
  int specular_texture_id;
  vec3 transmittance;
  float roughness;
  vec3 emission;
  int ior;
  uint material_type;
};

#define MATERIAL_TYPE_LAMBERTIAN 0
#define MATERIAL_TYPE_SPECULAR 1
#define MATERIAL_TYPE_TRANSMISSIVE 2
#define MATERIAL_TYPE_PRINCIPLED 3
#define MATERIAL_TYPE_EMISSION 4
