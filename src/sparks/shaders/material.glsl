
struct Material {
  vec3 ambient;
  int ambient_texture_id;
  vec3 diffuse;
  int diffuse_texture_id;
  vec3 specular;
  int specular_texture_id;
  vec3 transmittance;
  int reserve0;
  vec3 emission;
  int reserve1;
  uint material_type;
};

#define MATERIAL_TYPE_LAMBERTIAN 0
#define MATERIAL_TYPE_SPECULAR 1
#define MATERIAL_TYPE_TRANSMISSIVE 2
#define MATERIAL_TYPE_PRINCIPLED 3
#define MATERIAL_TYPE_EMISSION 4
