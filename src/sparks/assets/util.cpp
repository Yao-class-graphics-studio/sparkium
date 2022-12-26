#include "sparks/assets/util.h"

#include "glm/gtc/matrix_transform.hpp"

#include "sstream"
#include "vector"
#include "unordered_map"
#include "iostream"

namespace sparks {
glm::vec3 DecomposeRotation(glm::mat3 R) {
  return {
      std::atan2(-R[2][1], std::sqrt(R[0][1] * R[0][1] + R[1][1] * R[1][1])),
      std::atan2(R[2][0], R[2][2]), std::atan2(R[0][1], R[1][1])};
}

glm::mat4 ComposeRotation(glm::vec3 pitch_yaw_roll) {
  return glm::rotate(glm::mat4{1.0f}, pitch_yaw_roll.y,
                     glm::vec3{0.0f, 1.0f, 0.0f}) *
         glm::rotate(glm::mat4{1.0f}, pitch_yaw_roll.x,
                     glm::vec3{1.0f, 0.0f, 0.0f}) *
         glm::rotate(glm::mat4{1.0f}, pitch_yaw_roll.z,
                     glm::vec3{0.0f, 0.0f, 1.0f});
}

glm::vec3 string2vec3(const std::string &s) {
  std::istringstream ss(s);
  std::vector<float> v;
  std::string word;
  while (ss >> word)
    v.push_back(std::stof(word));
  return glm::vec3(v[0], v[1], v[2]);
}

glm::vec4 string2vec4(const std::string &s) {
  std::istringstream ss(s);
  std::vector<float> v;
  std::string word;
  while (ss >> word)
    v.push_back(std::stof(word));
  return glm::vec4(v[0], v[1], v[2], v[3]);
}

glm::mat4 transformMatrix(tinyxml2::XMLElement *transformElement) {
  if (!transformElement)
    return glm::mat4{1.0f};
  std::string transformType = transformElement->FindAttribute("type")->Value();
  if (transformType == "lookat") {
    std::unordered_map<std::string, glm::vec3> vec3mp;
    for (tinyxml2::XMLElement *childElement =
             transformElement->FirstChildElement("vec3");
         childElement; childElement = childElement->NextSiblingElement("vec3"))
      vec3mp[childElement->FindAttribute("name")->Value()] =
          string2vec3(childElement->FindAttribute("value")->Value());

    return glm::inverse(
        glm::lookAt(vec3mp["eye"], vec3mp["center"], vec3mp["up"]));
  } else if (transformType == "translate") {
    glm::vec3 translation = string2vec3(
        transformElement->FirstChildElement()->FindAttribute("value")->Value());
    return glm::translate(glm::mat4{1.0f}, translation);
  } else {
    std::cout << "Unknown Transformation Type: " << transformType << std::endl;
    return glm::mat4{1.0f};
  }
}

std::string getBaseName(const std::string &s) {
  size_t a = s.find_last_of('/');
  if (a == std::string::npos)
    a = 0;
  else
    a++;
  size_t b = s.find_last_of('.');
  return s.substr(a, b - a);
}

}  // namespace sparks
