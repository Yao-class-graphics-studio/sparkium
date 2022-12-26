#pragma once
#include "glm/glm.hpp"
#include "sparks/util/tinyxml2.h"
#include "string"

namespace sparks {
glm::vec3 DecomposeRotation(glm::mat3 R);

glm::mat4 ComposeRotation(glm::vec3 pitch_yaw_roll);

glm::vec3 string2vec3(const std::string &s);

glm::vec4 string2vec4(const std::string &s);

glm::mat4 transformMatrix(tinyxml2::XMLElement *transformElement);

std::string getBaseName(const std::string &s);
}  // namespace sparks
