float_attr_list = ["metallic", "eta", "roughness", "specularTint", "anisotropic", "sheen", "sheenTint", "clearcoat", "clearcoatGloss", "specTrans", "flatness", "diffTrans"]
vec3_attr_list = ["scatterDistance"]

for attr in float_attr_list:
    print("""  child_element = material_element->FirstChildElement("{attr}");
  if (child_element) {{
    {attr} = std::stof(child_element->FindAttribute("value")->Value());
  }}
""".format(attr=attr))
    
for attr in vec3_attr_list:
    print("""  child_element = material_element->FirstChildElement("{attr}");
  if (child_element) {{
    {attr} = StringToVec3(child_element->FindAttribute("value")->Value());
  }}
""".format(attr=attr))