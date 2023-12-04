#include "sparks/assets/subsurface.h"

namespace sparks {

BSSRDF* SubsurfaceMaterial::ComputeBSSRDF(const HitRecord &hit, const glm::vec3 direction, const Scene* scene) const {
    BSSRDF *newBSSRDF = new TabulatedBSSRDF(hit.position, direction, eta, BSSRDF_TABULATED, hit.geometry_normal, 
                                            this, sigma_a, sigma_s, table);
    return newBSSRDF;
}

}