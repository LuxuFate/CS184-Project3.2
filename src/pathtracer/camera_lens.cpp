#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "CGL/misc.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;

namespace CGL {

using Collada::CameraInfo;

Ray Camera::generate_ray_for_thin_lens(double x, double y, double rndR, double rndTheta) const {

  // TODO Project 3-2: Part 4
  // compute position and direction of ray from the input sensor sample coordinate.
  // Note: use rndR and rndTheta to uniformly sample a unit disk.
    double RhFov = tan(radians(hFov) * 0.5);
    double RvFov = tan(radians(vFov) * 0.5);
    Vector3D bottomLeft = Vector3D(-RhFov, -RvFov, -1.0);
    Vector3D topRight = Vector3D(RhFov, RvFov, -1.0);
    double Cx = bottomLeft.x + (topRight.x - bottomLeft.x) * x;
    double Cy = bottomLeft.y + (topRight.y - bottomLeft.y) * y;
    Vector3D sensor = Vector3D(Cx, Cy, -1.0);
    
    Vector3D pLens = Vector3D(lensRadius * sqrt(rndR) * cos(rndTheta),                                                lensRadius * sqrt(rndR) * sin(rndTheta), 0);
    Vector3D pFocus = (sensor * focalDistance) - pLens;
    
    Vector3D dir = c2w * pFocus;
    dir.normalize();
    
    Ray ray = Ray(pos + (c2w * pLens), dir);
    ray.min_t = nClip;
    ray.max_t = fClip;
    return ray;
}


} // namespace CGL
