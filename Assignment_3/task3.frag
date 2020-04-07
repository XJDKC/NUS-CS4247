//============================================================
// STUDENT NAME: Xing Rulin
// STUDENT NO.: A0214251W
// NUS EMAIL ADDRESS: E0518535@u.nus.edu
// COMMENTS TO GRADER: None
//
//============================================================

// FRAGMENT SHADER FOR SHADERTOY
// Run this at https://www.shadertoy.com/new
// See documentation at https://www.shadertoy.com/howto

// Your browser must support WebGL 2.0.
// Check your browser at https://webglreport.com/?v=2

//============================================================================
// Constants.
//============================================================================

const float PI = 3.1415926536;

const vec3 BACKGROUND_COLOR = vec3(0.1, 0.2, 0.6);

// Vertical field-of-view angle of camera. In degrees.
const float FOVY = 50.0;

// Use this for avoiding the "epsilon problem" or the shadow acne problem.
const float DEFAULT_TMIN = 10.0e-4;

// Use this for tmax for non-shadow ray intersection test.
const float DEFAULT_TMAX = 10.0e6;

// Equivalent to number of recursion levels (0 means ray-casting only).
// We are using iterations to replace recursions.
const int NUM_ITERATIONS = 2;

// Constants for the scene objects.
const int NUM_LIGHTS = 2;
const int NUM_MATERIALS = 3;
const int NUM_PLANES = 2;
const int NUM_SPHERES = 2;

//============================================================================
// Define new struct types.
//============================================================================
struct Ray_t {
  vec3 o; // Ray Origin.
  vec3 d; // Ray Direction. A unit vector.
};

struct Plane_t {
  // The plane equation is Ax + By + Cz + D = 0.
  float A, B, C, D;
  int materialID;
};

struct Sphere_t {
  vec3 center;
  float radius;
  int materialID;
};

struct Light_t {
  vec3 position; // Point light 3D position.
  vec3 I_a;      // For Ambient.
  vec3 I_source; // For Diffuse and Specular.
};

struct Material_t {
  vec3 k_a;  // Ambient coefficient.
  vec3 k_d;  // Diffuse coefficient.
  vec3 k_r;  // Reflected specular coefficient.
  vec3 k_rg; // Global reflection coefficient.
  float n;   // The specular reflection exponent. Ranges from 0.0 to 128.0.
};

//----------------------------------------------------------------------------
// The lighting model used here is similar to that on Slides 8 and 12 of
// Lecture Topic 9 (Ray Tracing). Here it is computed as
//
//     I_local = SUM_OVER_ALL_LIGHTS {
//                   I_a * k_a +
//                   k_shadow * I_source * [ k_d * (N.L) + k_r * (R.V)^n ]
//               }
// and
//     I = I_local  +  k_rg * I_reflected
//----------------------------------------------------------------------------

//============================================================================
// Global scene data.
//============================================================================
Plane_t Plane[NUM_PLANES];
Sphere_t Sphere[NUM_SPHERES];
Light_t Light[NUM_LIGHTS];
Material_t Material[NUM_MATERIALS];

/////////////////////////////////////////////////////////////////////////////
// Initializes the scene.
/////////////////////////////////////////////////////////////////////////////
void InitScene() {
  // Horizontal plane.
  Plane[0].A = 0.0;
  Plane[0].B = 1.0;
  Plane[0].C = 0.0;
  Plane[0].D = 0.0;
  Plane[0].materialID = 0;

  // Vertical plane.
  Plane[1].A = 0.0;
  Plane[1].B = 0.0;
  Plane[1].C = 1.0;
  Plane[1].D = 3.5;
  Plane[1].materialID = 0;

  // Center bouncing sphere.
  Sphere[0].center = vec3(0.0, abs(sin(2.0 * iTime)) + 0.7, 0.0);
  Sphere[0].radius = 0.7;
  Sphere[0].materialID = 1;

  // Circling sphere.
  Sphere[1].center = vec3(1.5 * cos(iTime), 0.5, 1.5 * sin(iTime));
  Sphere[1].radius = 0.5;
  Sphere[1].materialID = 2;

  // Silver material.
  Material[0].k_d = vec3(0.5, 0.5, 0.5);
  Material[0].k_a = 0.2 * Material[0].k_d;
  Material[0].k_r = 2.0 * Material[0].k_d;
  Material[0].k_rg = 0.5 * Material[0].k_r;
  Material[0].n = 64.0;

  // Gold material.
  Material[1].k_d = vec3(0.8, 0.7, 0.1);
  Material[1].k_a = 0.2 * Material[1].k_d;
  Material[1].k_r = 2.0 * Material[1].k_d;
  Material[1].k_rg = 0.5 * Material[1].k_r;
  Material[1].n = 64.0;

  // Green plastic material.
  Material[2].k_d = vec3(0.0, 0.8, 0.0);
  Material[2].k_a = 0.2 * Material[2].k_d;
  Material[2].k_r = vec3(1.0, 1.0, 1.0);
  Material[2].k_rg = 0.5 * Material[2].k_r;
  Material[2].n = 128.0;

  // Light 0.
  Light[0].position = vec3(4.0, 8.0, -3.0);
  Light[0].I_a = vec3(0.1, 0.1, 0.1);
  Light[0].I_source = vec3(1.0, 1.0, 1.0);

  // Light 1.
  Light[1].position = vec3(-4.0, 8.0, 0.0);
  Light[1].I_a = vec3(0.1, 0.1, 0.1);
  Light[1].I_source = vec3(1.0, 1.0, 1.0);
}

/////////////////////////////////////////////////////////////////////////////
// Returns a random number between 0 and 1.
//
// This pseudorandom number generator is based on the 32-bit combined LFSR
// generator proposed in the paper "Tables of Maximally-Equidistributed
// Combined LFSR Generators" by Pierre L'Ecuyer.
// (http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.43.3639)
/////////////////////////////////////////////////////////////////////////////

// VERY IMPORTANT: The initial seeds rand_z1, rand_z2, rand_z3, rand_z4
// must be larger than 1, 7, 15, and 127 respectively.
const uint CONST_RAND_SEED = 987654321U;
uint rand_z1 = uint(CONST_RAND_SEED + 2U);
uint rand_z2 = uint(CONST_RAND_SEED + 8U);
uint rand_z3 = uint(CONST_RAND_SEED + 16U);
uint rand_z4 = uint(CONST_RAND_SEED + 128U);

float rand(void) {
  uint b = ((rand_z1 << 6) ^ rand_z1) >> 13;
  rand_z1 = ((rand_z1 & 4294967294U) << 18) ^ b;
  b = ((rand_z2 << 2) ^ rand_z2) >> 27;
  rand_z2 = ((rand_z2 & 4294967288U) << 2) ^ b;
  b = ((rand_z3 << 13) ^ rand_z3) >> 21;
  rand_z3 = ((rand_z3 & 4294967280U) << 7) ^ b;
  b = ((rand_z4 << 3) ^ rand_z4) >> 12;
  rand_z4 = ((rand_z4 & 4294967168U) << 13) ^ b;
  return float(rand_z1 ^ rand_z2 ^ rand_z3 ^ rand_z4) * 2.3283064365386963e-10;
}

/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a plane and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is such an intersection, outputs the value of t, the position
// of the intersection (hitPos) and the normal vector at the intersection
// (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectPlane(in Plane_t pln, in Ray_t ray, in float tmin, in float tmax,
                    out float t, out vec3 hitPos, out vec3 hitNormal) {
  vec3 N = vec3(pln.A, pln.B, pln.C);
  float NRd = dot(N, ray.d);
  float NRo = dot(N, ray.o);
  float t0 = (-pln.D - NRo) / NRd;
  if (t0 < tmin || t0 > tmax)
    return false;

  // We have a hit -- output results.
  t = t0;
  hitPos = ray.o + t0 * ray.d;
  hitNormal = normalize(N);
  return true;
}

/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a plane and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectPlane(in Plane_t pln, in Ray_t ray, in float tmin,
                    in float tmax) {
  vec3 N = vec3(pln.A, pln.B, pln.C);
  float NRd = dot(N, ray.d);
  float NRo = dot(N, ray.o);
  float t0 = (-pln.D - NRo) / NRd;
  if (t0 < tmin || t0 > tmax)
    return false;
  return true;
}

/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a sphere and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is one or two such intersections, outputs the value of the
// smaller t, the position of the intersection (hitPos) and the normal
// vector at the intersection (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectSphere(in Sphere_t sph, in Ray_t ray, in float tmin,
                     in float tmax, out float t, out vec3 hitPos,
                     out vec3 hitNormal) {
  ///////////////////////////////////
  // TASK 1: WRITE YOUR CODE HERE. //
  ///////////////////////////////////

  vec3 bcRo = ray.o - sph.center;

  float a = 1.0;
  float b = 2.0 * dot(bcRo, ray.d);
  float c = dot(bcRo, bcRo) - sph.radius * sph.radius;
  float d = b * b - 4.0 * a * c;
  float t0 = -b / (2.0 * a);

  if (d > 0.0) {
    d = sqrt(d) / (2.0 * a);
    float t1 = t0 - d;
    float t2 = t0 + d;

    if (tmin <= t1 && t1 <= tmax) {
      t0 = t1;
    } else if (tmin <= t2 && t2 <= tmax) {
      t0 = t2;
    } else {
      return false;
    }
  } else if (d < 0.0 || t0 < tmin || t0 > tmax) {
    return false;
  }

  // We have a hit -- output results.
  t = t0;
  hitPos = bcRo + t0 * ray.d;
  hitNormal = normalize(hitPos);
  hitPos += sph.center;

  return true;
}

/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a sphere and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectSphere(in Sphere_t sph, in Ray_t ray, in float tmin,
                     in float tmax) {
  ///////////////////////////////////
  // TASK 1: WRITE YOUR CODE HERE. //
  ///////////////////////////////////

  vec3 bcRo = ray.o - sph.center;

  float a = 1.0;
  float b = 2.0 * dot(bcRo, ray.d);
  float c = dot(bcRo, bcRo) - sph.radius * sph.radius;
  float d = b * b - 4.0 * a * c;
  float t0 = -b / (2.0 * a);

  if (d > 0.0) {
    d = sqrt(d) / (2.0 * a);
    float t1 = t0 - d;
    float t2 = t0 + d;

    if (tmin <= t1 && t1 <= tmax) {
      t0 = t1;
    } else if (tmin <= t2 && t2 <= tmax) {
      t0 = t2;
    } else {
      return false;
    }
  } else if (d < 0.0 || t0 < tmin || t0 > tmax) {
    return false;
  }

  return true;
}

/////////////////////////////////////////////////////////////////////////////
// Computes (I_a * k_a) + k_shadow * I_source * [ k_d * (N.L) + k_r * (R.V)^n ].
// Input vectors L, N and V are pointing AWAY from surface point.
// Assume all vectors L, N and V are unit vectors.
/////////////////////////////////////////////////////////////////////////////
vec3 PhongLighting(in vec3 L, in vec3 N, in vec3 V, in bool inShadow,
                   in Material_t mat, in Light_t light) {
  if (inShadow) {
    return light.I_a * mat.k_a;
  } else {
    vec3 R = reflect(-L, N);
    float N_dot_L = max(0.0, dot(N, L));
    float R_dot_V = max(0.0, dot(R, V));
    float R_dot_V_pow_n = (R_dot_V == 0.0) ? 0.0 : pow(R_dot_V, mat.n);

    return light.I_a * mat.k_a +
           light.I_source * (mat.k_d * N_dot_L + mat.k_r * R_dot_V_pow_n);
  }
}

/////////////////////////////////////////////////////////////////////////////
// Casts a ray into the scene and returns color computed at the nearest
// intersection point. The color is the sum of light from all light sources,
// each computed using Phong Lighting Model, with consideration of
// whether the interesection point is being shadowed from the light.
// If there is no interesection, returns the background color, and outputs
// hasHit as false.
// If there is intersection, returns the computed color, and outputs
// hasHit as true, the 3D position of the intersection (hitPos), the
// normal vector at the intersection (hitNormal), and the k_rg value
// of the material of the intersected object.
/////////////////////////////////////////////////////////////////////////////
vec3 CastRay(in Ray_t ray, out bool hasHit, out vec3 hitPos, out vec3 hitNormal,
             out vec3 k_rg) {
  // Find whether and where the ray hits some object.
  // Take the nearest hit point.

  bool hasHitSomething = false;
  float nearest_t =
      DEFAULT_TMAX;       // The ray parameter t at the nearest hit point.
  vec3 nearest_hitPos;    // 3D position of the nearest hit point.
  vec3 nearest_hitNormal; // Normal vector at the nearest hit point.
  int nearest_hitMatID;   // MaterialID of the object at the nearest hit point.

  float temp_t;
  vec3 temp_hitPos;
  vec3 temp_hitNormal;
  bool temp_hasHit;

  ///////////////////////////////////////////////////////////////////////////
  // TASK 1:
  // * Try interesecting input ray with all the planes and spheres,
  //   and record the front-most (nearest) interesection.
  // * If there is interesection, need to record hasHitSomething,
  //   nearest_t, nearest_hitPos, nearest_hitNormal, nearest_hitMatID.
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////
  // TASK 1: WRITE YOUR CODE HERE. //
  ///////////////////////////////////

  for (int i = 0; i < NUM_PLANES; ++i) {
    temp_hasHit = IntersectPlane(Plane[i], ray, DEFAULT_TMIN, DEFAULT_TMAX,
                                 temp_t, temp_hitPos, temp_hitNormal);
    if (temp_hasHit) {
      hasHitSomething = true;
      if (temp_t < nearest_t) {
        nearest_t = temp_t;
        nearest_hitPos = temp_hitPos;
        nearest_hitNormal = temp_hitNormal;
        nearest_hitMatID = Plane[i].materialID;
      }
    }
  }

  for (int i = 0; i < NUM_SPHERES; ++i) {
    temp_hasHit = IntersectSphere(Sphere[i], ray, DEFAULT_TMIN, DEFAULT_TMAX,
                                  temp_t, temp_hitPos, temp_hitNormal);
    if (temp_hasHit) {
      hasHitSomething = true;
      if (temp_t < nearest_t) {
        nearest_t = temp_t;
        nearest_hitPos = temp_hitPos;
        nearest_hitNormal = temp_hitNormal;
        nearest_hitMatID = Sphere[i].materialID;
      }
    }
  }

  // One of the output results.
  hasHit = hasHitSomething;
  if (!hasHitSomething)
    return BACKGROUND_COLOR;

  vec3 I_local = vec3(0.0); // Result color will be accumulated here.

  ///////////////////////////////////////////////////////////////////////////
  // TASK 1:
  // * Accumulate lighting from each light source on the nearest hit point.
  //   They are all accumulated into I_local.
  // * For each light source, make a shadow ray, and check if the shadow ray
  //   intersects any of the objects (the planes and spheres) between the
  //   nearest hit point and the light source.
  // * Then, call PhongLighting() to compute lighting for this light source.
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////
  // TASK 1: WRITE YOUR CODE HERE. //
  ///////////////////////////////////

  // Populate output results.
  hitPos = nearest_hitPos;
  hitNormal = nearest_hitNormal;
  k_rg = Material[nearest_hitMatID].k_rg;

  Ray_t shadow_ray;
  shadow_ray.o = nearest_hitPos;
  vec3 V = normalize(ray.o - nearest_hitPos);
  for (int i = 0; i < NUM_LIGHTS; ++i) {
    vec3 L = normalize(Light[i].position - nearest_hitPos);
    float tmax = distance(Light[i].position, nearest_hitPos);
    tmax = min(tmax, DEFAULT_TMAX);

    shadow_ray.d = L;

    bool inShadow = false;
    for (int j = 0; j < NUM_PLANES && !inShadow; ++j) {
      if (IntersectPlane(Plane[j], shadow_ray, DEFAULT_TMIN, tmax)) {
        inShadow = true;
      }
    }

    for (int j = 0; j < NUM_SPHERES && !inShadow; ++j) {
      if (IntersectSphere(Sphere[j], shadow_ray, DEFAULT_TMIN, tmax)) {
        inShadow = true;
      }
    }

    I_local += PhongLighting(L, nearest_hitNormal, V, inShadow,
                             Material[nearest_hitMatID], Light[i]);
  }

  return I_local;
}

/////////////////////////////////////////////////////////////////////////////
// Execution of fragment shader starts here.
// 1. Initializes the scene.
// 2. Compute a primary ray for the current pixel (fragment).
// 3. Trace ray into the scene with NUM_ITERATIONS recursion levels.
/////////////////////////////////////////////////////////////////////////////
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  // Initialize random number generator before the first call to rand().
  uint RAND_SEED = uint((mod(iTime * 100.0, 100.0) + 101.01) *
                        (fragCoord.x + 17.0) * (fragCoord.y + 23.0));
  rand_z1 = uint(RAND_SEED + 2U);
  rand_z2 = uint(RAND_SEED + 8U);
  rand_z3 = uint(RAND_SEED + 16U);
  rand_z4 = uint(RAND_SEED + 128U);

  InitScene();

  // Camera position and orientation in world space.
  vec3 cam_pos = vec3(2.5, 1.0, 2.5);
  vec3 cam_lookat = vec3(0.25, 1.0, 0.0);
  vec3 cam_up_vec = vec3(0.0, 1.0, 0.0);

  // Camera coordinate frame in world space.
  vec3 cam_z_axis = normalize(cam_pos - cam_lookat);
  vec3 cam_x_axis = normalize(cross(cam_up_vec, cam_z_axis));
  vec3 cam_y_axis = normalize(cross(cam_z_axis, cam_x_axis));

  // Vertical field-of-view angle of camera. In radians.
  float cam_FOVY = FOVY * PI / 180.0;

  // Perpendicular distance of the image rectangle from the camera.
  // If implementing depth-of-field, the plane of the image rectangle
  // is the plane of focus.
  float image_dist = distance(cam_pos, vec3(0.0, 0.7, 0.0));

  float image_height = 2.0 * image_dist * tan(cam_FOVY / 2.0);
  float image_width = image_height * iResolution.x / iResolution.y;
  float pixel_width = image_width / iResolution.x;

  // Image rectangle origin (bottom-leftmost corner) position in camera space.
  vec3 image_origin =
      vec3(-image_width / 2.0, -image_height / 2.0, -image_dist);

  /////////////////////////////////////////////////////////////////////////
  // TASK 2:
  // * Trace multiple (SPP) random primary rays per pixel to produce
  //   depth-of-field effect and for image anti-aliasing (reduce jaggies).
  // * Each primary ray starts from a random position on the lens and
  //   points towards a random position inside the current pixel.
  // * The lens is assumed to have a square-shaped aperture of size
  //   aperture_width x aperture_width. The lens is centered at the
  //   the origin of the camera frame, and parallel to the x-y plane
  //   of the camera frame.
  // * The final color of the current pixel is the average color of
  //   all the primary rays.
  /////////////////////////////////////////////////////////////////////////

  //=======================================================================
  // These constants are used for distribution ray tracing to produce
  // depth-of-field effect and for image anti-aliasing (reduce jaggies).
  //=======================================================================
  // Number of samples (random primary rays) per pixel.
  const int SPP = 32;

  // Lens aperture width. Assume square aperture.
  const float aperture_width = 0.3;
  //=======================================================================

  ////////////////////////////////////
  // TASK 2: MODIFY THE CODE BELOW. //
  ////////////////////////////////////

  vec3 averageColor = vec3(0.0);
  for (int i = 0; i < SPP; ++i) {
    vec2 rectPoint = vec2(rand() - 0.5, rand() - 0.5) + fragCoord;

    // Current pixel 3D position in camera space.
    vec3 pixel_pos = image_origin + vec3(pixel_width * rectPoint, 0);

    // Create primary ray.
    Ray_t pRay;
    pRay.o = cam_pos + aperture_width * (cam_x_axis * (rand() - 0.5) +
                                         cam_y_axis * (rand() - 0.5));
    pRay.d =
        normalize(cam_pos + pixel_pos.x * cam_x_axis +
                  pixel_pos.y * cam_y_axis + pixel_pos.z * cam_z_axis - pRay.o);

    // Start Ray Tracing.
    // Use iterations to emulate the recursion.

    vec3 I_result = vec3(0.0);
    vec3 compounded_k_rg = vec3(1.0);
    Ray_t nextRay = pRay;

    for (int level = 0; level <= NUM_ITERATIONS; level++) {
      bool hasHit;
      vec3 hitPos, hitNormal, k_rg;

      vec3 I_local = CastRay(nextRay, hasHit, hitPos, hitNormal, k_rg);

      I_result += compounded_k_rg * I_local;

      if (!hasHit)
        break;

      compounded_k_rg *= k_rg;

      nextRay = Ray_t(hitPos, normalize(reflect(nextRay.d, hitNormal)));
    }

    averageColor += I_result / float(SPP);
  }

  fragColor = vec4(averageColor, 1.0);
}
