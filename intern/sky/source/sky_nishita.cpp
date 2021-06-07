/*
 * Copyright 2011-2020 Blender Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "sky_float3.h"
#include "sky_model.h"

#include <cstdint>
#include <cstring>

/* Constants */
static const float rayleigh_scale = 8e3f;        // Rayleigh scale height (m)
static const float mie_scale = 1.2e3f;           // Mie scale height (m)
static const float mie_coeff = 2e-5f;            // Mie scattering coefficient (m^-1)
static const float mie_G = 0.76f;                // aerosols anisotropy
static const float sqr_G = mie_G * mie_G;        // squared aerosols anisotropy
static const float earth_radius = 6360e3f;       // radius of Earth (m)
static const float atmosphere_radius = 6420e3f;  // radius of atmosphere (m)
static const int steps = 32;                     // segments of primary ray
static const int steps_light = 16;               // segments of sun connection ray
static const int samples = 10;                   // samples for indirect scattering
static const int bounce_depth = 5;               // recursion depth of indirect scattering
static const int num_wavelengths = 21;           // number of wavelengths
static const int min_wavelength = 380;           // lowest sampled wavelength (nm)
static const int max_wavelength = 780;           // highest sampled wavelength (nm)
// step between each sampled wavelength (nm)
static const float step_lambda = (max_wavelength - min_wavelength) / (num_wavelengths - 1);
/* Sun irradiance on top of the atmosphere (W*m^-2*nm^-1) */
static const float irradiance[] = {
    1.45756829855592995315f, 1.56596305559738380175f, 1.65148449067670455293f,
    1.71496242737209314555f, 1.75797983805020541226f, 1.78256407885924539336f,
    1.79095108475838560302f, 1.78541550133410664714f, 1.76815554864306845317f,
    1.74122069647250410362f, 1.70647127164943679389f, 1.66556087452739887134f,
    1.61993437242451854274f, 1.57083597368892080581f, 1.51932335059305478886f,
    1.46628494965214395407f, 1.41245852740172450623f, 1.35844961970384092709f,
    1.30474913844739281998f, 1.25174963272610817455f, 1.19975998755420620867f};
/* Rayleigh scattering coefficient (m^-1) */
static const float rayleigh_coeff[] = {
    0.00005424820087636473f, 0.00004418549866505454f, 0.00003635151910165377f,
    0.00003017929012024763f, 0.00002526320226989157f, 0.00002130859310621843f,
    0.00001809838025320633f, 0.00001547057129129042f, 0.00001330284977336850f,
    0.00001150184784075764f, 0.00000999557429990163f, 0.00000872799973630707f,
    0.00000765513700977967f, 0.00000674217203751443f, 0.00000596134125832052f,
    0.00000529034598065810f, 0.00000471115687557433f, 0.00000420910481110487f,
    0.00000377218381260133f, 0.00000339051255477280f, 0.00000305591531679811f};
/* Ozone absorption coefficient (m^-1) */
static const float ozone_coeff[] = {
    0.00000000325126849861f, 0.00000000585395365047f, 0.00000001977191155085f,
    0.00000007309568762914f, 0.00000020084561514287f, 0.00000040383958096161f,
    0.00000063551335912363f, 0.00000096707041180970f, 0.00000154797400424410f,
    0.00000209038647223331f, 0.00000246128056164565f, 0.00000273551299461512f,
    0.00000215125863128643f, 0.00000159051840791988f, 0.00000112356197979857f,
    0.00000073527551487574f, 0.00000046450130357806f, 0.00000033096079921048f,
    0.00000022512612292678f, 0.00000014879129266490f, 0.00000016828623364192f};
/* CIE XYZ color matching functions */
static const float cmf_xyz[][3] = {{0.00136800000f, 0.00003900000f, 0.00645000100f},
                                   {0.01431000000f, 0.00039600000f, 0.06785001000f},
                                   {0.13438000000f, 0.00400000000f, 0.64560000000f},
                                   {0.34828000000f, 0.02300000000f, 1.74706000000f},
                                   {0.29080000000f, 0.06000000000f, 1.66920000000f},
                                   {0.09564000000f, 0.13902000000f, 0.81295010000f},
                                   {0.00490000000f, 0.32300000000f, 0.27200000000f},
                                   {0.06327000000f, 0.71000000000f, 0.07824999000f},
                                   {0.29040000000f, 0.95400000000f, 0.02030000000f},
                                   {0.59450000000f, 0.99500000000f, 0.00390000000f},
                                   {0.91630000000f, 0.87000000000f, 0.00165000100f},
                                   {1.06220000000f, 0.63100000000f, 0.00080000000f},
                                   {0.85444990000f, 0.38100000000f, 0.00019000000f},
                                   {0.44790000000f, 0.17500000000f, 0.00002000000f},
                                   {0.16490000000f, 0.06100000000f, 0.00000000000f},
                                   {0.04677000000f, 0.01700000000f, 0.00000000000f},
                                   {0.01135916000f, 0.00410200000f, 0.00000000000f},
                                   {0.00289932700f, 0.00104700000f, 0.00000000000f},
                                   {0.00069007860f, 0.00024920000f, 0.00000000000f},
                                   {0.00016615050f, 0.00006000000f, 0.00000000000f},
                                   {0.00004150994f, 0.00001499000f, 0.00000000000f}};

static float3 geographical_to_direction(float lat, float lon)
{
  return make_float3(cosf(lat) * cosf(lon), cosf(lat) * sinf(lon), sinf(lat));
}

static float3 spec_to_xyz(const float *spectrum)
{
  float3 xyz = make_float3(0.0f, 0.0f, 0.0f);
  for (int i = 0; i < num_wavelengths; i++) {
    xyz.x += cmf_xyz[i][0] * spectrum[i];
    xyz.y += cmf_xyz[i][1] * spectrum[i];
    xyz.z += cmf_xyz[i][2] * spectrum[i];
  }
  return xyz * step_lambda;
}

/* Atmosphere volume models */
static float density_rayleigh(float height)
{
  return expf(-height / rayleigh_scale);
}

static float density_mie(float height)
{
  return expf(-height / mie_scale);
}

static float density_ozone(float height)
{
  float den = 0.0f;
  if (height >= 10000.0f && height < 25000.0f) {
    den = 1.0f / 15000.0f * height - 2.0f / 3.0f;
  }
  else if (height >= 25000 && height < 40000) {
    den = -(1.0f / 15000.0f * height - 8.0f / 3.0f);
  }
  return den;
}

static float phase_rayleigh(float mu)
{
  return 3.0f / (16.0f * M_PI_F) * (1.0f + sqr(mu));
}

static float phase_mie(float mu)
{
  return (3.0f * (1.0f - sqr_G) * (1.0f + sqr(mu))) /
         (8.0f * M_PI_F * (2.0f + sqr_G) * powf((1.0f + sqr_G - 2.0f * mie_G * mu), 1.5));
}

/* Old intersection helpers */
static bool old_surface_intersection(float3 pos, float3 dir)
{
  if (dir.z >= 0) {
    return false;
  }
  float b = -2.0f * dot(dir, -pos);
  float c = len_squared(pos) - sqr(earth_radius);
  float t = b * b - 4.0f * c;
  if (t >= 0.0f) {
    return true;
  }
  return false;
}

static float3 old_atmosphere_intersection(float3 pos, float3 dir)
{
  float b = -2.0f * dot(dir, -pos);
  float c = len_squared(pos) - sqr(atmosphere_radius);
  float t = (-b + sqrtf(b * b - 4.0f * c)) / 2.0f;
  return make_float3(pos.x + dir.x * t, pos.y + dir.y * t, pos.z + dir.z * t);
}

static float3 old_ray_optical_depth(float3 ray_origin, float3 ray_dir)
{
  /* this code computes the optical depth along a ray through the atmosphere */
  float3 ray_end = old_atmosphere_intersection(ray_origin, ray_dir);
  float ray_length = distance(ray_origin, ray_end);

  /* to compute the optical depth, we step along the ray in segments and
   * accumulate the optical depth along each segment */
  float segment_length = ray_length / steps_light;
  float3 segment = segment_length * ray_dir;

  /* instead of tracking the transmission spectrum across all wavelengths directly,
   * we use the fact that the density always has the same spectrum for each type of
   * scattering, so we split the density into a constant spectrum and a factor and
   * only track the factors */
  float3 optical_depth = make_float3(0.0f, 0.0f, 0.0f);

  /* the density of each segment is evaluated at its middle */
  float3 P = ray_origin + 0.5f * segment;

  for (int i = 0; i < steps_light; i++) {
    /* height above sea level */
    float height = len(P) - earth_radius;

    /* accumulate optical depth of this segment (density is assumed to be constant along it) */
    float3 density = make_float3(
        density_rayleigh(height), density_mie(height), density_ozone(height));
    optical_depth += density;

    /* advance along ray */
    P += segment;
  }

  return optical_depth * segment_length;
}

/* Intersection helpers */
static float surface_intersection(float3 pos, float3 dir)
{
  float u = dot(dir, pos);
  if (u > 0.0f)
    return 1e10f;  // Ray points away from sphere

  float v = len_squared(pos) - sqr(earth_radius);
  float d = u * u - v;

  if (d < 0.0f) {
    // Miss
    return 1e10f;
  }
  else {
    // Assume ray doesn't start inside earth => smaller of the two solutions is nearest hit.
    return -u - sqrtf(d);
  }
}

static float atmosphere_intersection(float3 pos, float3 dir)
{
  float u = dot(dir, pos);
  float v = len_squared(pos) - sqr(atmosphere_radius);
  float d = u * u - v;
  // Assume that the ray starts inside atmosphere => 2 guaranteed hits,
  // smaller one will be behind origin so use larger one.
  return -u + sqrtf(d);
}

static float3 ray_optical_depth(float3 ray_origin, float3 ray_dir, float ray_length)
{
  // Instead of ray marching, this code uses Gauss-Laguerre quadrature with four nodes.
  // This is faster and gives more precise results.
  float rescale = 0.1f;

  float nodes[4] = {0.322548f, 1.74576f, 4.53662f, 9.39507f};
  float weights[4] = {0.603154f, 0.357419f, 0.0388879f, 0.000539295f};

  float3 t_to_pos = ray_dir * ray_length * rescale;

  float3 v = make_float3(0.0f, 0.0f, 0.0f);
  for (int i = 0; i < 4; i++) {
    float t = nodes[i];
    float3 P = ray_origin + t_to_pos * t;
    float height = len(P) - earth_radius;
    float3 density = make_float3(
        density_rayleigh(height), density_mie(height), density_ozone(height));
    v += density * weights[i] * expf(t);
  }

  v = v * (rescale * ray_length);
  return v;
}

// RNG helpers
static inline uint32_t rotl(const uint32_t x, int k)
{
  return (x << k) | (x >> (32 - k));
}

static float rng_get(uint32_t *rng)
{
  const uint32_t result = rng[0] + rng[3];
  const uint32_t t = rng[1] << 9;

  rng[2] ^= rng[0];
  rng[3] ^= rng[1];
  rng[1] ^= rng[2];
  rng[0] ^= rng[3];
  rng[2] ^= t;
  rng[3] = rotl(rng[3], 11);

  return float(result) / float(0xffffffff);
}

static void rng_seed(uint32_t *rng, uint32_t seed1, uint32_t seed2)
{
  uint64_t v[2];
  v[0] = ((uint64_t)seed1) << 32 | seed2;

  uint64_t z = (v[0] += 0x9e3779b97f4a7c15);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  v[0] = z ^ (z >> 31);

  z = (v[0] += 0x9e3779b97f4a7c15);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  v[1] = z ^ (z >> 31);

  memcpy(rng, v, sizeof(uint32_t) * 4);
}

static float3 sample_uniform_sphere(float u1, float u2)
{
  float z = 1.0f - 2.0f * u1;
  float r = sqrtf(fmaxf(0.0f, 1.0f - z * z));
  float phi = M_2PI_F * u2;
  float x = r * cosf(phi);
  float y = r * sinf(phi);

  return make_float3(x, y, z);
}

static void single_scattering(float3 ray_dir,
                              float3 sun_dir,
                              float3 ray_origin,
                              float air_density,
                              float dust_density,
                              float ozone_density,
                              float *r_spectrum)
{
  /* this code computes single-inscattering along a ray through the atmosphere */
  float3 ray_end = old_atmosphere_intersection(ray_origin, ray_dir);
  float ray_length = distance(ray_origin, ray_end);

  /* to compute the inscattering, we step along the ray in segments and accumulate
   * the inscattering as well as the optical depth along each segment */
  float segment_length = ray_length / steps;
  float3 segment = segment_length * ray_dir;

  /* instead of tracking the transmission spectrum across all wavelengths directly,
   * we use the fact that the density always has the same spectrum for each type of
   * scattering, so we split the density into a constant spectrum and a factor and
   * only track the factors */
  float3 optical_depth = make_float3(0.0f, 0.0f, 0.0f);

  /* zero out light accumulation */
  for (int wl = 0; wl < num_wavelengths; wl++) {
    r_spectrum[wl] = 0.0f;
  }

  /* phase function for scattering and the density scale factor */
  float mu = dot(ray_dir, sun_dir);
  float3 phase_function = make_float3(phase_rayleigh(mu), phase_mie(mu), 0.0f);
  float3 density_scale = make_float3(air_density, dust_density, ozone_density);

  /* the density and in-scattering of each segment is evaluated at its middle */
  float3 P = ray_origin + 0.5f * segment;

  for (int i = 0; i < steps; i++) {
    /* height above sea level */
    float height = len(P) - earth_radius;

    /* evaluate and accumulate optical depth along the ray */
    float3 density = density_scale * make_float3(density_rayleigh(height),
                                                 density_mie(height),
                                                 density_ozone(height));
    optical_depth += segment_length * density;

    /* if the Earth isn't in the way, evaluate inscattering from the sun */
    if (!old_surface_intersection(P, sun_dir)) {
      float3 light_optical_depth = density_scale * old_ray_optical_depth(P, sun_dir);
      float3 total_optical_depth = optical_depth + light_optical_depth;

      /* attenuation of light */
      for (int wl = 0; wl < num_wavelengths; wl++) {
        float3 extinction_density = total_optical_depth * make_float3(rayleigh_coeff[wl],
                                                                      1.11f * mie_coeff,
                                                                      ozone_coeff[wl]);
        float attenuation = expf(-reduce_add(extinction_density));

        float3 scattering_density = density * make_float3(rayleigh_coeff[wl], mie_coeff, 0.0f);

        /* the total inscattered radiance from one segment is:
         * Tr(A<->B) * Tr(B<->C) * sigma_s * phase * L * segment_length
         *
         * These terms are:
         * Tr(A<->B): Transmission from start to scattering position (tracked in optical_depth)
         * Tr(B<->C): Transmission from scattering position to light (computed in
         * old_ray_optical_depth) sigma_s: Scattering density phase: Phase function of the
         * scattering type (Rayleigh or Mie) L: Radiance coming from the light source
         * segment_length: The length of the segment
         *
         * The code here is just that, with a bit of additional optimization to not store full
         * spectra for the optical depth
         */
        r_spectrum[wl] += attenuation * reduce_add(phase_function * scattering_density) *
                          irradiance[wl] * segment_length;
      }
    }

    /* advance along ray */
    P += segment;
  }
}

static void multiple_scattering(float3 ray_dir,
                                float3 sun_dir,
                                float3 ray_origin,
                                float air_density,
                                float dust_density,
                                float ozone_density,
                                float *r_spectrum,
                                uint32_t *rng,
                                int depth)
{
  float ray_length = std::min(atmosphere_intersection(ray_origin, ray_dir),
                              surface_intersection(ray_origin, ray_dir));

  if (ray_length < 1e-2f)
    return;

  float segment_length = ray_length / steps;
  float3 segment = segment_length * ray_dir;

  float3 optical_depth = make_float3(0.0f, 0.0f, 0.0f);

  float mu = dot(ray_dir, sun_dir);
  float3 phase_function = make_float3(phase_rayleigh(mu), phase_mie(mu), 0.0f);
  float3 density_scale = make_float3(air_density, dust_density, ozone_density);

  float3 P = ray_origin + 0.5f * segment;  // Note: We evaluate in middle, but use full segment OD?

  // Note: We split the in-scattered light into light coming directly from the sun
  // and light coming from a second (and maybe third etc.) scattering event.
  // Direct light from the sun is evaluated with ray marching, while indirect light
  // is evaluated using Monte-Carlo integration.
  for (int i = 0; i < steps; i++) {
    float height = len(P) - earth_radius;
    float3 density = density_scale * make_float3(density_rayleigh(height),
                                                 density_mie(height),
                                                 density_ozone(height));
    optical_depth += segment_length * density;

    // Integrate single scattering
    if (surface_intersection(P, sun_dir) > 1e9f) {
      float sun_distance = atmosphere_intersection(P, sun_dir);
      float3 light_optical_depth = density_scale * ray_optical_depth(P, sun_dir, sun_distance);
      float3 total_optical_depth = optical_depth + light_optical_depth;

      for (int wl = 0; wl < num_wavelengths; wl++) {
        float3 extinction_density = total_optical_depth * make_float3(rayleigh_coeff[wl],
                                                                      1.11f * mie_coeff,
                                                                      ozone_coeff[wl]);
        float attenuation = expf(-reduce_add(extinction_density));

        float3 scattering_density = density * make_float3(rayleigh_coeff[wl], mie_coeff, 0.0f);

        r_spectrum[wl] += attenuation * reduce_add(phase_function * scattering_density) *
                          irradiance[wl] * segment_length;
      }
    }

    P += segment;
  }

  // Decide how many samples to use for evaluating indirectly scattered light.
  int multiScatterSamples;
  if (depth == 0) {
    // If this ray is coming from the camera, take many samples.
    multiScatterSamples = samples;
  }
  else if (depth < bounce_depth) {
    // This is already a scattered ray, so only take one sample to avoid exponential growth.
    multiScatterSamples = 1;
  }
  else {
    // Stop iteration.
    multiScatterSamples = 0;
  }
  for (int sample = 0; sample < multiScatterSamples; sample++) {
    // Step 1: Sample random position along ray
    // The logic used here is just a guess from my side - it's designed so that the
    // sample density falls off exponentially along the ray. Since throughput decreases
    // exponentially, this places most effort in areas that end up contributing a lot.
#if 1
    /* pdf(t) = c*exp(-4t/d), where d is the ray length.
     * The probability pdf(t) is highest at the start and falls off to exp(-4) (~0.01) at the end.
     * cdf(t) = c*d/4 * (1 - exp(-4t/d))
     * t(zeta) = d/4 * log(c*d / (c*d - 4zeta))
     * Normalization: c = 4exp(4)/((exp(4)-1) * d)
     * Obvious optimization: Cancel d in c and t(zeta) */
    const float cd = ((4 * expf(4)) / (expf(4) - 1));
    float localT = 0.25f * logf(cd / (cd - 4 * rng_get(rng)));
    float scatterT = ray_length * localT;
    float pdfT = (cd / ray_length) * expf(-4 * localT);
#else  // Alternative: Just pick uniformly along the path. Gives much more noise.
    float scatterT = rng_get(rng) * ray_length;
    float pdfT = 1.0f / ray_length;
#endif
    float3 scatterP = ray_origin + scatterT * ray_dir;

    // Step 2: Sample random scattered direction
    // This code currently samples uniformly on the sphere.
    // Ideally, this should sample proportional to the Rayleigh or Mie phase function
    // so that phase function and pdfD cancel out for lower noise.
    // The two can be combined with MIS.
    float3 scatterD = sample_uniform_sphere(rng_get(rng), rng_get(rng));
    float mu = dot(ray_dir, scatterD);
    float3 phase_function = make_float3(phase_rayleigh(mu), phase_mie(mu), 0.0f);
    float pdfD = 0.25f / M_PI_F;

    // Step 3: Evaluate inscattering
    float scatterSpectrum[num_wavelengths] = {0};
    multiple_scattering(scatterD,
                        sun_dir,
                        scatterP,
                        air_density,
                        dust_density,
                        ozone_density,
                        scatterSpectrum,
                        rng,
                        depth + 1);

    // Step 4: Add contribution to spectrum
    float3 optical_density = ray_optical_depth(ray_origin, ray_dir, scatterT);
    float height = len(scatterP) - earth_radius;
    float3 density = density_scale * make_float3(density_rayleigh(height),
                                                 density_mie(height),
                                                 density_ozone(height));
    for (int wl = 0; wl < num_wavelengths; wl++) {
      float3 extinction_density = optical_density * make_float3(rayleigh_coeff[wl],
                                                                1.11f * mie_coeff,
                                                                ozone_coeff[wl]);
      float attenuation = expf(-reduce_add(extinction_density));
      float3 scattering_density = density * make_float3(rayleigh_coeff[wl], mie_coeff, 0.0f);

      r_spectrum[wl] += attenuation * reduce_add(phase_function * scattering_density) *
                        scatterSpectrum[wl] / (pdfD * pdfT * multiScatterSamples);
    }
  }
}

void SKY_nishita_skymodel_precompute_texture(float *pixels,
                                             int stride,
                                             int start_y,
                                             int end_y,
                                             int width,
                                             int height,
                                             float sun_elevation,
                                             float altitude,
                                             float air_density,
                                             float dust_density,
                                             float ozone_density,
                                             bool multi_scattering)
{
  /* calculate texture pixels */
  int half_width = width / 2;
  float3 cam_pos = make_float3(0, 0, earth_radius + altitude);
  float3 sun_dir = geographical_to_direction(sun_elevation, 0.0f);

  float latitude_step = M_PI_2_F / height;
  float longitude_step = M_2PI_F / width;
  float half_lat_step = latitude_step / 2.0f;

  // Todo: Use better RNG (e.g. LD sequence)
  uint32_t rng[4];
  rng_seed(rng, start_y, width * 0x9a43b471 + height);

  for (int y = start_y; y < end_y; y++) {
    /* sample more pixels toward the horizon */
    float latitude = (M_PI_2_F + half_lat_step) * sqr((float)y / height);

    float *pixel_row = pixels + (y * width * stride);
    for (int x = 0; x < half_width; x++) {
      float longitude = longitude_step * x - M_PI_F;

      float3 dir = geographical_to_direction(latitude, longitude);
      float spectrum[num_wavelengths] = {0};
      if (multi_scattering)
        multiple_scattering(
            dir, sun_dir, cam_pos, air_density, dust_density, ozone_density, spectrum, rng, 0);
      else
        single_scattering(
            dir, sun_dir, cam_pos, air_density, dust_density, ozone_density, spectrum);
      float3 xyz = spec_to_xyz(spectrum);

      /* store pixels */
      int pos_x = x * stride;
      pixel_row[pos_x] = xyz.x;
      pixel_row[pos_x + 1] = xyz.y;
      pixel_row[pos_x + 2] = xyz.z;
      /* mirror sky */
      int mirror_x = (width - x - 1) * stride;
      pixel_row[mirror_x] = xyz.x;
      pixel_row[mirror_x + 1] = xyz.y;
      pixel_row[mirror_x + 2] = xyz.z;
    }
  }
}

/*********** Sun ***********/
static void sun_radiation(float3 cam_dir,
                          float altitude,
                          float air_density,
                          float dust_density,
                          float solid_angle,
                          float *r_spectrum)
{
  float3 cam_pos = make_float3(0, 0, earth_radius + altitude);
  float sun_distance = atmosphere_intersection(cam_pos, cam_dir);
  float3 optical_depth = ray_optical_depth(cam_pos, cam_dir, sun_distance);

  /* compute final spectrum */
  for (int i = 0; i < num_wavelengths; i++) {
    /* combine spectra and the optical depth into transmittance */
    float transmittance = rayleigh_coeff[i] * optical_depth.x * air_density +
                          1.11f * mie_coeff * optical_depth.y * dust_density;
    r_spectrum[i] = irradiance[i] * expf(-transmittance) / solid_angle;
  }
}

void SKY_nishita_skymodel_precompute_sun(float sun_elevation,
                                         float angular_diameter,
                                         float altitude,
                                         float air_density,
                                         float dust_density,
                                         float *r_pixel_bottom,
                                         float *r_pixel_top)
{
  /* definitions */
  float half_angular = angular_diameter / 2.0f;
  float solid_angle = M_2PI_F * (1.0f - cosf(half_angular));
  float spectrum[num_wavelengths];
  float bottom = sun_elevation - half_angular;
  float top = sun_elevation + half_angular;
  float elevation_bottom, elevation_top;
  float3 pix_bottom, pix_top, sun_dir;

  /* compute 2 pixels for sun disc */
  elevation_bottom = (bottom > 0.0f) ? bottom : 0.0f;
  elevation_top = (top > 0.0f) ? top : 0.0f;
  sun_dir = geographical_to_direction(elevation_bottom, 0.0f);
  sun_radiation(sun_dir, altitude, air_density, dust_density, solid_angle, spectrum);
  pix_bottom = spec_to_xyz(spectrum);
  sun_dir = geographical_to_direction(elevation_top, 0.0f);
  sun_radiation(sun_dir, altitude, air_density, dust_density, solid_angle, spectrum);
  pix_top = spec_to_xyz(spectrum);

  /* store pixels */
  r_pixel_bottom[0] = pix_bottom.x;
  r_pixel_bottom[1] = pix_bottom.y;
  r_pixel_bottom[2] = pix_bottom.z;
  r_pixel_top[0] = pix_top.x;
  r_pixel_top[1] = pix_top.y;
  r_pixel_top[2] = pix_top.z;
}
