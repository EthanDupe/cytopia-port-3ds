#ifndef CYTOPIA_HEADER_ONLY_NOISE_H
#define CYTOPIA_HEADER_ONLY_NOISE_H

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace noise
{
namespace module
{

class Module
{
public:
  virtual ~Module() = default;
  virtual double GetValue(double x, double y, double z) const = 0;
};

namespace detail
{
inline double fade(double t) { return t * t * t * (t * (t * 6.0 - 15.0) + 10.0); }
inline double lerp(double t, double a, double b) { return a + t * (b - a); }

inline double grad(int hash, double x, double y, double z)
{
  int h = hash & 15;
  double u = h < 8 ? x : y;
  double v = h < 4 ? y : (h == 12 || h == 14 ? x : z);
  return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
}

inline int hashCoords(int x, int y, int z, int seed)
{
  uint32_t h = static_cast<uint32_t>(seed) * 374761393U;
  h += static_cast<uint32_t>(x) * 3266489917U;
  h = (h << 17) | (h >> 15);
  h += static_cast<uint32_t>(y) * 668265263U;
  h = (h << 15) | (h >> 17);
  h += static_cast<uint32_t>(z) * 2246822519U;
  h ^= h >> 15;
  h *= 2246822519U;
  h ^= h >> 13;
  h *= 3266489917U;
  h ^= h >> 16;
  return static_cast<int>(h & 255U);
}

inline double gradientCoherentNoise3D(double x, double y, double z, int seed)
{
  int x0 = static_cast<int>(std::floor(x));
  int y0 = static_cast<int>(std::floor(y));
  int z0 = static_cast<int>(std::floor(z));
  int x1 = x0 + 1;
  int y1 = y0 + 1;
  int z1 = z0 + 1;

  double xs = x - static_cast<double>(x0);
  double ys = y - static_cast<double>(y0);
  double zs = z - static_cast<double>(z0);

  double u = fade(xs);
  double v = fade(ys);
  double w = fade(zs);

  double n000 = grad(hashCoords(x0, y0, z0, seed), xs, ys, zs);
  double n100 = grad(hashCoords(x1, y0, z0, seed), xs - 1.0, ys, zs);
  double n010 = grad(hashCoords(x0, y1, z0, seed), xs, ys - 1.0, zs);
  double n110 = grad(hashCoords(x1, y1, z0, seed), xs - 1.0, ys - 1.0, zs);
  double n001 = grad(hashCoords(x0, y0, z1, seed), xs, ys, zs - 1.0);
  double n101 = grad(hashCoords(x1, y0, z1, seed), xs - 1.0, ys, zs - 1.0);
  double n011 = grad(hashCoords(x0, y1, z1, seed), xs, ys - 1.0, zs - 1.0);
  double n111 = grad(hashCoords(x1, y1, z1, seed), xs - 1.0, ys - 1.0, zs - 1.0);

  double nx00 = lerp(u, n000, n100);
  double nx10 = lerp(u, n010, n110);
  double nx01 = lerp(u, n001, n101);
  double nx11 = lerp(u, n011, n111);

  double nxy0 = lerp(v, nx00, nx10);
  double nxy1 = lerp(v, nx01, nx11);

  return lerp(w, nxy0, nxy1);
}
} // namespace detail

class Perlin : public Module
{
public:
  void SetSeed(int seed) { m_seed = seed; }
  void SetFrequency(double frequency) { m_frequency = frequency; }
  void SetLacunarity(double lacunarity) { m_lacunarity = lacunarity; }
  void SetOctaveCount(int octaveCount) { m_octaveCount = std::clamp(octaveCount, 1, 30); }
  void SetPersistence(double persistence) { m_persistence = persistence; }

  double GetValue(double x, double y, double z) const override
  {
    double value = 0.0;
    double curPersistence = 1.0;
    x *= m_frequency;
    y *= m_frequency;
    z *= m_frequency;

    for (int curOctave = 0; curOctave < m_octaveCount; ++curOctave)
    {
      int seed = (m_seed + curOctave) & 0x7fffffff;
      double signal = detail::gradientCoherentNoise3D(x, y, z, seed);
      value += signal * curPersistence;
      x *= m_lacunarity;
      y *= m_lacunarity;
      z *= m_lacunarity;
      curPersistence *= m_persistence;
    }
    return value;
  }

private:
  int m_seed = 0;
  double m_frequency = 1.0;
  double m_lacunarity = 2.0;
  int m_octaveCount = 6;
  double m_persistence = 0.5;
};

class RidgedMulti : public Module
{
public:
  void SetSeed(int seed) { m_seed = seed; }
  void SetFrequency(double frequency) { m_frequency = frequency; }
  void SetLacunarity(double lacunarity) { m_lacunarity = lacunarity; }
  void SetOctaveCount(int octaveCount) { m_octaveCount = std::clamp(octaveCount, 1, 30); }

  double GetValue(double x, double y, double z) const override
  {
    x *= m_frequency;
    y *= m_frequency;
    z *= m_frequency;

    double value = 0.0;
    double weight = 1.0;
    double offset = 1.0;
    double gain = 2.0;
    double spectralWeight = 1.0;

    for (int curOctave = 0; curOctave < m_octaveCount; ++curOctave)
    {
      int seed = (m_seed + curOctave) & 0x7fffffff;
      double signal = detail::gradientCoherentNoise3D(x, y, z, seed);
      signal = std::abs(signal);
      signal = offset - signal;
      signal *= signal;
      signal *= weight;
      weight = std::clamp(signal * gain, 0.0, 1.0);
      value += (signal * spectralWeight);
      x *= m_lacunarity;
      y *= m_lacunarity;
      z *= m_lacunarity;
      spectralWeight *= std::pow(m_lacunarity, -1.0);
    }
    return (value * 1.25) - 1.0;
  }

private:
  int m_seed = 0;
  double m_frequency = 1.0;
  double m_lacunarity = 2.0;
  int m_octaveCount = 6;
};

class ScaleBias : public Module
{
public:
  void SetSourceModule(int /*index*/, const Module &source) { m_source = &source; }
  void SetScale(double scale) { m_scale = scale; }
  void SetBias(double bias) { m_bias = bias; }

  double GetValue(double x, double y, double z) const override
  {
    return m_source ? (m_source->GetValue(x, y, z) * m_scale + m_bias) : m_bias;
  }

private:
  const Module *m_source = nullptr;
  double m_scale = 1.0;
  double m_bias = 0.0;
};

class Clamp : public Module
{
public:
  void SetSourceModule(int /*index*/, const Module &source) { m_source = &source; }
  void SetBounds(double lower, double upper)
  {
    m_lower = lower;
    m_upper = upper;
  }

  double GetValue(double x, double y, double z) const override
  {
    return m_source ? std::clamp(m_source->GetValue(x, y, z), m_lower, m_upper) : m_lower;
  }

private:
  const Module *m_source = nullptr;
  double m_lower = -1.0;
  double m_upper = 1.0;
};

class Blend : public Module
{
public:
  void SetSourceModule(int index, const Module &source)
  {
    if (index == 0)
      m_source0 = &source;
    else if (index == 1)
      m_source1 = &source;
  }
  void SetControlModule(const Module &control) { m_control = &control; }

  double GetValue(double x, double y, double z) const override
  {
    if (!m_source0 || !m_source1 || !m_control)
      return 0.0;
    double v0 = m_source0->GetValue(x, y, z);
    double v1 = m_source1->GetValue(x, y, z);
    double alpha = std::clamp(m_control->GetValue(x, y, z), 0.0, 1.0);
    return detail::lerp(alpha, v0, v1);
  }

private:
  const Module *m_source0 = nullptr;
  const Module *m_source1 = nullptr;
  const Module *m_control = nullptr;
};

} // namespace module
} // namespace noise

#endif // CYTOPIA_HEADER_ONLY_NOISE_H
