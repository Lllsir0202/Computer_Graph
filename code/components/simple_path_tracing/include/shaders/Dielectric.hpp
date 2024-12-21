#pragma once
#ifndef __DIELECTRIC_HPP__
#define __DIELECTRIC_HPP__

#include "Shader.hpp"

namespace SimplePathTracer
{
    class Dielectric : public Shader
    {
      private:
        Vec3  absorbed;
        float ior;

      public:
        Dielectric(Material& material, vector<Texture>& textures);
        Scattered shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal) const;
    };
}  // namespace SimplePathTracer

#endif