#pragma once
#ifndef __CONDUCTOR_HPP__
#define __CONDUCTOR_HPP__

#include "Shader.hpp"

namespace PhotonMapping
{
    class Conductor : public Shader
    {
      private:
        Vec3 reflect;

      public:
        Conductor(Material& material, vector<Texture>& textures);
        Scattered shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal) const;
    };
}  // namespace PhotonMapping

#endif