#pragma once
#ifndef __SHADER_CREATOR_HPP__
#define __SHADER_CREATOR_HPP__

#include "Shader.hpp"
#include "Conductor.hpp"
#include "Dielectric.hpp"
#include "Lambertian.hpp"

namespace PhotonMapping
{
    class ShaderCreator
    {
      public:
        ShaderCreator() = default;
        SharedShader create(Material& material, vector<Texture>& t)
        {
            SharedShader shader{nullptr};
            switch (material.type)
            {
                case 0: shader = make_shared<Lambertian>(material, t); break;
                case 2: shader = make_shared<Dielectric>(material, t); break;
                case 3: shader = make_shared<Conductor>(material, t); break;
                default: shader = make_shared<Lambertian>(material, t); break;
            }
            return shader;
        }
    };
}  // namespace PhotonMapping

#endif