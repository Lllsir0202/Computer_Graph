#include "shaders/Conductor.hpp"
#include "samplers/SamplerInstance.hpp"

#include "Onb.hpp"

namespace SimplePathTracer
{
    Conductor::Conductor(Material& material, vector<Texture>& textures) : Shader(material, textures)
    {
        auto reflect_ptr = material.getProperty<Property::Wrapper::RGBType>("reflect");
        if (reflect_ptr)
            reflect = (*reflect_ptr).value;
        else
            reflect = {1, 1, 1};
    }

    static float GeometrySchlickGGX(float NdotV, float roughness)
    {
        float k = (roughness * roughness) * 0.5f;
        return NdotV / (NdotV * (1.0f - k) + k);
    }

    static float TrowbridgeReitzGGX(float NdotH, float roughness)
    {
        float a2    = roughness * roughness;
        float denom = (NdotH * NdotH) * (a2 - 1.0f) + 1.0f;
        return a2 / (PI * denom * denom);
    }

    Scattered Conductor::shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal) const
    {
        Vec3  Wr    = glm::normalize(ray.direction - 2 * (glm::dot(ray.direction, normal)) * normal);
        float NdotH = glm::dot(normal, glm::normalize(-ray.direction + Wr));
        float D     = TrowbridgeReitzGGX(NdotH, 0.5f);
        Vec3  f     = reflect + (Vec3(1.0, 1.0, 1.0) - reflect) * pow(1.0f - glm::dot(normal, Wr), 5.0f);
        float NdotV = glm::abs(glm::dot(normal, -ray.direction));
        float G     = GeometrySchlickGGX(NdotV, 0.5f);
        Vec3  attenuation =
            f * D * G / (4.f * glm::abs(glm::dot(Wr, normal)) * glm::abs(glm::dot(ray.direction, normal)));
        return {Ray{hitPoint, Wr}, attenuation, Vec3{0}, 1.0f};
    }
}  // namespace SimplePathTracer