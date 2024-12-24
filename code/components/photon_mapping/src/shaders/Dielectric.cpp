#include "shaders/Dielectric.hpp"
#include "samplers/SamplerInstance.hpp"

#include "Onb.hpp"

namespace PhotonMapping
{
    Dielectric::Dielectric(Material& material, vector<Texture>& textures) : Shader(material, textures)
    {
        auto absorbed_ptr = material.getProperty<Property::Wrapper::RGBType>("absorbed");
        if (absorbed_ptr)
            absorbed = (*absorbed_ptr).value;
        else
            absorbed = {1, 1, 1};

        auto ior_ptr = material.getProperty<Property::Wrapper::FloatType>("ior");
        if (ior_ptr)
            ior = (*ior_ptr).value;
        else
            ior = 1;
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

    Scattered Dielectric::shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal) const
    {
        Vec3  n           = glm::dot(normal, ray.direction) > 0 ? -normal : normal;
        Vec3  Wr          = glm::normalize(ray.direction - 2 * (glm::dot(ray.direction, n)) * n);
        float NdotH       = glm::dot(n, glm::normalize(-ray.direction + Wr));
        float D           = TrowbridgeReitzGGX(NdotH, 0.5f);
        Vec3  f           = absorbed + (Vec3(1.0, 1.0, 1.0) - absorbed) * pow(1.0f - glm::dot(n, Wr), 5.0f);
        float NdotV       = glm::abs(glm::dot(n, -ray.direction));
        float G           = GeometrySchlickGGX(NdotV, 0.5f);
        Vec3  attenuation = f * D * G / (4.f * glm::abs(glm::dot(Wr, n)) * glm::abs(glm::dot(ray.direction, n)));

        float ni_nt    = glm::dot(normal, -ray.direction) > 0 ? (1.0f / ior) : (ior / 1.0f);
        float cos_val  = abs(glm::dot(n, -ray.direction));
        float sin_val  = max(0.f, 1.f - cos_val * cos_val);
        float sin_val2 = ni_nt * ni_nt * sin_val;
        if (sin_val2 >= 1.0f - 0.00001f) return {Ray{hitPoint, Wr}, attenuation, Vec3{0}, 1.0f};

        float cos_T = sqrt(1 - sin_val2);
        Vec3  Wt    = glm::normalize(ni_nt * ray.direction + (ni_nt * cos_val - cos_T) * n);
        f           = absorbed + (Vec3(1.0, 1.0, 1.0) - absorbed) * pow(1.0f - abs(glm::dot(n, -ray.direction)), 5.0f);
        Vec3 r_attenuation = (1.f / (ni_nt * ni_nt)) * (Vec3{1.0, 1.0, 1.0} - f) / abs(cos_val);

        return {Ray{hitPoint, Wr}, attenuation, Vec3{0}, 1.0f, true, Ray{hitPoint, Wt}, r_attenuation, 1.0f};
    }
}  // namespace PhotonMapping