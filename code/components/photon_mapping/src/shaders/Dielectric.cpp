#include "shaders/Dielectric.hpp"
#include "samplers/SamplerInstance.hpp"

#include "Onb.hpp"

namespace PhotonMapping
{
    Dielectric::Dielectric(Material& material, vector<Texture>& textures)
        : Shader(material, textures)
    {
        auto absorbed_ptr = material.getProperty<Property::Wrapper::RGBType>("absorbed");
        if (absorbed_ptr) absorbed = (*absorbed_ptr).value;
        else absorbed = { 1, 1, 1 };

        auto ior_ptr = material.getProperty<Property::Wrapper::FloatType>("ior");
        if (ior_ptr) ior = (*ior_ptr).value;
        else ior = 1;
    }

    Scattered Dielectric::shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal) const
    {
        Vec3 n = glm::dot(normal, ray.direction) > 0 ? -normal : normal;
        Vec3 Wr = glm::normalize(ray.direction - 2 * (glm::dot(ray.direction, n)) * n);
        Vec3 f = absorbed + (Vec3(1.0, 1.0, 1.0) - absorbed) * pow(1.0f - abs(glm::dot(normal, Wr)), 5.0f);
        Vec3 attenuation = f / abs(glm::dot(ray.direction, n));
        Vec3 Wi = -ray.direction;

        float ni_nt = glm::dot(normal, Wi) > 0 ? (1.0f / ior) : (ior / 1.0f);
        float cos_val = abs(glm::dot(n, Wi));
        float sin_val = max(0.f, 1.f - cos_val * cos_val);
        float sin_val2 = ni_nt * ni_nt * sin_val;
        if (sin_val2 >= 1.0f - 0.00001f)
            return {
                Ray{hitPoint, Wr},
                attenuation,
                Vec3{0},
                1.0f
            };

        float cos_T = sqrt(1 - sin_val2);
        Vec3 Wt = ni_nt * (-Wi) + (ni_nt * cos_val - cos_T) * n;
        Wt = glm::normalize(Wt);
        f = absorbed + (Vec3(1.0, 1.0, 1.0) - absorbed) * pow(1.0f - abs(glm::dot(n, Wi)), 5.0f);
        Vec3 r_attenuation = (1.f / (ni_nt * ni_nt)) * (Vec3{ 1.0,1.0,1.0 } - f) / abs(cos_val);

        return {
            Ray{hitPoint, Wr},
            attenuation,
            Vec3{0},
            1.0f,
            true,
            Ray{hitPoint,Wt},
            r_attenuation,
            1.0f
        };
    }
}