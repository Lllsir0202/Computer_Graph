#include "shaders/Conductor.hpp"
#include "samplers/SamplerInstance.hpp"

#include "Onb.hpp"

namespace PhotonMapping
{
    Conductor::Conductor(Material& material, vector<Texture>& textures)
        : Shader(material, textures)
    {
        auto reflect_ptr = material.getProperty<Property::Wrapper::RGBType>("reflect");
        if (reflect_ptr) reflect = (*reflect_ptr).value;
        else reflect = { 1, 1, 1 };
    }

    Scattered Conductor::shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal) const 
    {
        Vec3 Wr = glm::normalize(ray.direction - 2 * (glm::dot(ray.direction, normal)) * normal);
        Vec3 f = reflect + (Vec3(1.0, 1.0, 1.0) - reflect) * pow(1.0f - glm::dot(normal, Wr), 5.0f);
        Vec3 attenuation = f / abs(glm::dot(ray.direction, normal));
        return {
            Ray{hitPoint, Wr},
            attenuation,
            Vec3{0},
            1.0f
        };
    }
}