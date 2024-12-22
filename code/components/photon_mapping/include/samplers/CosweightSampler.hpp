#pragma once
#ifndef __COSWEIGHT_HPP__
#define __COSWEIGHT_HPP__

#include "Sampler3d.hpp"
#include <ctime>

namespace PhotonMapping
{
    using namespace std;
    class CosWeightSphere : public Sampler3d
    {
    private:
        constexpr static float C_PI = 3.14159265358979323846264338327950288f;

        default_random_engine e;
        uniform_real_distribution<float> u;
    public:
        CosWeightSphere()
            : e((unsigned int)time(0) + insideSeed())
            , u(0, 1)
        {}

        Vec3 sample3d() override {
            float rand1 = u(e);
            float rand2 = u(e);
            // 余弦加权，rand1控制余弦的极角
            float theta = acos(sqrt(1 - rand1));

            // 方向角，均匀分布
            float phi = 2 * C_PI * rand2;

            float x = sin(theta) * cos(phi);
            float y = sin(theta) * sin(phi);
            float z = cos(theta);
            return { x, y, z };
        }
    };
}

#endif