#pragma once
#ifndef __SCATTERED_HPP__
#define __SCATTERED_HPP__

#include "Ray.hpp"

namespace PhotonMapping
{
    struct Scattered
    {
        // 散射光线
        Ray ray = {};
        // 衰减值->或者应该说是衰减因子
        Vec3 attenuation = {};
        // 发射的光强，我们目前的应该都是0
        Vec3 emitted = {};
        // 散射光线的概率密度函数
        float pdf = {0.f};

        // 以下表示折射

        // 表示是否发生折射
        bool has_refraction = false;
        // 表示折射光线
        Ray r_ray = {};
        // 表示折射衰减
        Vec3 r_attenuation = {};
        // 表示折射概率密度函数
        float r_pdf = {0.f};
    };

}  // namespace PhotonMapping

#endif