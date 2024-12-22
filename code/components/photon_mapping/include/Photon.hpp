#pragma once
#ifndef __PHOTON_HPP__
#define __PHOTON_HPP__

#include "scene/Scene.hpp"
#include "Ray.hpp"

// 定义photon类，需要记录光子位置、强度、方向
namespace PhotonMapping
{
    using namespace std;
    using namespace NRenderer;
    class photon
    {
      private:
        // 位置
        Vec3 position;
        RGB  power;
        Ray  input;
        Ray  output;

      public:
        /*photon(Vec3 POS) :position(POS) {};
        photon() {};*/
        photon(Vec3 POS, RGB POW, Ray In, Ray Out) : position(POS), power(POW), input(In), output(Out) {};
        const auto& GetPosition() { return position; }
        const auto& GetPosition() const { return position; }
        const auto& GetPower() { return power; }
        const auto& GetPower() const { return power; }
        const auto& GetInput() { return input; }
        const auto& GetInput() const { return input; }
        const auto& GetOutput() { return output; }
        const auto& GetOutput() const { return output; }

      public:
        float  operator[](int index) const { return position[index]; }
        float& operator[](int index) { return position[index]; }
    };
}  // namespace PhotonMapping

#endif