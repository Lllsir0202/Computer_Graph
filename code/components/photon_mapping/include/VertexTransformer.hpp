#pragma once
#ifndef __VERTEX_TRANSFORM_HPP__
#define __VERTEX_TRANSFORM_HPP__

#include "scene/Scene.hpp"

namespace PhotonMapping
{
    using namespace NRenderer;
    // 由局部坐标转换为世界坐标
    class VertexTransformer
    {
      private:
      public:
        void exec(SharedScene spScene);
    };
}  // namespace PhotonMapping

#endif