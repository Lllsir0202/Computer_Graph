#pragma once
#ifndef __NR_RENDER_SETTINGS_MANAGER_HPP__
#define __NR_RENDER_SETTINGS_MANAGER_HPP__

#include "scene/Camera.hpp"

namespace NRenderer
{
    struct RenderSettings
    {
        unsigned int width;
        unsigned int height;
        unsigned int depth;
        unsigned int samplesPerPixel;
        unsigned int PhotonNums;
        unsigned int PhotonIters;
        bool ifcaustic;
        RenderSettings()
            : width(500), height(500), depth(4), samplesPerPixel(16), PhotonNums(100000), PhotonIters(10), ifcaustic(false)
        {}
    };
    struct AmbientSettings
    {
        enum Type
        {
            CONSTANT,
            ENVIROMENT_MAP
        };
        Type   type       = Type::CONSTANT;
        RGB    ambient    = {0, 0, 0};
        Handle mapTexture = {};
    };
    struct RenderSettingsManager
    {
        Camera          camera          = {};
        RenderSettings  renderSettings  = {};
        AmbientSettings ambientSettings = {};
    };

}  // namespace NRenderer

#endif