#pragma once
#ifndef __PHOTONMAPPING_HPP__
#define __PHOTONMAPPING_HPP__

#include "scene/Scene.hpp"
#include "Ray.hpp"
#include "Camera.hpp"
#include "intersections/HitRecord.hpp"

#include "shaders/ShaderCreator.hpp"

// 引入光子
#include "Photon.hpp"

#include <tuple>
namespace PhotonMapping
{
    using namespace NRenderer;
    using namespace std;

    class PhotonMappingRenderer
    {
    public:
    private:
        SharedScene spScene;
        Scene& scene;

        unsigned int width;
        unsigned int height;
        unsigned int depth;
        unsigned int samples;

        // 下面是记录光子数目
        unsigned int photonnum;

        // 用vector先存储记录光子
        vector<photon> Photons;

        // 记录下迭代次数
        unsigned int photoniters;

        using SCam = PhotonMapping::Camera;
        SCam camera;

        vector<SharedShader> shaderPrograms;
    public:
        PhotonMappingRenderer(SharedScene spScene)
            : spScene               (spScene)
            , scene                 (*spScene)
            , camera                (spScene->camera)
        {
            width = scene.renderOption.width;
            height = scene.renderOption.height;
            depth = scene.renderOption.depth;
            samples = scene.renderOption.samplesPerPixel;

            // 添加光子数目
            photonnum = scene.renderOption.photonnum;
            photoniters = scene.renderOption.photoniters;
        }
        ~PhotonMappingRenderer() = default;

        using RenderResult = tuple<RGBA*, unsigned int, unsigned int>;
        RenderResult render();
        void release(const RenderResult& r);

    private:
        void renderTask(RGBA* pixels, int width, int height, int off, int step);

        RGB gamma(const RGB& rgb);
        RGB trace(const Ray& ray, int currDepth);
        HitRecord closestHitObject(const Ray& r);
        tuple<float, Vec3> closestHitLight(const Ray& r);

        // 添加随机发射光子
        void RandomPhoton();

        // 添加光子追踪
        void TracePhoton(const Ray& r, const RGB& power, unsigned depth);
    };
}

#endif