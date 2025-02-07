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
// 引入KDTree
#include "KDTree.hpp"
// 引入BVH
#include "BVH.hpp"

#include <tuple>
namespace PhotonMapping
{
    using namespace NRenderer;
    using namespace std;

    struct DirectLightingRes
    {
        Vec3  radiance;
        float pdf;
    };

    class PhotonMappingRenderer
    {
      public:
      private:
        SharedScene spScene;
        Scene&      scene;

        unsigned int width;
        unsigned int height;
        unsigned int depth;
        unsigned int samples;

        // 添加一个选项设定是否考虑焦散
        bool ifcaustic;

        // 下面是记录光子数目
        unsigned int photonnum;

        // 用vector先存储记录光子
        vector<photon> Photons;

        // 记录焦散高强度光子
        vector<photon> Caustics_Photons;

        // 记录下迭代次数 maybe unused
        unsigned int photoniters;
        // 模板初始化方式有点固定，但不太想多改动
        // 因此收集到光子后再初始化kdtree
        KDTree<photon>* kdtree;

        // 添加焦散现象的考虑
        // 这里记录下所有反射一次以上的，且采用轮盘赌，仅记录深度较浅的较高质量光子
        KDTree<photon>* caustics_kdtree;
        AABBBVH         bvh;

        using SCam = PhotonMapping::Camera;
        SCam camera;

        vector<SharedShader> shaderPrograms;

      public:
        PhotonMappingRenderer(SharedScene spScene) : spScene(spScene), scene(*spScene), camera(spScene->camera)
        {
            width   = scene.renderOption.width;
            height  = scene.renderOption.height;
            depth   = scene.renderOption.depth;
            samples = scene.renderOption.samplesPerPixel;

            // 添加光子数目
            photonnum   = scene.renderOption.photonnum;
            photoniters = scene.renderOption.photoniters;

            // 是否考虑焦散
            ifcaustic = scene.renderOption.ifcaustic;

            kdtree = nullptr;
            caustics_kdtree = nullptr;
        }
        ~PhotonMappingRenderer()
        {
            if (kdtree) delete kdtree;
            if (caustics_kdtree) delete caustics_kdtree;
        }

        using RenderResult = tuple<RGBA*, unsigned int, unsigned int>;
        RenderResult render();
        void         release(const RenderResult& r);

      private:
        void renderTask(RGBA* pixels, int width, int height, int off, int step);

        RGB                gamma(const RGB& rgb);
        RGB                trace(const Ray& ray, int currDepth);
        HitRecord          closestHitObject(const Ray& r);
        tuple<float, Vec3> closestHitLight(const Ray& r);

        // 添加随机发射光子
        void RandomPhoton();

        // 添加光子追踪
        void TracePhoton(const Ray& r, const RGB& power, unsigned depth);

        // 用于处理焦散
        void TraceCaustics(const Ray& r, const RGB& power, unsigned depth);
        void RandomCaustics();


        // 用于估计光子提供的间接光照强度
        // 根据传入的hitpoint，查询到一定范围内的光子，并返回间接光照强度
        //RGB EstimateIndirectRadiance(const Ray& r, const HitRecord& Hit);

        tuple<Vec3, Vec3> sampleOnLight(const AreaLight& light);
        DirectLightingRes sampleDirectLighting(const HitRecord& hit, const AreaLight& light);

      private:
        // 这个函数是简单的暴力搜索，仅用于性能对比
        std::vector<photon> getkNearestPhotons(const Vec3& pos, size_t k);
    };
}  // namespace PhotonMapping

#endif