#include "server/Server.hpp"

#include "PhotonMapping.hpp"

#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"

#include "glm/gtc/matrix_transform.hpp"
#include "Onb.hpp"

#include <random>

namespace PhotonMapping
{
    static Vec3 RandomPhotonPositionGenerater(const AreaLight& area)
    {
        // 生成随机数
        static thread_local std::mt19937 rng(std::random_device{}());

        // 生成0-1的随机数
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);

        float r1 = dist(rng);
        float r2 = dist(rng);

        // 返回随机位置
        // 对于一个面光源，u、v分别表示两个方向，随机取值加上位置，可以得到任意的位置
        return area.position + r1 * area.u + r2 * area.v;
    }

    RGB PhotonMappingRenderer::gamma(const RGB& rgb) {
        return glm::sqrt(rgb);
    }

    void PhotonMappingRenderer::renderTask(RGBA* pixels, int width, int height, int off, int step) {
        for(int i=off; i<height; i+=step) {
            for (int j=0; j<width; j++) {
                Vec3 color{0, 0, 0};
                for (int k=0; k < samples; k++) {
                    auto r = defaultSamplerInstance<UniformInSquare>().sample2d();
                    float rx = r.x;
                    float ry = r.y;
                    float x = (float(j)+rx)/float(width);
                    float y = (float(i)+ry)/float(height);
                    auto ray = camera.shoot(x, y);
                    color += trace(ray, 0);
                }
                color /= samples;
                color = gamma(color);
                pixels[(height-i-1)*width+j] = {color, 1};
            }
        }
    }

    auto PhotonMappingRenderer::render() -> RenderResult {
        // shaders
        shaderPrograms.clear();
        ShaderCreator shaderCreator{};
        for (auto& m : scene.materials) {
            shaderPrograms.push_back(shaderCreator.create(m, scene.textures));
        }

        RGBA* pixels = new RGBA[width*height]{};

        // 局部坐标转换成世界坐标
        VertexTransformer vertexTransformer{};
        vertexTransformer.exec(spScene);

        const auto taskNums = 8;
        thread t[taskNums];
        for (int i=0; i < taskNums; i++) {
            t[i] = thread(&PhotonMappingRenderer::renderTask,
                this, pixels, width, height, i, taskNums);
        }
        for(int i=0; i < taskNums; i++) {
            t[i].join();
        }
        getServer().logger.log("Done...");
        return {pixels, width, height};
    }

    void PhotonMappingRenderer::release(const RenderResult& r) {
        auto [p, w, h] = r;
        delete[] p;
    }

    HitRecord PhotonMappingRenderer::closestHitObject(const Ray& r) {
        HitRecord closestHit = nullopt;
        float closest = FLOAT_INF;
        for (auto& s : scene.sphereBuffer) {
            auto hitRecord = Intersection::xSphere(r, s, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        for (auto& t : scene.triangleBuffer) {
            auto hitRecord = Intersection::xTriangle(r, t, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        for (auto& p : scene.planeBuffer) {
            auto hitRecord = Intersection::xPlane(r, p, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        return closestHit; 
    }
    
    tuple<float, Vec3> PhotonMappingRenderer::closestHitLight(const Ray& r) {
        Vec3 v = {};
        HitRecord closest = getHitRecord(FLOAT_INF, {}, {}, {});
        for (auto& a : scene.areaLightBuffer) {
            auto hitRecord = Intersection::xAreaLight(r, a, 0.000001, closest->t);
            if (hitRecord && closest->t > hitRecord->t) {
                closest = hitRecord;
                v = a.radiance;
            }
        }
        return { closest->t, v };
    }

    RGB PhotonMappingRenderer::trace(const Ray& r, int currDepth) {
        if (currDepth == depth) return scene.ambient.constant;
        auto hitObject = closestHitObject(r);
        auto [ t, emitted ] = closestHitLight(r);
        // hit object
        if (hitObject && hitObject->t < t) {
            auto mtlHandle = hitObject->material;
            auto scattered = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
            auto scatteredRay = scattered.ray;
            auto attenuation = scattered.attenuation;
            auto emitted = scattered.emitted;
            auto next = trace(scatteredRay, currDepth+1);
            float n_dot_in = glm::dot(hitObject->normal, scatteredRay.direction);
            float pdf = scattered.pdf;
            /**
             * emitted      - Le(p, w_0)
             * next         - Li(p, w_i)
             * n_dot_in     - cos<n, w_i>
             * atteunation  - BRDF
             * pdf          - p(w)
             **/

            RGB refrac_res{ 0.0f,0.0f,0.0f };
            if (scattered.has_refraction && currDepth <= 2) 
            {
                auto refraction_ray = scattered.r_ray;
                auto refraction_next = trace(refraction_ray, currDepth + 1);
                float r_n_dot_in = glm::dot(hitObject->normal, refraction_ray.direction);
                refrac_res += refraction_next * abs(r_n_dot_in) * scattered.r_attenuation / scattered.r_pdf;
            }

            return emitted + attenuation * next * n_dot_in / pdf + refrac_res;
        }
        // 
        else if (t != FLOAT_INF) {
            return emitted;
        }
        else {
            return Vec3{0};
        }
    }

    void PhotonMappingRenderer::RandomPhoton()
    {
        getServer().logger.log("Start to Random Shoot Photon...");
        // 接下来处理采样，首先对于随机的光子数目，产生随机位置和方向
        for (int i = 0; i < photonnum; i++)
        {
            for (auto& area : scene.areaLightBuffer)
            {
                // 得到面光源表面的随机位置
                Vec3 Position = RandomPhotonPositionGenerater(area);

                // 接下来需要计算面光源发出光线的随机方向
                // 这里先使用实现好的半球，后面考虑重新实现一个余弦加权随机
                Vec3 Local_dir = defaultSamplerInstance<HemiSphere>().sample3d();

                // 转换之前我们需要知道面光源的法向量
                Vec3 AreaLight_normal = glm::normalize(glm::cross(area.u, area.v));
                Vec3 World_dir = glm::normalize(Onb(AreaLight_normal).local(Local_dir));
                
                // 得到光强
                // 首先计算光源面积
                float AreaLight_Area = glm::length(glm::cross(area.u, area.v));
                // 由于我们需要计算的是
                RGB Power = (area.radiance * AreaLight_Area) / (static_cast<float>(photonnum));

                // 得到光
                Ray r(Position, World_dir);

                // 接下来对光子进行追踪
                TracePhoton(r, Power, 0);
            }
        }
    }

    // 用于追踪随机发射的光子
    void PhotonMappingRenderer::TracePhoton(const Ray& r, const RGB& power, unsigned depth)
    {
        getServer().logger.log("Current Depth is " + to_string(depth) + "/" + to_string(this->depth) + "\n");
        if (depth > this->depth)
        {
            return;
        }

        // 找到最近的hitobject
        HitRecord hitrecord = closestHitObject(r);
        // 如果没有hit，那么返回
        if (!hitrecord)
        {
            return;
        }

        // 分别得到hit点、法向量、材质
        const auto& hitpoint = hitrecord->hitPoint;
        const auto& normal = hitrecord->normal;
        const auto& material = hitrecord->material;

        // 得到打在hitobject上后的光线
        const auto& scatter = shaderPrograms[material.index()]->shade(r, hitpoint, normal);

        // 接下来根据材质进行判断
        if (spScene->materials[material.index()].type == 0) // 表示漫反射
        {

        }
    }
}

