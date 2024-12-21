#include "server/Server.hpp"

#include "PhotonMapping.hpp"

#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"

#include "glm/gtc/matrix_transform.hpp"
#include "Onb.hpp"

#include <random>
#include <chrono>

namespace
{
    static float random_double()
    {
        static thread_local std::mt19937             rng(std::random_device{}());
        static std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        return dist(rng);
    }
}  // namespace

namespace PhotonMapping
{
    static Vec3 RandomPhotonPositionGenerater(const AreaLight& area)
    {
        // 生成随机数
        static thread_local std::mt19937 rng(std::random_device{}());

        // 生成0-1的随机数
        // 此处可设为static
        static std::uniform_real_distribution<float> dist(0.0f, 1.0f);

        float r1 = dist(rng);
        float r2 = dist(rng);

        // 返回随机位置
        // 对于一个面光源，u、v分别表示两个方向，随机取值加上位置，可以得到任意的位置
        return area.position + r1 * area.u + r2 * area.v;
    }

    // 用于俄罗斯轮盘赌决定是否继续发射
    static bool Russian_Roulette(const float p)
    {
        // 生成随机数
        static thread_local std::mt19937 rng(std::random_device{}());

        // 生成0-1的随机数
        // 此处可设为static
        static std::uniform_real_distribution<float> dist(0.0f, 1.0f);

        // 生成随机数
        float rand = dist(rng);
        return rand > p;
    }

    RGB PhotonMappingRenderer::gamma(const RGB& rgb) { return glm::sqrt(rgb); }

    void PhotonMappingRenderer::renderTask(RGBA* pixels, int width, int height, int off, int step)
    {
        for (int i = off; i < height; i += step)
        {
            for (int j = 0; j < width; j++)
            {
                Vec3 color{0, 0, 0};
                for (int k = 0; k < samples; k++)
                {
                    auto  r   = defaultSamplerInstance<UniformInSquare>().sample2d();
                    float rx  = r.x;
                    float ry  = r.y;
                    float x   = (float(j) + rx) / float(width);
                    float y   = (float(i) + ry) / float(height);
                    auto  ray = camera.shoot(x, y);
                    color += trace(ray, 0);
                }
                color /= samples;
                color                                = gamma(color);
                pixels[(height - i - 1) * width + j] = {color, 1};
            }
        }
    }

    auto PhotonMappingRenderer::render() -> RenderResult
    {
        // shaders
        shaderPrograms.clear();
        ShaderCreator shaderCreator{};
        for (auto& m : scene.materials) { shaderPrograms.push_back(shaderCreator.create(m, scene.textures)); }
        RandomPhoton();
        getServer().logger.log("Finish Random Emit photons...");

        // 建立kdtree
        getServer().logger.log("Start to build KD_tree...");
        cout << "Start to build KD_tree..." << endl;

        kdtree = new KDTree<photon>(Photons.begin(), Photons.end());

        // 确认建立完成
        getServer().logger.log("Finish building KD_tree...");
        cout << "Finish building KD_tree..." << endl;

        RGBA* pixels = new RGBA[width * height]{};

        // 局部坐标转换成世界坐标
        VertexTransformer vertexTransformer{};
        vertexTransformer.exec(spScene);

        /*getServer().logger.log("Current photon num is " + to_string(this->photonnum) + "\n");
        getServer().logger.log("Current width is " + to_string(this->width) + "\n");*/

        const auto taskNums = 8;
        thread     t[taskNums];
        for (int i = 0; i < taskNums; i++)
        {
            t[i] = thread(&PhotonMappingRenderer::renderTask, this, pixels, width, height, i, taskNums);
        }
        for (int i = 0; i < taskNums; i++) { t[i].join(); }
        getServer().logger.log("Done...");

        // getServer().logger.log("Test kd-tree");

        // mt19937 rng(42);
        // uniform_real_distribution<float> dist(-1000.0f, 1000.0f);

        // size_t num_photons = 100000;
        // vector<photon> photons;
        // photons.reserve(num_photons);
        // for (size_t i = 0; i < num_photons; ++i)
        // {
        //     Vec3 pos(dist(rng), dist(rng), dist(rng));
        //     photons.emplace_back(pos, Vec3(0.f, 0.f, 0.f), Ray(), Ray());
        // }

        // auto start_build = chrono::high_resolution_clock::now();
        // kdtree = new KDTree<photon, 3>(photons.begin(), photons.end());
        // auto end_build = chrono::high_resolution_clock::now();
        // chrono::duration<double> build_duration = end_build - start_build;
        // getServer().logger.log("KDTree built in " + to_string(build_duration.count()) + " seconds.\n");

        // size_t num_queries = 1000;
        // int k = 10;

        // vector<Vec3> query_points;
        // query_points.reserve(num_queries);
        // for (size_t i = 0; i < num_queries; ++i)
        // {
        //     Vec3 q_pos(dist(rng), dist(rng), dist(rng));
        //     query_points.emplace_back(q_pos);
        // }

        // vector<vector<photon>> kdtree_results;
        // kdtree_results.reserve(num_queries);

        // auto start_kdtree = std::chrono::high_resolution_clock::now();
        // for (const auto& q : query_points)
        // {
        //     kdtree_results.emplace_back(kdtree->kNearest(q, k));
        // }
        // auto end_kdtree = std::chrono::high_resolution_clock::now();
        // chrono::duration<double> kdtree_duration = end_kdtree - start_kdtree;
        // getServer().logger.log("KDTree query done in " + to_string(kdtree_duration.count()) + " seconds.\n");

        return {pixels, width, height};
    }

    void PhotonMappingRenderer::release(const RenderResult& r)
    {
        auto [p, w, h] = r;
        delete[] p;
        if (kdtree)
        {
            delete kdtree;
            kdtree = nullptr;
        }
    }

    HitRecord PhotonMappingRenderer::closestHitObject(const Ray& r)
    {
        HitRecord closestHit = nullopt;
        float     closest    = FLOAT_INF;
        for (auto& s : scene.sphereBuffer)
        {
            auto hitRecord = Intersection::xSphere(r, s, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest)
            {
                closest    = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        for (auto& t : scene.triangleBuffer)
        {
            auto hitRecord = Intersection::xTriangle(r, t, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest)
            {
                closest    = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        for (auto& p : scene.planeBuffer)
        {
            auto hitRecord = Intersection::xPlane(r, p, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest)
            {
                closest    = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        return closestHit;
    }

    tuple<float, Vec3> PhotonMappingRenderer::closestHitLight(const Ray& r)
    {
        Vec3      v       = {};
        HitRecord closest = getHitRecord(FLOAT_INF, {}, {}, {});
        for (auto& a : scene.areaLightBuffer)
        {
            auto hitRecord = Intersection::xAreaLight(r, a, 0.000001, closest->t);
            if (hitRecord && closest->t > hitRecord->t)
            {
                closest = hitRecord;
                v       = a.radiance;
            }
        }
        return {closest->t, v};
    }

    RGB PhotonMappingRenderer::trace(const Ray& r, int currDepth)
    {
        // cout << currDepth << endl;
        if (currDepth == depth) return scene.ambient.constant;
        auto hitObject    = closestHitObject(r);
        auto [t, emitted] = closestHitLight(r);
        // hit object

        if (hitObject && hitObject->t < t)
        {
            auto  mtlHandle    = hitObject->material;
            auto  scattered    = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
            auto  scatteredRay = scattered.ray;
            auto  attenuation  = scattered.attenuation;
            auto  emitted      = scattered.emitted;
            auto  next         = trace(scatteredRay, currDepth + 1);
            float n_dot_in     = glm::dot(hitObject->normal, scatteredRay.direction);
            float pdf          = scattered.pdf;
            /**
             * emitted      - Le(p, w_0)
             * next         - Li(p, w_i)
             * n_dot_in     - cos<n, w_i>
             * atteunation  - BRDF
             * pdf          - p(w)
             **/

            Vec3 directLighting(0.0f);
            for (const auto& light : scene.areaLightBuffer)
            {
                auto [radiance, pdf] = sampleDirectLighting(hitObject, light);
                directLighting += scattered.attenuation * radiance;
            }

            RGB              indirectLighting(0.0f);
            const float      P_RR = 0.8f;
            static const int k    = 50;
            if (random_double() < P_RR)
            {
                auto nearPhotons = kdtree->kNearest(hitObject->hitPoint, k);

                if (!nearPhotons.empty())
                {
                    float maxDist = 0.0f;
                    for (const auto& photon : nearPhotons)
                        maxDist = std::max(maxDist, glm::distance(photon.GetPosition(), hitObject->hitPoint));

                    for (const auto& photon : nearPhotons)
                    {
                        float dist   = glm::distance(photon.GetPosition(), hitObject->hitPoint);
                        float weight = 1.0f - (dist * dist) / (maxDist * maxDist);

                        float cos_theta = glm::dot(hitObject->normal, -photon.GetInput().direction);
                        if (cos_theta <= 0.0f) continue;
                        indirectLighting +=
                            photon.GetPower() * weight * scattered.attenuation * cos_theta / (PI * maxDist * maxDist);
                    }
                    indirectLighting /= P_RR;
                }
            }

            return scattered.emitted + directLighting + indirectLighting;
        }
        else if (t != FLOAT_INF) { return emitted; }
        return Vec3(0.0f);
    }

    void PhotonMappingRenderer::RandomPhoton()
    {
        getServer().logger.log("Start to Random Shoot Photon...");

        const float P_RR = 0.8f;
        for (int i = 0; i < photonnum; ++i)
        {
            for (auto& area : scene.areaLightBuffer)
            {
                Vec3 Position  = RandomPhotonPositionGenerater(area);
                Vec3 Local_dir = defaultSamplerInstance<HemiSphere>().sample3d();
                Vec3 World_dir = glm::normalize(Onb(glm::normalize(glm::cross(area.u, area.v))).local(Local_dir));

                float area_size = glm::length(glm::cross(area.u, area.v));
                RGB   Power     = (area.radiance * area_size) / (static_cast<float>(photonnum) * PI);

                Ray r(Position, World_dir);
                TracePhoton(r, Power, 0);
            }
        }

        getServer().logger.log("Photon mapping complete. Total photons: " + std::to_string(Photons.size()));
    }

    // 用于追踪随机发射的光子
    void PhotonMappingRenderer::TracePhoton(const Ray& r, const RGB& power, unsigned depth)
    {
        // getServer().logger.log("Current Depth is " + to_string(depth) + "/" + to_string(this->depth) + "\n");
        if (depth > this->depth) { return; }

        // 找到最近的hitobject   ->  后续可以改进对是否命中的判断，添加包围
        HitRecord hitrecord = closestHitObject(r);
        // 如果没有hit，那么返回
        if (!hitrecord) { return; }
        /*if (depth > 1)
        {
            cout << "her" << endl;
        }*/
        // 分别得到hit点、法向量、材质
        const auto& hitpoint = hitrecord->hitPoint;
        const auto& normal   = hitrecord->normal;
        const auto& material = hitrecord->material;

        // 得到打在hitobject上后的光线
        const auto& scatter = shaderPrograms[material.index()]->shade(r, hitpoint, normal);

        // 使用亮度作为轮盘赌的参数
        // float L = 0.2126 * power.r + 0.7152 * power.g + 0.0722 * power.b;
        // 不知道为何，这里的RGB并不是0-1(?)
        // 所以使用亮度不太合适了
        // 这里考虑结合法线角度和路径长度来作为依据
        // 再加上一个depth / maxdepth
        // auto ndoti = glm::dot(hitrecord->normal, r.direction);
        // float p = 1.f - 0.5 * (ndoti > 0.0f ? ndoti : 0) - 0.5 * static_cast<float>(depth) / this->depth;
        float p = 0.8;
        // cout << p << endl;
        //  接下来根据材质进行判断
        if (spScene->materials[material.index()].type == 0)  // 表示漫反射
        {
            // 漫反射需要记录光子
            photon Photon(hitpoint, power, r, scatter.ray);
            Photons.push_back(Photon);

            // 根据轮盘赌策略计算是否继续
            if (Russian_Roulette(p))
            {
                // 计算新的ray和power
                // 漫反射的光强与表面法线相关
                auto cos_thera = (glm::abs(glm::dot(hitrecord->normal, -r.direction)));
                RGB  newpower  = power * scatter.attenuation * cos_thera / (scatter.pdf * p);
                TracePhoton(scatter.ray, newpower, depth + 1);
            }
        }

        // 如果是镜面的话，那么只有反射没有散射
        if (spScene->materials[material.index()].type == 2)  // 表示镜面反射(即导体)
        {
            if (Russian_Roulette(p))
            {
                // 同样的计算新的ray和power
                RGB newpower = power * scatter.attenuation / (scatter.pdf * p);

                // 对于导体来说，应该是需要计算反射方向的
                // glm提供了直接计算反射的()
                Vec3 reflectDir = glm::reflect(r.direction, hitrecord->normal);
                Ray  reflectRay(hitpoint, reflectDir);
                TracePhoton(reflectRay, newpower, depth + 1);
            }
        }

        // 如果是电介质的话，可能会同时存在折射和散射
        if (spScene->materials[material.index()].type == 3)  // 表示电介质
        {
            // 采用俄罗斯轮盘赌决定处理散射还是折射
            if (Russian_Roulette(p))
            {
                // 处理折射
                RGB newpower = power * scatter.r_attenuation / scatter.r_pdf;

                TracePhoton(scatter.r_ray, newpower, depth + 1);
            }
            else
            {
                // 处理反射
                RGB newpower = power * scatter.attenuation / scatter.pdf;

                TracePhoton(scatter.ray, newpower, depth + 1);
            }
        }

        getServer().logger.log("Current size of Photons is " + to_string(Photons.size()) + "\n");
    }

    // 用于计算hitpoint提供的间接光强
    RGB PhotonMappingRenderer::EstimateIndirectRadiance(const Ray& r, const HitRecord& Hit)
    {
        RGB         IndirectPower = {0.f, 0.f, 0.f};
        const auto& Scatter       = shaderPrograms[Hit->material.index()]->shade(r, Hit->hitPoint, Hit->normal);
        float       radius        = 10.0;
        // const auto& photons = kdtree->withinRadius(Hit->hitPoint, radius);
        const auto& photons = kdtree->kNearest(Hit->hitPoint, radius);

        cout << photons.size() << endl;
        for (auto& photon : photons)
        {
            // 首先计算光子位置到Hitpoint的方向向量
            const auto& pos_vec = glm::normalize(photon.GetPosition() - r.direction);
            // 计算光子和法向量的cos
            const auto& cos_thera = glm::dot(pos_vec, Hit->normal);

            // 表示在背面，所以继续
            if (cos_thera < 0.f) { continue; }

            IndirectPower += photon.GetPower() * Scatter.attenuation * cos_thera;
        }
        return IndirectPower;
    }

    tuple<Vec3, Vec3> PhotonMappingRenderer::sampleOnLight(const AreaLight& light_source)
    {
        const Vec3& light_pos = light_source.position;
        const Vec3& light_u   = light_source.u;
        const Vec3& light_v   = light_source.v;

        Vec3 light_normal = glm::normalize(glm::cross(light_u, light_v));
        Vec3 sample_point = light_pos + light_u * random_double() + light_v * random_double();

        return {sample_point, light_normal};
    }

    DirectLightingRes PhotonMappingRenderer::sampleDirectLighting(
        const HitRecord& hit_record, const AreaLight& light_source)
    {
        auto [sample_point, light_normal] = sampleOnLight(light_source);

        Vec3 shadow_ray_dir = glm::normalize(sample_point - hit_record->hitPoint);
        Ray  shadow_ray{hit_record->hitPoint, shadow_ray_dir};

        float dist_to_light = glm::length(sample_point - hit_record->hitPoint);
        float cos_theta     = glm::dot(-shadow_ray_dir, light_normal);
        float pdf_light     = 1.0f / (glm::length(light_source.u) * glm::length(light_source.v));

        auto shadow_hit = closestHitObject(shadow_ray);
        if (shadow_hit && shadow_hit->t < dist_to_light || cos_theta < 0.0001f) return {Vec3(0.0f), pdf_light};

        float surface_cos = glm::dot(hit_record->normal, shadow_ray_dir);
        Vec3  radiance = light_source.radiance * surface_cos * cos_theta / (dist_to_light * dist_to_light * pdf_light);

        return {radiance, pdf_light};
    }
}  // namespace PhotonMapping