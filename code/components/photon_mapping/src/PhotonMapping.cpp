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
    static float RandomFloat()
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
        return rand < p;
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

        // 局部坐标转换成世界坐标
        VertexTransformer vertexTransformer{};
        vertexTransformer.exec(spScene);
        bvh.build(spScene);

        RandomPhoton();

        getServer().logger.log("Finish Random Emit photons...");

        // 建立kdtree
        getServer().logger.log("Start to build KD_tree...");
        cout << "Start to build KD_tree..." << endl;

        kdtree = new KDTree<photon>(Photons.begin(), Photons.end());
        Photons.clear();

        // 确认建立完成
        getServer().logger.log("Finish building KD_tree...");
        cout << "Finish building KD_tree..." << endl;

        // 建立焦散
        if(ifcaustic)
        {
            RandomCaustics();
            getServer().logger.log("Start to build Caustic KD_tree...");
            cout << "Start to build Caustic KD_tree..." << endl;
            caustics_kdtree = new KDTree<photon>(Caustics_Photons.begin(), Caustics_Photons.end());
            Caustics_Photons.clear();
            getServer().logger.log("Finish building Caustic KD_tree...");
            cout << "Finish building Caustic KD_tree..." << endl;
        }


        RGBA* pixels = new RGBA[width * height]{};

        /*getServer().logger.log("Current photon num is " + to_string(this->photonnum) + "\n");
        getServer().logger.log("Current width is " + to_string(this->width) + "\n");*/

        const auto taskNums = 16;
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
        bvh.destory();
        if (caustics_kdtree)
        {
            delete caustics_kdtree;
            caustics_kdtree = nullptr;
        }
    }

    HitRecord PhotonMappingRenderer::closestHitObject(const Ray& r)
    {
        HitRecord closestHit = nullopt;
        float     closest    = FLOAT_INF;
        return bvh.intersect(r, 0.000001, closest);
        
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
        if (currDepth >= depth) return scene.ambient.constant;

        auto hitObject    = closestHitObject(r);
        auto [t, emitted] = closestHitLight(r);
        if (hitObject && hitObject->t < t)
        {
            auto type = spScene->materials[hitObject->material.index()].type;  // 0->漫反射 2->电介质 3->镜面

            if (type == 0) // 表示处理漫反射
            {
                auto  mtlHandle = hitObject->material;
                auto  scattered = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
                auto  scatteredRay      = scattered.ray;
                auto  attenuation       = scattered.attenuation;
                auto  scattered_emitted = scattered.emitted;
                auto  next              = trace(scatteredRay, currDepth + 1);
                float n_dot_in          = glm::dot(hitObject->normal, scatteredRay.direction);
                float pdf               = scattered.pdf;

                // 直接照明
                Vec3 directLighting(0.0f);
                for (const auto& light : scene.areaLightBuffer)
                {
                    auto [radiance, pdf_light] = sampleDirectLighting(hitObject, light);
                    directLighting += scattered.attenuation * radiance;
                }

                // 间接照明
                RGB   indirectLighting(0.0f);
                auto  nearPhotons = kdtree->kNearest(hitObject->hitPoint, photoniters);
            // auto  nearPhotons = getkNearestPhotons(hitObject->hitPoint, photoniters);
                float size        = nearPhotons.size();

                if (!nearPhotons.empty())
                {
                    float maxDist = 0.0f;
                    for (const auto& photon : nearPhotons)
                        maxDist = std::max(maxDist, glm::distance(photon.GetPosition(), hitObject->hitPoint));
                    for (const auto& photon : nearPhotons)
                    {
                        float dist   = glm::distance(photon.GetPosition(), hitObject->hitPoint);
                        float weight = 1.0f - (dist * dist) / (maxDist * maxDist);  // 或使用高斯衰减

                        float cos_theta = glm::dot(hitObject->normal, -photon.GetInput().direction);
                        if (cos_theta <= 0.0f) continue;
                        indirectLighting +=
                            photon.GetPower() * weight * scattered.attenuation * cos_theta / (PI * maxDist * maxDist);
                    }
                }

                // 焦散因素
                RGB CausticLighting(0.0f);
                if(caustics_kdtree && ifcaustic)
                {
                    auto CausticPhotons = caustics_kdtree->kNearest(hitObject->hitPoint, 25);
                    if (!CausticPhotons.empty())
                    {
                        float radius  = 50;
                        float maxDist = 0.0f;
                        for (const auto& photon : CausticPhotons)
                            maxDist = std::max(maxDist, glm::distance(photon.GetPosition(), hitObject->hitPoint));
                        if (maxDist > 1e-4)
                        {
                            for (const auto& photon : CausticPhotons)
                            {
                                float dist   = glm::distance(photon.GetPosition(), hitObject->hitPoint);
                                float weight = 1.0f - (dist * dist) / (maxDist * maxDist);  // 或使用高斯衰减

                                float cos_theta = glm::dot(hitObject->normal, -photon.GetInput().direction);
                                if (cos_theta <= 0.0f) continue;
                                CausticLighting += photon.GetPower() * weight * scattered.attenuation * cos_theta /
                                                   (PI * maxDist * maxDist);
                            }
                        }
                        
                    }
                }

                // 如果漫反射到光源，则不能计算
                auto NextHitobject = closestHitObject(scatteredRay);
                if (NextHitobject && spScene->materials[NextHitobject->material.index()].type == 0)
                {
                    indirectLighting += gamma(next) / size;  // 加入后续光照的处理
                }

                RGB finalColor = scattered_emitted + directLighting + indirectLighting + attenuation * next;
                return glm::clamp(finalColor, Vec3(0.0f), Vec3(1.0f));
            }
            else if (type == 2) // 电介质
            {
                auto  mtlHandle = hitObject->material;
                auto  scattered = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
                auto  scatteredRay      = scattered.ray;
                auto  attenuation       = scattered.attenuation;
                auto  scattered_emitted = scattered.emitted;
                auto  next              = trace(scatteredRay, currDepth + 1);
                float n_dot_in          = glm::dot(hitObject->normal, scatteredRay.direction);
                float pdf               = scattered.pdf;

                // 直接照明
                Vec3 directLighting(0.0f);

                if(scattered.has_refraction)
                {
                    // 如果存在折射，那么进行计算
                    auto  r_ray         = scattered.r_ray;
                    auto  r_attenuation = scattered.r_attenuation;
                    auto  r_pdf         = scattered.r_pdf;
                    float r_n_dot_in    = glm::dot(hitObject->normal, r_ray.direction);
                    auto  r_next        = trace(r_ray, currDepth + 1);

                    // 计算折射的光
                    auto r_Lighting = r_next * glm::abs(r_n_dot_in) * r_attenuation / r_pdf;
                    directLighting += r_Lighting;
                }

               
                directLighting += scattered_emitted + attenuation * next * n_dot_in + spScene->ambient.constant;
                return glm::clamp(directLighting, Vec3(0.0f), Vec3(1.0f));
            }
            else if (type == 3) // 镜面
            {
                auto  mtlHandle = hitObject->material;
                auto  scattered = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
                auto  scatteredRay      = scattered.ray;
                auto  attenuation       = scattered.attenuation;
                auto  scattered_emitted = scattered.emitted;
                auto  next              = trace(scatteredRay, currDepth + 1);
                float n_dot_in          = glm::dot(hitObject->normal, scatteredRay.direction);
                float pdf               = scattered.pdf;

                // 直接照明
                Vec3 directLighting(0.0f);
                directLighting =
                    scattered_emitted + attenuation * next * glm::abs(n_dot_in) + spScene->ambient.constant;
                return glm::clamp(directLighting, Vec3(0.0f), Vec3(1.0f));
            }
        }
        else if (t != FLOAT_INF) { return glm::clamp(emitted, Vec3(0.0f), Vec3(1.0f)); }
        return spScene->ambient.constant;
        return Vec3(0.0f);
    }

    void PhotonMappingRenderer::RandomPhoton()
    {
        getServer().logger.log("Start to Random Shoot Photon...");

        // 对面光源上进行光子数目个采样
        // for(int j = 0; j < 5 ; j++)
        {
            for (auto& area : scene.areaLightBuffer)
            {
                for (int i = 0; i < photonnum; i++)
                {
                    // 得到面光源表面的随机位置
                    Vec3 Position = RandomPhotonPositionGenerater(area);

                    // 接下来需要计算面光源发出光线的随机方向
                    // 这里先使用实现好的半球，后面考虑重新实现一个余弦加权随机
                    // Vec3 Local_dir = defaultSamplerInstance<HemiSphere>().sample3d();
                    // 采用余弦加权采样
                    Vec3 Local_dir = defaultSamplerInstance<CosWeightSphere>().sample3d();

                    // 转换之前我们需要知道面光源的法向量
                    Vec3 AreaLight_normal = glm::normalize(glm::cross(area.u, area.v));
                    // cout << "u is " << area.u << " v is " << area.v  << AreaLight_normal << endl;
                    Vec3 World_dir = glm::normalize(Onb(AreaLight_normal).local(Local_dir));
                    // cout << "world dir is " << World_dir << endl;
                    //  得到光强
                    //  首先计算光源面积
                    float AreaLight_Area = glm::length(glm::cross(area.u, area.v));
                    // 由于我们需要计算的是
                    // 由于面光源的radiance是光源在单位立体角和单位面积上的发射功率
                    // 所以我们还需要再除以一个PI
                    RGB Power = (area.radiance * AreaLight_Area) / ((static_cast<float>(photonnum) * PI));

                    // 这里的rgb是没有问题的
                    // cout << "power is : R " << Power.r << " G " << Power.g << " B " << Power.b << endl;

                    Ray r(Position, World_dir);
                    TracePhoton(r, Power, 0);
                }
            }
        }

        getServer().logger.log("Photon mapping complete. Total photons: " + std::to_string(Photons.size()));
    }

    // 用于处理焦散的
    void PhotonMappingRenderer::RandomCaustics()
    {
        getServer().logger.log("Start to Random Shoot Photon for caustic...");

        // 对面光源上进行光子数目个采样
        // for(int j = 0; j < 5 ; j++)
        {
            for (auto& area : scene.areaLightBuffer)
            {
                for (int i = 0; i < photonnum; i++)
                {
                    // 得到面光源表面的随机位置
                    Vec3 Position = RandomPhotonPositionGenerater(area);

                    // 接下来需要计算面光源发出光线的随机方向
                    // 这里先使用实现好的半球，后面考虑重新实现一个余弦加权随机
                    // Vec3 Local_dir = defaultSamplerInstance<HemiSphere>().sample3d();
                    // 采用余弦加权采样
                    Vec3 Local_dir = defaultSamplerInstance<CosWeightSphere>().sample3d();

                    // 转换之前我们需要知道面光源的法向量
                    Vec3 AreaLight_normal = glm::normalize(glm::cross(area.u, area.v));
                    // cout << "u is " << area.u << " v is " << area.v  << AreaLight_normal << endl;
                    Vec3 World_dir = glm::normalize(Onb(AreaLight_normal).local(Local_dir));
                    // cout << "world dir is " << World_dir << endl;
                    //  得到光强
                    //  首先计算光源面积
                    float AreaLight_Area = glm::length(glm::cross(area.u, area.v));
                    // 由于我们需要计算的是
                    // 由于面光源的radiance是光源在单位立体角和单位面积上的发射功率
                    // 所以我们还需要再除以一个PI
                    RGB Power = (area.radiance * AreaLight_Area) / ((static_cast<float>(photonnum) * PI));

                    // 这里的rgb是没有问题的
                    // cout << "power is : R " << Power.r << " G " << Power.g << " B " << Power.b << endl;

                    Ray r(Position, World_dir);
                    TraceCaustics(r, Position, 0);
                }
            }
        }

        getServer().logger.log("Photon mapping complete. Total photons: " + std::to_string(Caustics_Photons.size()));
    }

    // static int num = 0;
    // static int num1 = 0;
     //static int num2 = 0;
    // static int num3 = 0;
    //  用于追踪随机发射的光子
    void PhotonMappingRenderer::TracePhoton(const Ray& r, const RGB& power, unsigned depth)
    {
        // getServer().logger.log("Current Depth is " + to_string(depth) + "/" + to_string(this->depth) + "\n");
        if (depth > this->depth) { return; }
        // cout << "Total is " << num++ << endl;
        //  找到最近的hitobject   ->  后续可以改进对是否命中的判断，添加包围
        HitRecord hitrecord = closestHitObject(r);
        // 如果没有hit，那么返回
        if (!hitrecord)
        {
            /*cout << " Current Ray is " << r.direction << endl;
            cout << "Not Hit !!!???" << num << endl;*/
            return;
        }
        /*if (depth > 1)
        {
            cout << "her" << endl;
        }*/
        // 分别得到hit点、法向量、材质
        const auto& hitpoint = hitrecord->hitPoint;
        const auto& normal   = hitrecord->normal;
        const auto& material = hitrecord->material;

        // 得到打在hitobject上后的光线
        const auto scatter = shaderPrograms[material.index()]->shade(r, hitpoint, normal);

        // 使用亮度作为轮盘赌的参数
        // float L = 0.2126 * power.r + 0.7152 * power.g + 0.0722 * power.b;
        // 不知道为何，这里的RGB并不是0-1(?)
        // 所以使用亮度不太合适了
        // 这里考虑结合法线角度和路径长度来作为依据
        // 再加上一个depth / maxdepth
        auto  ndoti = glm::dot(hitrecord->normal, r.direction);
        float p     = 1.f - 0.5 * (ndoti > 0.0f ? ndoti : 0) - 0.5 * static_cast<float>(depth) / this->depth;
        // float p = 0.9;
        //  cout << p << endl;
        /*cout << spScene->materials[material.index()].type << endl;*/
        //  接下来根据材质进行判断
        if (spScene->materials[material.index()].type == 0)  // 表示漫反射
        {
            // cout << " Lab is " << num1++ << endl;
            //  漫反射需要记录光子
            photon Photon(hitpoint, power, r, scatter.ray);
            // cout << hitpoint << endl;
            Photons.push_back(Photon);

            // 根据轮盘赌策略计算是否继续
            if (Russian_Roulette(p))
            {
                // 计算新的ray和power
                // 漫反射的光强与表面法线相关
                auto cos_thera = (glm::abs(glm::dot(hitrecord->normal, -r.direction)));
                RGB  newpower  = power * scatter.attenuation * cos_thera / (scatter.pdf * p);

                auto newray = scatter.ray;
                // cout << "attenuation is " << scatter.attenuation << endl;
                // cout << "cos_thera is " << cos_thera << endl;
                // cout << "pdf is " << scatter.pdf << endl;
                // cout << "newpower is : R " << newpower.r << " G " << newpower.g << " B " << newpower.b << endl;
                TracePhoton(newray, newpower, depth + 1);
            }
        }
        
        // 如果是镜面的话，那么只有反射没有散射
        if (spScene->materials[material.index()].type == 3)  // 表示镜面反射(即导体)
        {
             //cout << num2++ << endl;

            if (Russian_Roulette(p))
            {
                // 同样的计算新的ray和power
                auto cos_thera = (glm::abs(glm::dot(hitrecord->normal, -r.direction)));

                RGB newpower = power * scatter.attenuation * cos_thera / (scatter.pdf * p);

                // 对于导体来说，应该是需要计算反射方向的
                // glm提供了直接计算反射的()
                Vec3 reflectDir = glm::reflect(r.direction, hitrecord->normal);
                Ray  reflectRay(hitpoint, reflectDir);
                TracePhoton(reflectRay, newpower, depth + 1);
            }
        }

        // 如果是电介质的话，可能会同时存在折射和散射
        if (spScene->materials[material.index()].type == 2)  // 表示电介质->玻璃
        {
            // cout << num3++ << endl;

            // 采用俄罗斯轮盘赌决定处理散射还是折射
            if (Russian_Roulette(p))
            {
                if (scatter.has_refraction)
                {
                    // 处理折射
                    auto r_cos_thera = (glm::abs(glm::dot(hitrecord->normal, -scatter.r_ray.direction)));

                    RGB r_newpower = power * scatter.r_attenuation * r_cos_thera / scatter.r_pdf;
                    //cout << r_newpower << endl;
                    TracePhoton(scatter.r_ray, r_newpower, depth + 1);
                }


                // 处理散射
                auto cos_thera = (glm::abs(glm::dot(hitrecord->normal, -r.direction)));
                auto newpower       = power * scatter.attenuation * cos_thera / scatter.pdf;
                //cout << newpower << endl;

                TracePhoton(scatter.ray, newpower, depth + 1);
            }
        }

        // getServer().logger.log("Current size of Photons is " + to_string(Photons.size()) + "\n");
    }

    tuple<Vec3, Vec3> PhotonMappingRenderer::sampleOnLight(const AreaLight& light_source)
    {
        const Vec3& light_pos = light_source.position;
        const Vec3& light_u   = light_source.u;
        const Vec3& light_v   = light_source.v;

        Vec3 light_normal = glm::normalize(glm::cross(light_u, light_v));
        Vec3 sample_point = light_pos + light_u * RandomFloat() + light_v * RandomFloat();

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

    void PhotonMappingRenderer::TraceCaustics(const Ray& r, const RGB& power, unsigned depth)
    {
        // getServer().logger.log("Current Depth is " + to_string(depth) + "/" + to_string(this->depth) + "\n");
        if (depth > this->depth) { return; }
        // cout << "Total is " << num++ << endl;
        //  找到最近的hitobject   ->  后续可以改进对是否命中的判断，添加包围
        HitRecord hitrecord = closestHitObject(r);
        // 如果没有hit，那么返回
        if (!hitrecord)
        {
            /*cout << " Current Ray is " << r.direction << endl;
            cout << "Not Hit !!!???" << num << endl;*/
            return;
        }
        /*if (depth > 1)
        {
            cout << "her" << endl;
        }*/
        // 分别得到hit点、法向量、材质
        const auto& hitpoint = hitrecord->hitPoint;
        const auto& normal   = hitrecord->normal;
        const auto& material = hitrecord->material;

        // 得到打在hitobject上后的光线
        const auto scatter = shaderPrograms[material.index()]->shade(r, hitpoint, normal);

        // 使用亮度作为轮盘赌的参数
        // float L = 0.2126 * power.r + 0.7152 * power.g + 0.0722 * power.b;
        // 不知道为何，这里的RGB并不是0-1(?)
        // 所以使用亮度不太合适了
        // 这里考虑结合法线角度和路径长度来作为依据
        // 再加上一个depth / maxdepth
        auto  ndoti = glm::dot(hitrecord->normal, r.direction);
        float p     = 1.f - 0.5 * (ndoti > 0.0f ? ndoti : 0) - 0.5 * static_cast<float>(depth) / this->depth;
        // float p = 0.9;
        //  cout << p << endl;
        /*cout << spScene->materials[material.index()].type << endl;*/
        //  接下来根据材质进行判断
        if (spScene->materials[material.index()].type == 0)  // 表示漫反射
        {
            //  漫反射需要记录光子
            // 在焦散处理时，我们只需要记录高质量的光子
            // 记录反射一次后的
            if(Russian_Roulette(0.1))
            {
                if (depth > 0)
                {
                    photon Photon(hitpoint, power, r, scatter.ray);
                    // cout << hitpoint << endl;
                    Caustics_Photons.push_back(Photon);
                    return;
                }
            }
        }

        // 如果是镜面的话，那么只有反射没有折射
        if (spScene->materials[material.index()].type == 3)  // 表示镜面反射(即导体)
        {
            // cout << num2++ << endl;

            if (Russian_Roulette(p))
            {
                // 同样的计算新的ray和power
                auto cos_thera = (glm::abs(glm::dot(hitrecord->normal, -r.direction)));

                RGB newpower = power * scatter.attenuation * cos_thera / (scatter.pdf * p);

                // 对于导体来说，应该是需要计算反射方向的
                // glm提供了直接计算反射的()
                Vec3 reflectDir = glm::reflect(r.direction, hitrecord->normal);
                Ray  reflectRay(hitpoint, reflectDir);
                TraceCaustics(reflectRay, newpower, depth + 1);
            }
        }

        // 如果是电介质的话，可能会同时存在折射和散射
        if (spScene->materials[material.index()].type == 2)  // 表示电介质->玻璃
        {
            // cout << num3++ << endl;

            // 采用俄罗斯轮盘赌决定处理散射还是折射
            if (Russian_Roulette(p))
            {
                if (scatter.has_refraction)
                {
                    // 处理折射
                    auto r_cos_thera = (glm::abs(glm::dot(hitrecord->normal, -scatter.r_ray.direction)));

                    RGB r_newpower = power * scatter.r_attenuation * r_cos_thera / scatter.r_pdf;
                    // cout << r_newpower << endl;
                    TraceCaustics(scatter.r_ray, r_newpower, depth + 1);
                }

                // 处理散射
                auto cos_thera = (glm::abs(glm::dot(hitrecord->normal, -r.direction)));
                auto newpower  = power * scatter.attenuation * cos_thera / scatter.pdf;
                // cout << newpower << endl;

                TraceCaustics(scatter.ray, newpower, depth + 1);
            }
        }

        // getServer().logger.log("Current size of Photons is " + to_string(Photons.size()) + "\n");
    }
    
    // 这个函数是简单的暴力搜索，仅用于性能对比
    vector<photon> PhotonMappingRenderer::getkNearestPhotons(const Vec3& pos, size_t k)
    {
        vector<pair<float, photon>> distPhotons;
        distPhotons.reserve(Photons.size());

        for (const auto& p : Photons)
        {
            float dist = glm::distance(pos, p.GetPosition());
            distPhotons.emplace_back(dist, p);
        }

        partial_sort(distPhotons.begin(),
            distPhotons.begin() + min(k, distPhotons.size()),
            distPhotons.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

        vector<photon> result;
        result.reserve(min(k, distPhotons.size()));

        for (size_t i = 0; i < min(k, distPhotons.size()); i++) result.push_back(distPhotons[i].second);

        return result;
    }
}  // namespace PhotonMapping