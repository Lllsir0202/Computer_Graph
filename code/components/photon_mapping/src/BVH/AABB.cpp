#include "BVH.hpp"
#include <iostream>
#include <limits>
#include <algorithm>
#include <stack>
using namespace std;

namespace
{
    constexpr float fInf  = numeric_limits<float>::infinity();
    constexpr float nfInf = -fInf;

    constexpr NRenderer::Vec3 infVec(fInf, fInf, fInf);
    constexpr NRenderer::Vec3 nInfVec(nfInf, nfInf, nfInf);
}  // namespace

namespace PhotonMapping
{
    AABBBVH::AABBBVH() : root(nullptr), objects() {}
    AABBBVH::~AABBBVH() { destory(); }

    HitRecord AABBBVH::_intersect(AABB* node, const Ray& ray, float tMin, float tMax) const
    {
        if (!node) return getMissRecord();

        auto boxHit = node->intersect(ray);
        if (!boxHit) return getMissRecord();
        float nearT = boxHit->t;
        if (nearT > tMax) return getMissRecord();
        tMin = max(tMin, nearT);

        if (node->isLeaf)
        {
            switch (node->type)
            {
                case AABB::T::SPHERE: return Intersection::xSphere(ray, *((Sphere*)node->entity), tMin, tMax);
                case AABB::T::PLANE: return Intersection ::xPlane(ray, *((Plane*)node->entity), tMin, tMax);
                case AABB::T::TRIANGLE: return Intersection ::xTriangle(ray, *((Triangle*)node->entity), tMin, tMax);
                default:
                {
                    cerr << "not here\n";
                    return getMissRecord();
                }
            }
        }

        auto leftHit = _intersect(node->left, ray, tMin, tMax);
        if (leftHit && leftHit->t < tMax) tMax = leftHit->t;
        auto rightHit = _intersect(node->right, ray, tMin, tMax);

        if (!leftHit) return rightHit;
        if (!rightHit) return leftHit;
        return (leftHit->t < rightHit->t) ? leftHit : rightHit;
    }

    bool AABBBVH::build(const SharedScene& s)
    {
        for (const auto& node : s->nodes)
        {
            switch (node.type)
            {
                case AABB::T::SPHERE:
                {
                    if (node.entity >= s->sphereBuffer.size())
                    {
                        cerr << "Error: SPHERE entity index out of range: " << node.entity << "\n";
                        break;
                    }
                    objects.emplace_back(new AABB(&s->sphereBuffer[node.entity]));
                    break;
                }
                case AABB::T::TRIANGLE:
                {
                    if (node.entity >= s->triangleBuffer.size())
                    {
                        cerr << "Error: TRIANGLE entity index out of range: " << node.entity << "\n";
                        break;
                    }
                    objects.emplace_back(new AABB(&s->triangleBuffer[node.entity]));
                    break;
                }
                case AABB::T::PLANE:
                {
                    if (node.entity >= s->planeBuffer.size())
                    {
                        cerr << "Error: PLANE entity index out of range: " << node.entity << "\n";
                        break;
                    }
                    objects.emplace_back(new AABB(&s->planeBuffer[node.entity]));
                    break;
                }
                case AABB::T::MESH:
                {
                    cerr << "MESH type not yet supported in AABB BVH\n";
                    break;
                }
                default:
                {
                    cerr << "Meet a unknown model type with id " << static_cast<int>(node.type) << "\n";
                }
            }
        }

        size_t cnt = objects.size();
        if (cnt == 0) return false;
        cout << "Start build BVH with size " << cnt << "\n";

        struct BuildTask
        {
            size_t begin;
            size_t end;
            AABB** target;
            BuildTask(size_t b, size_t e, AABB** t) : begin(b), end(e), target(t) {}
        };

        stack<BuildTask> buildStack;
        buildStack.push(BuildTask(0, cnt - 1, &root));

        while (!buildStack.empty())
        {
            auto task = buildStack.top();
            buildStack.pop();

            if (task.begin > task.end)
            {
                *task.target = nullptr;
                continue;
            }
            if (task.begin == task.end)
            {
                *task.target = objects[task.begin];
                continue;
            }

            AABB* node   = new AABB();
            *task.target = node;

            Vec3 minp = infVec;
            Vec3 maxp = nInfVec;
            for (size_t i = task.begin; i <= task.end; ++i)
            {
                if (!objects[i]) continue;
                minp = Vec3(min(minp[0], objects[i]->minp[0]),
                    min(minp[1], objects[i]->minp[1]),
                    min(minp[2], objects[i]->minp[2]));
                maxp = Vec3(max(maxp[0], objects[i]->maxp[0]),
                    max(maxp[1], objects[i]->maxp[1]),
                    max(maxp[2], objects[i]->maxp[2]));
            }

            node->minp   = minp;
            node->maxp   = maxp;
            node->isLeaf = false;

            Vec3  s      = maxp - minp;
            int   axis   = (s[1] > s[0]) ? ((s[2] > s[1]) ? 2 : 1) : ((s[2] > s[0]) ? 2 : 0);
            float midVal = (minp[axis] + maxp[axis]) / 2.0f;

            size_t idx = task.begin;
            AABB*  tmp = nullptr;
            for (size_t i = task.begin; i <= task.end; ++i)
            {
                float center = (objects[i]->minp[axis] + objects[i]->maxp[axis]) / 2.0f;
                if (center >= midVal) continue;
                if (i != idx)
                {
                    tmp          = objects[i];
                    objects[i]   = objects[idx];
                    objects[idx] = tmp;
                }
                ++idx;
            }

            if (idx == task.begin || idx == task.end + 1) idx = task.begin + (task.end - task.begin + 1) / 2;

            buildStack.push(BuildTask(idx, task.end, &node->right));
            buildStack.push(BuildTask(task.begin, idx - 1, &node->left));
        }

        if (!root)
        {
            cerr << "BVH root is nullptr after build\n";
            return false;
        }
        cout << "BVH build done with size " << cnt << "\n";
        return true;
    }

    void AABBBVH::destory()
    {
        if (root)
        {
            delete root;
            root = nullptr;
        }
        objects.clear();
    }

    HitRecord AABBBVH::intersect(const Ray& ray, float tMin, float tMax) const
    {
        if (!root) return getMissRecord();
        return _intersect(root, ray, tMin, tMax);
    }
}  // namespace PhotonMapping