#include "BoundingVolume.hpp"
using namespace std;

namespace PhotonMapping
{
    AABB::AABB()
        : entity(nullptr), minp(Vec3()), maxp(Vec3()), isLeaf(false), type(T::SPHERE), left(nullptr), right(nullptr)
    {}
    AABB::~AABB()
    {
        if (left)
        {
            delete left;
            left = nullptr;
        }
        if (right)
        {
            delete right;
            right = nullptr;
        }
    }
    AABB::AABB(AABB* left, AABB* right)
        : entity(nullptr),
          minp{
              min(left->minp[0], right->minp[0]),
              min(left->minp[1], right->minp[1]),
              min(left->minp[2], right->minp[2]),
          },
          maxp{
              max(left->maxp[0], right->maxp[0]),
              max(left->maxp[1], right->maxp[1]),
              max(left->maxp[2], right->maxp[2]),
          },
          isLeaf(false),
          type(T::SPHERE),
          left(left),
          right(right)
    {}
    AABB::AABB(Sphere* sp) : entity(sp), isLeaf(true), type(T::SPHERE), left(nullptr), right(nullptr)
    {
        float radius = sp->radius;
        if (radius < 0) radius = -radius;
        if (radius == 0.f) radius = 0.1f;
        Vec3 r(radius, radius, radius);
        minp = sp->position - r;
        maxp = sp->position + r;
    }
    AABB::AABB(Triangle* tr) : entity(tr), isLeaf(true), type(T::TRIANGLE), left(nullptr), right(nullptr)
    {
        Vec3& v1 = tr->v[0];
        Vec3& v2 = tr->v[1];
        Vec3& v3 = tr->v[2];

        minp = {
            min(v1[0], min(v2[0], v3[0])),
            min(v1[1], min(v2[1], v3[1])),
            min(v1[2], min(v2[2], v3[2])),
        };
        maxp = {
            max(v1[0], max(v2[0], v3[0])),
            max(v1[1], max(v2[1], v3[1])),
            max(v1[2], max(v2[2], v3[2])),
        };

        if (minp[0] == maxp[0]) maxp[0] += 0.1f;
        if (minp[1] == maxp[1]) maxp[1] += 0.1f;
        if (minp[1] == maxp[1]) maxp[1] += 0.1f;
    }
    AABB::AABB(Plane* pl) : entity(pl), isLeaf(true), type(T::PLANE), left(nullptr), right(nullptr)
    {
        static constexpr float h = 0.1f;

        // position at corner
        Vec3 c1 = pl->position;
        Vec3 c2 = pl->position + pl->u;
        Vec3 c3 = c2 + pl->v;
        Vec3 c4 = pl->position + pl->v;

        Vec3 tmp = h * pl->normal;
        c1 -= tmp;
        c2 -= tmp;
        c3 += tmp;
        c4 += tmp;

        minp = {min(c1[0], min(c2[0], min(c3[0], c4[0]))),
            min(c1[1], min(c2[1], min(c3[1], c4[1]))),
            min(c1[2], min(c2[2], min(c3[2], c4[2])))};
        maxp = {max(c1[0], max(c2[0], max(c3[0], c4[0]))),
            max(c1[1], max(c2[1], max(c3[1], c4[1]))),
            max(c1[2], max(c2[2], max(c3[2], c4[2])))};
    }

    HitRecord AABB::intersect(const Ray& ray) const
    {
        // slabs
        Vec3  invDir = 1.0f / ray.direction;
        Vec3  tMin   = (minp - ray.origin) * invDir;
        Vec3  tMax   = (maxp - ray.origin) * invDir;
        Vec3  tNear  = min(tMax, tMin);
        Vec3  tFar   = max(tMax, tMin);
        float nearT  = glm::max(tNear[0], glm::max(tNear[1], tNear[2]));
        float farT   = glm::min(tFar[0], glm::min(tFar[1], tFar[2]));
        if (farT - nearT >= 1e-6f && farT >= 0.f) { return getHitRecord(nearT, ray.at(nearT), {}, {}); }
        return getMissRecord();
    }

    Vec3 AABB::center() const
    {
        return Vec3((minp[0] + maxp[0]) / 2, (minp[1] + maxp[1]) / 2, (minp[2] + maxp[2]) / 2);
    }
}  // namespace PhotonMapping