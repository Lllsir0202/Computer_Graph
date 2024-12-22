#pragma once
#ifndef __BOUNDINGVOLUME_HPP__
#define __BOUNDINGVOLUME_HPP__

#include "geometry/vec.hpp"
#include "scene/Model.hpp"
#include "Ray.hpp"

namespace PhotonMapping
{
    class BoundingVolume
    {
      protected:
        using T = Node::Type;

      public:
        virtual ~BoundingVolume() = default;

      public:
        virtual bool intersect(const Ray&) const = 0;
        virtual Vec3 center() const                            = 0;
    };

    class AABB : public BoundingVolume
    {
      private:
        Entity* entity;

        Vec3 minp;
        Vec3 maxp;
        bool isLeaf;
        T    type;

        AABB* left;
        AABB* right;

      public:
        AABB();
        AABB(AABB*, AABB*);
        AABB(Sphere*);
        AABB(Triangle*);
        AABB(Plane*);
        ~AABB();

      public:
        virtual bool intersect(const Ray&) const;
        virtual Vec3 center() const;
    };
}  // namespace PhotonMapping

#endif