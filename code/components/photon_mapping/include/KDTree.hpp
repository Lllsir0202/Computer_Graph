#pragma once
#ifndef __KDTREE_HPP__
#define __KDTREE_HPP__

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <queue>
#include <type_traits>
#include <vector>

namespace PhotonMapping
{
    template <typename Point, int K = 3>
    class KDTree
    {
      public:
        struct Node
        {
            Point point;
            int   axis;
            Node* left;
            Node* right;

            Node(const Point& pt, int ax) : point(pt), axis(ax), left(nullptr), right(nullptr) {}
            ~Node()
            {
                delete left;
                delete right;
            }
        };

      private:
        struct HeapItem
        {
            float distance;
            Point point;

            bool operator<(const HeapItem& other) const { return distance < other.distance; }
        };

      private:
        Node* root;

      public:
        KDTree() : root(nullptr) {}

        template <typename Iterator,
            typename = typename std::enable_if<std::is_convertible<
                typename std::iterator_traits<Iterator>::iterator_category, std::input_iterator_tag>::value>::type>
        KDTree(Iterator begin, Iterator end);

        KDTree(const KDTree&)            = delete;
        KDTree& operator=(const KDTree&) = delete;
        KDTree(KDTree&&)                 = delete;
        KDTree& operator=(KDTree&&)      = delete;

        ~KDTree() { delete root; }

      private:
        template <typename Target>
        float distance(const Point& a, const Target& b) const;

        float variance(const std::vector<Point>& points, int dim) const;
        int   chooseSplit(const std::vector<Point>& points) const;

      private:
        Node* buildTree(std::vector<Point>& points, int depth);
        void  insertRecursive(Node* node, const Point& point);

        template <typename Target>
        void searchKNearest(const Node* node, const Target& target, int k, std::priority_queue<HeapItem>& heap) const;

        template <typename Target>
        void searchWithinRadius(const Node* node, const Target& target, float r, std::vector<Point>& result) const;

      public:
        void insert(const Point& point);

        template <typename Target>
        std::vector<Point> kNearest(const Target& target, int k) const;

        template <typename Target>
        std::vector<Point> withinRadius(const Target& target, float r) const;

        template <typename Iterator>
        void insert(Iterator begin, Iterator end);
    };
}  // namespace PhotonMapping

#include "KDTree.tpp"

#endif