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
    // 此部分是之前别的项目中用过的KDTree模板，此处直接搬过来用
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
                if (left) delete left;
                if (right) delete right;
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

        KDTree(const KDTree&) = delete;
        KDTree& operator=(const KDTree&) = delete;
        KDTree(KDTree&&) = delete;
        KDTree& operator=(KDTree&&) = delete;
        ~KDTree()
        {
            if (root) delete root;
        }

    private:
        float distance(const Point& a, const Point& b) const;
        float variance(const std::vector<Point>& points, int dim) const;
        int   chooseSplit(const std::vector<Point>& points) const;

    private:
        Node* buildTree(std::vector<Point>& points, int depth);
        void  insertRecursive(Node* node, const Point& point);
        void  searchKNearest(const Node* node, const Point& target, int k, std::priority_queue<HeapItem>& heap) const;

    public:
        void               insert(const Point& point);
        std::vector<Point> kNearest(const Point& target, int k) const;

        template <typename Iterator>
        void insert(Iterator begin, Iterator end);
    };
}

#include "KDTree.tpp"

#endif