#ifndef __KDTREE_TPP__
#define __KDTREE_TPP__

namespace PhotonMapping
{
    template <typename Point, int K>
    template <typename Iterator, typename>
    KDTree<Point, K>::KDTree(Iterator begin, Iterator end)
    {
        std::vector<Point> pts(begin, end);
        if (!pts.empty()) root = buildTree(pts, 0);
    }

    template <typename Point, int K>
    template <typename Target>
    float KDTree<Point, K>::distance(const Point& a, const Target& b) const
    {
        float dist_sq = 0.0f;
        for (int i = 0; i < K; ++i)
        {
            float d = a[i] - b[i];
            dist_sq += d * d;
        }
        return std::sqrt(dist_sq);
    }

    template <typename Point, int K>
    float KDTree<Point, K>::variance(const std::vector<Point>& points, int dim) const
    {
        if (points.empty()) return 0.0f;

        float sum    = 0.0f;
        float sum_sq = 0.0f;
        for (const auto& point : points)
        {
            float val = point[dim];
            sum += val;
            sum_sq += val * val;
        }
        float mean = sum / points.size();
        float var  = (sum_sq / points.size()) - (mean * mean);
        return var;
    }

    template <typename Point, int K>
    int KDTree<Point, K>::chooseSplit(const std::vector<Point>& points) const
    {
        float max_variance = -1.0f;
        int   best_axis    = 0;
        for (int i = 0; i < K; ++i)
        {
            float var = variance(points, i);
            if (var > max_variance)
            {
                max_variance = var;
                best_axis    = i;
            }
        }
        return best_axis;
    }

    template <typename Point, int K>
    typename KDTree<Point, K>::Node* KDTree<Point, K>::buildTree(std::vector<Point>& points, int depth)
    {
        if (points.empty()) return nullptr;
        int axis = chooseSplit(points);

        std::nth_element(
            points.begin(), points.begin() + points.size() / 2, points.end(), [axis](const Point& a, const Point& b) {
                return a[axis] < b[axis];
            });

        size_t median_idx = points.size() / 2;
        Node*  node       = new Node(points[median_idx], axis);

        std::vector<Point> left_points(points.begin(), points.begin() + median_idx);
        std::vector<Point> right_points(points.begin() + median_idx + 1, points.end());

        node->left  = buildTree(left_points, depth + 1);
        node->right = buildTree(right_points, depth + 1);

        return node;
    }

    template <typename Point, int K>
    void KDTree<Point, K>::insertRecursive(Node* node, const Point& point)
    {
        if (!node) return;

        int axis = node->axis;

        if (point[axis] < node->point[axis])
        {
            if (node->left)
                insertRecursive(node->left, point);
            else
                node->left = new Node(point, (axis + 1) % K);
        }
        else
        {
            if (node->right)
                insertRecursive(node->right, point);
            else
                node->right = new Node(point, (axis + 1) % K);
        }
    }

    template <typename Point, int K>
    template <typename Target>
    void KDTree<Point, K>::searchKNearest(
        const Node* node, const Target& target, int k, std::priority_queue<HeapItem>& heap) const
    {
        if (!node) return;

        float dist = distance(node->point, target);
        if (heap.size() < static_cast<size_t>(k))
            heap.push(HeapItem{dist, node->point});
        else if (dist < heap.top().distance)
        {
            heap.pop();
            heap.push(HeapItem{dist, node->point});
        }

        int         axis = node->axis;
        float       diff = target[axis] - node->point[axis];
        const Node* near = diff < 0 ? node->left : node->right;
        const Node* far  = diff < 0 ? node->right : node->left;

        searchKNearest(near, target, k, heap);

        if (heap.size() < static_cast<size_t>(k) || std::abs(diff) < heap.top().distance)
            searchKNearest(far, target, k, heap);
    }

    template <typename Point, int K>
    template <typename Target>
    void KDTree<Point, K>::searchWithinRadius(
        const Node* node, const Target& target, float r, std::vector<Point>& result) const
    {
        if (!node) return;

        float dist = distance(node->point, target);
        if (dist <= r) { result.emplace_back(node->point); }

        int         axis = node->axis;
        float       diff = target[axis] - node->point[axis];
        const Node* near = diff < 0 ? node->left : node->right;
        const Node* far  = diff < 0 ? node->right : node->left;

        searchWithinRadius(near, target, r, result);

        if (std::abs(diff) <= r) { searchWithinRadius(far, target, r, result); }
    }

    template <typename Point, int K>
    void KDTree<Point, K>::insert(const Point& point)
    {
        if (!root)
        {
            root = new Node(point, 0);
            return;
        }

        insertRecursive(root, point);
    }

    template <typename Point, int K>
    template <typename Target>
    std::vector<Point> KDTree<Point, K>::kNearest(const Target& target, int k) const
    {
        std::priority_queue<HeapItem> heap;
        searchKNearest(root, target, k, heap);

        std::vector<Point> result;
        result.reserve(heap.size());
        while (!heap.empty())
        {
            result.emplace_back(heap.top().point);
            heap.pop();
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

    template <typename Point, int K>
    template <typename Target>
    std::vector<Point> KDTree<Point, K>::withinRadius(const Target& target, float r) const
    {
        std::vector<Point> result;
        searchWithinRadius(root, target, r, result);
        return result;
    }

    template <typename Point, int K>
    template <typename Iterator>
    void KDTree<Point, K>::insert(Iterator begin, Iterator end)
    {
        for (auto it = begin; it != end; ++it) insert(*it);
    }
}  // namespace PhotonMapping

#endif