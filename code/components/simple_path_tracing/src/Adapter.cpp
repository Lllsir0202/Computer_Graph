#include "server/Server.hpp"
#include "scene/Scene.hpp"
#include "component/RenderComponent.hpp"
#include "Camera.hpp"

#include "SimplePathTracer.hpp"

using namespace std;
using namespace NRenderer;

namespace SimplePathTracer
{
    class Adapter : public RenderComponent
    {
        void render(SharedScene spScene)
        {
            SimplePathTracerRenderer renderer{spScene};
            auto                     renderResult = renderer.render();
            auto [pixels, width, height]          = renderResult;
            getServer().screen.set(pixels, width, height);
            renderer.release(renderResult);
        }
    };
}  // namespace SimplePathTracer

const static string description = "A Simple Path Tracer. "
                                  "Optimized by BVH and support conduct and dieletronic"
                                  "\n";

REGISTER_RENDERER(OptPathTracer, description, SimplePathTracer::Adapter);