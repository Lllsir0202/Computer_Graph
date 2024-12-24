#include "server/Server.hpp"
#include "scene/Scene.hpp"
#include "component/RenderComponent.hpp"
#include "Camera.hpp"

#include "PhotonMapping.hpp"

using namespace std;
using namespace NRenderer;

namespace PhotonMapping
{
    class Adapter : public RenderComponent
    {
        void render(SharedScene spScene)
        {
            PhotonMappingRenderer renderer{spScene};
            auto                  renderResult = renderer.render();
            auto [pixels, width, height]       = renderResult;
            getServer().screen.set(pixels, width, height);
            renderer.release(renderResult);
        }
    };
}  // namespace PhotonMapping

const static string description = "A Simple PhotonMapper. "
                                  "Support conduct and dieletric as well as lambertian"
                                  "\n";

REGISTER_RENDERER(PhotonMap, description, PhotonMapping::Adapter);