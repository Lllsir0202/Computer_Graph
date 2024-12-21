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
        void render(SharedScene spScene) {
            PhotonMappingRenderer renderer{ spScene };
            auto renderResult = renderer.render();
            auto [pixels, width, height] = renderResult;
            getServer().screen.set(pixels, width, height);
            renderer.release(renderResult);
        }
    };
}

const static string description =
"A Simple PhotonMapping. "
"Only some simple primitives and materials(Lambertian) are supported."
"\nPlease use scene file : cornel_area_light.scn";

REGISTER_RENDERER(Test, description, PhotonMapping::Adapter);