#include <stdexcept>
#include "raytracescene.h"
#include "utils/sceneparser.h"

RayTraceScene::RayTraceScene(int width, int height, const RenderData &metaData):
    m_width(width),
    m_height(height),
    m_globalData(metaData.globalData),
    scenelights(metaData.lights),
    sceneshapes(metaData.shapes),
    camera(metaData.cameraData)
    {}

const int& RayTraceScene::width() const {
    // Optional TODO: implement the getter or make your own design
    return m_width;
}

const int& RayTraceScene::height() const {
    // Optional TODO: implement the getter or make your own design
    return m_height;
}

const SceneGlobalData& RayTraceScene::getGlobalData() const {
    // Optional TODO: implement the getter or make your own design
    return m_globalData;
}

const std::vector<SceneLightData>& RayTraceScene::getsceneLights() const {
    // Optional TODO: implement the getter or make your own design
    return scenelights;
}

const std::vector<RenderShapeData>& RayTraceScene::getsceneShapes() const {
    return sceneshapes;
}

const Camera& RayTraceScene::getCamera() const {
    // Optional TODO: implement the getter or make your own design
    return camera;
}
