#pragma once

// #include "utils/scenedata.h"
#include "utils/sceneparser.h"
#include "camera/camera.h"

// A class representing a scene to be ray-traced


class RayTraceScene
{
public:
    RayTraceScene(int width, int height, const RenderData &metaData);

    // The getter of the width of the scene
    const int& width() const;

    // The getter of the height of the scene
    const int& height() const;

    // The getter of the global data of the scene
    const SceneGlobalData& getGlobalData() const;

    const std::vector<SceneLightData>& getsceneLights() const;

    const std::vector<RenderShapeData>& getsceneShapes() const;

    // The getter of the shared pointer to the camera instance of the scene
    const Camera& getCamera() const;

private:
    int m_width;
    int m_height;
    SceneGlobalData m_globalData;
    std::vector<SceneLightData> scenelights;
    std::vector<RenderShapeData> sceneshapes;
    Camera camera;
};
