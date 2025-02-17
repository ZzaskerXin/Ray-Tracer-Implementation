#include "sceneparser.h"
#include "scenefilereader.h"
#include <glm/gtx/transform.hpp>

// #include <chrono>
// #include <iostream>

void traverseScene(SceneNode* node, glm::mat4 ctm, RenderData &renderData) {
    for (SceneTransformation* transform : node->transformations) {
        switch (transform->type) {
        case TransformationType::TRANSFORMATION_TRANSLATE:
            ctm = glm::translate(ctm, transform->translate);
            break;
        case TransformationType::TRANSFORMATION_SCALE:
            ctm = glm::scale(ctm, transform->scale);
            break;
        case TransformationType::TRANSFORMATION_ROTATE:
            ctm = glm::rotate(ctm, transform->angle, transform->rotate);
            break;
        case TransformationType::TRANSFORMATION_MATRIX:
            ctm = ctm * transform->matrix;
            break;
        }
    }

    for (ScenePrimitive* primitive : node->primitives) {
        RenderShapeData shapeData;
        shapeData.primitive = *primitive;
        shapeData.ctm = ctm;
        if (primitive->material.textureMap.isUsed) {
            shapeData.texture = loadImageFromFile(primitive->material.textureMap.filename);
            std::cout << shapeData.texture->height << std::endl;
            std::cout << shapeData.texture->width << std::endl;
            std::cout << shapeData.primitive.material.textureMap.filename << std::endl;
        }
        renderData.shapes.push_back(shapeData);
    }

    for (SceneLight* light : node->lights) {
        SceneLightData lightData;
        lightData.id = light->id;
        lightData.type = light->type;
        lightData.color = light->color;
        lightData.function = light->function;

        lightData.pos = ctm * glm::vec4(0, 0, 0, 1);
        lightData.dir = glm::normalize(ctm * light->dir);

        lightData.penumbra = light->penumbra;
        lightData.angle = light->angle;
        renderData.lights.push_back(lightData);
    }

    for (SceneNode* child : node->children) {
        traverseScene(child, ctm, renderData);
    }
}

bool SceneParser::parse(std::string filepath, RenderData &renderData) {
    ScenefileReader fileReader = ScenefileReader(filepath);
    bool success = fileReader.readJSON();
    if (!success) {
        return false;
    }

    SceneGlobalData global = fileReader.getGlobalData();
    SceneCameraData camera = fileReader.getCameraData();

    renderData.globalData = global;
    renderData.cameraData = camera;

    renderData.shapes.clear();

    SceneNode* rootNode = fileReader.getRootNode();

    glm::mat4 initctm = glm::mat4(1.0f);
    traverseScene(rootNode, initctm, renderData);

    return true;
}
