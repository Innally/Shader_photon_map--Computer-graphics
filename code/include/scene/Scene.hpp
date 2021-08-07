#pragma once
#ifndef __NR_SCENE_HPP__
#define __NR_SCENE_HPP__

#include "Texture.hpp"
#include "Material.hpp"
#include "Model.hpp"
#include "Light.hpp"
#include "Camera.hpp"

namespace NRenderer
{
    struct RenderOption
    {
        unsigned int width; // 宽度
        unsigned int height; // 高度
        unsigned int depth; // 深度
        unsigned int samplesPerPixel;  // 采样数量，以上都可以在nrenderer app里面调整
        RenderOption()
            : width             (500)
            , height            (500)
            , depth             (4)
            , samplesPerPixel   (16)
        {}
    };

    struct Ambient // 环境光结构，分成constant和environment map两种，可以在app里选择
    {
        enum class Type
        {
            CONSTANT, ENVIROMENT_MAP
        };
        Type type;
        Vec3 constant = {};
        Handle environmentMap = {};
    };

    struct Scene
    {
        Camera camera; // 相机

        RenderOption renderOption;

        Ambient ambient; //环境光

        // buffers 材料、纹理、模型、nodes的缓冲区
        vector<Material> materials;
        vector<Texture> textures;

        vector<Model> models;
        vector<Node> nodes;
        // object buffer 场景内物体的缓冲区
        vector<Sphere> sphereBuffer;  //球
        vector<Triangle> triangleBuffer; // 三角
        vector<Plane> planeBuffer; // 或许是平面
        vector<Mesh> meshBuffer; // 网格

        vector<Light> lights;
        // light buffer 光缓冲区
        vector<PointLight> pointLightBuffer;
        vector<AreaLight> areaLightBuffer;
        vector<DirectionalLight> directionalLightBuffer;
        vector<SpotLight> spotLightBuffer;
    };
    using SharedScene = shared_ptr<Scene>;
} // namespace NRenderer


#endif