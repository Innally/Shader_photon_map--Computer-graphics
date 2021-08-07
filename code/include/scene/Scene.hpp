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
        unsigned int width; // ���
        unsigned int height; // �߶�
        unsigned int depth; // ���
        unsigned int samplesPerPixel;  // �������������϶�������nrenderer app�������
        RenderOption()
            : width             (500)
            , height            (500)
            , depth             (4)
            , samplesPerPixel   (16)
        {}
    };

    struct Ambient // ������ṹ���ֳ�constant��environment map���֣�������app��ѡ��
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
        Camera camera; // ���

        RenderOption renderOption;

        Ambient ambient; //������

        // buffers ���ϡ�����ģ�͡�nodes�Ļ�����
        vector<Material> materials;
        vector<Texture> textures;

        vector<Model> models;
        vector<Node> nodes;
        // object buffer ����������Ļ�����
        vector<Sphere> sphereBuffer;  //��
        vector<Triangle> triangleBuffer; // ����
        vector<Plane> planeBuffer; // ������ƽ��
        vector<Mesh> meshBuffer; // ����

        vector<Light> lights;
        // light buffer �⻺����
        vector<PointLight> pointLightBuffer;
        vector<AreaLight> areaLightBuffer;
        vector<DirectionalLight> directionalLightBuffer;
        vector<SpotLight> spotLightBuffer;
    };
    using SharedScene = shared_ptr<Scene>;
} // namespace NRenderer


#endif