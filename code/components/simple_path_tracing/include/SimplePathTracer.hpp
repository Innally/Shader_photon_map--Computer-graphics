#pragma once
#ifndef __SIMPLE_PATH_TRACER_HPP__
#define __SIMPLE_PATH_TRACER_HPP__

#include "scene/Scene.hpp"
#include "Ray.hpp"
#include "Camera.hpp"
#include "intersections/HitRecord.hpp"

#include "shaders/ShaderCreator.hpp"
#include<map>
#include <tuple>
namespace SimplePathTracer
{
    enum reflectionType{
        diffusion = 0,
        specular,
        absorb
    };
    
    class kdtree {
    public:
        kdtree *left, *right,*parent;
        Photon val;
        kdtree() {
            left = nullptr;
            right = nullptr;
            parent = nullptr;
        }
        kdtree(kdtree* p) :parent(p) {
            left = nullptr;
            right = nullptr;
            parent = nullptr;   
        }
        ~kdtree() { if (left != nullptr) delete left; if (right != nullptr)delete right; }
    };
    

    class photonMap {
    public:
        int p_num;//the number of photon
        int max;
        Photon* p_list;//the array of photon
        kdtree *kd;

    public:
        // functions
        photonMap(int maxnum) :max(maxnum) { p_list = new Photon[maxnum]; p_num = 0; kd = new kdtree; };
        photonMap() { max = 5; p_list = new Photon[max]; p_num = 0; kd = new kdtree; };
        void add_p(Photon* p);// to add photon into my photon list
        //void splitByMed(Photon* t,int start,int end,int med,int axis);
        void balance(Photon* tempPhoton, int start, int end,kdtree* root);
        int calMedian(int, int);
        int findDim(Photon*,int,int);
        void heapSort(Photon* plist, int start, int end, int axis);
        ~photonMap() { delete[]p_list; delete kd; };

    };

    class nearestNeighbor
    {
    public:
        Vec3 pos;
        int max_photons, found, sample; //max number, number already found
        float l2dis; // l2 distance, don't use sqrt to save time
        //Photon photons[100];// to store 
        float diskey[50];// to sort and store the distances
        map<float, Photon> dict;
        nearestNeighbor() { };
        nearestNeighbor(Vec3 v)
        {
            sample = 6;
            pos = v;
            found = 0;
            max_photons =50;
            l2dis = 500;//500
        }

        ~nearestNeighbor() {
        }
        void findNeighbor(nearestNeighbor* n, kdtree* kd);
        void getNearest(photonMap* phtmap);  //使用kdtree查找
        void getNeighborwithoutKDtree(photonMap *phtmap);// 不使用kdttee查找
    };


    using namespace NRenderer;
    using namespace std;

    class SimplePathTracerRenderer
    {
    private:
        default_random_engine eng;
    private:
        SharedScene spScene;
        Scene& scene;

        unsigned int width;
        unsigned int height;
        unsigned int depth;
        unsigned int samples;

        photonMap* phtmap;

        using SCam = SimplePathTracer::Camera;
        SCam camera;

        vector<SharedShader> shaderPrograms;
    public:
        SimplePathTracerRenderer(SharedScene spScene)
            : spScene               (spScene)
            , scene                 (*spScene)
            , camera                (spScene->camera)
        {
            width = scene.renderOption.width;
            height = scene.renderOption.height;
            depth = scene.renderOption.depth;
            samples = scene.renderOption.samplesPerPixel;
            phtmap = new photonMap(10000);//这里是初始化photonMap
        }
        ~SimplePathTracerRenderer() = default;

        using RenderResult = tuple<RGBA*, unsigned int, unsigned int>;
        RenderResult render();
        void release(const RenderResult& r);

    private:
        void renderTask(RGBA* pixels, int width, int height, int off, int step);

        RGB gamma(const RGB& rgb);
        RGB trace(const Ray& ray, int currDepth);
        RGB PMtrace(const Ray& r, int currDepth);
        HitRecord closestHitObject(const Ray& r);
        HitRecord closestHitObject(const Photon& p);
        tuple<float, Vec3> closestHitLight(const Ray& r);
        void photontrace();
    public:
        void generatePhoton(AreaLight a, Vec3& ori, Vec3& dir);
        void generatePhotonMap();
        void tracePhoton(Photon& p ,int depth);
        Vec3 getIrradiance(Ray r);
        int russian(float dif,float spec);

    };
}

#endif