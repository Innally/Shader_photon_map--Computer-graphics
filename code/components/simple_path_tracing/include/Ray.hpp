#pragma once
#ifndef __RAY_HPP__
#define __RAY_HPP__

#include "geometry/vec.hpp"

#include <limits>

#define FLOAT_INF numeric_limits<float>::infinity()
namespace SimplePathTracer
{
    using namespace NRenderer;
    using namespace std;


    struct Photon
    {
        Vec3 origin;
        // keep it as a unit vector
        Vec3 direction;
        Vec3 power;
        int axis;// this indicates the kdtree dimention

        void setOrigin(const Vec3& v) {
            origin = v;
        }

        void setDirection(const Vec3& v) {
            direction = glm::normalize(v);
        }

        void setPow(const Vec3& v) {
            power = v;
        }

        inline
            Vec3 at(float t) const {
            return origin + t * direction;
        }

        Photon(const Vec3& origin, const Vec3& direction,const Vec3& power)
            : origin(origin)
            , direction(direction)
            ,power(power)
        {
            axis = -1;
        }

        Photon()
            : origin{}
            , direction{}
            , power{}
        {
            axis = -1;
        }

        float getAxis(int axis) {
            switch (axis) {
            case 0:
                return origin.x; break;
            case 1:
                return origin.y; break;
            case 2:
                return origin.z; break;
            default:
                return origin.x; break;
            }
        }

    };



    struct Ray
    {
        Vec3 origin;
        // keep it as a unit vector
        Vec3 direction;

        void setOrigin(const Vec3& v) {
            origin = v;
        }

        void setDirection(const Vec3& v) {
            direction = glm::normalize(v);
        }


        inline
        Vec3 at(float t) const {
            return origin + t*direction;
        }


      /*  Ray(const Vec3& origin, const Vec3& direction, const Vec3& pow)
        {
            origin = origin;
            pow = pow;
            direction = direction;
        }*/


        Ray(const Vec3& origin, const Vec3& direction)
            : origin                (origin)
            , direction             (direction)
        {}
    
        Ray()
            : origin        {}
            , direction     {}
        {}
    };
}

#endif