#ifndef PARTICLES_H
#define PARTICLES_H

#include <cmath>

template<typename scalar=double>
struct vec3 {
    scalar x, y, z;
    struct vec3<scalar> operator-(struct vec3<scalar> v)
    {
        return {this->x - v.x, this->y - v.y, this->z - v.z};
    }
    scalar length()
    {
        return sqrt(x*x + y*y + z*z);
    }
    void operator*=(scalar m)
    {
        x *= m;
        y *= m;
        z *= m;
        return;
    }
    vec3<scalar> operator*(double m)
    {
        vec3<scalar> res = *this;
        res *= m;
        return res;
    }
    vec3<scalar> operator*(vec3<scalar> v)
    {
        vec3<scalar> res = *this;
        res.x *= v.x;
        res.y *= v.y;
        res.z *= v.z;
        return res;
    }
    vec3<scalar> operator+(vec3<scalar> v)
    {
        vec3<scalar> res = *this;
        res += v;
        return res;
    }
    void operator+=(vec3<scalar> v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return;
    }
    void operator-=(vec3<scalar> v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return;
    }
    vec3<scalar> operator/(vec3<scalar> v)
    {
        vec3<scalar> res = *this;
        res.x /= v.x;
        res.y /= v.y;
        res.z /= v.z;
    }
    vec3<scalar> operator/(double d)
    {
        vec3<scalar> res = *this;
        res.x /= d;
        res.y /= d;
        res.z /= d;
        return res;
    }
};

template<typename scalar=double>
struct particle {
    vec3<scalar> pos;
    vec3<scalar> vel;
};

#endif