#pragma once
#include <string>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include "vector.cpp"

class Ray
{
public:
    Vector origin;
    Vector direction;

    explicit Ray(Vector O, Vector d)
    {
        origin = O;
        direction = normalization(d);
    }
};

class Intersection
{
public:
    bool is_intersection = 0;
    Vector position;
    Vector normal;
    double distance;
    int index;

    Intersection() {}

    Intersection(bool b, Vector P, Vector N, double d, int i)
    {
        is_intersection = b;
        position = P;
        normal = N;
        distance = d;
        index = i;
    }
};

class Geometry
{
public:
    Vector albedo = Vector(1, 1, 1);
    int index;
    bool mirror = false;
    bool transparent = false;
    double refract_index;
    virtual Intersection intersect(const Ray &r) const = 0;
};
