#pragma once
#include <cmath>

class Vector
{
public:
    explicit Vector(double x = 0., double y = 0.)
    {
        coords[0] = x;
        coords[1] = y;
    };
    Vector &operator+=(const Vector &b)
    {
        coords[0] += b[0];
        coords[1] += b[1];
        return *this;
    }
    Vector &operator-=(const Vector &b)
    {
        coords[0] -= b[0];
        coords[1] -= b[1];
        return *this;
    }
    const double &operator[](int i) const { return coords[i]; }
    double &operator[](int i) { return coords[i]; }

private:
    double coords[2];
};

bool operator==(const Vector &a, const Vector &b)
{
    if ((a[0] == b[0]) && (a[1] == b[1]))
    {
        return true;
    }
    return false;
}

Vector operator+(const Vector &a, const Vector &b)
{
    return Vector(a[0] + b[0], a[1] + b[1]);
}

Vector operator-(const Vector &a, const Vector &b)
{
    return Vector(a[0] - b[0], a[1] - b[1]);
}

Vector operator-(const Vector &b)
{
    return Vector(-b[0], -b[1]);
}

double dot(const Vector &a, const Vector &b)
{
    return a[0] * b[0] + a[1] * b[1];
}

double cross(const Vector &a, const Vector &b)
{
    return a[0] * b[1] - a[1] * b[0];
}

double norm(const Vector &a)
{
    return sqrt(pow(a[0], 2) + pow(a[1], 2));
}

Vector operator*(const double b, const Vector &a)
{
    return Vector(a[0] * b, a[1] * b);
}
Vector operator*(const Vector &a, const double b)
{
    return Vector(a[0] * b, a[1] * b);
}
Vector operator*(const Vector &a, const Vector &b)
{
    return Vector(a[0] * b[0], a[1] * b[1]);
}

Vector operator/(const Vector &a, const double t)
{
    return Vector(a[0] / t, a[1] / t);
}

Vector normalization(const Vector &a)
{
    return a / norm(a);
}

Vector intersect(Vector A, Vector B, std::pair<Vector, Vector> line)
{
    // returns the point of intersection between the Edge [A,B] and line (u,v)
    Vector u = line.first;
    Vector v = line.second;
    Vector N = Vector(v[1] - u[1], u[0] - v[0]);
    double t = dot(u - A, N) / dot(B - A, N);
    if (t < 0 || t > 1)
    { // no intersection
        return Vector(INFINITY, INFINITY);
    }
    return A + t * (B - A);
};

bool inside(Vector P, std::pair<Vector, Vector> clipEdge)
{
    Vector u = clipEdge.first;
    Vector v = clipEdge.second;
    Vector N = Vector(v[1] - u[1], u[0] - v[0]); // outwards normal to clipEdge
    if (dot(P - u, N) <= 0)
        return true;
    return false;
}
