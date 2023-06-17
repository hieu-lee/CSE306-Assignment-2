#include <vector>
#include "vector.cpp"

class Polygon
{
protected:
    std::vector<Vector> the_vertices;

public:
    Polygon(){};
    Polygon(std::vector<Vector> &a_vertices)
    {
        the_vertices = a_vertices;
    }

    Vector get_centroid() const
    {
        Vector my_centroid(0, 0);
        const int N = the_vertices.size();
        for (int i = 0; i < N; i++)
        {
            my_centroid += the_vertices[i];
        }
        return my_centroid / N;
    }

    size_t get_n_vertices() const
    {
        return the_vertices.size();
    }

    Vector get_vertice(size_t a_index) const
    {
        return the_vertices[a_index];
    }

    void clip(Polygon a_clip_poly)
    {
        Polygon my_out_poly;
        for (int i = 0; i < a_clip_poly.the_vertices.size(); i++)
        {
            std::pair<Vector, Vector> my_clip_edge(a_clip_poly.the_vertices[i], a_clip_poly.the_vertices[(i != 0) ? i - 1 : a_clip_poly.the_vertices.size() - 1]);
            my_out_poly = Polygon();
            for (size_t i = 0; i < the_vertices.size(); i++)
            {
                Vector my_cur_vertex = the_vertices[i];
                Vector my_prev_vertex = the_vertices[(i > 0) ? (i - 1) : the_vertices.size() - 1];
                Vector my_intersection = intersect(my_prev_vertex, my_cur_vertex, my_clip_edge);
                if (inside(my_cur_vertex, my_clip_edge))
                {
                    if (!inside(my_prev_vertex, my_clip_edge))
                        my_out_poly.the_vertices.push_back(my_intersection);
                    my_out_poly.the_vertices.push_back(my_cur_vertex);
                }
                else if (inside(my_prev_vertex, my_clip_edge))
                    my_out_poly.the_vertices.push_back(my_intersection);
            }
            *this = my_out_poly;
        }
    };

    double get_area() const
    {
        double my_area = 0;
        if (the_vertices.size() > 0)
        {
            Vector o = the_vertices[0];
            for (int i = 1; i + 1 < the_vertices.size(); i++)
            {
                Vector a = the_vertices[i] - o;
                Vector b = the_vertices[i + 1] - o;
                my_area += cross(a, b);
            }
        }

        return my_area / 2;
    }
};

class Triangle : public Polygon
{
public:
    explicit Triangle(std::vector<Vector> a_vertices = {}) : Polygon(a_vertices) {}

    bool in_circle(Vector Pi) const
    {
        Vector u = the_vertices[1] - the_vertices[0];
        Vector v = the_vertices[2] - the_vertices[0];
        Vector M = (the_vertices[0] + the_vertices[1]) / 2;
        Vector N = (the_vertices[0] + the_vertices[2]) / 2;
        Vector up = Vector(u[1], -u[0]);
        Vector vp = Vector(v[1], -v[0]);

        Vector K = (dot(u, M) * vp - up * dot(v, N)) / (u[0] * v[1] - u[1] * v[0]);
        double r = abs(sqrt((pow(K[0] - the_vertices[0][0], 2)) + (pow(K[1] - the_vertices[0][1], 2))));

        return (abs(sqrt(pow(Pi[0] - K[0], 2) + (pow(Pi[1] - K[1], 2)))) <= r);
    }
};