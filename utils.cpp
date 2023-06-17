#include "svg.cpp"

std::vector<Polygon> compute_vor(std::vector<Vector> &a_points, std::vector<double> &a_weights)
{
    if (a_weights.size() == 0)
    {
        a_weights = std::vector<double>(a_points.size());
    }
    std::vector<Polygon> my_polygons;
    for (int i = 0; i < a_points.size(); i++)
    {
        std::vector<Vector> my_vertices = {Vector(1, 0), Vector(1, 1), Vector(0, 1), Vector(0, 0)};
        Polygon a_polygon = Polygon(my_vertices);
        for (int j = 0; j < a_points.size(); j++)
        {
            if (j != i)
            {
                Vector Pi = a_points[i];
                Vector Pj = a_points[j];
                Vector o = (Pi + Pj) / 2 + (a_weights[i] - a_weights[j]) / (2 * pow(norm(Pj - Pi), 2)) * (Pj - Pi);
                Vector d = normalization(Pi - Pj);
                Vector N = Vector(d[1], -d[0]);

                Vector A = o + N * 10;
                Vector D = o - N * 10;
                Vector B = A + d * 20;
                Vector my_centroid = D + d * 20;

                my_vertices = {A, D, my_centroid, B};
                Polygon my_square = Polygon(my_vertices);

                a_polygon.clip(my_square);
            }
        }
        my_polygons.push_back(a_polygon);
    }
    return my_polygons;
}

std::vector<double> lbfgs(std::vector<Vector> &a_points, std::vector<double> &a_lambdas)
{
    int n_points = a_points.size();
    std::vector<double> my_weights(n_points);

    const double eps = .1;
    for (int i = 0; i < 20; i++)
    {
        auto my_polygons = compute_vor(a_points, my_weights);
        for (int j = 0; j < n_points; j++)
        {
            my_weights[j] += eps * (-my_polygons[j].get_area() + a_lambdas[j]);
        }
    }

    return my_weights;
}

std::pair<std::vector<Vector>, std::vector<Vector>> gallouet_method(
    std::vector<Vector> &X, std::vector<Vector> &v, std::vector<double> &m)
{
    auto f = [](Vector a_vector)
    {
        return Vector(std::min(1., std::max(0., a_vector[0])), std::min(1., std::max(0., a_vector[1])));
    };
    const int N = X.size();
    std::vector<double> Uniform(N, 1. / N);
    std::vector<double> Vw = lbfgs(X, Uniform);
    std::vector<Polygon> my_polygons = compute_vor(X, Vw);
    std::vector<Vector> vnext(N);
    std::vector<Vector> Xnext(N);
    const double eps = 0.004;
    const double dt = 0.002;
    for (int i = 0; i < N; i++)
    {
        Vector my_centroid = (my_polygons[i].get_n_vertices() > 0) ? my_polygons[i].get_centroid() : f(X[i]);

        Vector Fspring = 1. / pow(eps, 2) * (my_centroid - X[i]);
        Vector F = Fspring + Vector(0, -9.81);
        vnext[i] = v[i] + dt * F / m[i];
        Xnext[i] = X[i] + dt * v[i];
    }
    return {Xnext, vnext};
}