#include <iostream>
#include "utils.cpp"
#include <random>

std::random_device my_random_device;
std::default_random_engine my_random_engine(my_random_device());
std::uniform_real_distribution<float> my_distribution(0, 1);

double generate()
{
    return my_distribution(my_random_engine);
}

void generate_diagram()
{
    const int n_points = 500;
    std::vector<Vector> my_points(n_points);
    for (int i = 0; i < n_points; i++)
    {
        my_points.push_back(Vector(generate(), generate()));
    }

    std::vector<double> my_lambda(n_points);
    const Vector C = Vector(0.5, 0.5);
    for (int i = 0; i < my_points.size(); i++)
    {
        my_lambda[i] = exp(-norm(my_points[i] - C) / 0.02);
    }

    auto w = lbfgs(my_points, my_lambda);
    std::vector<Polygon> my_polygons = compute_vor(my_points, w);
    save_svg(my_polygons, "image.svg");
}

void generate_fluid()
{
    const int WATER_MASS = 200;
    const int AIR_MASS = 100;
    int n_water = 50;
    int n_air = 200;
    int n = n_water + n_air;
    std::vector<Vector> my_points;
    std::vector<double> my_mass;
    for (int i = 0; i < n_water; i++)
    {
        my_points.push_back(Vector(generate(), generate()));
        my_mass.push_back(WATER_MASS);
    }
    for (int i = 0; i < n_air; i++)
    {
        my_points.push_back(Vector(generate(), generate()));
        my_mass.push_back(AIR_MASS);
    }
    std::vector<Vector> my_velocity(n);
    std::vector<double> my_lambda(n);

    const int n_frames = 60;
    for (int a_frame = 0; a_frame < n_frames; a_frame++)
    {
        std::cout << a_frame << std::endl;
        auto [x, v] = gallouet_method(my_points, my_velocity, my_mass);
        my_points = x;
        my_velocity = v;

        auto weights = lbfgs(my_points, my_lambda);
        auto vor = compute_vor(my_points, weights);
        save_svg_animated(vor, "water.svg", a_frame, n_frames, n_water);
    }
}

int main(int argc, char **argv)
{
    // generate_diagram();
    generate_fluid();
}