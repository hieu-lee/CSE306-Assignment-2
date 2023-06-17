#include <string>
#include <vector>
#include <algorithm>
#include <utility>

#include "vector.cpp"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

Vector random_direction()
{
    double r1 = ((double)rand() / (RAND_MAX));
    double r2 = ((double)rand() / (RAND_MAX));
    double x = cos(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    double y = sin(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    double z = 1 - 2 * r2;
    return Vector(x, y, z);
}

int main()
{
    const char *my_input_img_filename = "input.png";
    const char *my_model_img_filename = "model.png";
    const size_t n_iters = 1000;

    typedef unsigned char byte;

    int input_W, input_H, my_input_channels;
    byte *input_image = stbi_load(my_input_img_filename, &input_W, &input_H, &my_input_channels, 0);

    int model_W, model_H, my_model_channels;
    byte *model_image = stbi_load(my_model_img_filename, &model_W, &model_H, &my_model_channels, 0);

    size_t n_pixels = input_W * input_H;
    std::vector<std::pair<int, int>> projection_I(n_pixels);
    std::vector<std::pair<int, int>> projection_M(n_pixels);
    Vector my_pixel, my_model_pixel, v;

    for (size_t iter = 0; iter < n_iters; iter++)
    {
        v = random_direction();

        for (size_t i = 0; i < n_pixels; i++)
        {
            byte *I = input_image + my_input_channels * i;
            byte *M = model_image + my_model_channels * i;
            my_pixel = Vector(*I, *(I + 1), *(I + 2));
            my_model_pixel = Vector(*M, *(M + 1), *(M + 2));
            projection_I[i] = std::pair<int, int>(dot(my_pixel, v), i);
            projection_M[i] = std::pair<int, int>(dot(my_model_pixel, v), i);
        }

        std::sort(projection_I.begin(), projection_I.end());
        std::sort(projection_M.begin(), projection_M.end());

        for (size_t i = 0; i < n_pixels; i++)
        {
            int my_perm_index = projection_I[i].second;
            byte *I = input_image + my_input_channels * my_perm_index;
            my_pixel = Vector(*I, *(I + 1), *(I + 2)) + (projection_M[i].first - projection_I[i].first) * v;
            *I = my_pixel[0];
            *(I + 1) = my_pixel[1];
            *(I + 2) = my_pixel[2];
        }
    }

    stbi_write_png("output.png", input_W, input_H, my_input_channels, &input_image[0], 0);

    return 0;
}