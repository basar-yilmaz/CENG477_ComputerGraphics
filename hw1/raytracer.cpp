#include <iostream>
#include <math.h>
#include "parser.h"
#include "ppm.h"

using namespace parser;
using namespace std;

typedef unsigned char RGB[3];

// class Ray {
//     private:
//         Vec3f origin;
//         Vec3f direction;
// }

typedef struct
{
    Vec3f origin;
    Vec3f direction;

} Ray;

/*
 * Generate a ray from the camera
 * @param camera: the camera object (position(e), gaze(-w), up(u), near_plane(Vec4f), near_distance, image_width, image_height, image_name)
 * @param i: the i-th pixel
 * @param j: the j-th pixel
 * @return: a ray
 */
Ray generateRay(Camera &camera, int i, int j)
{
    Ray ray;

    // all the abbreviations are from the lecture slides s.t. l, r, s...

    // Calculate the camera's local coordinate system (unit vectors)
    Vec3f gaze = normalizeVector(camera.gaze); // -w
    Vec3f up = normalizeVector(camera.up);
    Vec3f right = normalizeVector(crossProduct(up, gaze)); // gives v of len 1

    float l = camera.near_plane.x; // left
    float r = camera.near_plane.y; // right
    float b = camera.near_plane.z; // bottom
    float t = camera.near_plane.w; // top

    Vec3f m = camera.position scalarMulti(gaze, camera.near_distance); // m = e + w * d
    Vec3f q = m + scalarMulti(right, l) + scalarMulti(up, t);            // q = m + r * u + t * v
    

    // Set the ray's origin and direction
    ray.origin = camera.position;
    ray.direction = normalizeVector(ray_direction);

    return ray;
}

/** TODO
 * can be written a non returning function
 */
float dotProduct(const Vec3f &a, const Vec3f &b)
{
    float result = 0;
    result += a.x * b.x;
    result += a.y * b.y;
    result += a.z * b.z;
    return result;
}

/** TODO
 * can be written a non returning function
 */
Vec3f scalarMulti(const Vec3f &a, const float k)
{
    Vec3f result;
    result.x *= k;
    result.y *= k;
    result.z *= k;
    return result;
}

Vec3f add(const Vec3f &a, const Vec3f &b)
{
    Vec3f result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}

Vec3f subtract(const Vec3f &a, const Vec3f &b)
{
    Vec3f result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;

    // Check for values close to zero and round them to zero
    const float epsilon = 1e-6;
    if (std::abs(result.x) < epsilon)
        result.x = 0.0f;
    if (std::abs(result.y) < epsilon)
        result.y = 0.0f;
    if (std::abs(result.z) < epsilon)
        result.z = 0.0f;

    return result;
}

// get cross product of two vectors s.t. A X B
Vec3f crossProduct(const Vec3f &a, const Vec3f &b)
{
    Vec3f result;
    result.x = (a.y * b.z) - (a.z * b.y);
    result.y = (a.z * b.x) - (a.x * b.z);
    result.z = (a.x * b.y) - (a.y * b.x);
    return result;
}

float length(const Vec3f &a)
{
    return sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2));
}

// get the normalized vector w/ length = 1
// return a new vector
Vec3f normalizeVector(const Vec3f &a)
{
    Vec3f normalized_vec;
    float len = length(a);
    normalized_vec.x = a.x / len;
    normalized_vec.y = a.y / len;
    normalized_vec.z = a.z / len;
    return normalized_vec;
}

// modify the current vector (more efficient)
void normalizeVector2(Vec3f &a)
{
    float len = length(a);
    if (len != 0.0f)
    {
        float invLen = 1.0f / len;
        a.x *= invLen;
        a.y *= invLen;
        a.z *= invLen;
    }
}

// Find the distance between two vectors
float distance(const Vec3f &a, const Vec3f &b)
{
    float dx = b.x - a.x;
    float dy = b.y - a.y;
    float dz = b.z - a.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

// calculate the determinant of a 3x3 matrix
float determinant(const float coef_matrix[3][3])
{
    float res;
    res =
        (coef_matrix[0][0] * (coef_matrix[1][1] * coef_matrix[2][2] - coef_matrix[1][2] * coef_matrix[2][1])) -
        (coef_matrix[0][1] * (coef_matrix[1][0] * coef_matrix[2][2] - coef_matrix[1][2] * coef_matrix[2][0])) +
        (coef_matrix[0][2] * (coef_matrix[1][0] * coef_matrix[2][1] - coef_matrix[1][1] * coef_matrix[2][0]));
    return res;
}

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_scene_file.xml>" << std::endl;
        return 1;
    }

    // Sample usage for reading an XML scene file
    parser::Scene scene;

    try
    {
        scene.loadFromXml(argv[1]);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error loading the scene from XML: " << e.what() << std::endl;
        return 1;
    }

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

    const RGB BAR_COLOR[8] =
        {
            {255, 255, 255}, // 100% White
            {255, 255, 0},   // Yellow
            {0, 255, 255},   // Cyan
            {0, 255, 0},     // Green
            {255, 0, 255},   // Magenta
            {255, 0, 0},     // Red
            {0, 0, 255},     // Blue
            {0, 0, 0},       // Black
        };

    int width = 640, height = 480;
    int columnWidth = width / 8;

    unsigned char *image = new unsigned char[width * height * 3];

    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int colIdx = x / columnWidth;
            image[i++] = BAR_COLOR[colIdx][0];
            image[i++] = BAR_COLOR[colIdx][1];
            image[i++] = BAR_COLOR[colIdx][2];
        }
    }

    write_ppm("test.ppm", image, width, height);

    return 0;
}
