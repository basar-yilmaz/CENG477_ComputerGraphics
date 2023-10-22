#include <iostream>
#include <math.h>
#include "parser.h"
#include "ppm.h"

using namespace parser;
using namespace std;

typedef unsigned char RGB[3];

// created obj enum to make it easier to check the type of the object
enum obj
{
    TRIANGLE,
    SPHERE,
    MESH
};

typedef struct
{
    int objId;
    obj objType;                // 0 for triangle, 1 for sphere, 2 for mesh
    Vec3f normal;               // normal of the object
    Vec3f intersectionPoint;    // intersection point of the ray and the object
    bool isIntersected = false; // if the ray intersects with the object
    int material_id;            // material id of the object we get this from xml file
    double t = -1.0f;           // distance between the ray's origin and the intersection point
} Intersection;

typedef struct
{
    Vec3f origin;
    Vec3f direction;

} Ray;

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

// get the normalized vector with length = 1
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

// modify the current vector
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

    // all the abbreviations are from the lecture slides such that l, r, s...

    // Calculate the camera's local coordinate system (unit vectors)
    // we can optimize this part to avoid calculating the same thing again and again
    Vec3f gaze = normalizeVector(camera.gaze); // -w
    Vec3f u = normalizeVector(camera.up);
    Vec3f v = normalizeVector(crossProduct(u, gaze)); // gives v of len 1

    float l = camera.near_plane.x; // left
    float r = camera.near_plane.y; // right
    float b = camera.near_plane.z; // bottom
    float t = camera.near_plane.w; // top

    Vec3f m = add(camera.position, scalarMulti(gaze, camera.near_distance)); // m = e - w * d or e + gaze * d
    Vec3f q = add(m, add(scalarMulti(u, r), scalarMulti(v, t)));             // q = m + r * u + t * v

    float su = (i + 0.5) * (r - l) / camera.image_width;  // su = (i+0.5) * (r-l) / image_width
    float sv = (j + 0.5) * (t - b) / camera.image_height; // sv = (j+0.5) * (t-b) / image_height

    Vec3f s = add(q, add(scalarMulti(u, su), scalarMulti(v, -sv))); // s = q + su * u - sv * v

    // Set the ray's origin and direction
    ray.origin = camera.position;
    ray.direction = normalizeVector(subtract(s, camera.position));

    return ray;
}

/*
 * Calculate the intersection point of a ray and a plane
 * @param ray: the ray
 * @param t: the distance between the ray's origin and the intersection point
 * @return: the intersection point
 */
Vec3f getIntersectionPoint(Ray &ray, double t)
{
    Vec3f intersectionPoint;
    intersectionPoint = add(ray.origin, scalarMulti(ray.direction, t));
    return intersectionPoint;
}

Intersection rayTriangleIntersection(const Scene &scene, const Ray &ray, const Triangle &triangle)
{
}

Intersection raySphereIntersection(const Scene &scene, const Ray &ray, const Sphere &sphere)
{
}

Intersection rayMeshIntersection(const Scene &scene, const Ray &ray, const Mesh &mesh)
{
}

Intersection rayObjectIntersection(const Scene &scene, const Ray &ray)
{
}

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_scene_file.xml>" << std::endl;
        return 1;
    }

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

    for (int cameraIndex = 0; cameraIndex < scene.cameras.size(); cameraIndex++)
    {
        int imageWidth = scene.cameras[cameraIndex].image_width;
        int imageHeight = scene.cameras[cameraIndex].image_height;

        unsigned char *image = new unsigned char[imageWidth * imageHeight * 3];

        // at this point we set up our camera and created space for our image
        // now we need to loop through each pixel and generate a ray
        for (int i = 0; i < imageHeight; i++)
        {
            for (int j = 0; j < imageWidth; j++)
            {
                Ray myRay = generateRay(scene.cameras[cameraIndex], i, j);

                // Perform ray-object intersection tests here

                // Calculate shading and color for the intersection point

                // Write the color to the image buffer
            }
        }

        // read the image name from the camera object and write the image to that file
        write_ppm(scene.cameras[cameraIndex].image_name.c_str(), image, imageWidth, imageHeight);

        delete[] image;
    }

    return 0;
}
