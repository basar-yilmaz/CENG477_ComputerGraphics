#include <iostream>
#include <math.h>
#include "parser.h"
#include "ppm.h"
#include "raytracer.h"

using namespace parser;
using namespace std;

typedef unsigned char RGB[3];

#define EPSILON 0.0001

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
    bool isIntersected = false; // if the ray intersects with the object assume it is false
    int material_id;            // material id of the object we get this from xml file, probably will be used for reflection part
    float t = -1.0f;            // distance between the ray's origin and the intersection point
} Intersection;

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

    // all the abbreviations are from the lecture slides such that l, r, s...

    // Calculate the camera's local coordinate system (unit vectors)
    // we can optimize this part to avoid calculating the same thing again and again
    Vec3f gaze = normalizeVector(camera.gaze); // -w
    Vec3f u = normalizeVector(crossProduct(gaze, camera.up));
    Vec3f v = crossProduct(u, gaze); // as gaze and u has length 1, v will also have length 1

    float l = camera.near_plane.x; // left
    float r = camera.near_plane.y; // right
    float b = camera.near_plane.z; // bottom
    float t = camera.near_plane.w; // top

    float su = (j + 0.5) * (r - l) / camera.image_width;  // su = (i+0.5) * (r-l) / image_width
    float sv = (i + 0.5) * (t - b) / camera.image_height; // sv = (j+0.5) * (t-b) / image_height

    Vec3f m, q, s;

    m = add(camera.position, scalarMulti(gaze, camera.near_distance)); // m = e - w * d or e + gaze * d

    q = add(m, add(scalarMulti(u, l), scalarMulti(v, t))); // q = m + l * u + t * v

    s = add(q, subtract(scalarMulti(u, su), scalarMulti(v, sv))); // s = q + su * u - sv * v

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
Vec3f getIntersectionPoint(const Ray &ray, float t)
{
    Vec3f intersectionPoint;
    intersectionPoint = add(ray.origin, scalarMulti(ray.direction, t));
    return intersectionPoint;
}

/*
 * Calculate the intersection point of a ray and a triangle
 * @param scene: the scene object
 * @param ray: the ray
 * @param triangle: the triangle object
 * @return: the intersection point
 */
Intersection rayTriangleIntersection(const Scene &scene, const Ray &ray, const Triangle &triangle)
{
    // we will use barycentric coordinates to calculate the intersection point
    Intersection point;

    // we subtract 1 because the vertex ids start from 1
    Vec3f a, b, c;

    a = scene.vertex_data[triangle.indices.v0_id - 1];
    b = scene.vertex_data[triangle.indices.v1_id - 1];
    c = scene.vertex_data[triangle.indices.v2_id - 1];

    // after this point I followed the lecture slides (cramer's rule)

    float matrix[3][3] = {
        {a.x - b.x, a.x - c.x, ray.direction.x},
        {a.y - b.y, a.y - c.y, ray.direction.y},
        {a.z - b.z, a.z - c.z, ray.direction.z}};

    float A = determinant(matrix);

    float beta[3][3] = {
        {a.x - ray.origin.x, a.x - c.x, ray.direction.x},
        {a.y - ray.origin.y, a.y - c.y, ray.direction.y},
        {a.z - ray.origin.z, a.z - c.z, ray.direction.z}};

    float gamma[3][3] = {
        {a.x - b.x, a.x - ray.origin.x, ray.direction.x},
        {a.y - b.y, a.y - ray.origin.y, ray.direction.y},
        {a.z - b.z, a.z - ray.origin.z, ray.direction.z}};

    float t[3][3] = {
        {a.x - b.x, a.x - c.x, a.x - ray.origin.x},
        {a.y - b.y, a.y - c.y, a.y - ray.origin.y},
        {a.z - b.z, a.z - c.z, a.z - ray.origin.z}};

    float beta_val = determinant(beta) / A;
    float gamma_val = determinant(gamma) / A;
    float t_val = determinant(t) / A;

    if (t_val < 0)
        return point;

    // check if the ray intersects with the triangle
    if (beta_val >= 0 && gamma_val >= 0 && beta_val + gamma_val <= 1)
    {
        point.isIntersected = true;
        point.intersectionPoint = getIntersectionPoint(ray, t_val);
        point.t = t_val;
        point.normal = normalizeVector(crossProduct(subtract(b, a), subtract(c, a)));
        point.objId = triangle.indices.v0_id;
        point.objType = TRIANGLE;
        point.material_id = triangle.material_id;
    }

    return point;
}

/*
 * Calculate the intersection point of a ray and a sphere
 * @param scene: the scene object
 * @param ray: the ray
 * @param sphere: the sphere object
 * @return: the intersection point
 */
Intersection raySphereIntersection(const Scene &scene, const Ray &ray, const Sphere &sphere)
{

    // we subtract 1 because the vertex ids start from 1
    Vec3f sphereCenter = scene.vertex_data[sphere.center_vertex_id - 1];

    // formula from the lecture slides ,all same
    // trying to get At^2 + Bt + C = 0

    // since it repeats at least 3 times make it constant
    const Vec3f centerToEye = subtract(ray.origin, sphereCenter); // e-c

    float C = dotProduct(centerToEye, centerToEye) - (sphere.radius * sphere.radius); // C = (e-c).(e-c) - r^2

    float B = 2 * dotProduct(ray.direction, centerToEye); // B = 2d.(e-c)

    float A = dotProduct(ray.direction, ray.direction); // A = d.d

    float discriminant = B * B - 4 * A * C;

    // we got all we need to calculate the intersection point
    Intersection point;

    // there is no real solution
    if (discriminant < 0)
    {
        point.isIntersected = false;
        return point;
    }

    // add this in order to avoid unnecessary calculations
    // not sure if I need to implement this case
    else if (discriminant == 0)
    {
        float sol = -B / (2 * A);
        point.isIntersected = true;
        point.intersectionPoint = getIntersectionPoint(ray, sol);
        point.t = sol;
        point.normal = normalizeVector(subtract(point.intersectionPoint, sphereCenter));
        point.objId = sphere.center_vertex_id;
        point.objType = SPHERE;
        point.material_id = sphere.material_id;
    }

    // two solutions
    else
    {
        float sol1 = (-B + sqrtf(discriminant)) / (2 * A);
        float sol2 = (-B - sqrtf(discriminant)) / (2 * A);

        // we need to check which one is closer to the camera
        float valid_sol = sol1 < sol2 ? sol1 : sol2;

        point.isIntersected = true;
        point.intersectionPoint = getIntersectionPoint(ray, valid_sol);
        point.t = valid_sol;
        point.normal = normalizeVector(subtract(point.intersectionPoint, sphereCenter));
        point.objId = sphere.center_vertex_id;
        point.objType = SPHERE;
        point.material_id = sphere.material_id;
    }

    return point;
}

/*
 * Calculate the intersection point of a ray and a mesh
 * @param scene: the scene object
 * @param ray: the ray
 * @param mesh: the mesh object
 * @return: the intersection point
 */
Intersection rayMeshIntersection(const Scene &scene, const Ray &ray, const Mesh &mesh)
{
    Intersection closestIntersection, triangleIntersection;

    closestIntersection.t = 1000000; // TODO: set this to the max value of float or some very high value temporarily

    for (const Face &face : mesh.faces)
    {
        // Get the vertices of the current triangle, subtract 1 because the vertex ids start from 1
        Vec3f a = scene.vertex_data[face.v0_id - 1];
        Vec3f b = scene.vertex_data[face.v1_id - 1];
        Vec3f c = scene.vertex_data[face.v2_id - 1];

        // Check for intersection with the current triangle
        // declare a new triangle with the material id and face
        triangleIntersection = rayTriangleIntersection(scene, ray, {mesh.material_id, face});

        if (triangleIntersection.isIntersected && triangleIntersection.t < closestIntersection.t)
        {
            closestIntersection = triangleIntersection;
        }
    }

    // Return the closest intersection found within the mesh
    return closestIntersection;
}

/*
 * Calculate the optimal intersection point
 * @param scene: the scene object
 * @param ray: the ray
 * @return: the optimal intersection point
 */

Intersection rayObjectIntersection(const Scene &scene, const Ray &ray)
{

    // we need to loop through all the objects and find the closest intersection point
    Intersection closestIntersection;
    closestIntersection.t = 1000000; // TODO: set this to the max value of float or some very high value temporarily
    int count = 0;
    // sphere intersection
    for (const Sphere &sphere : scene.spheres)
    {
        Intersection sphereIntersection = raySphereIntersection(scene, ray, sphere);

        if (sphereIntersection.isIntersected && sphereIntersection.t < closestIntersection.t)
        {
            closestIntersection = sphereIntersection;
        }
    }

    // triangle intersection
    for (const Triangle &triangle : scene.triangles)
    {
        Intersection triangleIntersection = rayTriangleIntersection(scene, ray, triangle);
        if (triangleIntersection.isIntersected && triangleIntersection.t < closestIntersection.t)
        {
            closestIntersection = triangleIntersection;
        }
    }

    // mesh intersection
    for (const Mesh &mesh : scene.meshes)
    {
        Intersection meshIntersection = rayMeshIntersection(scene, ray, mesh);

        if (meshIntersection.isIntersected && meshIntersection.t < closestIntersection.t)
        {
            closestIntersection = meshIntersection;
        }
    }

    return closestIntersection;
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

    std::cout << "Loaded scene from XML successfully." << std::endl;

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
                Intersection intersection = rayObjectIntersection(scene, myRay);

                // TODO If the ray intersects with an object, calculate the color of the pixel
                // implement lighting shading etc. in color function
                // if (intersection.isIntersected)
                // {
                // }

                // this part is for testing purposes
                if (intersection.isIntersected)
                {
                    // Set the color to blue if there's an intersection
                    int pixelIndex = (i * imageWidth + j) * 3;
                    image[pixelIndex] = 0;       // Red component
                    image[pixelIndex + 1] = 0;   // Green component
                    image[pixelIndex + 2] = 255; // Blue component
                }
                else
                {
                    // Set a different color for pixels without intersection
                    int pixelIndex = (i * imageWidth + j) * 3;
                    image[pixelIndex] = (unsigned char)scene.background_color.x;
                    image[pixelIndex + 1] = (unsigned char)scene.background_color.y;
                    image[pixelIndex + 2] = (unsigned char)scene.background_color.z;
                }
            }
        }

        write_ppm(scene.cameras[cameraIndex].image_name.c_str(), image, imageWidth, imageHeight);

        delete[] image;
    }

    return 0;
}
