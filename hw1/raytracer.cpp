#include <iostream>
#include <math.h>
#include "parser.h"
#include "ppm.h"
#include "raytracer.h"
#include <limits>

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

    if (fabs(A) < EPSILON) {
        return point;
    }

    if (t_val < 0.0) {
        return point;
    }
    if (gamma_val < 0 || gamma_val > 1) {
        return point;
    }

    if (beta_val < 0 || beta_val >(1 - gamma_val)) {
        return point;
    }
    // check if the ray intersects with the triangle

    point.isIntersected = true;
    point.intersectionPoint = getIntersectionPoint(ray, t_val);
    point.t = t_val;
    point.normal = normalizeVector(crossProduct(subtract(b, a), subtract(c, a)));
    point.objId = triangle.indices.v0_id;
    point.objType = TRIANGLE;
    point.material_id = triangle.material_id;

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
        if (valid_sol < 0) {
            point.isIntersected = false;
            return point;
        }
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

    closestIntersection.t = numeric_limits<float>::max();

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
    closestIntersection.t = numeric_limits<float>::max(); // TODO: set this to the max value of float or some very high value temporarily
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

/*
 * Calculate the ambient shading
 * @param scene: the scene object
 * @param inter: the intersection point
 * @param material: the material object
 * @return: the ambient shading radiance
 */
Vec3f ambientShading(const Scene &scene, const Intersection &inter, const Material &material)
{
    // all the abrevation from slides
    // outgoing radiance = k_a(coefficient) * I_a(ambient light intensity)

    Vec3f ambientCoef = material.ambient;
    Vec3f ambientLight = scene.ambient_light;

    Vec3f radiance = multiplyVectors(ambientCoef, ambientLight);

    return radiance;
}

Vec3f diffuseShading(const Scene &scene, const PointLight &light, const Intersection &inter, const Material &material)
{
    // all the abrevation from slides
    // outgoing radiance = k_d(coefficient) * max(0, w_i.n)(cos approx) * I(light intensity) / r^2(light distance) (inverse square law)
    Vec3f lightDirection, normal, diffuseCoef, radiance;

    lightDirection = subtract(light.position, inter.intersectionPoint); // w_i
    normal = inter.normal;
    diffuseCoef = material.diffuse;

    double r = length(lightDirection);

    // Normalize the light direction vector
    lightDirection = normalizeVector(lightDirection);

    radiance = scalarMulti(diffuseCoef, max(dotProduct(lightDirection, normal), 0.0f));
    Vec3f intensity = light.intensity;
    radiance = multiplyVectors(radiance, intensity);
    radiance = scalarDivision(radiance, r * r);

    return radiance;
}

/*
 * Calculate the specular shading
 * @param scene: the scene object
 * @param light: the light object
 * @param inter: the intersection point
 * @param material: the material object
 * @param camera: the camera object (differs from ambient and diffuse shading)
 * as the specular shading depends on the view direction
 * @return: the specular shading radiance
 */
Vec3f specularShading(const Scene& scene, const PointLight& light, const Intersection& inter, const Material& material, const Camera& camera)
{
    Vec3f specularCoef, lightDirection, cameraDirection, normal, halfVector, receivedRadiance, radiance;

    specularCoef = material.specular;
    lightDirection = subtract(light.position, inter.intersectionPoint);
    cameraDirection = normalizeVector(subtract(camera.position, inter.intersectionPoint)); // Fixed this
    normal = inter.normal;

    double r = length(lightDirection);

    // Calculate received radiance considering the inverse square law
    receivedRadiance = scalarDivision(light.intensity, r * r); // Adjusted for inverse square law

    lightDirection = normalizeVector(lightDirection);

    halfVector = normalizeVector(add(lightDirection, cameraDirection));

    radiance = scalarMulti(specularCoef, pow(max(dotProduct(normal, halfVector), 0.0f), material.phong_exponent));

    radiance = multiplyVectors(radiance, receivedRadiance);

    return radiance;
}

/*
 * Calculate the color of the pixel
 * @param scene: the scene object
 * @param intersection: the intersection point
 * @param recDepth: the recursion depth
 * @return: the color of the pixel
 */
Vec3f coloring(const Scene& scene, const Intersection& intersection, int recDepth, const Camera& camera)
{
    Vec3f color = { 0, 0, 0 };

    if (intersection.isIntersected)
    {

        int materialId = intersection.material_id;

        color = ambientShading(scene, intersection, scene.materials[materialId - 1]); // color so far we get from only ambient light

        for (PointLight light : scene.point_lights)
        {
            int isShadow = 0;

            // calculate the shadowray and check if it intersects with any object
            Vec3f shadowRayDirection = normalizeVector(subtract(light.position, intersection.intersectionPoint));
            Vec3f offset = scalarMulti(intersection.normal, scene.shadow_ray_epsilon);
            Ray shadowRay = { add(intersection.intersectionPoint, offset), shadowRayDirection };
            Intersection shadowRayIntersection = rayObjectIntersection(scene, shadowRay);

            // if the shadow ray intersects with an object, we need to check if the intersection point is behind the light source
            if (shadowRayIntersection.isIntersected)
            {
                Vec3f shadowRayIntersectionPoint = shadowRayIntersection.intersectionPoint;
                Vec3f lightToIntersection = subtract(shadowRayIntersectionPoint, light.position);
                float distanceToLight = length(subtract(light.position, intersection.intersectionPoint));
                //float distanceToIntersection = length(subtract(shadowRayIntersectionPoint, intersection.intersectionPoint));
                if (shadowRayIntersection.t < distanceToLight)
                {
                    isShadow = 1;
                }
            }

            if (!isShadow)
            {
                /*Vec3f lightDirection = normalizeVector(subtract(light.position, intersection.intersectionPoint));
                Vec3f normal = intersection.normal;*/

                // Compute and add diffuse shading
                Vec3f diffuseColor = diffuseShading(scene, light, intersection, scene.materials[materialId - 1]);
                color = add(color, diffuseColor);

                // Compute and add specular shading
                Vec3f specularColor = specularShading(scene, light, intersection, scene.materials[materialId - 1], camera);
                color = add(color, specularColor);
            }
        }

        Vec3f reflection;
        reflection.x = 0;
        reflection.y = 0;
        reflection.z = 0;
        bool isMirror = (scene.materials[materialId - 1].mirror.x > 0 || scene.materials[materialId - 1].mirror.y > 0 || scene.materials[materialId - 1].mirror.z > 0);

        // at this point we processed the light yet there might be reflections so we need to check the recursion depth
        if (recDepth > 0  && isMirror)
        {
            // Compute the incoming direction from the intersection point back to the camera
            Vec3f incomingDirection = normalizeVector(subtract(intersection.intersectionPoint, camera.position));
            float wi = -2 * dotProduct(incomingDirection, intersection.normal);
            Vec3f normal_ref;
            normal_ref.x = intersection.normal.x * wi + incomingDirection.x;
            normal_ref.y = intersection.normal.y * wi + incomingDirection.y;
            normal_ref.z = intersection.normal.z * wi + incomingDirection.z;

            normal_ref = normalizeVector(normal_ref);

            Vec3f ref_epsilon;
            ref_epsilon.x = normal_ref.x * scene.shadow_ray_epsilon;
            ref_epsilon.y = normal_ref.y * scene.shadow_ray_epsilon;
            ref_epsilon.z = normal_ref.z * scene.shadow_ray_epsilon;

            Ray reflectionRay = { add(intersection.intersectionPoint, ref_epsilon), normal_ref };
            Intersection intRes = rayObjectIntersection(scene, reflectionRay);

            if (!(intRes.objId == intersection.objId && intRes.objType == intersection.objType))
            {
                reflection = coloring(scene, intRes, (recDepth -1), camera);
            }

            color.x += reflection.x * scene.materials[materialId - 1].mirror.x;
            color.y += reflection.y * scene.materials[materialId - 1].mirror.y;
            color.z += reflection.z * scene.materials[materialId - 1].mirror.z;
        }

    }

    else
    {
        color.x = min(scene.background_color.x, 255);
        color.y = min(scene.background_color.y, 255);
        color.z = min(scene.background_color.z, 255);
    }
    color.x = min(color.x, 255.0f);
    color.y = min(color.y, 255.0f);
    color.z = min(color.z, 255.0f);
    return color;
}

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_scene_file.xml>" << std::endl;
        return 1;
    }

    parser::Scene scene;
    //int numberOfCameras = scene.cameras.size();

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
        int pixelIndex = 0;
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
                Vec3f color = coloring(scene, intersection, scene.max_recursion_depth, scene.cameras[cameraIndex]);
                if (color.x > 255)
                    image[pixelIndex] = 255;
                else
                    image[pixelIndex] = round(color.x);
                if (color.y > 255)
                    image[pixelIndex + 1] = 255;
				else
					image[pixelIndex + 1] = round(color.y);
                if (color.z > 255)
					image[pixelIndex + 2] = 255;
                else
                    image[pixelIndex + 2] = round(color.z);
                pixelIndex+= 3;
                // this part is for testing purposes
                /*if (intersection.isIntersected)
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
                }*/
            }
        }

        write_ppm(scene.cameras[cameraIndex].image_name.c_str(), image, imageWidth, imageHeight);

        delete[] image;
    }

    return 0;
}
