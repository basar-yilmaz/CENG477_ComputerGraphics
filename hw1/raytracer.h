#pragma once

/*
 * This file contains the functions for vector operations
 * such as dot product, cross product, etc.
 * There might be various ways to implement these functions more efficiently!
 * Even some operants can be overloaded to make the code more readable.
 */

#include <math.h>
#include "parser.h"

using parser::Vec3f;

float dotProduct(const Vec3f &a, const Vec3f &b)
{
    float result = 0;
    result += a.x * b.x;
    result += a.y * b.y;
    result += a.z * b.z;

    printf("dot product: %f\n", result);
    return result;
}

Vec3f scalarMulti(const Vec3f &a, const float k)
{
    Vec3f result;
    result.x = a.x * k;
    result.y = a.y * k;
    result.z = a.z * k;
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
    // some floating point precautions
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