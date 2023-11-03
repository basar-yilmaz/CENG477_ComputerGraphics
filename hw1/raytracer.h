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
// #define EPSILON 1e-6 // TODO find an optimal epsilon value

/**
 * Calculates the dot product of two 3D vectors.
 *
 * @param a The first vector.
 * @param b The second vector.
 *
 * @return The dot product of the two vectors.
 */
float dotProduct(const Vec3f &a, const Vec3f &b)
{
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

/**
 * Divides each component of a Vec3f by a scalar value.
 *
 * @param a the Vec3f to be divided
 * @param k the scalar value to divide by
 *
 * @return a new Vec3f with each component divided by the scalar value
 *
 * @throws None
 */
Vec3f scalarDivision(const Vec3f &a, const float k)
{
    return {a.x / k, a.y / k, a.z / k};
}

/**
 * Multiplies each component of a Vec3f by a scalar value.
 *
 * @param a The Vec3f to be multiplied.
 * @param k The scalar value to multiply the Vec3f by.
 *
 * @return A new Vec3f with each component multiplied by the scalar value.
 *
 * @throws None.
 */
Vec3f scalarMulti(const Vec3f &a, const float k)
{
    return {a.x * k, a.y * k, a.z * k};
}

/**
 * Adds two Vec3f objects together.
 *
 * @param a The first Vec3f object.
 * @param b The second Vec3f object.
 *
 * @return The sum of the two Vec3f objects.
 *
 * @throws None.
 */
Vec3f add(const Vec3f &a, const Vec3f &b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

/**
 * Subtracts two Vec3f objects.
 *
 * @param a the first Vec3f object
 * @param b the second Vec3f object
 *
 * @return a new Vec3f object that is the result of subtracting b from a
 *
 * @throws None
 */
Vec3f subtract(const Vec3f &a, const Vec3f &b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

// get cross product of two vectors s.t. A X B
Vec3f crossProduct(const Vec3f &a, const Vec3f &b)
{
    return {
        (a.y * b.z) - (a.z * b.y),
        (a.z * b.x) - (a.x * b.z),
        (a.x * b.y) - (a.y * b.x)};
}

float length(const Vec3f &a)
{
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

// get the normalized vector with length = 1
// return a new vector
Vec3f normalizeVector(const Vec3f &a)
{
    float len = length(a);

    return {a.x / len, a.y / len, a.z / len};
}

// Calculate the distance between two vectors using the Euclidean distance formula
float distance(const Vec3f &a, const Vec3f &b)
{
    return std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y) + (b.z - a.z) * (b.z - a.z));
}

// calculate the determinant of a 3x3 matrix
float determinant(const float coef_matrix[3][3])
{
    return coef_matrix[0][0] * (coef_matrix[1][1] * coef_matrix[2][2] - coef_matrix[1][2] * coef_matrix[2][1]) -
           coef_matrix[0][1] * (coef_matrix[1][0] * coef_matrix[2][2] - coef_matrix[1][2] * coef_matrix[2][0]) +
           coef_matrix[0][2] * (coef_matrix[1][0] * coef_matrix[2][1] - coef_matrix[1][1] * coef_matrix[2][0]);
}

// Calculates the determinant of a 3x3 matrix.
// Arguments:
//   v0, v1, v2: Vec3f objects representing the matrix rows.
// Returns:
//   A float representing the determinant of the matrix.
float determinantA(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2)
{
    float det = v0.x * (v1.y * v2.z - v2.y * v1.z) +
                v0.y * (v2.x * v1.z - v1.x * v2.z) +
                v0.z * (v1.x * v2.y - v1.y * v2.x);
    return det;
}

Vec3f multiplyVectors(const Vec3f &a, const Vec3f &b)
{
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}
