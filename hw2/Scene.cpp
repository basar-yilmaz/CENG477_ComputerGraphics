#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;
/*
	Get modelling transformations
*/
Matrix4 Scene::modellingTransformations(Mesh *mesh)
{
	Matrix4 modellingMatrix = getIdentityMatrix();
	for (int i = 0; i < (mesh->numberOfTransformations); i++)
	{
		// rotation check this one again
		if (mesh->transformationTypes[i] == 'r')
		{
			Rotation *rotate1 = rotations[mesh->transformationIds[i] - 1];

			double angle = (rotate1->angle * M_PI) / 180;

			// final rotation matrix
			double rotationMatrix[4][4] =
				{
					{cos(angle) + (rotate1->ux * rotate1->ux) * (1. - cos(angle)), rotate1->ux * rotate1->uy * (1. - cos(angle)) - rotate1->uz * sin(angle), rotate1->ux * rotate1->uz * (1. - cos(angle)) + rotate1->uy * sin(angle), 0},
					{rotate1->uy * rotate1->ux * (1. - cos(angle)) + rotate1->uz * sin(angle), cos(angle) + (rotate1->uy * rotate1->uy) * (1. - cos(angle)), rotate1->uy * rotate1->uz * (1. - cos(angle)) - rotate1->ux * sin(angle), 0},
					{rotate1->uz * rotate1->ux * (1. - cos(angle)) - rotate1->uy * sin(angle), rotate1->uz * rotate1->uy * (1. - cos(angle)) + rotate1->ux * sin(angle), cos(angle) + (rotate1->uz * rotate1->uz) * (1. - cos(angle)), 0},
					{0, 0, 0, 1}};
			modellingMatrix = multiplyMatrixWithMatrix(modellingMatrix, Matrix4(rotationMatrix));
		}

		// translation
		else if (mesh->transformationTypes[i] == 't')
		{
			Translation *translation1 = translations[mesh->transformationIds[i] - 1];
			double translationMatrix[4][4] = {{1, 0, 0, translation1->tx},
											  {0, 1, 0, translation1->ty},
											  {0, 0, 1, translation1->tz},
											  {0, 0, 0, 1}};
			modellingMatrix = multiplyMatrixWithMatrix(modellingMatrix, Matrix4(translationMatrix));
		}

		// scaling
		else if (mesh->transformationTypes[i] == 's')
		{
			Scaling *scale1 = scalings[mesh->transformationIds[i] - 1];
			double scalingMatrix[4][4] = {{scale1->sx, 0, 0, 0},
										  {0, scale1->sy, 0, 0},
										  {0, 0, scale1->sz, 0},
										  {0, 0, 0, 1}};
			modellingMatrix = multiplyMatrixWithMatrix(modellingMatrix, Matrix4(scalingMatrix));
		}
	}
	return modellingMatrix;
}

/*
	Get orthographic projection transformation
*/
Matrix4 Scene::orthTransformation(Camera *camera)
{
	double temp[4][4] = {{2 / (camera->right - camera->left), 0, 0, -(camera->right + camera->left) / (camera->right - camera->left)},
						 {0, 2 / (camera->top - camera->bottom), 0, -(camera->top + camera->bottom) / (camera->top - camera->bottom)},
						 {0, 0, -2 / (camera->far - camera->near), -(camera->far + camera->near) / (camera->far - camera->near)},
						 {0, 0, 0, 1}};
	return Matrix4(temp);
}

/*
	Get perspective projection transformation
*/
Matrix4 Scene::perspectiveTransformation(Camera *camera)
{
	double temp[4][4] = {{2 * camera->near / (camera->right - camera->left), 0, (camera->right + camera->left) / (camera->right - camera->left), 0},
						 {0, 2 * camera->near / (camera->top - camera->bottom), (camera->top + camera->bottom) / (camera->top - camera->bottom), 0},
						 {0, 0, -(camera->far + camera->near) / (camera->far - camera->near), -2 * camera->far * camera->near / (camera->far - camera->near)},
						 {0, 0, -1, 0}};
	return Matrix4(temp);
}

/*
	perspective division
*/
Vec4 Scene::perspectiveDivide(Vec4 &v)
{
	Vec4 res;
	res.x = v.x / v.t;
	res.y = v.y / v.t;
	res.z = v.z / v.t;
	res.t = 1;
	res.colorId = v.colorId;
	res.color = v.color;
	return res;
}

/*
	Get camera transformation
*/
Matrix4 Scene::camTransformation(Camera *camera)
{
	// formula from viewving transformation slides page 6
	double tMatrix[4][4] = {{1, 0, 0, -(camera->position.x)},
							{0, 1, 0, -(camera->position.y)},
							{0, 0, 1, -(camera->position.z)},
							{0, 0, 0, 1}};

	// formula from viewving transformation slides page 7
	double rMatrix[4][4] = {{camera->u.x, camera->u.y, camera->u.z, 0},
							{camera->v.x, camera->v.y, camera->v.z, 0},
							{camera->w.x, camera->w.y, camera->w.z, 0},
							{0, 0, 0, 1}};

	// formula from viewing transformation slides page 8
	return multiplyMatrixWithMatrix(Matrix4(rMatrix), Matrix4(tMatrix));
}

/*
	Get viewport transformation
*/
Matrix4 Scene::viewTransformation(Camera *camera)
{
	double temp[4][4] = {{camera->horRes / (double)2, 0, 0, (camera->horRes - 1) / (double)2},
						 {0, camera->verRes / (double)2, 0, (camera->verRes - 1) / (double)2},
						 {0, 0, (double)0.5, (double)0.5},
						 {0, 0, 0, 0}}; // not sure maybe {0, 0, 0, 1}
	return Matrix4(temp);
}
/*
	Parses XML file
*/

Color colorDifference(Color &c1, Color &c2)
{
	Color result;
	result.r = c1.r - c2.r;
	result.g = c1.g - c2.g;
	result.b = c1.b - c2.b;
	return result;
}

Color colorDivision(Color &c1, int divisor)
{
	Color result;
	result.r = c1.r / divisor;
	result.g = c1.g / divisor;
	result.b = c1.b / divisor;
	return result;
}

Color colorAddition(Color &c1, Color &c2)
{
	Color result;
	result.r = c1.r + c2.r;
	result.g = c1.g + c2.g;
	result.b = c1.b + c2.b;
	return result;
}

/*
	Checks if line is visible helper function for Liang-Barsky algorithm
	Algorithm from slides page 46
*/
bool visible(double den, double num, double &enteringT, double &leavingT)
{
	if (den > 0) // potentially entering
	{
		double t = num / den;
		if (t > leavingT)
		{
			return false;
		}
		if (t > enteringT)
		{
			enteringT = t;
		}
	}
	else if (den < 0) // potentially leaving
	{
		double t = num / den;
		if (t < enteringT)
		{
			return false;
		}
		if (t < leavingT)
		{
			leavingT = t;
		}
	}
	else if (num > 0) // line parallel to edge
	{
		return false;
	}

	return true;
}
/*
	Liang-Barsky algorithm
	Algorithm from slides page 47
*/
bool Scene::lineClipping(LineVec4 &vertices)
{
	// clipping slides page 46-47
	double enteringT = 0;
	double leavingT = 1;
	bool visibility = 0;

	// calculate dx dy dz
	double dx = vertices.vertex2.x - vertices.vertex1.x;
	double dy = vertices.vertex2.y - vertices.vertex1.y;
	double dz = vertices.vertex2.z - vertices.vertex1.z;

	// we will be around -1 to 1
	int x_min = -1, x_max = 1;
	int y_min = -1, y_max = 1;
	int z_min = -1, z_max = 1;

	// calculate color difference
	Color colorDifference;
	colorDifference.r = (vertices.vertex2.color.r - vertices.vertex1.color.r);
	colorDifference.g = (vertices.vertex2.color.g - vertices.vertex1.color.g);
	colorDifference.b = (vertices.vertex2.color.b - vertices.vertex1.color.b);

	if (visible(dx, (x_min - vertices.vertex1.x), enteringT, leavingT)) // left
	{
		if (visible(-dx, (vertices.vertex1.x - x_max), enteringT, leavingT)) // right
		{
			if (visible(dy, (y_min - vertices.vertex1.y), enteringT, leavingT)) // bottom
			{
				if (visible(-dy, (vertices.vertex1.y - y_max), enteringT, leavingT)) // top
				{
					if (visible(dz, (z_min - vertices.vertex1.z), enteringT, leavingT)) // front
					{
						if (visible(-dz, (vertices.vertex1.z - z_max), enteringT, leavingT)) // back
						{
							if (leavingT < 1)
							{
								vertices.vertex2.x = vertices.vertex1.x + leavingT * dx;
								vertices.vertex2.y = vertices.vertex1.y + leavingT * dy;
								vertices.vertex2.z = vertices.vertex1.z + leavingT * dz;
								vertices.vertex2.color.r = vertices.vertex1.color.r + leavingT * colorDifference.r;
								vertices.vertex2.color.g = vertices.vertex1.color.g + leavingT * colorDifference.g;
								vertices.vertex2.color.b = vertices.vertex1.color.b + leavingT * colorDifference.b;
							}
							if (enteringT > 0)
							{
								vertices.vertex1.x = vertices.vertex1.x + enteringT * dx;
								vertices.vertex1.y = vertices.vertex1.y + enteringT * dy;
								vertices.vertex1.z = vertices.vertex1.z + enteringT * dz;
								vertices.vertex1.color.r = vertices.vertex1.color.r + enteringT * colorDifference.r;
								vertices.vertex1.color.g = vertices.vertex1.color.g + enteringT * colorDifference.g;
								vertices.vertex1.color.b = vertices.vertex1.color.b + enteringT * colorDifference.b;
							}
							visibility = 1;
						}
					}
				}
			}
		}
	}

	return visibility;
}
/*
	Transformations, clipping, culling, rasterization are done here.
*/

bool Scene::isCulled(Vec4 vertices[3])
{
	Vec3 vector1 = subtractVec3(Vec3(vertices[1].x, vertices[1].y, vertices[1].z, vertices[1].colorId), Vec3(vertices[0].x, vertices[0].y, vertices[0].z, vertices[0].colorId));
	Vec3 vector2 = subtractVec3(Vec3(vertices[2].x, vertices[2].y, vertices[2].z, vertices[2].colorId), Vec3(vertices[0].x, vertices[0].y, vertices[0].z, vertices[0].colorId));

	Vec3 normal = crossProductVec3(vector1, vector2);
	normal = normalizeVec3(normal);
	double dotProductOfViewVsOrigin = dotProductVec3(normal, Vec3(vertices[1].x, vertices[1].y, vertices[1].z, vertices[1].colorId));
	if (dotProductOfViewVsOrigin < 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Scene::lineRasterization(std::pair<Vec3, Vec3> &vertices)
{
	double dx = vertices.second.x - vertices.first.x;
	double dy = vertices.second.y - vertices.first.y;
	int movementInImage = 1; //
	Color colorChange;

	// we need to implement modified rasterization algorithm
	// since the slope of the line may be greater than 1
	// implementation similar to rasterization slides page 22

	if (std::abs(dy) > std::abs(dx))
	{
		// the case where slope is greater than 1
		// difference is that x and y values should be swapped

		// if x2 < x1 swap them
		if (vertices.second.y < vertices.first.y)
		{
			Vec3 tmp = vertices.first;
			vertices.first = vertices.second;
			vertices.second = tmp;
		}

		movementInImage = -1 ? vertices.second.x < vertices.first.x : 1;

		int x0 = vertices.first.x;
		Color c0 = vertices.first.color;
		int d = (vertices.first.x - vertices.second.x) + (movementInImage * 0.5 * (vertices.second.y - vertices.first.y));
		Color tempDiff = colorDifference(vertices.second.color, vertices.first.color);
		colorChange = colorDivision(tempDiff, vertices.second.y - vertices.first.y); // skip alpha value by directly computing color increment

		// we got our constants now iterate through whole y values
		for (int i = vertices.first.y; i <= vertices.second.y; i++)
		{
			Color rounded_c0;
			rounded_c0.r = (int)(c0.r + 0.5);
			rounded_c0.g = (int)(c0.g + 0.5);
			rounded_c0.b = (int)(c0.b + 0.5);

			image[x0][i] = rounded_c0;

			// choose between x0 and x0+1 (NE or E)

			if (d * movementInImage <= 0)
			{
				d += (vertices.first.y - vertices.second.y); // move horizontally only
			}
			else
			{ // move diagonally (NE)
				d += (vertices.first.y - vertices.second.y) + (movementInImage * (vertices.second.x - vertices.first.x));
				x0 += movementInImage;
			}
			c0 = colorAddition(c0, colorChange); // interpolate color
		}
	}
	else
	{
		// the case where slope is in between 0 and <=1

		// if x2 < x1 swap them
		if (vertices.second.x < vertices.first.x)
		{
			Vec3 tmp = vertices.first;
			vertices.first = vertices.second;
			vertices.second = tmp;
		}

		movementInImage = -1 ? vertices.second.y < vertices.first.y : 1;

		int y0 = vertices.first.y;
		Color c0 = vertices.first.color;
		int d = (vertices.first.y - vertices.second.y) + (movementInImage * 0.5 * (vertices.second.x - vertices.first.x));
		Color tempDiff = colorDifference(vertices.second.color, vertices.first.color);
		colorChange = colorDivision(tempDiff, vertices.second.x - vertices.first.x);

		// we got our constants now iterate through whole x values
		for (int i = vertices.first.x; i <= vertices.second.x; i++)
		{
			Color rounded_c0;
			rounded_c0.r = (int)(c0.r + 0.5);
			rounded_c0.g = (int)(c0.g + 0.5);
			rounded_c0.b = (int)(c0.b + 0.5);

			image[i][y0] = rounded_c0;

			// choose between y0 and y0+1 (NE or E)

			if (d * movementInImage <= 0)
			{
				d += (vertices.first.y - vertices.second.y); // move horizontally only
			}
			else
			{ // move diagonally (NE)
				d += (vertices.first.y - vertices.second.y) + (movementInImage * (vertices.second.x - vertices.first.x));
				y0 += movementInImage;
			}
			c0 = colorAddition(c0, colorChange); // interpolate color
		}
	}
}

double Scene::f_alpha(vector<Vec3> &vertices, double currentX, double currentY)
{
	double temp1 = currentX * (vertices[1].y - vertices[2].y) + currentY * (vertices[2].x - vertices[1].x) + (vertices[1].x * vertices[2].y) - (vertices[1].y * vertices[2].x);
	double temp2 = vertices[0].x * (vertices[1].y - vertices[2].y) + vertices[0].y * (vertices[2].x - vertices[1].x) + (vertices[1].x * vertices[2].y) - (vertices[1].y * vertices[2].x);
	return temp1 / temp2;
}

double Scene::f_beta(vector<Vec3> &vertices, double currentX, double currentY)
{
	double temp1 = currentX * (vertices[2].y - vertices[0].y) + currentY * (vertices[0].x - vertices[2].x) + (vertices[2].x * vertices[0].y) - (vertices[2].y * vertices[0].x);
	double temp2 = vertices[1].x * (vertices[2].y - vertices[0].y) + vertices[1].y * (vertices[0].x - vertices[2].x) + (vertices[2].x * vertices[0].y) - (vertices[2].y * vertices[0].x);
	return temp1 / temp2;
}

double Scene::f_gamma(vector<Vec3> &vertices, double currentX, double currentY)
{
	double temp1 = currentX * (vertices[0].y - vertices[1].y) + currentY * (vertices[1].x - vertices[0].x) + (vertices[0].x * vertices[1].y) - (vertices[0].y * vertices[1].x);
	double temp2 = vertices[2].x * (vertices[0].y - vertices[1].y) + vertices[2].y * (vertices[1].x - vertices[0].x) + (vertices[0].x * vertices[1].y) - (vertices[0].y * vertices[1].x);
	return temp1 / temp2;
}

void Scene::triangleRasterization(vector<Vec3> &vertices, int maxFrameX, int maxFrameY)
{
	// implementation from rasterization slides page 30

	double x_min = std::min(vertices[0].x, std::min(vertices[1].x, vertices[2].x));
	if (x_min < 0)
	{
		x_min = 0;
	}
	// if (x_min > maxFrameX - 1)
	// {
	// 	x_min = maxFrameX - 1;
	// }
	double x_max = std::max(vertices[0].x, std::max(vertices[1].x, vertices[2].x));
	if (x_max < 0)
	{
		x_max = 0;
	}
	if (x_max > maxFrameX - 1)
	{
		x_max = maxFrameX - 1;
	}
	double y_min = std::min(vertices[0].y, std::min(vertices[1].y, vertices[2].y));
	if (y_min < 0)
	{
		y_min = 0;
	}
	// if (y_min > maxFrameY - 1)
	// {
	// 	y_min = maxFrameY - 1;
	// }
	double y_max = std::max(vertices[0].y, std::max(vertices[1].y, vertices[2].y));
	if (y_max < 0)
	{
		y_max = 0;
	}
	if (y_max > maxFrameY - 1)
	{
		y_max = maxFrameY - 1;
	}

	for (int y = y_min; y <= y_max; y++)
	{
		for (int x = x_min; x <= x_max; x++)
		{
			// TODO check this
			double alpha, beta, gamma;
			alpha = f_alpha(vertices, x, y);
			beta = f_beta(vertices, x, y);
			gamma = f_gamma(vertices, x, y);

			if (alpha >= 0 && beta >= 0 && gamma >= 0)
			{
				Color color;
				color.r = alpha * vertices[0].color.r + beta * vertices[1].color.r + gamma * vertices[2].color.r;
				color.g = alpha * vertices[0].color.g + beta * vertices[1].color.g + gamma * vertices[2].color.g;
				color.b = alpha * vertices[0].color.b + beta * vertices[1].color.b + gamma * vertices[2].color.b;

				Color rounded_color;
				rounded_color.r = (int)(color.r + 0.5);
				rounded_color.g = (int)(color.g + 0.5);
				rounded_color.b = (int)(color.b + 0.5);

				image[x][y] = rounded_color;
			}
		}
	}
}

void Scene::forwardRenderingPipeline(Camera *camera)
{
	// TODO: Implement z-depth buffering
	// TODO: some pictures does not show up correctly

	bool cullingFlag = this->cullingEnabled;

	Matrix4 viewingTransformationMatrix = viewTransformation(camera);

	for (const auto &currentMesh : meshes)
	{
		for (const auto &currentTriangle : currentMesh->triangles)
		{
			// Vertex 0
			Vec4 vertex0 = Vec4(vertices[currentTriangle.vertexIds[0] - 1]->x,
								vertices[currentTriangle.vertexIds[0] - 1]->y,
								vertices[currentTriangle.vertexIds[0] - 1]->z,
								1,
								*(colorsOfVertices[vertices[currentTriangle.vertexIds[0] - 1]->colorId - 1]));

			// Vertex 1
			Vec4 vertex1 = Vec4(vertices[currentTriangle.vertexIds[1] - 1]->x,
								vertices[currentTriangle.vertexIds[1] - 1]->y,
								vertices[currentTriangle.vertexIds[1] - 1]->z,
								1,
								*(colorsOfVertices[vertices[currentTriangle.vertexIds[1] - 1]->colorId - 1]));

			// Vertex 2
			Vec4 vertex2 = Vec4(vertices[currentTriangle.vertexIds[2] - 1]->x,
								vertices[currentTriangle.vertexIds[2] - 1]->y,
								vertices[currentTriangle.vertexIds[2] - 1]->z,
								1,
								*(colorsOfVertices[vertices[currentTriangle.vertexIds[2] - 1]->colorId - 1]));

			// apply modelling transformations
			Matrix4 modelMatrix = modellingTransformations(currentMesh);

			// apply camera transformations
			Matrix4 cameraMatrix = camTransformation(camera);

			// apply projection transformations depends on type
			Matrix4 projectionMatrix = (camera->projectionType == 0) ? orthTransformation(camera) : perspectiveTransformation(camera);

			Matrix4 wholeTransformationsMatrix = multiplyMatrixWithMatrix(projectionMatrix, multiplyMatrixWithMatrix(cameraMatrix, modelMatrix));

			Vec4 transformedVertex0 = multiplyMatrixWithVec4(wholeTransformationsMatrix, vertex0);
			Vec4 transformedVertex1 = multiplyMatrixWithVec4(wholeTransformationsMatrix, vertex1);
			Vec4 transformedVertex2 = multiplyMatrixWithVec4(wholeTransformationsMatrix, vertex2);

			transformedVertex0.color = vertex0.color;
			transformedVertex1.color = vertex1.color;
			transformedVertex2.color = vertex2.color;

			transformedVertex0 = camera->projectionType ? perspectiveDivide(transformedVertex0) : transformedVertex0;
			transformedVertex1 = camera->projectionType ? perspectiveDivide(transformedVertex1) : transformedVertex1;
			transformedVertex2 = camera->projectionType ? perspectiveDivide(transformedVertex2) : transformedVertex2;

			Vec4 verticesArray[3];
			verticesArray[0] = transformedVertex0;
			verticesArray[1] = transformedVertex1;
			verticesArray[2] = transformedVertex2;

			if (cullingFlag && isCulled(verticesArray))
				continue;

			if (currentMesh->type)
			{
				// solid mode

				vector<Vec3> viewPortedVertices;

				for (int i = 0; i < 3; i++)
				{
					Vec4 temporaryTransferMatrix = multiplyMatrixWithVec4(viewingTransformationMatrix, verticesArray[i]);
					Vec3 currentTransformedVertex = Vec3(temporaryTransferMatrix.x, temporaryTransferMatrix.y, temporaryTransferMatrix.z, verticesArray[i].color);
					viewPortedVertices.push_back(currentTransformedVertex);
				}

				triangleRasterization(viewPortedVertices, camera->horRes, camera->verRes);
			}
			else
			{
				// wireframe mode
				LineVec4 v0 = LineVec4(verticesArray[0], verticesArray[1]);
				LineVec4 v1 = LineVec4(verticesArray[1], verticesArray[2]);
				LineVec4 v2 = LineVec4(verticesArray[2], verticesArray[0]);

				// check clipping here
				bool firstLineVisible = lineClipping(v0);
				bool secondLineVisible = lineClipping(v1);
				bool thirdLineVisible = lineClipping(v2);

				// check if it is totally clipped out or not draw it
				if (firstLineVisible)
				{
					Vec4 tempMatrix1 = multiplyMatrixWithVec4(viewingTransformationMatrix, v0.vertex1);
					Vec3 viewportedVertex0 = Vec3(tempMatrix1.x, tempMatrix1.y, tempMatrix1.z, v0.vertex1.color);

					Vec4 tempMatrix2 = multiplyMatrixWithVec4(viewingTransformationMatrix, v0.vertex1);
					Vec3 viewportedVertex1 = Vec3(tempMatrix2.x, tempMatrix2.y, tempMatrix2.z, v1.vertex2.color);

					std::pair<Vec3, Vec3> viewportedVertices = std::make_pair(viewportedVertex0, viewportedVertex1);

					lineRasterization(viewportedVertices);
				}

				if (secondLineVisible)
				{

					Vec4 tempMatrix1 = multiplyMatrixWithVec4(viewingTransformationMatrix, v1.vertex1);
					Vec3 viewportedVertex0 = Vec3(tempMatrix1.x, tempMatrix1.y, tempMatrix1.z, v1.vertex1.color);

					Vec4 tempMatrix2 = multiplyMatrixWithVec4(viewingTransformationMatrix, v1.vertex1);
					Vec3 viewportedVertex1 = Vec3(tempMatrix2.x, tempMatrix2.y, tempMatrix2.z, v1.vertex2.color);

					std::pair<Vec3, Vec3> viewportedVertices = std::make_pair(viewportedVertex0, viewportedVertex1);

					lineRasterization(viewportedVertices);
				}

				if (thirdLineVisible)
				{

					Vec4 tempMatrix1 = multiplyMatrixWithVec4(viewingTransformationMatrix, v2.vertex1);
					Vec3 viewportedVertex0 = Vec3(tempMatrix1.x, tempMatrix1.y, tempMatrix1.z, v2.vertex1.color);

					Vec4 tempMatrix2 = multiplyMatrixWithVec4(viewingTransformationMatrix, v2.vertex1);
					Vec3 viewportedVertex1 = Vec3(tempMatrix2.x, tempMatrix2.y, tempMatrix2.z, v2.vertex2.color);

					std::pair<Vec3, Vec3> viewportedVertices = std::make_pair(viewportedVertex0, viewportedVertex1);

					lineRasterization(viewportedVertices);
				}
			}
		}
	}
}

Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	command = "convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}
