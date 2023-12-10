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
	Parses XML file
*/
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

// we will use this struct to store lines for Liang-Barsky algorithm
struct LineVec4
{
	Vec4 vec;
	Color color;

	LineVec4(Vec4 &v, Color &c) : vec(v), color(c) {}
};

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

void perspectiveDivide(Vec4 &v)
{
	v.x = v.x / v.t;
	v.y = v.y / v.t;
	v.z = v.z / v.t;
	v.t = 1;
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
bool lineClipping(LineVec4 &v1, LineVec4 &v2)
{
	// clipping slides page 46-47
	double enteringT = 0;
	double leavingT = 1;
	bool visibility = 0;

	// calculate dx dy dz
	double dx = v2.vec.x - v1.vec.x;
	double dy = v2.vec.y - v1.vec.y;
	double dz = v2.vec.z - v1.vec.z;

	// we will be around -1 to 1
	int x_min = -1, x_max = 1;
	int y_min = -1, y_max = 1;
	int z_min = -1, z_max = 1;

	// calculate color difference
	Color colorDifference;
	colorDifference.r = (v2.color.r - v1.color.r);
	colorDifference.g = (v2.color.g - v1.color.g);
	colorDifference.b = (v2.color.b - v1.color.b);

	if (visible(dx, (x_min - v1.vec.x), enteringT, leavingT)) // left
	{
		if (visible(-dx, (v1.vec.x - x_max), enteringT, leavingT)) // right
		{
			if (visible(dy, (y_min - v1.vec.y), enteringT, leavingT)) // bottom
			{
				if (visible(-dy, (v1.vec.y - y_max), enteringT, leavingT)) // top
				{
					if (visible(dz, (z_min - v1.vec.z), enteringT, leavingT)) // front
					{
						if (visible(-dz, (v1.vec.z - z_max), enteringT, leavingT)) // back
						{
							if (leavingT < 1)
							{
								v2.vec.x = v1.vec.x + leavingT * dx;
								v2.vec.y = v1.vec.y + leavingT * dy;
								v2.vec.z = v1.vec.z + leavingT * dz;
								v2.color.r = v1.color.r + leavingT * colorDifference.r;
								v2.color.g = v1.color.g + leavingT * colorDifference.g;
								v2.color.b = v1.color.b + leavingT * colorDifference.b;
							}
							if (enteringT > 0)
							{
								v1.vec.x = v1.vec.x + enteringT * dx;
								v1.vec.y = v1.vec.y + enteringT * dy;
								v1.vec.z = v1.vec.z + enteringT * dz;
								v1.color.r = v1.color.r + enteringT * colorDifference.r;
								v1.color.g = v1.color.g + enteringT * colorDifference.g;
								v1.color.b = v1.color.b + enteringT * colorDifference.b;
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

void lineRasterization(vector<vector<Color>> &image, LineVec4 &v1, LineVec4 &v2)
{
	double dx = v2.vec.x - v1.vec.x;
	double dy = v2.vec.y - v1.vec.y;
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
		if (v2.vec.y < v1.vec.y)
		{
			LineVec4 tmp = v1;
			v1 = v2;
			v2 = tmp;
		}

		movementInImage = -1 ? v2.vec.x < v1.vec.x : 1;

		int x0 = v1.vec.x;
		Color c0 = v1.color;
		int d = (v1.vec.x - v2.vec.x) + (movementInImage * 0.5 * (v2.vec.y - v1.vec.y));
		Color tempDiff = colorDifference(v2.color, v1.color);
		colorChange = colorDivision(tempDiff, v2.vec.y - v1.vec.y); // skip alpha value by directly computing color increment

		// we got our constants now iterate through whole y values
		for (int i = v1.vec.y; i <= v2.vec.y; i++)
		{
			Color rounded_c0;
			rounded_c0.r = (int)(c0.r + 0.5);
			rounded_c0.g = (int)(c0.g + 0.5);
			rounded_c0.b = (int)(c0.b + 0.5);

			image[x0][i] = rounded_c0;

			// choose between x0 and x0+1 (NE or E)

			if (d * movementInImage <= 0)
			{
				d += (v1.vec.y - v2.vec.y); // move horizontally only
			}
			else
			{ // move diagonally (NE)
				d += (v1.vec.y - v2.vec.y) + (movementInImage * (v2.vec.x - v1.vec.x));
				x0 += movementInImage;
			}
			c0 = colorAddition(c0, colorChange); // interpolate color
		}
	}
	else
	{
		// the case where slope is in between 0 and <=1

		// if x2 < x1 swap them
		if (v2.vec.x < v1.vec.x)
		{
			LineVec4 tmp = v1;
			v1 = v2;
			v2 = tmp;
		}

		movementInImage = -1 ? v2.vec.y < v1.vec.y : 1;

		int y0 = v1.vec.y;
		Color c0 = v1.color;
		int d = (v1.vec.y - v2.vec.y) + (movementInImage * 0.5 * (v2.vec.x - v1.vec.x));
		Color tempDiff = colorDifference(v2.color, v1.color);
		colorChange = colorDivision(tempDiff, v2.vec.x - v1.vec.x);

		// we got our constants now iterate through whole x values
		for (int i = v1.vec.x; i <= v2.vec.x; i++)
		{
			Color rounded_c0;
			rounded_c0.r = (int)(c0.r + 0.5);
			rounded_c0.g = (int)(c0.g + 0.5);
			rounded_c0.b = (int)(c0.b + 0.5);

			image[i][y0] = rounded_c0;

			// choose between y0 and y0+1 (NE or E)

			if (d * movementInImage <= 0)
			{
				d += (v1.vec.y - v2.vec.y); // move horizontally only
			}
			else
			{ // move diagonally (NE)
				d += (v1.vec.y - v2.vec.y) + (movementInImage * (v2.vec.x - v1.vec.x));
				y0 += movementInImage;
			}
			c0 = colorAddition(c0, colorChange); // interpolate color
		}
	}
}

double f01(double x, double y, LineVec4 &v0, LineVec4 &v1)
{
	return (v0.vec.y - v1.vec.y) * x + (v1.vec.x - v0.vec.x) * y + v0.vec.x * v1.vec.y - v0.vec.y * v1.vec.x; // x(y0 - y1) + y(x1 - x0) + x0y1 - y0x1
}

double f12(double x, double y, LineVec4 &v1, LineVec4 &v2)
{
	return (v1.vec.y - v2.vec.y) * x + (v2.vec.x - v1.vec.x) * y + v1.vec.x * v2.vec.y - v1.vec.y * v2.vec.x; // x(y1 - y2) + y(x2 - x1) + x1y2 - y1x2
}

double f20(double x, double y, LineVec4 &v2, LineVec4 &v0)
{
	return (v2.vec.y - v0.vec.y) * x + (v0.vec.x - v2.vec.x) * y + v2.vec.x * v0.vec.y - v2.vec.y * v0.vec.x; // x(y2 - y0) + y(x0 - x2) + x2y0 - y2x0
}

void triangleRasterization(vector<vector<Color>> &image, LineVec4 &v0, LineVec4 &v1, LineVec4 &v2, int &maxFrameX, int &maxFrameY)
{
	// implementation from rasterization slides page 30

	double x_min = std::min(v0.vec.x, std::min(v1.vec.x, v2.vec.x));
	if (x_min < 0)
	{
		x_min = 0;
	}
	if (x_min > maxFrameX - 1)
	{
		x_min = maxFrameX - 1;
	}
	double x_max = std::max(v0.vec.x, std::max(v1.vec.x, v2.vec.x));
	if (x_max < 0)
	{
		x_max = 0;
	}
	if (x_max > maxFrameX - 1)
	{
		x_max = maxFrameX - 1;
	}
	double y_min = std::min(v0.vec.y, std::min(v1.vec.y, v2.vec.y));
	if (y_min < 0)
	{
		y_min = 0;
	}
	if (y_min > maxFrameY - 1)
	{
		y_min = maxFrameY - 1;
	}
	double y_max = std::max(v0.vec.y, std::max(v1.vec.y, v2.vec.y));
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
			double alpha = f12(x, y, v1, v2) / f12(v0.vec.x, v0.vec.y, v1, v2);
			double beta = f20(x, y, v2, v0) / f20(v1.vec.x, v1.vec.y, v2, v0);
			double gamma = f01(x, y, v0, v1) / f01(v2.vec.x, v2.vec.y, v0, v1);

			if (alpha >= 0 && beta >= 0 && gamma >= 0)
			{
				Color color;
				color.r = alpha * v0.color.r + beta * v1.color.r + gamma * v2.color.r;
				color.g = alpha * v0.color.g + beta * v1.color.g + gamma * v2.color.g;
				color.b = alpha * v0.color.b + beta * v1.color.b + gamma * v2.color.b;

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
	// TODO: Implement this function

	bool cullingFlag = this->cullingEnabled;

	// 1- calculate camera transformation matrix
	Matrix4 camTransMatrix;
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

	// formula from viewving transformation slides page 8
	camTransMatrix = multiplyMatrixWithMatrix(rMatrix, tMatrix);

	// 2- calculate orthographic or perspective projection transformation matrix
	Matrix4 projTransMatrix;
	if (camera->projectionType)
	{
		// calculate perspective projection transformation matrix
		// formula from viewving transformation slides page 26
		// TODO not sure about implementation of perspective divide check later didn't implemented here I believe .s
		double temp[4][4] = {{2 * camera->near / (camera->right - camera->left), 0, (camera->right + camera->left) / (camera->right - camera->left), 0},
							 {0, 2 * camera->near / (camera->top - camera->bottom), (camera->top + camera->bottom) / (camera->top - camera->bottom), 0},
							 {0, 0, -(camera->far + camera->near) / (camera->far - camera->near), -2 * camera->far * camera->near / (camera->far - camera->near)},
							 {0, 0, -1, 0}};

		projTransMatrix = Matrix4(temp);
	}
	else
	{
		// calculate orthographic projection transformation matrix
		// formula from viewing transformation slides page 14
		double temp[4][4] = {{2 / (camera->right - camera->left), 0, 0, -(camera->right + camera->left) / (camera->right - camera->left)},
							 {0, 2 / (camera->top - camera->bottom), 0, -(camera->top + camera->bottom) / (camera->top - camera->bottom)},
							 {0, 0, -2 / (camera->far - camera->near), -(camera->far + camera->near) / (camera->far - camera->near)},
							 {0, 0, 0, 1}};
		projTransMatrix = Matrix4(temp);
	}

	// 3- calculate viewport transformation matrix
	double tmpMvp[4][4] = {{camera->horRes / (double)2, 0, 0, (camera->horRes - 1) / (double)2},
						   {0, camera->verRes / (double)2, 0, (camera->verRes - 1) / (double)2},
						   {0, 0, 0.5, 0.5},
						   {0, 0, 0, 1}};
	Matrix4 viewportMatrix = Matrix4(tmpMvp);

	// iterate through all meshes
	for (auto &currentMesh : this->meshes)
	{
		// 4- implement modelling transformation
		Matrix4 modellingMatrix = getIdentityMatrix();
		// iterate through all modelling transformations of current mesh
		for (int i = 0; i < currentMesh->numberOfTransformations; i++)
		{
			if (currentMesh->transformationTypes[i] == 'r')
			{
				Vec3 tmpVec = {this->rotations[currentMesh->transformationIds[i] - 1]->ux,
							   this->rotations[currentMesh->transformationIds[i] - 1]->uy,
							   this->rotations[currentMesh->transformationIds[i] - 1]->uz};
				tmpVec = normalizeVec3(tmpVec);

				// smallest component of v
				double min = std::min(abs(tmpVec.y), abs(tmpVec.z));
				if (abs(tmpVec.x) < min)
				{
					min = abs(tmpVec.x);
				}

				Vec3 v, w;

				if (min == abs(tmpVec.x))
				{
					v.x = 0;
					v.y = -tmpVec.z;
					v.z = tmpVec.y;
				}
				else if (min == abs(tmpVec.y))
				{
					v.x = -tmpVec.z;
					v.y = 0;
					v.z = tmpVec.x;
				}
				else
				{
					v.x = -tmpVec.y;
					v.y = tmpVec.x;
					v.z = 0;
				}

				v = normalizeVec3(v);
				w = crossProductVec3(tmpVec, v);

				double rotationMatrix[4][4] = {{tmpVec.x, tmpVec.y, tmpVec.z, 0},
											   {v.x, v.y, v.z, 0},
											   {w.x, w.y, w.z, 0},
											   {0, 0, 0, 1}};
				double inverseRotationMatrix[4][4] = {{tmpVec.x, v.x, w.x, 0},
													  {tmpVec.y, v.y, w.y, 0},
													  {tmpVec.z, v.z, w.z, 0},
													  {0, 0, 0, 1}};

				// now we are on x-axis no problem if we rotate around x-axis
				// TODO check correctness of this
				double rotationAngleTheta = rotationAngleTheta * M_PI / 180;
				double rotationMatrix2[4][4] = {{1, 0, 0, 0},
												{0, cos(rotationAngleTheta), -sin(rotationAngleTheta), 0},
												{0, sin(rotationAngleTheta), cos(rotationAngleTheta), 0},
												{0, 0, 0, 1}};

				Matrix4 firstRotationMatrix = multiplyMatrixWithMatrix(Matrix4(rotationMatrix2), Matrix4(rotationMatrix));
				Matrix4 secondRotationMatrix = multiplyMatrixWithMatrix(Matrix4(inverseRotationMatrix), firstRotationMatrix);
				modellingMatrix = multiplyMatrixWithMatrix(secondRotationMatrix, modellingMatrix);
			}

			if (currentMesh->transformationTypes[i] == 't')
			{
				double translationMatrix[4][4] = {{1, 0, 0, this->translations[currentMesh->transformationIds[i] - 1]->tx},
												  {0, 1, 0, this->translations[currentMesh->transformationIds[i] - 1]->ty},
												  {0, 0, 1, this->translations[currentMesh->transformationIds[i] - 1]->tz},
												  {0, 0, 0, 1}};
				modellingMatrix = multiplyMatrixWithMatrix(modellingMatrix, Matrix4(translationMatrix));
			}

			if (currentMesh->transformationTypes[i] == 's')
			{
				double scalingMatrix[4][4] = {{this->scalings[currentMesh->transformationIds[i] - 1]->sx, 0, 0, 0},
											  {0, this->scalings[currentMesh->transformationIds[i] - 1]->sy, 0, 0},
											  {0, 0, this->scalings[currentMesh->transformationIds[i] - 1]->sz, 0},
											  {0, 0, 0, 1}};
				modellingMatrix = multiplyMatrixWithMatrix(modellingMatrix, Matrix4(scalingMatrix));
			}
		}
		Matrix4 allMatrices = multiplyMatrixWithMatrix(camTransMatrix, modellingMatrix);
		allMatrices = multiplyMatrixWithMatrix(allMatrices, projTransMatrix); // projection matrix * camera matrix * modelling matrix
		for (auto &currentTriangle : currentMesh->triangles)
		{
			// create projected triangle
			Vec3 *vector0 = this->vertices[currentTriangle.vertexIds[0] - 1];
			Vec3 *vector1 = this->vertices[currentTriangle.vertexIds[1] - 1];
			Vec3 *vector2 = this->vertices[currentTriangle.vertexIds[2] - 1];

			// check culling if culling is enabled
			if (cullingFlag)
			{
				// check if mesh is backfacing the camera
				// in this case just skip this mesh
				// calculate normal vector of triangle
				Vec3 normal = crossProductVec3(subtractVec3(*vector1, *vector0),
											   subtractVec3(*vector2, *vector0));
				normal = normalizeVec3(normal);
				// TODO check this implementation
				// double resultingDotProduct = dotProductVec3(camera->position, normal);
				double resultingDotProduct = dotProductVec3(normal, *vector0);
				if (resultingDotProduct > 0)
				{
					continue; // skip the mesh
				}
			}

			// raw vectors before transformation
			Vec4 unprocessedVector0 = Vec4(vector0->x, vector0->y, vector0->z, 1, vector0->colorId);
			Vec4 unprocessedVector1 = Vec4(vector1->x, vector1->y, vector1->z, 1, vector1->colorId);
			Vec4 unprocessedVector2 = Vec4(vector2->x, vector2->y, vector2->z, 1, vector2->colorId);

			// transformed vectors after transformation we will use these vectors for clipping
			Vec4 transformedVector0 = multiplyMatrixWithVec4(allMatrices, unprocessedVector0);
			Vec4 transformedVector1 = multiplyMatrixWithVec4(allMatrices, unprocessedVector1);
			Vec4 transformedVector2 = multiplyMatrixWithVec4(allMatrices, unprocessedVector2);

			// we need to use clipping only if we are using wireframe mode
			// check wireframe => 0 for wireframe 1 for solid
			if (currentMesh->type)
			{
				// TODO implement solid mode

				// perspective divison time
				perspectiveDivide(transformedVector0);
				perspectiveDivide(transformedVector1);
				perspectiveDivide(transformedVector2);

				// implement vp transformations now (perspective divide already done)
				transformedVector0 = multiplyMatrixWithVec4(viewportMatrix, transformedVector0);
				transformedVector1 = multiplyMatrixWithVec4(viewportMatrix, transformedVector1);
				transformedVector2 = multiplyMatrixWithVec4(viewportMatrix, transformedVector2);

				// TODO implement rasterization
				LineVec4 line01 = LineVec4(transformedVector0, *this->colorsOfVertices[vector0->colorId - 1]);
				LineVec4 line12 = LineVec4(transformedVector1, *this->colorsOfVertices[vector1->colorId - 1]);
				LineVec4 line20 = LineVec4(transformedVector2, *this->colorsOfVertices[vector2->colorId - 1]);
				triangleRasterization(this->image, line01, line12, line20, camera->horRes, camera->verRes);
			}
			else
			{
				// TODO implement wireframe mode
				// Liang-Barsky algorithm (faster than Cohen-Sutherland)
				// we got 3 lines for each triangle

				// perspective divison time
				perspectiveDivide(transformedVector0);
				perspectiveDivide(transformedVector1);
				perspectiveDivide(transformedVector2);

				// Create line 0-1 and 1-2 and 2-0 (logical triangle)
				LineVec4 line01 = LineVec4(transformedVector0, *this->colorsOfVertices[vector0->colorId - 1]);
				LineVec4 line10 = LineVec4(transformedVector1, *this->colorsOfVertices[vector1->colorId - 1]);
				bool visibility01 = lineClipping(line01, line10);

				LineVec4 line12 = LineVec4(transformedVector1, *this->colorsOfVertices[vector1->colorId - 1]);
				LineVec4 line21 = LineVec4(transformedVector2, *this->colorsOfVertices[vector2->colorId - 1]);
				bool visibility12 = lineClipping(line12, line21);

				LineVec4 line20 = LineVec4(transformedVector2, *this->colorsOfVertices[vector2->colorId - 1]);
				LineVec4 line02 = LineVec4(transformedVector0, *this->colorsOfVertices[vector0->colorId - 1]);
				bool visibility20 = lineClipping(line20, line02);

				// implement vp transformations now (perspective divide already done)
				transformedVector0 = multiplyMatrixWithVec4(viewportMatrix, transformedVector0);
				transformedVector1 = multiplyMatrixWithVec4(viewportMatrix, transformedVector1);
				transformedVector2 = multiplyMatrixWithVec4(viewportMatrix, transformedVector2);

				// TODO implement rasterization
				if (visibility01)
				{
					lineRasterization(this->image, line01, line10);
				}
				if (visibility12)
				{
					lineRasterization(this->image, line12, line21);
				}
				if (visibility20)
				{
					lineRasterization(this->image, line20, line02);
				}
			}
			// TODO implement z-buffering have no idea how to do this
		}
	}
}