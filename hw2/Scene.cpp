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
	command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	// TODO: Implement this function

	bool cullingFlag = this->cullingEnabled;

	// 1- calculate camera transformation matrix
	Matrix4 camTransMatrix;
	// formula from viewving transformation slides page 6
	double tMatrix[4][4] = {{1, 0, 0, -camera->position.x},
							{0, 1, 0, -camera->position.y},
							{0, 0, 1, -camera->position.z},
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
		// formula from viewving transformation slides page 14
		double temp[4][4] = {{2 / (camera->right - camera->left), 0, 0, -(camera->right + camera->left) / (camera->right - camera->left)},
							 {0, 2 / (camera->top - camera->bottom), 0, -(camera->top + camera->bottom) / (camera->top - camera->bottom)},
							 {0, 0, -2 / (camera->far - camera->near), -(camera->far + camera->near) / (camera->far - camera->near)},
							 {0, 0, 0, 1}};
		projTransMatrix = Matrix4(temp);
	}

	// 3- calculate viewport transformation matrix
	double tmpMvp[4][4] = {{camera->horRes / 2, 0, 0, (camera->horRes - 1) / 2},
						   {0, camera->verRes / 2, 0, (camera->verRes - 1) / 2},
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
				// TODO implement rotation
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

			// raw vectors before transformation
			Vec4 unprocessedVector0 = Vec4(vector0->x, vector0->y, vector0->z, 1, vector0->colorId);
			Vec4 unprocessedVector1 = Vec4(vector1->x, vector1->y, vector1->z, 1, vector1->colorId);
			Vec4 unprocessedVector2 = Vec4(vector2->x, vector2->y, vector2->z, 1, vector2->colorId);

			// transformed vectors after transformation we will use these vectors for clipping
			Vec4 transformedVector0 = multiplyMatrixWithVec4(allMatrices, unprocessedVector0);
			Vec4 transformedVector1 = multiplyMatrixWithVec4(allMatrices, unprocessedVector1);
			Vec4 transformedVector2 = multiplyMatrixWithVec4(allMatrices, unprocessedVector2);

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

			// TODO implement clipping
			// we need to use clipping only if we are using wireframe mode
			// check wireframe => 0 for wireframe 1 for solid
			if (currentMesh->type)
			{
				// TODO implement solid mode
			}
			else
			{
				// TODO implement wireframe mode
				// Liang-Barsky algorithm (faster than Cohen-Sutherland)
				// we got 3 lines for each triangle

				Vec4 *line1 = new Vec4(transformedVector0);

				// TODO implement rasterization
			}
			// TODO implement z-buffering
		}
	}
}