#ifndef _SCENE_H_
#define _SCENE_H_
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Matrix4.h"

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color>> image;
	std::vector<std::vector<double>> depth;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	// will be used for liang-barsky line clipping algorithm
	struct LineVec4
	{
		Vec4 vertex1, vertex2;

		LineVec4(Vec4 &v, Vec4 &c) : vertex1(v), vertex2(c) {}
	};

	Scene(const char *xmlPath);

	void assignColorToPixel(int i, int j, Color c);
	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName);
	void forwardRenderingPipeline(Camera *camera);

	Matrix4 modellingTransformations(Mesh *mesh);
	Matrix4 orthTransformation(Camera *cam);
	Matrix4 perspectiveTransformation(Camera *cam);
	Vec4 perspectiveDivide(Vec4 &v);
	Matrix4 camTransformation(Camera *cam);
	Matrix4 viewTransformation(Camera *cam);
	bool lineClipping(LineVec4 &vertices);
	bool isCulled(Vec4 vertices[3]);
	void triangleRasterization(vector<Vec3> &vertices, int maxFrameX, int maxFrameY);
	double f_alpha(vector<Vec3> &vertices, double x, double y);
	double f_beta(vector<Vec3> &vertices, double x, double y);
	double f_gamma(vector<Vec3> &vertices, double x, double y);
	void lineRasterization(std::pair<Vec3, Vec3> &vertices);
};

#endif
