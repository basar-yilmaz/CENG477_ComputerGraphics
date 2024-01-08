#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <GL/glew.h>    // The GL Header File
#include <GL/gl.h>      // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <ft2build.h>
#include <ctime>
#include FT_FREETYPE_H
#include <string>

#define BUFFER_OFFSET(i) ((char *)NULL + (i))

using namespace std;

GLuint gProgram[6];
// GLint gIntensityLoc;
// float gIntensity = 1000;
int gWidth = 1024, gHeight = 800;

#include <vector>
#include <GL/glew.h>

static float angle = 0;
static float translateX = 0;
static float translateY = 0;
static int jumpFlag = 1;

enum Direction
{
    NORMAL = 0,
    RIGHT = 1,
    LEFT = 2
};

static Direction direction = NORMAL;

static int score = 0;
struct Vertex
{
    Vertex(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) {}
    GLfloat x, y, z;
};

struct Texture
{
    Texture(GLfloat inU, GLfloat inV) : u(inU), v(inV) {}
    GLfloat u, v;
};

struct Normal
{
    Normal(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) {}
    GLfloat x, y, z;
};

struct Face
{
    Face(int v[], int t[], int n[])
    {
        vIndex[0] = v[0];
        vIndex[1] = v[1];
        vIndex[2] = v[2];
        tIndex[0] = t[0];
        tIndex[1] = t[1];
        tIndex[2] = t[2];
        nIndex[0] = n[0];
        nIndex[1] = n[1];
        nIndex[2] = n[2];
    }
    GLuint vIndex[3], tIndex[3], nIndex[3];
};

// vector<Vertex> gVertices;
// vector<Texture> gTextures;
// vector<Normal> gNormals;
// vector<Face> gFaces;
// GLuint gVertexAttribBuffer, gIndexBuffer;

GLuint gTextVBO;
// int gVertexDataSizeInBytes, gNormalDataSizeInBytes;

/// Holds all state information relevant to a character as loaded using FreeType
struct Character
{
    GLuint TextureID;   // ID handle of the glyph texture
    glm::ivec2 Size;    // Size of glyph
    glm::ivec2 Bearing; // Offset from baseline to left/top of glyph
    GLuint Advance;     // Horizontal offset to advance to next glyph
};

std::map<GLchar, Character> Characters;

class Object;
bool ParseObj2(const string &fileName, Object &obj);
bool ReadDataFromFile(const string &fileName, string &data);
void createVS(GLuint &program, const string &filename);
void createFS(GLuint &program, const string &filename);

class Object
{

public:
    GLuint vertexBuffer;
    GLuint indexBuffer;
    GLuint program_id;
    int gVertexDataSizeInBytes;
    int gNormalDataSizeInBytes;

    int min_X, max_X;
    int min_Y, max_Y;
    int min_Z, max_Z;

    std::string name;

    std::vector<Vertex> vertices;
    std::vector<Normal> normals;
    std::vector<Texture> texture;
    std::vector<Face> faces;

    Object(GLuint programID, string name) : vertexBuffer(0), indexBuffer(0), program_id(programID), name(name) {}

    void initVBO()
    {
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);
        // cout << "program_id = " << this->program_id << endl;
        // cout << "program_name = " << this->name << endl;
        assert(glGetError() == GL_NONE);

        gluErrorString(glGetError());

        glGenBuffers(1, &vertexBuffer);
        glGenBuffers(1, &indexBuffer);

        gluErrorString(glGetError());

        assert(vertexBuffer > 0 && indexBuffer > 0);

        glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);

        unsigned verticesSize = (this->vertices).size();
        unsigned normalsSize = (this->normals).size();
        unsigned facesSize = (this->faces).size();

        gVertexDataSizeInBytes = verticesSize * 3 * sizeof(GLfloat);
        gNormalDataSizeInBytes = normalsSize * 3 * sizeof(GLfloat);
        int indexDataSizeInBytes = facesSize * 3 * sizeof(GLuint);
        GLfloat *vertexData = new GLfloat[verticesSize * 3];
        GLfloat *normalData = new GLfloat[normalsSize * 3];
        GLuint *indexData = new GLuint[facesSize * 3];

        float minX = 1e6, maxX = -1e6;
        float minY = 1e6, maxY = -1e6;
        float minZ = 1e6, maxZ = -1e6;

        for (int i = 0; i < verticesSize; ++i)
        {
            vertexData[3 * i] = (this->vertices)[i].x;
            vertexData[3 * i + 1] = (this->vertices)[i].y;
            vertexData[3 * i + 2] = (this->vertices)[i].z;

            minX = std::min(minX, (this->vertices)[i].x);
            maxX = std::max(maxX, (this->vertices)[i].x);
            minY = std::min(minY, (this->vertices)[i].y);
            maxY = std::max(maxY, (this->vertices)[i].y);
            minZ = std::min(minZ, (this->vertices)[i].z);
            maxZ = std::max(maxZ, (this->vertices)[i].z);
        }

        // std::cout << "minX = " << minX << std::endl;
        // std::cout << "maxX = " << maxX << std::endl;
        // std::cout << "minY = " << minY << std::endl;
        // std::cout << "maxY = " << maxY << std::endl;
        // std::cout << "minZ = " << minZ << std::endl;
        // std::cout << "maxZ = " << maxZ << std::endl;

        this->min_X = minX;
        this->max_X = maxX;
        this->min_Y = minY;
        this->max_Y = maxY;
        this->min_Z = minZ;
        this->max_Z = maxZ;

        for (int i = 0; i < normalsSize; ++i)
        {
            normalData[3 * i] = (this->normals)[i].x;
            normalData[3 * i + 1] = (this->normals)[i].y;
            normalData[3 * i + 2] = (this->normals)[i].z;
        }

        for (int i = 0; i < facesSize; ++i)
        {
            indexData[3 * i] = (this->faces)[i].vIndex[0];
            indexData[3 * i + 1] = (this->faces)[i].vIndex[1];
            indexData[3 * i + 2] = (this->faces)[i].vIndex[2];
        }

        glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, 0, GL_STATIC_DRAW);
        glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
        glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

        // done copying; can free now
        delete[] vertexData;
        delete[] normalData;
        delete[] indexData;

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
    }

    void initShaders()
    {
        gProgram[this->program_id] = glCreateProgram();

        if (this->program_id > 2)
        {
            createVS(gProgram[this->program_id], "cube_vert.glsl");
            createFS(gProgram[this->program_id], "cube_frag.glsl");
        }
        else
        {
            createVS(gProgram[this->program_id], name + "_vert.glsl");
            createFS(gProgram[this->program_id], name + "_frag.glsl");
        }

        glBindAttribLocation(gProgram[this->program_id], 0, "inVertex"); //
        glBindAttribLocation(gProgram[this->program_id], 1, "inNormal");

        glLinkProgram(gProgram[this->program_id]);
        glUseProgram(gProgram[program_id]);

        cout << "shader_program_id = " << this->program_id << endl;
        cout << "program_name = " << this->name << endl;
    }

    void draw()
    {
        glClear(GL_DEPTH_BUFFER_BIT);
        glUseProgram(gProgram[this->program_id]);

        // cout << "program_id = " << this->program_id << endl;
        glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));

        glDrawElements(GL_TRIANGLES, faces.size() * 3, GL_UNSIGNED_INT, 0);
    }
};

Object bunny(0, "bunny"); // programId = 0
Object quad(2, "quad");   // programId = 2
Object cube1(3, "cube1"); // programId = 3
Object cube2(4, "cube2"); // programId = 4
Object cube3(5, "cube3"); // programId = 5

/*
bool ParseObj(const string &fileName)
{
    fstream myfile;

    // Open the input
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            stringstream str(curLine);
            GLfloat c1, c2, c3;
            GLuint index[9];
            string tmp;

            if (curLine.length() >= 2)
            {
                if (curLine[0] == '#') // comment
                {
                    continue;
                }
                else if (curLine[0] == 'v')
                {
                    if (curLine[1] == 't') // texture
                    {
                        str >> tmp; // consume "vt"
                        str >> c1 >> c2;
                        gTextures.push_back(Texture(c1, c2));
                    }
                    else if (curLine[1] == 'n') // normal
                    {
                        str >> tmp; // consume "vn"
                        str >> c1 >> c2 >> c3;
                        gNormals.push_back(Normal(c1, c2, c3));
                    }
                    else // vertex
                    {
                        str >> tmp; // consume "v"
                        str >> c1 >> c2 >> c3;
                        gVertices.push_back(Vertex(c1, c2, c3));
                    }
                }
                else if (curLine[0] == 'f') // face
                {
                    str >> tmp; // consume "f"
                    char c;
                    int vIndex[3], nIndex[3], tIndex[3];
                    str >> vIndex[0];
                    str >> c >> c; // consume "//"
                    str >> nIndex[0];
                    str >> vIndex[1];
                    str >> c >> c; // consume "//"
                    str >> nIndex[1];
                    str >> vIndex[2];
                    str >> c >> c; // consume "//"
                    str >> nIndex[2];

                    assert(vIndex[0] == nIndex[0] &&
                           vIndex[1] == nIndex[1] &&
                           vIndex[2] == nIndex[2]); // a limitation for now

                    // make indices start from 0
                    for (int c = 0; c < 3; ++c)
                    {
                        vIndex[c] -= 1;
                        nIndex[c] -= 1;
                        tIndex[c] -= 1;
                    }

                    gFaces.push_back(Face(vIndex, tIndex, nIndex));
                }
                else
                {
                    cout << "Ignoring unidentified line in obj file: " << curLine << endl;
                }
            }

            // data += curLine;
            if (!myfile.eof())
            {
                // data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    assert(gVertices.size() == gNormals.size());

    return true;
}
*/
bool ParseObj2(const string &fileName, Object &obj)
{
    fstream myfile;

    // Open the input
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            stringstream str(curLine);
            GLfloat c1, c2, c3;
            GLuint index[9];
            string tmp;

            if (curLine.length() >= 2)
            {
                if (curLine[0] == '#') // comment
                {
                    continue;
                }
                else if (curLine[0] == 'v')
                {
                    if (curLine[1] == 't') // texture
                    {
                        str >> tmp; // consume "vt"
                        str >> c1 >> c2;
                        obj.texture.push_back(Texture(c1, c2));
                    }
                    else if (curLine[1] == 'n') // normal
                    {
                        str >> tmp; // consume "vn"
                        str >> c1 >> c2 >> c3;
                        obj.normals.push_back(Normal(c1, c2, c3));
                    }
                    else // vertex
                    {
                        str >> tmp; // consume "v"
                        str >> c1 >> c2 >> c3;
                        obj.vertices.push_back(Vertex(c1, c2, c3));
                    }
                }
                else if (curLine[0] == 'f') // face
                {
                    str >> tmp; // consume "f"
                    char c;
                    int vIndex[3], nIndex[3], tIndex[3];
                    str >> vIndex[0];
                    str >> c >> c; // consume "//"
                    str >> nIndex[0];
                    str >> vIndex[1];
                    str >> c >> c; // consume "//"
                    str >> nIndex[1];
                    str >> vIndex[2];
                    str >> c >> c; // consume "//"
                    str >> nIndex[2];

                    assert(vIndex[0] == nIndex[0] &&
                           vIndex[1] == nIndex[1] &&
                           vIndex[2] == nIndex[2]); // a limitation for now

                    // make indices start from 0
                    for (int c = 0; c < 3; ++c)
                    {
                        vIndex[c] -= 1;
                        nIndex[c] -= 1;
                        tIndex[c] -= 1;
                    }

                    obj.faces.push_back(Face(vIndex, tIndex, nIndex));
                }
                else
                {
                    cout << "Ignoring unidentified line in obj file: " << curLine << endl;
                }
            }

            // data += curLine;
            if (!myfile.eof())
            {
                // data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    assert(obj.vertices.size() == obj.normals.size());

    return true;
}

bool ReadDataFromFile(
    const string &fileName, ///< [in]  Name of the shader file
    string &data)           ///< [out] The contents of the file
{
    fstream myfile;

    // Open the input
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            data += curLine;
            if (!myfile.eof())
            {
                data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    return true;
}

void createVS(GLuint &program, const string &filename)
{
    string shaderSource;

    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar *shader = (const GLchar *)shaderSource.c_str();

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &shader, &length);
    glCompileShader(vs);

    char output[1024] = {0};
    glGetShaderInfoLog(vs, 1024, &length, output);
    printf("VS compile log: %s\n", output);

    glAttachShader(program, vs);
}

void createFS(GLuint &program, const string &filename)
{
    string shaderSource;

    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar *shader = (const GLchar *)shaderSource.c_str();

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &shader, &length);
    glCompileShader(fs);

    char output[1024] = {0};
    glGetShaderInfoLog(fs, 1024, &length, output);
    printf("FS compile log: %s\n", output);

    glAttachShader(program, fs);
}

void initTextShaders()
{
    gProgram[1] = glCreateProgram();

    createVS(gProgram[1], "vert_text.glsl");
    createFS(gProgram[1], "frag_text.glsl");

    glBindAttribLocation(gProgram[1], 2, "vertex"); // text

    glLinkProgram(gProgram[1]);
}

void initFonts(int windowWidth, int windowHeight)
{
    // Set OpenGL options
    // glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(windowWidth), 0.0f, static_cast<GLfloat>(windowHeight));
    glUseProgram(gProgram[1]);
    glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "projection"), 1, GL_FALSE, glm::value_ptr(projection));

    // FreeType
    FT_Library ft;
    // All functions return a value different than 0 whenever an error occurred
    if (FT_Init_FreeType(&ft))
    {
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
    }

    // Load font as face
    FT_Face face;
    if (FT_New_Face(ft, "/usr/share/fonts/truetype/ubuntu/Ubuntu-M.ttf", 0, &face))
    {
        std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
    }

    // Set size to load glyphs as
    FT_Set_Pixel_Sizes(face, 0, 48);

    // Disable byte-alignment restriction
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    // Load first 128 characters of ASCII set
    for (GLubyte c = 0; c < 128; c++)
    {
        // Load character glyph
        if (FT_Load_Char(face, c, FT_LOAD_RENDER))
        {
            std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
            continue;
        }
        // Generate texture
        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(
            GL_TEXTURE_2D,
            0,
            GL_RED,
            face->glyph->bitmap.width,
            face->glyph->bitmap.rows,
            0,
            GL_RED,
            GL_UNSIGNED_BYTE,
            face->glyph->bitmap.buffer);
        // Set texture options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // Now store character for later use
        Character character = {
            texture,
            glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            face->glyph->advance.x};
        Characters.insert(std::pair<GLchar, Character>(c, character));
    }

    glBindTexture(GL_TEXTURE_2D, 0);
    // Destroy FreeType once we're finished
    FT_Done_Face(face);
    FT_Done_FreeType(ft);

    //
    // Configure VBO for texture quads
    //
    glGenBuffers(1, &gTextVBO);
    glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void init()
{

    ParseObj2("bunny.obj", bunny);
    ParseObj2("quad.obj", quad);
    ParseObj2("cube.obj", cube1);
    ParseObj2("cube.obj", cube2);
    ParseObj2("cube.obj", cube3);
    gluErrorString(glGetError());

    glEnable(GL_DEPTH_TEST);

    initTextShaders();
    bunny.initShaders();
    quad.initShaders();
    cube1.initShaders();
    cube2.initShaders();
    cube3.initShaders();

    bunny.initVBO();
    quad.initVBO();
    cube1.initVBO();
    cube2.initVBO();
    cube3.initVBO();

    cout << cube1.faces.size() << endl;

    cout << quad.name << endl;
    cout << bunny.name << endl;
    cout << cube1.name << endl;

    initFonts(gWidth, gHeight);
}

void renderText(const std::string &text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color)
{
    // Activate corresponding render state
    glUseProgram(gProgram[1]);
    glUniform3f(glGetUniformLocation(gProgram[1], "textColor"), color.x, color.y, color.z);
    glActiveTexture(GL_TEXTURE0);

    // Iterate through all characters
    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++)
    {
        Character ch = Characters[*c];

        GLfloat xpos = x + ch.Bearing.x * scale;
        GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

        GLfloat w = ch.Size.x * scale;
        GLfloat h = ch.Size.y * scale;

        // Update VBO for each character
        GLfloat vertices[6][4] = {
            {xpos, ypos + h, 0.0, 0.0},
            {xpos, ypos, 0.0, 1.0},
            {xpos + w, ypos, 1.0, 1.0},

            {xpos, ypos + h, 0.0, 0.0},
            {xpos + w, ypos, 1.0, 1.0},
            {xpos + w, ypos + h, 1.0, 0.0}};

        // Render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.TextureID);

        // Update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

        // glBindBuffer(GL_ARRAY_BUFFER, 0);

        // Render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)

        x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }

    glBindTexture(GL_TEXTURE_2D, 0);
}

int findIndexOfSmallest(float values[])
{
    float smallest = values[0];
    int smallestIndex = 0;

    for (int i = 1; i < 3; ++i)
    {
        if (values[i] < smallest)
        {
            smallest = values[i];
            smallestIndex = i;
        }
    }

    return smallestIndex;
}

int generateAndPrintRandomFloats()
{
    // Seed the random number generator
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    // Generate 3 random floats
    float randomFloats[3];
    for (int i = 0; i < 3; ++i)
    {
        randomFloats[i] = static_cast<float>(std::rand()) / RAND_MAX;
    }

    // std::cout << randomFloats[0] << endl;
    // std::cout << randomFloats[1] << endl;
    // std::cout << randomFloats[2] << endl;

    // Find the index of the smallest value
    int smallestIndex = findIndexOfSmallest(randomFloats);

    return smallestIndex;
}

#define MAX_MOVE 9
#define BUNNY_SCALE 1.25

static glm::vec3 RED = glm::vec3(0.8, 0.1, 0.1);
static glm::vec3 YELLOW = glm::vec3(0.8, 0.8, 0.1);
static int yellow_index = 0;

static float quadX = 0;
static float quadY = 0;
static float quadZ = 0;

static float cubeX = 0;
static float cubeY = 0;
static float cubeZ = 0;

static float jump_multiplier = 0.1f;

static float accumulatior = 0.f;

static float speed_up = 0.2f;

static int draw_flag = -1;

static bool deadFlag = false;

static float bunnyRotationAngle = 0.f;

static bool scoreFlag = false;
static bool hasRotated = false;
static float rotationThreshold = 360.f;
void display()
{
    glClearColor(0, 0, 0, 1);
    glClearDepth(1.0f);
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    const float epsilon = 0.2f;
    if (direction == RIGHT)
    {
        if (cubeZ > 69.5f - epsilon && cubeZ < 70.f + epsilon)
        {
            draw_flag = 2;
            if (yellow_index == 1 && !deadFlag)
            {
                scoreFlag = true;
                score += 1000;
            }
            else
            {
                deadFlag = true;
                speed_up = 0;
                jumpFlag = 0;
            }
        }
    }
    else if (direction == LEFT)
    {
        if (cubeZ > 69.5f - epsilon && cubeZ < 70.f + epsilon)
        {
            draw_flag = 3;
            if (yellow_index == 2 && !deadFlag)
            {
                scoreFlag = true;
                score += 1000;
            }
            else
            {
                deadFlag = true;
                speed_up = 0;
                jumpFlag = 0;
            }
        }
    }
    else
    {
        // NORMAL
        if (cubeZ > 69.5f - epsilon && cubeZ < 70.f + epsilon)
        {
            draw_flag = 1;
            if (yellow_index == 0 && !deadFlag)
            {
                scoreFlag = true;
                score += 1000;
            }
            else
            {
                deadFlag = true;
                speed_up = 0;
                jumpFlag = 0;
            }
        }
    }

    // QUAD PART
    glUseProgram(gProgram[quad.program_id]);

    float scaleFrag = 0.2f; // Adjust this to control the size of the checker squares
    glUniform1f(glGetUniformLocation(gProgram[quad.program_id], "scale"), scaleFrag);
    glm::vec3 color1_quad = glm::vec3(0.0, 0, 0.0);   // Black
    glm::vec3 color2_quad = glm::vec3(1.0, 1.f, 1.0); // White
    glUniform3fv(glGetUniformLocation(gProgram[quad.program_id], "color1"), 1, glm::value_ptr(color1_quad));
    glUniform3fv(glGetUniformLocation(gProgram[quad.program_id], "color2"), 1, glm::value_ptr(color2_quad));

    glm::mat4 tTranslate = glm::translate(glm::mat4(1.0), glm::vec3(0.f, 0.f, -10.f));

    glm::mat4 tRotate = glm::rotate<float>(glm::mat4(1.0), glm::radians(-90.f), glm::vec3(1.0, 0.0, 0.f));
    glm::mat4 zScale = glm::scale(glm::mat4(1.0), glm::vec3(10.f, 100.f, 1.0f));

    glm::mat4 movingQuad = glm::translate(glm::mat4(1.0), glm::vec3(0, 0, quadZ));
    glm::mat4 modelMatQuad = movingQuad * tTranslate * tRotate * zScale;

    glm::mat4 viewMatrixQuad = glm::translate(glm::mat4(1.0), glm::vec3(0, -5, -quadZ));

    glm::mat4 modelMatInvQuad = glm::transpose(glm::inverse(modelMatQuad));

    glm::mat4 perspMatQuad = glm::perspective(glm::radians(90.0f), GLfloat(gWidth / gHeight), 0.1f, 200.0f);

    glUniformMatrix4fv(glGetUniformLocation(gProgram[quad.program_id], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMatQuad));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[quad.program_id], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatInvQuad));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[quad.program_id], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMatQuad));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[quad.program_id], "viewingMat"), 1, GL_FALSE, glm::value_ptr(viewMatrixQuad));

    quad.draw();

    // CUBE PART

    glUseProgram(gProgram[cube1.program_id]);

    glm::mat4 cubeTranslate = glm::translate(glm::mat4(1.0), glm::vec3(0.f, 0.f, -150.f));
    glm::mat4 matRyCube1 = glm::rotate<float>(glm::mat4(1.0), (-90. / 180.) * M_PI, glm::vec3(0.0, 1.0, 0.0));
    glm::mat4 matSCube1 = glm::scale(glm::mat4(1.0), glm::vec3(1.f, 3.f, 1.f));

    glm::mat4 tempMove1 = glm::translate(glm::mat4(1.0), glm::vec3(0, 0, cubeZ));

    glm::mat4 modelMatCube1 = tempMove1 * cubeTranslate * matRyCube1 * matSCube1;

    glm::mat4 viewMatrixCube1 = glm::translate(glm::mat4(1.0), glm::vec3(0, -5, cubeZ));
    glm::mat4 modelMatCube1Inv = glm::transpose(glm::inverse(modelMatCube1));
    glm::mat4 perspMatCube1 = glm::perspective(glm::radians(90.0f), GLfloat(gWidth / gHeight), 0.1f, 200.f);

    glm::vec3 color1;

    if (yellow_index == 0)
    {
        color1 = YELLOW;
    }
    else
    {
        color1 = RED;
    }

    glUniformMatrix4fv(glGetUniformLocation(gProgram[cube1.program_id], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMatCube1));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[cube1.program_id], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatCube1Inv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[cube1.program_id], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMatCube1));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[cube1.program_id], "viewingMat"), 1, GL_FALSE, glm::value_ptr(viewMatrixCube1));
    glUniform3fv(glGetUniformLocation(gProgram[cube1.program_id], "cubeColor"), 1, glm::value_ptr(color1));
    if (draw_flag != 1)
        cube1.draw();

    glUseProgram(gProgram[cube2.program_id]);

    glm::vec3 color2;
    if (yellow_index == 1)
    {
        color2 = YELLOW;
    }
    else
    {
        color2 = RED;
    }

    glm::mat4 cubeTranslate2 = glm::translate(glm::mat4(1.0), glm::vec3(10.f, 0.f, -150.f));

    glm::mat4 modelMatCube2 = tempMove1 * cubeTranslate2 * matRyCube1 * matSCube1;

    glm::mat4 viewMatrixCube2 = glm::translate(glm::mat4(1.0), glm::vec3(0, -5, cubeZ));
    glm::mat4 modelMatCube2Inv = glm::transpose(glm::inverse(modelMatCube1));
    glm::mat4 perspMatCube2 = glm::perspective(glm::radians(90.0f), GLfloat(gWidth / gHeight), 0.1f, 200.f);

    glUniformMatrix4fv(glGetUniformLocation(gProgram[cube2.program_id], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMatCube2));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[cube2.program_id], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatCube2Inv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[cube2.program_id], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMatCube2));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[cube2.program_id], "viewingMat"), 1, GL_FALSE, glm::value_ptr(viewMatrixCube2));
    glUniform3fv(glGetUniformLocation(gProgram[cube2.program_id], "cubeColor"), 1, glm::value_ptr(color2));
    if (draw_flag != 2)
        cube2.draw();

    glUseProgram(gProgram[cube3.program_id]);

    glm::vec3 color3;
    if (yellow_index == 2)
    {
        color3 = YELLOW;
    }
    else
    {
        color3 = RED;
    }

    glm::mat4 cubeTranslate3 = glm::translate(glm::mat4(1.0), glm::vec3(-10.f, 0.f, -150.f));

    glm::mat4 modelMatCube3 = tempMove1 * cubeTranslate3 * matRyCube1 * matSCube1;

    glm::mat4 viewMatrixCube3 = glm::translate(glm::mat4(1.0), glm::vec3(0, -5, cubeZ));

    glm::mat4 modelMatCube3Inv = glm::transpose(glm::inverse(modelMatCube3));
    glm::mat4 perspMatCube3 = glm::perspective(glm::radians(90.0f), GLfloat(gWidth / gHeight), 0.1f, 200.f);

    glUniformMatrix4fv(glGetUniformLocation(gProgram[cube3.program_id], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMatCube3));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[cube3.program_id], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatCube3Inv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[cube3.program_id], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMatCube3));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[cube3.program_id], "viewingMat"), 1, GL_FALSE, glm::value_ptr(viewMatrixCube3));
    glUniform3fv(glGetUniformLocation(gProgram[cube3.program_id], "cubeColor"), 1, glm::value_ptr(color3));
    if (draw_flag != 3)
        cube3.draw();

    // BUNNY PART
    glUseProgram(gProgram[bunny.program_id]);

    glm::mat4 initTranslationY = glm::translate(glm::mat4(1.0), glm::vec3(0.f, -6.f, -10.f));
    glm::mat4 matRy = glm::rotate<float>(glm::mat4(1.0), (-90. / 180.) * M_PI + glm::radians(bunnyRotationAngle), glm::vec3(0.0, 1.0, 0.0));

    glm::mat4 matS = glm::scale(glm::mat4(1.0), glm::vec3(BUNNY_SCALE, BUNNY_SCALE, BUNNY_SCALE));

    glm::mat4 deadRotation = glm::rotate<float>(glm::mat4(1.0), glm::radians(-90.f), glm::vec3(0.0, 0.0, 1.0f));

    glm::mat4 modelMat = initTranslationY;
    if (deadFlag)
    {
        modelMat *= deadRotation;
    }

    modelMat *= matRy * matS;

    if (scoreFlag && !deadFlag)
    {
        bunnyRotationAngle += 6.f;
        if (bunnyRotationAngle >= rotationThreshold)
        {
            bunnyRotationAngle = 0.f;
            scoreFlag = false;
        }
    }
    // at this point bunny is at bottom of the window

    // move left and right
    glm::mat4 moveRightLeft;
    if (direction == NORMAL)
    {
        if (translateX < -1e-2)
        {
            translateX += 1;
            moveRightLeft = glm::translate(glm::mat4(1.0), glm::vec3(translateX, 0.f, 0.f));
            modelMat = moveRightLeft * modelMat;
        }
        else if (translateX > 1e-2)
        {
            translateX -= 1;
            moveRightLeft = glm::translate(glm::mat4(1.0), glm::vec3(translateX, 0.f, 0.f));
            modelMat = moveRightLeft * modelMat;
        }
    }
    else if (direction == RIGHT)
    {
        translateX += 1;

        if (translateX >= MAX_MOVE)
        {
            translateX = MAX_MOVE;

            moveRightLeft = glm::translate(glm::mat4(1.0), glm::vec3(translateX, 0.f, 0.f));
        }
        else
        {
            moveRightLeft = glm::translate(glm::mat4(1.0), glm::vec3(translateX, 0.f, 0.f));
        }
        modelMat = moveRightLeft * modelMat;
    }
    else if (direction == LEFT)
    {
        translateX -= 1;
        if (translateX <= -MAX_MOVE)
        {
            translateX = -MAX_MOVE;
            moveRightLeft = glm::translate(glm::mat4(1.0), glm::vec3(translateX, 0.f, 0.f));
        }
        else
        {
            moveRightLeft = glm::translate(glm::mat4(1.0), glm::vec3(translateX, 0.f, 0.f));
        }
        modelMat = moveRightLeft * modelMat;
    }

    // jump up and down
    glm::mat4 jump = glm::translate(glm::mat4(1.0), glm::vec3(0.f, 0.f, 0.f));
    if (translateY < 2 && jumpFlag && !deadFlag)
    {
        translateY += jump_multiplier;
        if (translateY >= 2)
        {
            jumpFlag = 0;
        }
        jump = glm::translate(glm::mat4(1.0), glm::vec3(0.f, translateY, 0.f));
    }
    else if (translateY > -2 && !jumpFlag && !deadFlag)
    {
        translateY -= jump_multiplier;
        if (translateY <= -2)
        {
            jumpFlag = 1;
        }
        jump = glm::translate(glm::mat4(1.0), glm::vec3(0.f, translateY, 0.f));
    }

    if (!deadFlag)
    {
        jump_multiplier += 0.0002f;
        modelMat = jump * modelMat;
    }

    modelMat = glm::translate(glm::mat4(1.0), glm::vec3(0, 0, quadZ)) * modelMat;
    glm::mat4 viewMatrix = glm::translate(glm::mat4(1.0), glm::vec3(0, 0, -quadZ));

    glm::mat4 modelMatInv = glm::transpose(glm::inverse(modelMat));
    glm::mat4 perspMat = glm::perspective(glm::radians(90.0f), GLfloat(gWidth / gHeight), 0.1f, 200.f);

    glUniformMatrix4fv(glGetUniformLocation(gProgram[bunny.program_id], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[bunny.program_id], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatInv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[bunny.program_id], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[bunny.program_id], "viewingMat"), 1, GL_FALSE, glm::value_ptr(viewMatrix));
    bunny.draw();

    // it is hardcoding, but it works (change!)
    if (cubeZ > 76)
    {
        // pick yellow cube
        yellow_index = generateAndPrintRandomFloats();
        cubeZ = 0;
        cubeY = 0;
        cubeX = 0;
        draw_flag = -1;
    }

    if (!deadFlag)
    {
        speed_up += 0.0005f;
        quadZ -= speed_up;
        cubeZ += speed_up;
        score += 1 + jump_multiplier;
    }

    float scaleX = static_cast<float>(gWidth) / 1024;
    float scaleY = static_cast<float>(gHeight) / 800;
    float scale = 1 * std::min(scaleX, scaleY);

    //  TEXT PART
    std::string scoreString = "Score: " + std::to_string(score);

    float textHeight = 100 * scale; // Assuming the text height is around 48 pixels
    renderText(scoreString, 10, gHeight - textHeight - 10, scale, glm::vec3(0, 1, 1));

    gluErrorString(glGetError());
}

void reshape(GLFWwindow *window, int w, int h)
{
    w = w < 1 ? 1 : w;
    h = h < 1 ? 1 : h;

    gWidth = w;
    gHeight = h;

    glViewport(0, 0, w, h);
    glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(gWidth), 0.0f, static_cast<GLfloat>(gHeight));
    glUseProgram(gProgram[0]);
    gluErrorString(glGetError());
    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "projection"), 1, GL_FALSE, glm::value_ptr(projection));
}

void keyboard(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE || key == GLFW_KEY_Q && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }
    else if (!deadFlag && key == GLFW_KEY_A && action == GLFW_PRESS)
    {
        switch (direction)
        {
        case NORMAL:
            direction = LEFT;
            break;
        case RIGHT:
            direction = NORMAL;
            break;
        default:
            direction = LEFT;
            break;
        }
        // cout << "left" << endl;
        // Update the translation matrix to move the object left
        // translateX -= 5;
    }
    else if (!deadFlag && key == GLFW_KEY_D && action == GLFW_PRESS)
    {
        switch (direction)
        {
        case NORMAL:
            direction = RIGHT;
            break;
        case LEFT:
            direction = NORMAL;
            break;
        default:
            direction = RIGHT;
            break;
        }
        // cout << "right" << endl;
        // Update the translation matrix to move the object right
        // translateX += 5;
    }
    else if (key == GLFW_KEY_R && action == GLFW_PRESS)
    {
        direction = NORMAL;
        translateX = 0;
        translateY = 0;
        jumpFlag = 0;
        score = 0;
        speed_up = 0.2f;
        jump_multiplier = 0.1f;
        quadZ = 0;
        cubeZ = 0;
        cubeY = 0;
        cubeX = 0;
    }
    // else if (key == GLFW_KEY_P && action == GLFW_PRESS)
    // {
    //     cout << "cubeX = " << cubeX << endl;
    //     cout << "cubeY = " << cubeY << endl;
    //     cout << "cubeZ = " << cubeZ << endl;

    //     cout << "quadX = " << quadX << endl;
    //     cout << "quadY = " << quadY << endl;
    //     cout << "quadZ = " << quadZ << endl;
    // }
    // else if (key == GLFW_KEY_O && action == GLFW_PRESS)
    // {
    //     cubeY += 0.5;
    //     cout << "cubeY = " << cubeY << endl;
    // }
    // else if (key == GLFW_KEY_9 && action == GLFW_PRESS)
    // {
    //     cubeY -= 0.5;
    //     cout << "cubeY = " << cubeY << endl;
    // }
    // else if (key == GLFW_KEY_I && action == GLFW_PRESS)
    // {
    //     cubeZ += 0.5;
    //     cout << "cubeZ = " << cubeZ << endl;
    // }
    // else if (key == GLFW_KEY_8 && action == GLFW_PRESS)
    // {
    //     cubeZ -= 0.5;
    //     cout << "cubeZ = " << cubeZ << endl;
    // }
    // else if (key == GLFW_KEY_X && action == GLFW_PRESS)
    // {
    //     quadX += 0.1;
    //     cout << "quadX = " << quadX << endl;
    // }
    // else if (key == GLFW_KEY_3 && action == GLFW_PRESS)
    // {
    //     quadX -= 0.1;
    //     cout << "quadX = " << quadX << endl;
    // }
    // else if (key == GLFW_KEY_Y && action == GLFW_PRESS)
    // {
    //     quadY += 0.1;
    //     cout << "quadY = " << quadY << endl;
    // }
    // else if (key == GLFW_KEY_4 && action == GLFW_PRESS)
    // {
    //     quadY -= 0.1;
    //     cout << "quadY = " << quadY << endl;
    // }
    // else if (key == GLFW_KEY_Z && action == GLFW_PRESS)
    // {
    //     quadZ += 1;
    //     cout << "quadZ = " << quadZ << endl;
    // }
    // else if (key == GLFW_KEY_5 && action == GLFW_PRESS)
    // {
    //     quadZ -= 1;
    //     cout << "quadZ = " << quadZ << endl;
    // }
}

void mainLoop(GLFWwindow *window)
{
    while (!glfwWindowShouldClose(window))
    {
        display();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

int main(int argc, char **argv) // Create Main Function For Bringing It All Together
{
    GLFWwindow *window;
    if (!glfwInit())
    {
        exit(-1);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    // glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    window = glfwCreateWindow(gWidth, gHeight, "Simple Example", NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        exit(-1);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize GLEW to setup the OpenGL Function pointers
    if (GLEW_OK != glewInit())
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }

    char rendererInfo[512] = {0};
    strcpy(rendererInfo, (const char *)glGetString(GL_RENDERER));
    strcat(rendererInfo, " - ");
    strcat(rendererInfo, (const char *)glGetString(GL_VERSION));
    glfwSetWindowTitle(window, rendererInfo);

    init();

    glfwSetKeyCallback(window, keyboard);
    glfwSetWindowSizeCallback(window, reshape);

    reshape(window, gWidth, gHeight); // need to call this once ourselves
    mainLoop(window);                 // this does not return unless the window is closed

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
