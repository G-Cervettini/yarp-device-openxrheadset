///////////////////////////////////////////////////////////////////////////////
// Sphere.h
// ========
// Sphere for OpenGL with (radius, pan angle, tilt angle, degrees per triangle)
//
//  AUTHOR: Guglielmo Cervettini (guglielmo.cervettini@gmail.com)
// CREATED: 2022-09-29
// UPDATED: 2022-09-29
///////////////////////////////////////////////////////////////////////////////

#ifndef GEOMETRY_SPHERE_H
#define GEOMETRY_SPHERE_H

#include <vector>

class Sphere
{
public:
    // ctor/dtor
    Sphere(float radius=1.0f, int panAngle=180, int tiltAngle=90, unsigned int dpt=10, bool flip = false, bool smooth=true);
    ~Sphere() {}

    // getters/setters
    float getRadius() const                 { return radius; }
    int getPanAngle() const                 { return panAngle; }
    int getTiltAngle() const                { return tiltAngle; }
    unsigned int getDPT() const             { return dpt; }
    int getSectorCount() const              { return sectorCount; }
    int getStackCount() const               { return stackCount; }
    void set(float radius, int panAngle, int tiltAngle, unsigned int dpt, bool flip, bool smooth=true);
    void setRadius(float radius);
    void setPanAngle(int panAngle);
    void setTiltAngle(int tiltAngle);
    void setDPT(unsigned int dpt);
    void setFlip(bool flip);
    void setSmooth(bool smooth);

    // for vertex data
    unsigned int getVertexCount() const                         { return (unsigned int)vertices.size() / 3; }
    unsigned int getNormalCount() const                         { return (unsigned int)normals.size() / 3; }
    unsigned int getTexCoordCount() const                       { return (unsigned int)texCoords.size() / 2; }
    unsigned int getIndexCount() const                          { return (unsigned int)indices.size(); }
    unsigned int getLineIndexCount() const                      { return (unsigned int)lineIndices.size(); }
    unsigned int getTriangleCount() const                       { return getIndexCount() / 3; }
    unsigned int getVertexSize() const                          { return (unsigned int)vertices.size() * sizeof(float); }
    unsigned int getNormalSize() const                          { return (unsigned int)normals.size() * sizeof(float); }
    unsigned int getTexCoordSize() const                        { return (unsigned int)texCoords.size() * sizeof(float); }
    unsigned int getIndexSize() const                           { return (unsigned int)indices.size() * sizeof(unsigned int); }
    unsigned int getLineIndexSize() const                       { return (unsigned int)lineIndices.size() * sizeof(unsigned int); }
    const std::vector<float> getVertices() const                { return vertices; }
    const float* getNormals() const                             { return normals.data(); }
    const float* getTexCoords() const                           { return texCoords.data(); }
    const std::vector<unsigned int> getIndices() const          { return indices; }
    const std::vector<unsigned int> getLineIndices() const      { return lineIndices; }

    // for interleaved vertices: V/N/T
    unsigned int getInterleavedVertexCount() const              { return getVertexCount(); }    // # of vertices
    unsigned int getInterleavedVertexSize() const               { return (unsigned int)interleavedVertices.size() * sizeof(float); }    // # of bytes
    int getInterleavedStride() const                            { return interleavedStride; }   // should be 32 bytes
    const std::vector<float> getInterleavedVertices() const     { return interleavedVertices; }

    // draw in VertexArray mode
    void draw() const;                                  // draw surface
    void drawLines(const float lineColor[4]) const;     // draw lines only
    void drawWithLines(const float lineColor[4]) const; // draw surface and lines

    // debug
    void printSelf() const;

private:
    // member functions
    void buildVerticesSmooth(bool flip);
    void buildVerticesFlat(bool flip);
    void buildInterleavedVertices();
    void clearArrays();
    void addVertex(float x, float y, float z);
    void addNormal(float x, float y, float z);
    void addTexCoord(float s, float t);
    void addIndices(unsigned int i1, unsigned int i2, unsigned int i3);
    std::vector<float> computeFaceNormal(float x1, float y1, float z1,
                                         float x2, float y2, float z2,
                                         float x3, float y3, float z3);

    // memeber vars
    float radius;
    int panAngle;
    int tiltAngle;
    unsigned int dpt;                       // Degrees Per Triangle (in both longitude and latitude)
    bool smooth;
    bool flip { false };

    int sectorCount;                        // longitude, # of slices
    int stackCount;                         // latitude, # of stacks
    std::vector<float> vertices;
    std::vector<float> normals;
    std::vector<float> texCoords;
    std::vector<unsigned int> indices;
    std::vector<unsigned int> lineIndices;

    // interleaved
    std::vector<float> interleavedVertices;
    int interleavedStride;                  // # of bytes to hop to the next vertex (should be 32 bytes)

};

#endif
