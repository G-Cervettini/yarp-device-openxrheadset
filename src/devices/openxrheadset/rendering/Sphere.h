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
    float getRadius() const                 { return m_radius; }
    int getPanAngle() const                 { return m_panAngle; }
    int getTiltAngle() const                { return m_tiltAngle; }
    unsigned int getDPT() const             { return m_dpt; }
    int getSectorCount() const              { return m_sectorCount; }
    int getStackCount() const               { return m_stackCount; }
    void set(float radius, int panAngle, int tiltAngle, unsigned int dpt, bool flip, bool smooth=true);
    void setRadius(float radius);
    void setPanAngle(int panAngle);
    void setTiltAngle(int tiltAngle);
    void setDPT(unsigned int dpt);
    void setFlip(bool flip);
    void setSmooth(bool smooth);

    // for vertex data
    unsigned int getVertexCount() const                         { return (unsigned int)m_vertices.size() / 3; }
    unsigned int getNormalCount() const                         { return (unsigned int)m_normals.size() / 3; }
    unsigned int getTexCoordCount() const                       { return (unsigned int)m_texCoords.size() / 2; }
    unsigned int getIndexCount() const                          { return (unsigned int)m_indices.size(); }
    unsigned int getLineIndexCount() const                      { return (unsigned int)m_lineIndices.size(); }
    unsigned int getTriangleCount() const                       { return getIndexCount() / 3; }
    unsigned int getVertexSize() const                          { return (unsigned int)m_vertices.size() * sizeof(float); }
    unsigned int getNormalSize() const                          { return (unsigned int)m_normals.size() * sizeof(float); }
    unsigned int getTexCoordSize() const                        { return (unsigned int)m_texCoords.size() * sizeof(float); }
    unsigned int getIndexSize() const                           { return (unsigned int)m_indices.size() * sizeof(unsigned int); }
    unsigned int getLineIndexSize() const                       { return (unsigned int)m_lineIndices.size() * sizeof(unsigned int); }
    const std::vector<float> getVertices() const                { return m_vertices; }
    const float* getNormals() const                             { return m_normals.data(); }
    const float* getTexCoords() const                           { return m_texCoords.data(); }
    const std::vector<unsigned int> getIndices() const          { return m_indices; }
    const std::vector<unsigned int> getLineIndices() const      { return m_lineIndices; }

    // for interleaved vertices: V/N/T
    unsigned int getInterleavedVertexCount() const              { return getVertexCount(); }    // # of vertices
    unsigned int getInterleavedVertexSize() const               { return (unsigned int)m_interleavedVertices.size() * sizeof(float); }    // # of bytes
    int getInterleavedStride() const                            { return m_interleavedStride; }   // should be 32 bytes
    const std::vector<float> getInterleavedVertices() const     { return m_interleavedVertices; }

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
    float m_radius;
    int m_panAngle;
    int m_tiltAngle;
    unsigned int m_dpt;                                // Degrees Per Triangle (in both longitude and latitude)
    bool m_smooth;
    bool m_flip { false };

    int m_sectorCount;                        // longitude, # of slices
    int m_stackCount;                         // latitude, # of stacks
    std::vector<float> m_vertices;
    std::vector<float> m_normals;
    std::vector<float> m_texCoords;
    std::vector<unsigned int> m_indices;
    std::vector<unsigned int> m_lineIndices;

    // interleaved
    std::vector<float> m_interleavedVertices;
    int m_interleavedStride;                  // # of bytes to hop to the next vertex (should be 32 bytes)

};

#endif
