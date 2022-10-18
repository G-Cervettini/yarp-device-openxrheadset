///////////////////////////////////////////////////////////////////////////////
// Sphere.cpp
// ==========
// Sphere for OpenGL with (radius, pan angle, tilt angle, degrees per triangle)
//
//  AUTHOR: Guglielmo Cervettini (guglielmo.cervettini@gmail.com)
// CREATED: 2022-09-29
// UPDATED: 2022-09-29
///////////////////////////////////////////////////////////////////////////////


#ifdef _WIN32
#include <windows.h>    // include windows.h to avoid thousands of compile errors even though this class is not depending on Windows
#endif

#include <GL/gl.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Sphere.h"


// constants //////////////////////////////////////////////////////////////////
const int MIN_PAN_ANGLE = 1;
const int MAX_PAN_ANGLE = 360;
const int MIN_TILT_ANGLE = 1;
const int MAX_TILT_ANGLE = 180;



///////////////////////////////////////////////////////////////////////////////
// ctor
///////////////////////////////////////////////////////////////////////////////
Sphere::Sphere(float radius, int pan, int tilt, unsigned int dpt, bool flip, bool smooth) : m_interleavedStride(32)
{
    set(radius, pan, tilt, dpt, flip, smooth);
}



///////////////////////////////////////////////////////////////////////////////
// setters
///////////////////////////////////////////////////////////////////////////////
void Sphere::set(float radius, int pan, int tilt, unsigned int dpt, bool flip, bool smooth)
{
    this->m_radius = radius;

    this->m_flip = flip;

    this->m_panAngle = pan;
    if(flip)
    {
        if (pan < MIN_TILT_ANGLE)
            this->m_panAngle = MIN_TILT_ANGLE;
        if (pan > MAX_TILT_ANGLE)
            this->m_panAngle = MAX_TILT_ANGLE;
    }else
    {
        if (pan < MIN_PAN_ANGLE)
            this->m_panAngle = MIN_PAN_ANGLE;
        if (pan > MAX_PAN_ANGLE)
            this->m_panAngle = MAX_PAN_ANGLE;
    }

    this->m_tiltAngle = tilt;
    if(flip)
    {
        if (tilt < MIN_PAN_ANGLE)
            this->m_tiltAngle = MIN_PAN_ANGLE;
        if (tilt > MAX_PAN_ANGLE)
            this->m_tiltAngle = MAX_PAN_ANGLE;
    }else
    {
        if (tilt < MIN_TILT_ANGLE)
            this->m_tiltAngle = MIN_TILT_ANGLE;
        if (tilt > MAX_TILT_ANGLE)
            this->m_tiltAngle = MAX_TILT_ANGLE;
    }

    this->m_dpt = dpt;

    this->m_smooth = smooth;
    
    if(flip)
    {
        this->m_sectorCount = (this->m_tiltAngle) / dpt;
        this->m_stackCount = (this->m_panAngle) / dpt;
    }
    else
    {
        this->m_sectorCount = (this->m_panAngle) / dpt;
        this->m_stackCount = (this->m_tiltAngle) / dpt;
    }

    if(smooth)
        buildVerticesSmooth(this->m_flip);
    else
        buildVerticesFlat(this->m_flip);
}

void Sphere::setRadius(float radius)
{
    if(radius != this->m_radius)
        set(radius, m_panAngle, m_tiltAngle, m_dpt, m_flip, m_smooth);
}

void Sphere::setPanAngle(int pan)
{
    if(pan != this->m_panAngle)
        set(m_radius, pan, m_tiltAngle, m_dpt, m_flip, m_smooth);
}

void Sphere::setTiltAngle(int tilt)
{
    if(tilt != this->m_tiltAngle)
        set(m_radius, m_panAngle, tilt, m_dpt, m_flip, m_smooth);
}

void Sphere::setDPT(unsigned int dpt)
{
    if (dpt != this->m_dpt)
        set(m_radius, m_panAngle, m_tiltAngle, dpt, m_flip, m_smooth);
}

void Sphere::setFlip(bool flip)
{
    if (flip != this->m_flip)
        set(m_radius, m_panAngle, m_tiltAngle, m_dpt, flip, m_smooth);
}

void Sphere::setSmooth(bool smooth)
{
    if(this->m_smooth == smooth)
        return;

    this->m_smooth = smooth;

    if(smooth)
        buildVerticesSmooth(this->m_flip);
    else
        buildVerticesFlat(this->m_flip);
}



///////////////////////////////////////////////////////////////////////////////
// print itself
///////////////////////////////////////////////////////////////////////////////
void Sphere::printSelf() const
{
    std::cout << "=========== Sphere ===========\n"
              << "              Radius: " << m_radius << "\n"
              << "           Pan Angle: " << m_panAngle << " deg" << "\n"
              << "          Tilt Angle: " << m_tiltAngle << " deg" << "\n"
              << "Degrees Per Triangle: " << m_dpt << "\n"
              << "      Smooth Shading: " << (m_smooth ? "true" : "false") << "\n"
              << "        Sector Count: " << m_sectorCount << "\n"
              << "         Stack Count: " << m_stackCount << "\n"
              << "      Triangle Count: " << getTriangleCount() << "\n"
              << "         Index Count: " << getIndexCount() << "\n"
              << "        Vertex Count: " << getVertexCount() << "\n"
              << "        Normal Count: " << getNormalCount() << "\n"
              << "      TexCoord Count: " << getTexCoordCount() << std::endl;
}



///////////////////////////////////////////////////////////////////////////////
// draw a sphere in VertexArray mode
// OpenGL RC must be set before calling it
///////////////////////////////////////////////////////////////////////////////
void Sphere::draw() const
{
    // interleaved array
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glVertexPointer(3, GL_FLOAT, m_interleavedStride, &m_interleavedVertices[0]);
    glNormalPointer(GL_FLOAT, m_interleavedStride, &m_interleavedVertices[3]);
    glTexCoordPointer(2, GL_FLOAT, m_interleavedStride, &m_interleavedVertices[6]);

    glDrawElements(GL_TRIANGLES, (unsigned int)m_indices.size(), GL_UNSIGNED_INT, m_indices.data());

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
}



///////////////////////////////////////////////////////////////////////////////
// draw lines only
// the caller must set the line width before call this
///////////////////////////////////////////////////////////////////////////////
void Sphere::drawLines(const float lineColor[4]) const
{
    // set line colour
    glColor4fv(lineColor);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   lineColor);

    // draw lines with VA
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, m_vertices.data());

    glDrawElements(GL_LINES, (unsigned int)m_lineIndices.size(), GL_UNSIGNED_INT, m_lineIndices.data());

    glDisableClientState(GL_VERTEX_ARRAY);
    glEnable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);
}



///////////////////////////////////////////////////////////////////////////////
// draw a sphere surfaces and lines on top of it
// the caller must set the line width before call this
///////////////////////////////////////////////////////////////////////////////
void Sphere::drawWithLines(const float lineColor[4]) const
{
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f); // move polygon backward
    this->draw();
    glDisable(GL_POLYGON_OFFSET_FILL);

    // draw lines with VA
    drawLines(lineColor);
}



///////////////////////////////////////////////////////////////////////////////
// dealloc vectors
///////////////////////////////////////////////////////////////////////////////
void Sphere::clearArrays()
{
    std::vector<float>().swap(m_vertices);
    std::vector<float>().swap(m_normals);
    std::vector<float>().swap(m_texCoords);
    std::vector<unsigned int>().swap(m_indices);
    std::vector<unsigned int>().swap(m_lineIndices);
}



///////////////////////////////////////////////////////////////////////////////
// build vertices of sphere with smooth shading using parametric equation
// x = r * cos(u) * cos(v)
// y = r * cos(u) * sin(v)
// z = r * sin(u)
// where u: stack(latitude) angle (- tiltAngle/2 <= u <= tiltAngle/2)
//       v: sector(longitude) angle (0 <= v <= panAngle)
///////////////////////////////////////////////////////////////////////////////
void Sphere::buildVerticesSmooth(bool flip)
{
    const float PI = acos(-1);
    float PAN = 0.0f;
    float TILT = 0.0f;
    if(flip)
    {
        PAN = (float)m_tiltAngle * (PI / 180);
        TILT = (float)m_panAngle * (PI / 180);
    }else
    {
        PAN = (float)m_panAngle * (PI / 180);
        TILT = (float)m_tiltAngle * (PI / 180);
    }

    // clear memory of prev arrays
    clearArrays();

    float x, y, z, xy;                                   // vertex position
    float nx, ny, nz, lengthInv = 1.0f / m_radius;       // normal
    float s, t;                                          // texCoord
                                                         
    float sectorStep = PAN / m_sectorCount;
    float stackStep = TILT / m_stackCount;
    float sectorAngle, stackAngle;                       
                                                         
    for(int i = 0; i <= m_stackCount; ++i)
    {                                                    
        stackAngle = TILT / 2 - i * stackStep;           // starting from tiltAngle/2 to -tiltAngle/2
        xy = m_radius * cosf(stackAngle);                // r * cos(u)
        z = m_radius * sinf(stackAngle);                 // r * sin(u)

        // add (sectorCount+1) vertices per stack
        // the first and last vertices may have same position and normal, but different tex coords (if panAngle=360)
        for(int j = 0; j <= m_sectorCount; ++j)
        {
            sectorAngle = PAN / 2 - j * sectorStep;      // starting from panAngle/2 to -panAngle/2

            // vertex position
            x = xy * cosf(sectorAngle);                  // r * cos(u) * cos(v)
            y = xy * sinf(sectorAngle);                  // r * cos(u) * sin(v)

            // normalized vertex normal
            nx = - x * lengthInv;
            ny = - y * lengthInv;
            nz = - z * lengthInv;

            // vertex tex coord between [0, 1]
            s = (float)j / m_sectorCount;
            t = (float)i / m_stackCount;
            
            if (flip)
            {
                addVertex(-z, -y, -x);
                addNormal(-nz, -ny, -nx);
                addTexCoord(t, s);
            }else
            {
                addVertex(y, -z, -x);
                addNormal(ny, -nz, -nx);
                addTexCoord(1 - s, t);
            }
            
        }
    }

    // indices
    //  k1--k1+1
    //  |  / |
    //  | /  |
    //  k2--k2+1
    unsigned int k1, k2;
    for(int i = 0; i < m_stackCount; ++i)
    {
        k1 = i * (m_sectorCount + 1);                    // beginning of current stack
        k2 = k1 + m_sectorCount + 1;                     // beginning of next stack

        for(int j = 0; j < m_sectorCount; ++j, ++k1, ++k2)
        {
            // 2 triangles per sector
            addIndices(k1, k2, k1+1);                    // k1---k2---k1+1
            addIndices(k1+1, k2, k2+1);                  // k1+1---k2---k2+1
            
            // vertical lines
            m_lineIndices.push_back(k1);
            m_lineIndices.push_back(k2);

            if (j == (m_sectorCount - 1))
            {
                // adding last vertical line
                m_lineIndices.push_back(k1 + 1);
                m_lineIndices.push_back(k2 + 1);
            }

            // horizontal lines
            m_lineIndices.push_back(k1);
            m_lineIndices.push_back(k1 + 1);
            
            if (i == (m_stackCount - 1))
            {
                // adding last horizontal line
                m_lineIndices.push_back(k2);
                m_lineIndices.push_back(k2 +1);
            }
        }
    }

    // generate interleaved vertex array as well
    buildInterleavedVertices();
}



///////////////////////////////////////////////////////////////////////////////
// generate vertices with flat shading
// each triangle is independent (no shared vertices)
///////////////////////////////////////////////////////////////////////////////
void Sphere::buildVerticesFlat(bool flip)
{
    const float PI = acos(-1);
    float PAN = 0.0f;
    float TILT = 0.0f;
    if(flip)
    {
        PAN = (float)m_tiltAngle * (PI / 180);
        TILT = (float)m_panAngle * (PI / 180);
    }else
    {
        PAN = (float)m_panAngle * (PI / 180);
        TILT = (float)m_tiltAngle * (PI / 180);
    }

    // tmp vertex definition (x,y,z,s,t)
    struct Vertex
    {
        float x, y, z, s, t;
    };
    std::vector<Vertex> tmpVertices;

    float sectorStep = PAN / m_sectorCount;
    float stackStep = TILT / m_stackCount;
    float sectorAngle, stackAngle;

    // compute all vertices first, each vertex contains (x,y,z,s,t) except normal
    for (int i = 0; i <= m_stackCount; ++i)
    {
        stackAngle = TILT / 2 - i * stackStep;           // starting from tiltAngle/2 to -tiltAngle/2
        float xy = m_radius * cosf(stackAngle);          // r * cos(u)
        float z = m_radius * sinf(stackAngle);           // r * sin(u)

        // add (sectorCount+1) vertices per stack
        // the first and last vertices may have same position and normal, but different tex coords (if panAngle=360)
        for (int j = 0; j <= m_sectorCount; ++j)
        {
            sectorAngle = PAN / 2 - j * sectorStep;      // starting from 0 to panAngle

            Vertex vertex;
            if (flip)
            {
                vertex.x = -z;                              // x = - r * sin(u)
                vertex.y = -xy * sinf(sectorAngle);         // y = - r * cos(u) * sin(v)
                vertex.z = -xy * cosf(sectorAngle);         // z = - r * cos(u) * cos(v)
                vertex.s = (float)i / m_stackCount;         // t
                vertex.t = (float)j / m_sectorCount;        // s
            }else
            {
                vertex.x = +xy * sinf(sectorAngle);         // x = + r * cos(u) * sin(v)
                vertex.y = -z;                              // y = - r * sin(u)
                vertex.z = -xy * cosf(sectorAngle);         // z = - r * cos(u) * cos(v)
                vertex.s = 1 - ((float)j / m_sectorCount);  // s
                vertex.t = (float)i / m_stackCount;         // t
            }
            tmpVertices.push_back(vertex);
        }
    }

    // clear memory of prev arrays
    clearArrays();

    Vertex v1, v2, v3, v4;                               // 4 vertex positions and tex coords
    std::vector<float> n;                                // 1 face normal
                                                         
    int i, j, k, vi1, vi2;                               
    int index = 0;                                       // index for vertex
    for (i = 0; i < m_stackCount; ++i)
    {                                                    
        vi1 = i * (m_sectorCount + 1);                   // index of tmpVertices
        vi2 = (i + 1) * (m_sectorCount + 1);

        for (j = 0; j < m_sectorCount; ++j, ++vi1, ++vi2)
        {
            // get 4 vertices per sector
            //  v1--v3
            //  |    |
            //  v2--v4
            v1 = tmpVertices[vi1];
            v2 = tmpVertices[vi2];
            v3 = tmpVertices[vi1 + 1];
            v4 = tmpVertices[vi2 + 1];

            // store 2 triangles per sector (quad) =========================

            // put quad vertices: v1-v2-v3-v4
            addVertex(v1.x, v1.y, v1.z);
            addVertex(v2.x, v2.y, v2.z);
            addVertex(v3.x, v3.y, v3.z);
            addVertex(v4.x, v4.y, v4.z);

            // put tex coords of quad
            addTexCoord(v1.s, v1.t);
            addTexCoord(v2.s, v2.t);
            addTexCoord(v3.s, v3.t);
            addTexCoord(v4.s, v4.t);

            // put normal
            n = computeFaceNormal(v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v3.x, v3.y, v3.z);
            for (k = 0; k < 4; ++k)  // same normals for 4 vertices
            {
                addNormal(-n[0], -n[1], -n[2]);
            }

            // put indices of quad (2 triangles)
            addIndices(index, index + 1, index + 2);
            addIndices(index + 2, index + 1, index + 3);

            // indices for vertical lines
            m_lineIndices.push_back(index);
            m_lineIndices.push_back(index + 1);

            if (j == (m_sectorCount - 1))
            {
                // adding last vertical line
                m_lineIndices.push_back(index + 2);
                m_lineIndices.push_back(index + 3);
            }

            // indices for horizontal lines
            m_lineIndices.push_back(index);
            m_lineIndices.push_back(index + 2);

            if (i == (m_stackCount - 1))
            {
                // adding last horizontal line
                m_lineIndices.push_back(index + 1);
                m_lineIndices.push_back(index + 3);
            }

            index += 4;     // for next
        }
    }

    // generate interleaved vertex array as well
    buildInterleavedVertices();
}



///////////////////////////////////////////////////////////////////////////////
// generate interleaved vertices: V/N/T
// stride must be 32 bytes
///////////////////////////////////////////////////////////////////////////////
void Sphere::buildInterleavedVertices()
{
    std::vector<float>().swap(m_interleavedVertices);

    std::size_t i, j;
    std::size_t count = m_vertices.size();
    for(i = 0, j = 0; i < count; i += 3, j += 2)
    {
        m_interleavedVertices.push_back(m_vertices[i]);
        m_interleavedVertices.push_back(m_vertices[i+1]);
        m_interleavedVertices.push_back(m_vertices[i+2]);

        m_interleavedVertices.push_back(m_texCoords[j]);
        m_interleavedVertices.push_back(m_texCoords[j+1]);

        m_interleavedVertices.push_back(m_normals[i]);
        m_interleavedVertices.push_back(m_normals[i+1]);
        m_interleavedVertices.push_back(m_normals[i+2]);
    }
}



///////////////////////////////////////////////////////////////////////////////
// add single vertex to array
///////////////////////////////////////////////////////////////////////////////
void Sphere::addVertex(float x, float y, float z)
{
    m_vertices.push_back(x);
    m_vertices.push_back(y);
    m_vertices.push_back(z);
}



///////////////////////////////////////////////////////////////////////////////
// add single normal to array
///////////////////////////////////////////////////////////////////////////////
void Sphere::addNormal(float nx, float ny, float nz)
{
    m_normals.push_back(nx);
    m_normals.push_back(ny);
    m_normals.push_back(nz);
}



///////////////////////////////////////////////////////////////////////////////
// add single texture coord to array
///////////////////////////////////////////////////////////////////////////////
void Sphere::addTexCoord(float s, float t)
{
    m_texCoords.push_back(s);
    m_texCoords.push_back(t);
}



///////////////////////////////////////////////////////////////////////////////
// add 3 indices to array
///////////////////////////////////////////////////////////////////////////////
void Sphere::addIndices(unsigned int i1, unsigned int i2, unsigned int i3)
{
    m_indices.push_back(i1);
    m_indices.push_back(i2);
    m_indices.push_back(i3);
}



///////////////////////////////////////////////////////////////////////////////
// return face normal of a triangle v1-v2-v3
// if a triangle has no surface (normal length = 0), then return a zero vector
///////////////////////////////////////////////////////////////////////////////
std::vector<float> Sphere::computeFaceNormal(float x1, float y1, float z1,  // v1
                                             float x2, float y2, float z2,  // v2
                                             float x3, float y3, float z3)  // v3
{
    const float EPSILON = 0.000001f;

    std::vector<float> normal(3, 0.0f);     // default return value (0,0,0)
    float nx, ny, nz;

    // find 2 edge vectors: v1-v2, v1-v3
    float ex1 = x2 - x1;
    float ey1 = y2 - y1;
    float ez1 = z2 - z1;
    float ex2 = x3 - x1;
    float ey2 = y3 - y1;
    float ez2 = z3 - z1;

    // cross product: e1 x e2
    nx = ey1 * ez2 - ez1 * ey2;
    ny = ez1 * ex2 - ex1 * ez2;
    nz = ex1 * ey2 - ey1 * ex2;

    // normalize only if the length is > 0
    float length = sqrtf(nx * nx + ny * ny + nz * nz);
    if(length > EPSILON)
    {
        // normalize
        float lengthInv = 1.0f / length;
        normal[0] = nx * lengthInv;
        normal[1] = ny * lengthInv;
        normal[2] = nz * lengthInv;
    }

    return normal;
}
