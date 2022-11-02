/*
 * Copyright (C) 2022 Istituto Italiano di Tecnologia (IIT)
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms of the
 * BSD-2-Clause license. See the accompanying LICENSE file for details.
 */

#ifndef YARP_DEV_OPENGLSPHERELAYER_H
#define YARP_DEV_OPENGLSPHERELAYER_H

#include <vector>

#include <OpenXrConfig.h>
#include <OpenXrInterface.h>

#include <Renderer.h>
#include <VertexBuffer.h>
#include <VertexBufferLayout.h>
#include <IndexBuffer.h>
#include <VertexArray.h>
#include <Shader.h>
#include <Texture.h>
#include <FrameBuffer.h>
#include <Sphere.h>

class OpenGLSphereLayer : public IOpenXrSphereLayer
{
    Sphere m_sphere;                                                          // radius, pan angle (deg), tilt angle (deg), degrees per triangle, flip N-S poles to E-W poles, smooth shading (default: true)

    std::vector<float> m_positions;
    std::vector<unsigned int> m_indices;


    IOpenXrQuadLayer::Visibility m_visibility{ IOpenXrQuadLayer::Visibility::NONE };
    int32_t m_imageMaxWidth = 0;
    int32_t m_imageMaxHeight = 0;
    VertexArray m_va;
    VertexBuffer m_vb;
    VertexBufferLayout m_layout;
    IndexBuffer m_ib;
    Shader m_shader, m_shaderLine;
    Texture m_userTexture, m_internalTexture;
    FrameBuffer m_userBuffer;
    FrameBuffer m_internalBuffer;
    bool m_useAlpha{true};
    bool m_isEnabled{true};
    bool m_isReleased{false};
    GridVisibility m_isGridVisible{ GridVisibility::VISIBLE_GRID };

    glm::mat4 m_offsetTra = glm::mat4(1.0f);                                  // position of the Headset Frame WRT the Left or Right Screen Frame
    bool m_offsetIsSet{false};

    Eigen::Vector3f m_modelTraEig {0.0f, 0.0f, 0.0f};
    Eigen::Quaternionf m_modelRotEig {1.0f, 0.0f, 0.0f, 0.0f};
    glm::mat4 m_modelTra = glm::mat4(1.0f);
    glm::mat4 m_modelRot = glm::mat4(1.0f);

    glm::vec3 m_modelScale{1.0f, 1.0f, 1.0f};

    float m_fovY = glm::radians(60.0f);                                                      // Field Of View
    float m_zNear = 0.1f;
    float m_zFar = 100.0f;
    float m_aspectRatio = 1.0f;

public:

    OpenGLSphereLayer();

    ~OpenGLSphereLayer();

    OpenGLSphereLayer(const OpenGLSphereLayer&) = delete;

    OpenGLSphereLayer(OpenGLSphereLayer&&) = delete;

    OpenGLSphereLayer& operator=(const OpenGLSphereLayer&) = delete;

    OpenGLSphereLayer& operator=(OpenGLSphereLayer&&) = delete;

    bool initialize(int32_t imageMaxWidth, int32_t imageMaxHeight);
    
    void render();

    virtual void setSphereRadius(float radius) override;

    virtual void setViewAngles(int pan, int tilt) override;

    virtual void setGridResolution(unsigned int degreesPerTriangle) override;

    virtual void setGridPolesDirection(const GridPolesDirection& gridPolesDirection) override;

    virtual void setGridVisibility(const GridVisibility& gridVisibility) override;

    void setFOVs(float fovX, float fovY);

    void setDepthLimits(float zNear, float zFar);

    void setOffsetPosition(const Eigen::Vector3f& offset);

    bool offsetIsSet() const;

    Texture& getUserTexture();

    const  IOpenXrQuadLayer::Visibility& visibility() const;

    bool shouldRender() const;

    virtual void setPose(const Eigen::Vector3f& position,
                         const Eigen::Quaternionf &quaternion) override;

    virtual void setPosition(const Eigen::Vector3f& position) override;

    virtual void setQuaternion(const Eigen::Quaternionf &quaternion) override;

    void setAxisScale(float scaleX, float scaleY, float scaleZ);

    virtual void setDimensions(float widthInMeters, float heightInMeters) override;

    virtual void setVisibility(const Visibility& visibility) override;

    virtual void useAlphaChannel(bool useAlphaChannel = true) override;

    virtual bool getImage(uint32_t& glImage) override;

    virtual bool submitImage() override;

    virtual bool submitImage(int32_t xOffset, int32_t yOffset, int32_t imageWidth, int32_t imageHeight) override;

    virtual int32_t imageMaxHeight() const override;

    virtual int32_t imageMaxWidth() const override;

    virtual float layerWidth() const override;

    virtual float layerHeight() const override;

    virtual Eigen::Vector3f layerPosition() const override;

    virtual Eigen::Quaternionf layerQuaternion() const override;

    virtual void setEnabled(bool enabled) override;
};

#endif // YARP_DEV_OPENGLSPHERELAYER_H
