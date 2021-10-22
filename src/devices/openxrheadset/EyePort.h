/*
 * Copyright (C) 2021 Istituto Italiano di Tecnologia (IIT)
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms of the
 * BSD-2-Clause license. See the accompanying LICENSE file for details.
 */

#ifndef YARP_DEV_EYEPORT_H
#define YARP_DEV_EYEPORT_H

#include <yarp/dev/IFrameTransform.h>
#include <yarp/sig/Image.h>
#include <yarp/sig/Matrix.h>
#include <yarp/sig/Vector.h>
#include <yarp/os/Stamp.h>
#include <yarp/os/BufferedPort.h>
#include <OpenXrInterface.h>
#include <PortToQuadLayer.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

class EyePort
{
    yarp::sig::Matrix m_localPose;
    yarp::dev::IFrameTransform* m_tfPublisher{nullptr};
    std::string m_tfFrame;
    std::string m_rootFrame;
    float m_azimuthOffset {0.0};
    float m_elevationOffset {0.0};
    Eigen::Quaternionf m_desiredRotation;
    Eigen::Quaternionf m_rotationOffset;
    Eigen::Vector3f m_eyePosition;
    Eigen::Vector3f m_eyeRelativeImagePosition;
    PortToQuadLayer<yarp::sig::ImageOf<yarp::sig::PixelRgb>> m_layer;
    yarp::os::BufferedPort<yarp::sig::Vector> m_eyeAnglesPort;
    bool m_initialized{false};

public:

    bool open(std::shared_ptr<IOpenXrQuadLayer> quadLayer, const std::string& imagePortName, const std::string& anglesPortName,
              yarp::dev::IFrameTransform* tfPublisher, const std::string& tfFrame, const std::string& rootFrame);

    void close();

    void setEyePosition(const Eigen::Vector3f& position);

    void setEyeRotationOffset(double azimuth, double elevation);

    void setEyeRotation(double azimuth, double elevation);

    double azimuthOffset() const;

    double elevationOffset() const;

    void setEyeRelativeImagePosition(const Eigen::Vector3f& position);

    void setVisibility(const IOpenXrQuadLayer::Visibility& visibility);

    float layerWidth() const;

    float layerHeight() const;

    bool update();

    void publishEyeTransform();
};

#endif // YARP_DEV_EYEPORT_H
