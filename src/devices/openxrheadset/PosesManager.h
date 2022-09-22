/*
 * Copyright (C) 2022 Istituto Italiano di Tecnologia (IIT)
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms of the
 * BSD-2-Clause license. See the accompanying LICENSE file for details.
 */

#ifndef YARP_DEV_POSESMANAGER_H
#define YARP_DEV_POSESMANAGER_H

#include <OpenXrInterface.h>
#include <PosePublisher.h>
#include <yarp/dev/IFrameTransform.h>
#include <string>

class PosesManager
{
    std::shared_ptr<PosePublisherSettings> m_settings{nullptr};

    std::vector<OpenXrInterface::NamedPoseVelocity> m_posesInputList;

    std::unordered_map<std::string, PosePublisher> m_poses;

public:

    struct Label
    {
        std::string original;
        std::string modified;
    };

    void initialize(const std::vector<Label> &labels,
                    const PosePublisherSettings &settings);

    std::vector<OpenXrInterface::NamedPoseVelocity>& inputs();

    void publishFrames();

};


#endif // YARP_DEV_POSESMANAGER_H
