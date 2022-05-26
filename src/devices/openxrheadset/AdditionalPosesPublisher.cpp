/*
 * Copyright (C) 2022 Istituto Italiano di Tecnologia (IIT)
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms of the
 * BSD-2-Clause license. See the accompanying LICENSE file for details.
 */

#include <AdditionalPosesPublisher.h>
#include <OpenXrHeadsetLogComponent.h>
#include <OpenXrYarpUtilities.h>
#include <yarp/os/LogStream.h>

void AdditionalPosesPublisher::initialize(const std::vector<Label>& labels, const PosePublisherSettings &settings)
{
    m_settings = std::make_shared<PosePublisherSettings>(settings);

    for (auto& label : labels)
    {
        m_additionalPoses[label.original].setLabel(label.modified);
    }
}

std::vector<OpenXrInterface::NamedPoseVelocity> &AdditionalPosesPublisher::inputs()
{
    return m_additionalPosesInputList;
}

void AdditionalPosesPublisher::publishFrames()
{
    for (auto& inputPose : m_additionalPosesInputList)
    {
        PosePublisher& publisher = m_additionalPoses[inputPose.name];

        if (!publisher.configured())
        {
            publisher.configure(m_settings);
        }

        publisher.update(inputPose);
    }

    for (auto& poseIt : m_additionalPoses)
    {
        poseIt.second.publish();
    }
}
