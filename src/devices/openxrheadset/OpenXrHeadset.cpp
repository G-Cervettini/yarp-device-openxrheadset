/*
 * Copyright (C) 2021 Istituto Italiano di Tecnologia (IIT)
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms of the
 * BSD-2-Clause license. See the accompanying LICENSE file for details.
 */

#include <yarp/os/LogStream.h>

#include <map>
#include <algorithm>

#include <OpenXrHeadset.h>
#include <OpenXrHeadsetLogComponent.h>
#include <OpenXrYarpUtilities.h>

typedef bool(yarp::os::Value::*valueIsType)(void) const;

yarp::dev::OpenXrHeadset::OpenXrHeadset()
    : yarp::dev::DeviceDriver(),
      yarp::os::PeriodicThread(0.011, yarp::os::ShouldUseSystemClock::Yes), // ~90 fps
      m_stamp(0,0.0)
{}

yarp::dev::OpenXrHeadset::~OpenXrHeadset()
{
    this->stop();
}

bool yarp::dev::OpenXrHeadset::open(yarp::os::Searchable &cfg)
{
    std::string name = cfg.check("name", yarp::os::Value("OpenXrHeadset")).toString();
    if (name.front() != '/')
    {
        m_prefix = '/' + name;
    }
    else
    {
        m_prefix = name;
    }

    //checking the additional guis parameter in the configuration file..
    {
        constexpr unsigned int STRING = 0;
        constexpr unsigned int BOOL   = 1;
        constexpr unsigned int INT    = 2;
        constexpr unsigned int DOUBLE = 3;
        struct paramDescription
        {
            std::string name;
            int type;
        };

        std::map<int, std::string>                err_msgs;
        std::map<int, valueIsType>                isFunctionMap;
        std::vector<paramDescription> paramParser;

        err_msgs[STRING]      = "a string";
        err_msgs[BOOL]        = "a boolean type";
        err_msgs[INT]         = "an integer type";
        err_msgs[DOUBLE]      = "a real type";
        isFunctionMap[STRING] = &yarp::os::Value::isString;
        isFunctionMap[BOOL]   = &yarp::os::Value::isBool;
        isFunctionMap[INT]    = &yarp::os::Value::isInt32;
        isFunctionMap[DOUBLE] = &yarp::os::Value::isFloat64;

        int guiCount = cfg.find("gui_elements").asInt32();
        paramParser.clear();
        if (guiCount)
        {
            paramParser.push_back({"width",  DOUBLE});
            paramParser.push_back({"height", DOUBLE});
            paramParser.push_back({"x",      DOUBLE});
            paramParser.push_back({"y",      DOUBLE});
            paramParser.push_back({"z",      DOUBLE});

            for (int i = 0; i < guiCount; ++i)
            {
                std::string       groupName  = "GUI_" + std::to_string(i);
                yarp::os::Bottle& guip       = cfg.findGroup(groupName);

                if (guip.isNull())
                {
                    yCError(OPENXRHEADSET) << "group:" << groupName << "not found in configuration file..";
                    return false;
                }

                for (auto& p : paramParser)
                {
                    if (!guip.check(p.name) || !(guip.find(p.name).*isFunctionMap[p.type])())
                    {
                        std::string err_type = err_msgs.find(p.type) == err_msgs.end() ? "[unknow type]" : err_msgs[p.type];
                        yCError(OPENXRHEADSET) << "parameter" << p.name << "not found or not" << err_type << "in" << groupName << "group in configuration file";
                        return false;
                    }
                }

                m_huds.emplace_back();

                m_huds.back().width = guip.find("width").asFloat64();
                m_huds.back().height = guip.find("height").asFloat64();
                m_huds.back().x       = guip.find("x").asFloat64();
                m_huds.back().y       = guip.find("y").asFloat64();
                m_huds.back().z       = -std::max(0.01, std::abs(guip.find("z").asFloat64())); //make sure that z is negative and that is at least 0.01 in modulus
                std::transform(groupName.begin(), groupName.end(), groupName.begin(), ::tolower);
                m_huds.back().portName = m_prefix + "/" + groupName;
            }
        }

        paramParser.clear();
        int labelCount = cfg.find("labels").asInt32();
        if (labelCount)
        {
            paramParser.push_back({"width",  DOUBLE});
            paramParser.push_back({"height", DOUBLE});
            paramParser.push_back({"x",      DOUBLE});
            paramParser.push_back({"y",      DOUBLE});
            paramParser.push_back({"z",      DOUBLE});

            for (int i = 0; i < labelCount; ++i)
            {
                std::string       groupName  = "LABEL_" + std::to_string(i);
                yarp::os::Bottle& labelGroup = cfg.findGroup(groupName);

                if (labelGroup.isNull())
                {
                    yCError(OPENXRHEADSET) << "group:" << groupName << "not found in configuration file..";
                    return false;
                }

                for (auto& p : paramParser)
                {
                    if (!labelGroup.check(p.name) || !(labelGroup.find(p.name).*isFunctionMap[p.type])())
                    {
                        std::string err_type = err_msgs.find(p.type) == err_msgs.end() ? "[unknow type]" : err_msgs[p.type];
                        yCError(OPENXRHEADSET) << "parameter" << p.name << "not found or not" << err_type << "in" << groupName << "group in configuration file";
                        return false;
                    }
                }

                m_labels.emplace_back();
                LabelLayer& label = m_labels.back();

                label.width = labelGroup.find("width").asFloat64();
                label.height = labelGroup.find("height").asFloat64();
                label.x       = labelGroup.find("x").asFloat64();
                label.y       = labelGroup.find("y").asFloat64();
                label.z       = -std::max(0.01, std::abs(labelGroup.find("z").asFloat64())); //make sure that z is negative and that is at least 0.01 in modulus

                std::transform(groupName.begin(), groupName.end(), groupName.begin(), ::tolower);
                std::string portName = m_prefix + "/" + groupName;

                if (!label.options.parseFromConfigurationFile(portName, labelGroup))
                {
                    yCError(OPENXRHEADSET) << "Failed to parse" << groupName;
                    return false;
                }
            }
        }
    }

    m_getStickAsAxis = cfg.check("stick_as_axis", yarp::os::Value(false)).asBool();
    m_leftFrame = cfg.check("tf_left_hand_frame", yarp::os::Value("openxr_left_hand")).asString();
    m_rightFrame = cfg.check("tf_right_hand_frame", yarp::os::Value("openxr_right_hand")).asString();
    m_headFrame = cfg.check("tf_head_frame", yarp::os::Value("openxr_head")).asString();
    m_eyesManager.options().leftEyeFrame = cfg.check("tf_left_eye_frame", yarp::os::Value("openxr_left_eye")).asString();
    m_eyesManager.options().rightEyeFrame = cfg.check("tf_right_eye_frame", yarp::os::Value("openxr_right_eye")).asString();
    m_rootFrame = cfg.check("tf_root_frame", yarp::os::Value("openxr_origin")).asString();
    m_rootFrameRaw = m_rootFrame + "_raw";
    m_eyesManager.options().leftAzimuthOffset = cfg.check("left_azimuth_offset", yarp::os::Value(0.0)).asFloat64();
    m_eyesManager.options().leftElevationOffset = cfg.check("left_elevation_offset", yarp::os::Value(0.0)).asFloat64();
    m_eyesManager.options().eyeZPosition = -std::max(0.01, std::abs(cfg.check("eye_z_position", yarp::os::Value(-1.0)).asFloat64())); //make sure that z is negative and that is at least 0.01 in modulus
    m_eyesManager.options().interCameraDistance = std::abs(cfg.check("inter_camera_distance", yarp::os::Value(0.07)).asFloat64()); //Distance between the cameras of the iCub robot
    m_eyesManager.options().rightAzimuthOffset = cfg.check("right_azimuth_offset", yarp::os::Value(0.0)).asFloat64();
    m_eyesManager.options().rightElevationOffset = cfg.check("right_elevation_offset", yarp::os::Value(0.0)).asFloat64();
    m_eyesManager.options().splitEyes = cfg.check("split_eye_ports", yarp::os::Value(true)).asBool();
    m_eyesManager.options().portPrefix = m_prefix;

    m_rawRootFrameTransform.resize(4,4);
    m_rawRootFrameTransform.eye();

    std::vector<AdditionalPosesPublisher::Label> labels;
    yarp::os::Bottle& labelsGroup = cfg.findGroup("POSES_LABELS");
    for (size_t i = 1; i < labelsGroup.size(); ++i) //The first element is the name of the group itself
    {
        yarp::os::Value& labelElement = labelsGroup.get(i);
        if (!labelElement.isList() || labelElement.asList()->size() != 2)
        {
            yCError(OPENXRHEADSET) << "Each entry of the POSES_LABELS group is supposed to contain only two elements. The original name and the modified one. Cause: " << labelElement.toString();
            return false;
        }
        yarp::os::Bottle* labelList = labelElement.asList();
        labels.push_back({labelList->get(0).asString(), labelList->get(1).asString()});
    }


    //opening tf client
    yarp::os::Property tfClientCfg;
    tfClientCfg.put("device", cfg.check("tfDevice", yarp::os::Value("transformClient")).asString());
    tfClientCfg.put("local",  cfg.check("tfLocal", yarp::os::Value(m_prefix + "/tf")).asString());
    tfClientCfg.put("remote", cfg.check("tfRemote", yarp::os::Value("/transformServer")).asString());

    if (!m_driver.open(tfClientCfg))
    {
        yCError(OPENXRHEADSET) << "Unable to open polydriver with the following options:" << tfClientCfg.toString();
        return false;
    }

    if (!m_driver.view(m_tfPublisher) || m_tfPublisher == nullptr)
    {
        yCError(OPENXRHEADSET) << "Unable to view IFrameTransform interface.";
        return false;
    }
    yCInfo(OPENXRHEADSET) << "TransformCLient successfully opened at port: " << cfg.find("tfLocal").asString();

    if (!m_headFramePorts.open("Head", m_prefix + "/headpose", m_tfPublisher, m_headFrame, m_rootFrameRaw))
    {
        return false;
    }

    if (!m_leftHandFramePorts.open("Left Hand", m_prefix + "/left_hand", m_tfPublisher, m_leftFrame, m_rootFrameRaw))
    {
        return false;
    }

    if (!m_rightHandFramePorts.open("Right Hand", m_prefix + "/right_hand", m_tfPublisher, m_rightFrame, m_rootFrameRaw))
    {
        return false;
    }

    PosePublisherSettings posePublisherSettings;
    posePublisherSettings.tfPublisher = m_tfPublisher;
    posePublisherSettings.rootFrame = m_rootFrameRaw;
    posePublisherSettings.period = getPeriod();
    posePublisherSettings.checks.maxDistance = cfg.check("pose_check_max_distance", yarp::os::Value(0.1)).asFloat64();
    posePublisherSettings.checks.maxAngularDistanceInRad = cfg.check("pose_check_max_angle_rad", yarp::os::Value(0.5)).asFloat64();
    posePublisherSettings.checks.lastDataExpirationTime = cfg.check("pose_check_expiration_time", yarp::os::Value(5.0)).asFloat64();
    posePublisherSettings.checks.maxConvergenceTime = cfg.check("pose_check_max_convergence_time", yarp::os::Value(3.0)).asFloat64();
    posePublisherSettings.checks.convergenceRatio = cfg.check("pose_check_convergence_ratio", yarp::os::Value(0.05)).asFloat64();

    if (posePublisherSettings.checks.convergenceRatio < 0 || posePublisherSettings.checks.convergenceRatio > 1)
    {
        posePublisherSettings.checks.convergenceRatio = std::min(1.0, std::max(0.0, posePublisherSettings.checks.convergenceRatio));
        yCWarning(OPENXRHEADSET) << "pose_check_convergence_ratio is supposed to be in the range [0, 1]. Clamping to" << posePublisherSettings.checks.convergenceRatio << ".";
    }

    m_additionalPosesPublisher.initialize(labels, posePublisherSettings);

    // Start the thread
    if (!this->start()) {
        yCError(OPENXRHEADSET) << "Thread start failed, aborting.";
        this->close();
        return false;
    }

    return true;
}

bool yarp::dev::OpenXrHeadset::close()
{
    this->askToStop();
    return true;
}

bool yarp::dev::OpenXrHeadset::threadInit()
{
    {
        std::lock_guard<std::mutex> lock(m_mutex);

        if (!m_openXrInterface.initialize())
        {
            yCError(OPENXRHEADSET) << "Failed to initialize OpenXr interface.";
            return false;
        }

        m_eyesManager.options().leftEyeQuadLayer = m_openXrInterface.addHeadFixedQuadLayer();
        m_eyesManager.options().rightEyeQuadLayer = m_openXrInterface.addHeadFixedQuadLayer();

        if (!m_eyesManager.initialize(m_tfPublisher, m_headFrame))
        {
            yCError(OPENXRHEADSET) << "Failed to initialize eyes.";
        }

        for (GuiParam& gui : m_huds)
        {
            if (!gui.layer.initialize(m_openXrInterface.addHeadFixedQuadLayer(), gui.portName)) {
                yCError(OPENXRHEADSET) << "Cannot initialize" << gui.portName << "display texture.";
                return false;
            }
            gui.layer.setVisibility(IOpenXrQuadLayer::Visibility::BOTH_EYES);
            gui.layer.setDimensions(gui.width, gui.height);
            gui.layer.setPosition({gui.x, gui.y, gui.z});
        }

        for (LabelLayer& label : m_labels)
        {
            label.options.quadLayer = m_openXrInterface.addHeadFixedQuadLayer();

            if (!label.layer.initialize(label.options)) {
                yCError(OPENXRHEADSET) << "Cannot initialize" << label.options.portName << "label.";
                return false;
            }
            label.layer.setVisibility(IOpenXrQuadLayer::Visibility::BOTH_EYES);
            label.layer.setDimensions(label.width, label.height);
            label.layer.setPosition({label.x, label.y, label.z});
        }
    }

    for (size_t i = 0; i < 10 && m_openXrInterface.isRunning(); ++i)
    {
        run(); //dry run. This is to make sure that the number of buttons is correctly retrieved by the JoypadControlServer
        yarp::os::Time::delay(this->getPeriod());
    }

    {
        std::lock_guard<std::mutex> lock(m_mutex);

        this->yarp().attachAsServer(this->m_rpcPort);
        if(!m_rpcPort.open(m_prefix + "/rpc"))
        {
            yCError(OPENXRHEADSET) << "Could not open" << m_prefix + "/rpc" << " RPC port.";
            return false;
        }
    }

    return true;
}

void yarp::dev::OpenXrHeadset::threadRelease()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    if (m_closed)
        return;

    for (auto& hud : m_huds)
    {
        hud.layer.close();
    }

    for (auto& label : m_labels)
    {
        label.layer.close();
    }

    m_openXrInterface.close();

    if (m_tfPublisher)
    {
        m_driver.close();
        m_tfPublisher = nullptr;
    }

    m_headFramePorts.close();
    m_leftHandFramePorts.close();
    m_rightHandFramePorts.close();
    m_eyesManager.close();

    m_rpcPort.close();

    m_closed = true;
}

void yarp::dev::OpenXrHeadset::run()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    if (m_openXrInterface.isRunning())
    {
        if (!m_eyesManager.update()) {
            yCError(OPENXRHEADSET) << "Failed to update eyes.";
            return;
        }

        for (GuiParam& gui : m_huds)
        {
            if (!gui.layer.updateTexture()) {
                yCError(OPENXRHEADSET) << "Failed to update" << gui.portName << "display texture.";
                return;
            }
        }
        for (LabelLayer& label : m_labels)
        {
            if (!label.layer.updateTexture()) {
                yCError(OPENXRHEADSET) << "Failed to update" << label.options.portName << "display texture.";
                return;
            }
        }
        m_openXrInterface.draw();

        m_openXrInterface.getButtons(m_buttons);
        m_openXrInterface.getAxes(m_axes);
        m_openXrInterface.getThumbsticks(m_thumbsticks);
        m_openXrInterface.getAdditionalPoses(m_additionalPosesPublisher.inputs());

        m_stamp.update(m_openXrInterface.currentNanosecondsSinceEpoch() * 1e-9);

        if (m_openXrInterface.shouldResetLocalReferenceSpace())
        {
            //The local reference space has been changed by the user.
            m_rawRootFrameTransform.eye();
        }

        //Publish the transformation from the root frame to the OpenXR root frame
        if (!m_tfPublisher->setTransform(m_rootFrameRaw, m_rootFrame, m_rawRootFrameTransform))
        {
            yCWarning(OPENXRHEADSET) << "Failed to update the transformation of the raw root frame.";
        }

        m_headFramePorts.publishFrame(m_openXrInterface.headPose(), m_openXrInterface.headVelocity(), m_stamp);
        m_leftHandFramePorts.publishFrame(m_openXrInterface.leftHandPose(), m_openXrInterface.leftHandVelocity(), m_stamp);
        m_rightHandFramePorts.publishFrame(m_openXrInterface.rightHandPose(), m_openXrInterface.rightHandVelocity(), m_stamp);

        m_eyesManager.publishEyesTransforms();

        m_additionalPosesPublisher.publishFrames();
    }
    else
    {
        close();
        return;
    }
}

bool yarp::dev::OpenXrHeadset::startService()
{
    //To let the device driver knowing that it need to poll updateService continuosly
    return false;
}

bool yarp::dev::OpenXrHeadset::updateService()
{
    //To let the device driver that we are still alive
    return !m_closed;
}

bool yarp::dev::OpenXrHeadset::stopService()
{
    return this->close();
}

bool yarp::dev::OpenXrHeadset::getAxisCount(unsigned int &axis_count)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    axis_count = m_axes.size();

    if (m_getStickAsAxis)
    {
        axis_count += 2 * m_thumbsticks.size();
    }

    return true;
}

bool yarp::dev::OpenXrHeadset::getButtonCount(unsigned int &button_count)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    button_count = m_buttons.size();

    return true;
}

bool yarp::dev::OpenXrHeadset::getTrackballCount(unsigned int &trackball_count)
{
    trackball_count = 0;

    return true;
}

bool yarp::dev::OpenXrHeadset::getHatCount(unsigned int &hat_count)
{
    hat_count = 0; //These are handled as buttons in OpenXR

    return true;
}

bool yarp::dev::OpenXrHeadset::getTouchSurfaceCount(unsigned int &touch_count)
{
    touch_count = 0;

    return true;
}

bool yarp::dev::OpenXrHeadset::getStickCount(unsigned int &stick_count)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    if (m_getStickAsAxis)
    {
        stick_count = 0;
    }
    else
    {
        stick_count = m_thumbsticks.size();
    }

    return true;
}

bool yarp::dev::OpenXrHeadset::getStickDoF(unsigned int stick_id, unsigned int &dof)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    dof = 2; //Al thumbsticks have two degrees of freedom in OpenXR

    if (m_getStickAsAxis)
    {
        yCError(OPENXRHEADSET) << "The sticks are considered as axis, so there are none.";
        return false;
    }
    else if (stick_id >= m_thumbsticks.size())
    {
        yCError(OPENXRHEADSET) << "The stick_id" << stick_id << "is out of bound. Only" << m_thumbsticks.size() << "sticks are available." ;
        return false;
    }

    return true;
}

bool yarp::dev::OpenXrHeadset::getButton(unsigned int button_id, float &value)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    if (button_id < m_buttons.size())
    {
        value = m_buttons[button_id];
    }
    else
    {
        yCError(OPENXRHEADSET) << "Requested button with index" << button_id << ", but there are" << m_buttons.size() << "buttons.";
        return false;
    }

    return true;
}

bool yarp::dev::OpenXrHeadset::getTrackball(unsigned int /*trackball_id*/, yarp::sig::Vector &value)
{
    value.zero();
    yCError(OPENXRHEADSET) << "No trackball are considered in this device.";
    return false;
}

bool yarp::dev::OpenXrHeadset::getHat(unsigned int /*hat_id*/, unsigned char &value)
{
    value = 0;
    yCError(OPENXRHEADSET) << "No hats are considered in this device.";
    return false;
}

bool yarp::dev::OpenXrHeadset::getAxis(unsigned int axis_id, double &value)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    unsigned int inputId = axis_id;

    if (inputId < m_axes.size())
    {
        value = m_axes[inputId];
    }
    else
    {
        if (m_getStickAsAxis)
        {
            inputId -= m_axes.size();
            if (inputId < 2 * m_thumbsticks.size())
            {
                unsigned int thumbstickId = inputId / 2;

                value = m_thumbsticks[thumbstickId][inputId % 2]; //Each thumbstick counts as two axes
            }
            else
            {
                yCError(OPENXRHEADSET) << "The axis_id" << axis_id << "is out of bounds. There are"
                                       << m_axes.size() << "axes and" << m_thumbsticks.size()
                                       << "thumbsticks (counting as two axes each).";
                return false;
            }
        }
        else
        {
            yCError(OPENXRHEADSET) << "The axis_id" << axis_id << "is out of bounds. There are"
                                   << m_axes.size() << "axes.";
            return false;
        }
    }

    return true;
}

bool yarp::dev::OpenXrHeadset::getStick(unsigned int stick_id, yarp::sig::Vector &value,
                                        yarp::dev::IJoypadController::JoypadCtrl_coordinateMode coordinate_mode)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    if (m_getStickAsAxis)
    {
        yCError(OPENXRHEADSET) << "The sticks are considered axis, so there are none";
        return false;
    }

    if (stick_id < m_thumbsticks.size())
    {
        value.resize(2);
        if (coordinate_mode == JoypadCtrl_coordinateMode::JypCtrlcoord_POLAR)
        {
            value[0] = m_thumbsticks[stick_id].norm();
            value[1] = atan2(m_thumbsticks[stick_id][1], m_thumbsticks[stick_id][0]);
        }
        else
        {
            value[0] = m_thumbsticks[stick_id][0];
            value[1] = m_thumbsticks[stick_id][1];
        }
    }
    else
    {
        yCError(OPENXRHEADSET) << "The stick_id" << stick_id << "is out of bound. Only" << m_thumbsticks.size() << "sticks are available." ;
        return false;
    }

    return true;
}

bool yarp::dev::OpenXrHeadset::getTouch(unsigned int /*touch_id*/, yarp::sig::Vector &value)
{
    value.clear();
    yCError(OPENXRHEADSET) << "No touch devices are considered in this device.";
    return false;
}

std::string yarp::dev::OpenXrHeadset::getLeftHandInteractionProfile()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_openXrInterface.currentLeftHandInteractionProfile();
}

std::string yarp::dev::OpenXrHeadset::getRightHandInteractionProfile()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_openXrInterface.currentRightHandInteractionProfile();
}

std::vector<double> yarp::dev::OpenXrHeadset::getLeftImageDimensions()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_eyesManager.getLeftImageDimensions();
}

std::vector<double> yarp::dev::OpenXrHeadset::getRightImageDimensions()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_eyesManager.getRightImageDimensions();
}

std::vector<double> yarp::dev::OpenXrHeadset::getLeftImageAnglesOffsets()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_eyesManager.getLeftImageAnglesOffsets();
}

std::vector<double> yarp::dev::OpenXrHeadset::getRightImageAnglesOffsets()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_eyesManager.getRightImageAnglesOffsets();
}

bool yarp::dev::OpenXrHeadset::setLeftImageAnglesOffsets(const double azimuth, const double elevation)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_eyesManager.setLeftImageAnglesOffsets(azimuth, elevation);
}

bool yarp::dev::OpenXrHeadset::setRightImageAnglesOffsets(const double azimuth, const double elevation)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_eyesManager.setRightImageAnglesOffsets(azimuth, elevation);
}

bool yarp::dev::OpenXrHeadset::isLeftEyeActive()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_eyesManager.isLeftEyeActive();
}

bool yarp::dev::OpenXrHeadset::isRightEyeActive()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_eyesManager.isRightEyeActive();
}

double yarp::dev::OpenXrHeadset::getEyesZPosition()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_eyesManager.getEyesZPosition();
}

bool yarp::dev::OpenXrHeadset::setEyesZPosition(const double eyesZPosition)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_eyesManager.setEyesZPosition(eyesZPosition);
}

double yarp::dev::OpenXrHeadset::getInterCameraDistance()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_eyesManager.getInterCameraDistance();
}

bool yarp::dev::OpenXrHeadset::setInterCameraDistance(const double distance)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_eyesManager.setInterCameraDistance(distance);
}

std::string yarp::dev::OpenXrHeadset::getLeftImageControlPortName()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_eyesManager.getLeftImageControlPortName();
}

std::string yarp::dev::OpenXrHeadset::getRightImageControlPortName()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_eyesManager.getRightImageControlPortName();
}

bool yarp::dev::OpenXrHeadset::setLabelEnabled(const int32_t labelIndex, const bool enabled)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (labelIndex >= m_labels.size())
    {
        return false;
    }

    m_labels[labelIndex].layer.setEnabled(enabled);

    return true;
}

bool yarp::dev::OpenXrHeadset::alignRootFrameToHeadset()
{
    std::lock_guard<std::mutex> lock(m_mutex);

    OpenXrInterface::Pose headPose = m_openXrInterface.headPose();

    if (!headPose.positionValid || !headPose.rotationValid)
    {
        yCError(OPENXRHEADSET) << "Cannot align the root frame to the headset. The headset pose is not valid";
        return false;
    }

    const Eigen::Matrix3f &root_R_headset = headPose.rotation.matrix();
    // This code was taken from https://www.geometrictools.com/Documentation/EulerAngles.pdf
    // Section 2.2. It computes the XZY inverse kinematics, and we consider only the rotation around Y
    double gravityAngle = 0.0;
    if ((root_R_headset(0,1) < +1.0) && (root_R_headset(0, 1) > -1.0))
    {
        gravityAngle = std::atan2(root_R_headset(0, 2), root_R_headset(0, 0));
    }

    OpenXrInterface::Pose rootFramePose;

    rootFramePose.rotation = Eigen::AngleAxisf(-gravityAngle, Eigen::Vector3f::UnitY());
    rootFramePose.position = -(rootFramePose.rotation * headPose.position);

    poseToYarpMatrix(rootFramePose, m_rawRootFrameTransform);

    return true;
}
