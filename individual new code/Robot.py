__author__ = "Xuan Lin"
__email__ = "xuanlin1991@gmail.com"
__copyright__ = "Copyright 2017 RoMeLa"
__date__ = "July 1, 2017"

__version__ = "0.1.0"
__status__ = "Production"

"""
Robot class for hexapod code
"""

import numpy as np
import sys
sys.path.insert(0, 'dxl_manager')
import dxl_manager.dcm_controller as dcmc
import math

Run = 1

# The body class of the robot
class pointXY:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class pointXYZ:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        #This function also enables adding pointXYZ with pointXYZR
        total_x = self.x + other.x
        total_y = self.y + other.y
        total_z = self.z + other.z
        return pointXYZ(total_x, total_y, total_z)

    def __sub__(self, other):
        #This function also enables adding pointXYZ with pointXYZR
        total_x = self.x - other.x
        total_y = self.y - other.y
        total_z = self.z - other.z
        return pointXYZ(total_x, total_y, total_z)

    def arr(self):
        return np.array([[self.x], [self.y], [self.z]])

class pointXYZR:
    def __init__(self, x, y, z, r):
        self.x = x
        self.y = y
        self.z = z
        self.r = r

class motor123:
    def __init__(self, coxa, femur, tibia):
        self.coxa = coxa
        self.femur = femur
        self.tibia = tibia

    def __add__(self, other):
        #This function also enables adding pointXYZ with pointXYZR
        total_coxa = self.coxa + other.coxa
        total_femur = self.femur + other.femur
        total_tibia = self.tibia + other.tibia
        return motor123(total_coxa, total_femur, total_tibia)

    def __sub__(self, other):
        #This function also enables adding pointXYZ with pointXYZR
        total_coxa = self.coxa - other.coxa
        total_femur = self.femur - other.femur
        total_tibia = self.tibia - other.tibia
        return motor123(total_coxa, total_femur, total_tibia)

    def __div__(self, other):
        total_coxa = self.coxa/other
        total_femur = self.femur/other
        total_tibia = self.tibia/other
        return motor123(total_coxa, total_femur, total_tibia)

    # def min(self, other):
    #     if self.coxa > other.coxa:
    #         c = other.coxa
    #     else:
    #         c = self.coxa
    #     if self.femur > other.femur:
    #         f = other.femur
    #     else:
    #         f = self.femur
    #
    #     if self.tibia > other.tibia:
    #         t = other.tibia
    #     else:
    #         t = self.tibia
    #     return motor123(c,f,t)


def xyz_RotationMatrix(rotX, rotY, rotZ):
    Rx_bw = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(rotX), np.sin(rotX)], [0.0, -np.sin(rotX), np.cos(rotX)]])

    Ry_bw = np.array([[np.cos(rotY), 0.0, -np.sin(rotY)], [0.0, 1.0, 0.0], [np.sin(rotY), 0.0, np.cos(rotY)]])

    Rz_bw = np.array([[np.cos(rotZ), np.sin(rotZ), 0.0], [-np.sin(rotZ), np.cos(rotZ), 0], [0.0, 0.0, 1.0]])

    R_bw = np.dot(np.dot(Rz_bw, Ry_bw), Rx_bw)

    return R_bw

class Body:
    #Body dimension, common for all instances

    Y_COXA = 194  # Distance between front and back legs /2, in mm
    X_COXA = 111  # Distance between two front(or back) legs /2, in mm
    M_COXA = 224  # Distance between two middle legs /2, in mm
    #Body vertecies
    bodyVertices = [];
    #RF
    RF = pointXYZ(X_COXA, Y_COXA, 0.0)
    bodyVertices.append(RF)
    #RM
    RM = pointXYZ(M_COXA, 0, 0.0)
    bodyVertices.append(RM)
    #RR
    RR = pointXYZ(X_COXA, -Y_COXA, 0.0)
    bodyVertices.append(RR)
    #LR
    LR = pointXYZ(-X_COXA, -Y_COXA, 0.0)
    bodyVertices.append(LR)
    #LM
    LM = pointXYZ(-M_COXA, 0, 0.0)
    bodyVertices.append(LM)
    #LF
    LF = pointXYZ(-X_COXA, Y_COXA, 0.0)
    bodyVertices.append(LF)

class Leg:
    #Leg length, common for all instances
    # unit:mm

    L_coxa = 57
    L_femur = 199
    L_tibia = 381

    def __init__(self, installation_0):
        #For each leg, installazion zero is different
        self.theta_coxa0 = installation_0.coxa
        self.theta_femur0 = installation_0.femur
        self.theta_tibia0 = installation_0.tibia

#FK member function
def leg_FK(motorPos, installation_0):
    "This function takes three joint angles and compute the transformation\
        matrix from origin to endeffector, returning T0_1, T1_2, T2_3, T3_ef will be enough information"
    ## Angle Definition

    theta_coxa = motorPos.coxa
    theta_femur = motorPos.femur
    theta_tibia = motorPos.tibia

    #The installation zero here is not correct!!!!!!!
    theta_coxa_tot = theta_coxa + installation_0.coxa
    theta_femur_tot = theta_femur + installation_0.femur
    theta_tibia_tot = theta_tibia + installation_0.tibia

    ## DH Parameters
    alpha = [0, np.pi / 2, 0, np.pi / 2]
    a = [0, Leg.L_coxa, Leg.L_femur, 0]
    d = [0, 0, 0, Leg.L_tibia]
    theta = [theta_coxa_tot, theta_femur_tot, theta_tibia_tot, 0]

    i = 0
    T0_1 = np.array([[np.cos(theta[i]), -np.sin(theta[i]), 0, a[i]],
                     [np.sin(theta[i]) * np.cos(alpha[i]), np.cos(theta[i]) * np.cos(alpha[i]), -np.sin(alpha[i]),
                      -np.sin(alpha[i]) * d[i]],
                     [np.sin(theta[i]) * np.sin(alpha[i]), np.cos(theta[i]) * np.sin(alpha[i]), np.cos(alpha[i]),
                      np.cos(alpha[i]) * d[i]],
                     [0, 0, 0, 1]])

    i = 1
    T1_2 = np.array([[np.cos(theta[i]), -np.sin(theta[i]), 0, a[i]],
                     [np.sin(theta[i]) * np.cos(alpha[i]), np.cos(theta[i]) * np.cos(alpha[i]), -np.sin(alpha[i]),
                      -np.sin(alpha[i]) * d[i]],
                     [np.sin(theta[i]) * np.sin(alpha[i]), np.cos(theta[i]) * np.sin(alpha[i]), np.cos(alpha[i]),
                      np.cos(alpha[i]) * d[i]],
                     [0, 0, 0, 1]])

    i = 2
    T2_3 = np.array([[np.cos(theta[i]), -np.sin(theta[i]), 0, a[i]],
                     [np.sin(theta[i]) * np.cos(alpha[i]), np.cos(theta[i]) * np.cos(alpha[i]), -np.sin(alpha[i]),
                      -np.sin(alpha[i]) * d[i]],
                     [np.sin(theta[i]) * np.sin(alpha[i]), np.cos(theta[i]) * np.sin(alpha[i]), np.cos(alpha[i]),
                      np.cos(alpha[i]) * d[i]],
                     [0, 0, 0, 1]])

    i = 3
    T3_ef = np.array([[np.cos(theta[i]), -np.sin(theta[i]), 0, a[i]],
                      [np.sin(theta[i]) * np.cos(alpha[i]), np.cos(theta[i]) * np.cos(alpha[i]), -np.sin(alpha[i]),
                       -np.sin(alpha[i]) * d[i]],
                      [np.sin(theta[i]) * np.sin(alpha[i]), np.cos(theta[i]) * np.sin(alpha[i]), np.cos(alpha[i]),
                       np.cos(alpha[i]) * d[i]],
                      [0, 0, 0, 1]])

    return [T0_1, T1_2, T2_3, T3_ef]

#IK member function
def leg_IK(toeXYZ, installation_0):
    " Solve Inverse Kinematics, give three joint angles \
       Note here that X, Y, Z are values w.r.t the coxa zero point (not body center of mass)\
       In the axis base frame (X0-Y0-Z0)"
    # The actual motor angle is angle_bodyFrame + installation_0

    x = toeXYZ.x
    y = toeXYZ.y
    z = toeXYZ.z

    offset = 9.3
    ###Add in triangular error check
    error = 0

    # first, make this a 2DOF problem... by solving coxa
    coxa = np.arctan2(y, x) + 1.0*installation_0.coxa/180*np.pi  # The grammar of atan2 is atan2(Y,X)

    if np.imag(coxa) != 0:
        error = 11

    trueX = np.sqrt(x ** 2 + y ** 2) - Leg.L_coxa
    im = np.sqrt((trueX ** 2) + (z ** 2))  # length of imaginary leg

    # get femur angle above horizon...
    q1 = np.arctan2(z, trueX)  # Angle between tibia and ground (negative)
    d1 = Leg.L_femur ** 2 + im ** 2 - Leg.L_tibia ** 2
    d2 = 2 * Leg.L_femur * im
    q2 = np.arccos(d1 / d2)  # acos always gives a positive angle
    femur = - (q1 + q2) + 1.0*installation_0.femur/180*np.pi - offset/180*np.pi  # This angle is supposed to be positive
    # It is large positive(q2) + small negative(q1) =
    # positive

    # if (imag(femur)~=0) || (femur>(100/180*pi)) || (femur<-pi/2)  #This range
    # doesn't work well...
    if np.imag(femur) != 0:
        error = 12

    # and tibia angle from femur...

    d1 = Leg.L_femur ** 2 + Leg.L_tibia ** 2 - im ** 2
    d2 = 2 * Leg.L_tibia * Leg.L_femur
    tibia = np.arccos(d1 / d2) - np.pi / 2 + 1.0*installation_0.tibia/180*np.pi - offset/180*np.pi
    # if (imag(tibia)~=0) || (tibia>pi*1/2) || (tibia<-pi*1/2)
    if np.imag(tibia) != 0:
        error = 13

    if coxa < 0:
        coxa += 2.0*np.pi
    elif coxa > 2.0*np.pi:
        coxa -= 2.0*np.pi

    if femur < 0:
        femur += 2.0*np.pi
    elif femur > 2.0*np.pi:
        femur -= 2.0*np.pi

    if tibia < 0:
        tibia += 2.0*np.pi
    elif tibia > 2.0*np.pi:
        tibia -= 2.0*np.pi

    ####Should return coxa+self.theta_coxa0, add in installation zeros!
    #Throw exception instread of give back error variable?
    return motor123(coxa-np.pi, femur-np.pi , tibia-np.pi)

#Jacobian member function
def leg_motorTorque(coxa, femur, tibia, f, R0_ef):
    "This function takes the joint angle values and calculate motor torque\
    by calculate the Jacobian"
    # The algebra is verified

    torqueCoxa_Lim = 20 * 1000  # unit: N*mm
    torqueFemur_Lim = 20 * 1000
    torqueTibia_Lim = 20 * 1000

    torqueCoxa_Lim = 20 * 1000  # unit:N*mm
    torqueFemur_Lim = 20 * 1000
    torqueTibia_Lim = 20 * 1000

    M_2motorAssembly = 0.33  # unit:kg

    # Jacobian in end effector frame
    Jacobian_EF = np.array([[0, -Leg.L_coxa - Leg.L_tibia * np.sin(femur + tibia) - Leg.L_femur * np.cos(femur), 0,
                             np.sin(femur + tibia), 0, -np.cos(femur + tibia)],
                            [Leg.L_tibia + Leg.L_femur * np.sin(tibia), 0, -Leg.L_femur * np.cos(tibia), 0, 1, 0],
                            [Leg.L_tibia, 0, 0, 0, 1, 0]])
    # Gravity vector
    Tgravity = np.array([[0], [Leg.L_femur * M_2motorAssembly * np.cos(femur)], [0]])

    T66_EF_0 = np.bmat([[R0_ef.T, np.zeros((3, 3))], [np.zeros((3, 3)), R0_ef.T]])
    Wrench = np.bmat([[f], [np.zeros((3, 1))]])
    torque = Jacobian_EF.dot(T66_EF_0).dot(Wrench) + Tgravity

    return torque

class Robot:
    #Define constants

    #Leg numbers
    RIGHT_FRONT = 0
    RIGHT_MIDDLE = 1
    RIGHT_REAR = 2
    LEFT_REAR = 3
    LEFT_MIDDLE = 4
    LEFT_FRONT = 5
    legNo = [RIGHT_FRONT, RIGHT_MIDDLE, RIGHT_REAR, LEFT_REAR, LEFT_MIDDLE, LEFT_FRONT]

    #Gaits
    Ripple = 0
    Amble = 1
    Tripod = 2
    AmbleClimb = 3

    def gaitParameterGenerator(self):
        #Each time the gait is changed, should run this
        if self.gait == Robot.Ripple:
            self.legGaitNo = [0,4,8,12,16,20]
            self.stepsInCycle = 24
            self.pushSteps = 20
        elif self.gait == Robot.Amble:
            self.legGaitNo = [0, 4, 8, 0, 4, 8]
            self.stepsInCycle = 12
            self.pushSteps = 8
        elif self.gait == Robot.AmbleClimb:
            self.legGaitNo = [0, 4, 2, 2, 4, 0]
            self.stepsInCycle = 12
            self.pushSteps = 6
        elif self.gait == Robot.Tripod:
            self.legGaitNo = [0, 4, 0, 4, 0, 4]
            self.stepsInCycle = 8
            self.pushSteps = 4


    Body = Body()  # Import body class into robot class

    def __init__(self):

        if Run == 1:
            dcmc.main_import()  # Turn on controller
            for i in range(18):       # Enable torque
                dcmc.set_torque_enable((i+1, 1))
                dcmc.set_p_gain((i+1, 32))
                dcmc.set_d_gain((i+1, 10))

        self.legInstallation_0 = []   #Leg installation zeros, in degs
        # RF
        self.legInstallation_0.append(motor123(120, 180, 180))
        # RM
        self.legInstallation_0.append(motor123(180, 180, 180))
        # RR
        self.legInstallation_0.append(motor123(240, 180, 180))
        # LR
        self.legInstallation_0.append(motor123(-60, 180, 180))
        # LM
        self.legInstallation_0.append(motor123(0, 180, 180))
        # LF
        self.legInstallation_0.append(motor123(60, 180, 180))

        self.gait = Robot.Tripod
        self.gaitParameterGenerator()

        self.step = self.stepsInCycle-1
        self.hardWareTime = 100.0   #btw key frames, unit:millisecends
        self.angleSprawl = 60.0
        self.angleFemur = 70.0   #The angle that femur joint makes w.r.t. horizontal plane, in degrees
        self.angleTibia = 75.0   #The angle that tibia joint makes w.r.t. horizontal plane, in degrees
        self.setStaticToePos()

        self.liftHeight = 50.0
        self.walkTime = 0.2

        # Initialize key frame, key frames are sto1-6 as class motor123red by leg
        # There are two types of key frame: XYZ and motor positions
        # Only dynamic gait variable "gait" will have r information
        # make this a member of robot or a member of robot.leg?

        self.toeKeyFrame = []  # 1-6 with structure pointXYZ
        self.motorKeyFrame = []  # 1-6 motor positions
        self.motorCurrentPos = []  # 1-6 motor current positions
        self.motorSpeed = []   #1-6 motor speed
        p0 = pointXYZ(0.0, 0.0, 0.0)
        p1 = motor123(0, 0, 0)
        for i in self.legNo:
            self.toeKeyFrame.append(p0)
            self.motorKeyFrame.append(p1)
            self.motorCurrentPos.append(p1)
            self.motorSpeed.append(p1)

        self.dynamicToe = []  # Dynamic toe posistion variable used by gaitGenerator
        p0 = pointXYZR(0.0, 0.0, 0.0, 0.0)
        for i in self.legNo:
            self.dynamicToe.append(p0)

        #Get speeds and rotation angles from controller class?
        self.xSpeed = 0.0
        self.ySpeed = 0.0
        self.rSpeed = 0.0
        self.bodyShiftX = 0.0
        self.bodyShiftY = 0.0
        self.bodyShiftZ = 0.0
        self.bodyRotX = 0.0
        self.bodyRotY = 0.0
        self.bodyRotZ = 0.0

    def setStaticToePos(self):
        #staticToe is the static toe position w.r.t. shoulder
        r = Leg.L_coxa + Leg.L_femur * np.cos(self.angleFemur * (np.pi) / 180) \
            + Leg.L_tibia * np.cos(self.angleTibia * (np.pi) / 180)

        h = Leg.L_tibia * np.sin(self.angleTibia * (np.pi) / 180) \
            - Leg.L_femur * np.sin(self.angleFemur * (np.pi) / 180)

        zHeight = -h

        xFrontback = r * np.cos(self.angleSprawl * (np.pi) / 180)
        yFrontback = r * np.sin(self.angleSprawl * (np.pi) / 180)
        xMid = r
        yMid = 0

        self.staticToe = []

        # Right front
        RF = pointXYZ(xFrontback, yFrontback, zHeight)
        self.staticToe.append(RF)

        # Right middle
        RM = pointXYZ(xMid, yMid, zHeight)
        self.staticToe.append(RM)

        # Right rear
        RR = pointXYZ(xFrontback, -yFrontback, zHeight)
        self.staticToe.append(RR)

        # Left rear
        LR = pointXYZ(-xFrontback, -yFrontback, zHeight)
        self.staticToe.append(LR)

        # Left middle
        LM = pointXYZ(-xMid, yMid, zHeight)
        self.staticToe.append(LM)

        # Left front
        LF = pointXYZ(-xFrontback, yFrontback, zHeight)
        self.staticToe.append(LF)

    def generateGait(self):
        #Provide dynamicToePos as XYZR information

        ######################This part is the normal walking gait generator############################################
        for whichLeg in self.legNo:
            if self.step == self.legGaitNo[whichLeg]:
                x = self.dynamicToe[whichLeg].x/2
                y = self.dynamicToe[whichLeg].y/2
                z = self.liftHeight/2
                r = self.dynamicToe[whichLeg].r/2
                PP = pointXYZR(x, y, z, r)
                self.dynamicToe[whichLeg] = PP
            elif (self.step - self.legGaitNo[whichLeg]) == 1:
                x = 0.0
                y = 0.0
                z = self.liftHeight
                r = 0.0
                PP = pointXYZR(x, y, z, r)
                self.dynamicToe[whichLeg] = PP
            elif (self.step - self.legGaitNo[whichLeg]) == 2:
                x = self.xSpeed*self.pushSteps*self.hardWareTime/4000.0
                y = self.ySpeed*self.pushSteps*self.hardWareTime/4000.0
                z = self.liftHeight/2
                r = self.rSpeed*self.pushSteps*self.hardWareTime/4000.0
                PP = pointXYZR(x, y, z, r)
                self.dynamicToe[whichLeg] = PP
            elif (self.step - self.legGaitNo[whichLeg]) == 3:
                x = self.xSpeed*self.pushSteps*self.hardWareTime/2000.0
                y = self.ySpeed*self.pushSteps*self.hardWareTime/2000.0
                z = 0.0
                r = self.rSpeed*self.pushSteps*self.hardWareTime/2000.0
                PP = pointXYZR(x, y, z, r)
                self.dynamicToe[whichLeg] = PP
            else:
                x = self.dynamicToe[whichLeg].x - self.xSpeed*self.hardWareTime/1000.0
                y = self.dynamicToe[whichLeg].y - self.ySpeed*self.hardWareTime/1000.0
                z = 0.0
                r = self.dynamicToe[whichLeg].r - self.rSpeed*self.hardWareTime/1000.0
                PP = pointXYZR(x, y, z, r)
                self.dynamicToe[whichLeg] = PP

        ######################This part is used for climbing, pushing body motion with all toes on the wall#############
        # for whichLeg in self.legNo:
        #     if self.step == self.legGaitNo[whichLeg]:
        #         x = 0.0
        #         y = 0.0
        #         z = self.liftHeight
        #         r = 0.0
        #         PP = pointXYZR(x, y, z, r)
        #         self.dynamicToe[whichLeg] = PP
        #     elif (self.step - self.legGaitNo[whichLeg]) == 1:
        #         x = self.xSpeed*self.pushSteps*self.hardWareTime/2000.0
        #         y = self.ySpeed*self.pushSteps*self.hardWareTime/2000.0
        #         z = 0.0
        #         r = self.rSpeed*self.pushSteps*self.hardWareTime/2000.0
        #         PP = pointXYZR(x, y, z, r)
        #         self.dynamicToe[whichLeg] = PP
        #     elif (self.step<=5):
        #         #toe position stays at the same place
        #         x = self.dynamicToe[whichLeg].x
        #         y = self.dynamicToe[whichLeg].y
        #         z = self.dynamicToe[whichLeg].z
        #         r = self.dynamicToe[whichLeg].r
        #         PP = pointXYZR(x, y, z, r)
        #         self.dynamicToe[whichLeg] = PP
        #     else:
        #         #Push body forward
        #         x = self.dynamicToe[whichLeg].x - self.xSpeed*self.hardWareTime/1000.0
        #         y = self.dynamicToe[whichLeg].y - self.ySpeed*self.hardWareTime/1000.0
        #         z = 0.0
        #         r = self.dynamicToe[whichLeg].r - self.rSpeed*self.hardWareTime/1000.0
        #         PP = pointXYZR(x, y, z, r)
        #         self.dynamicToe[whichLeg] = PP

    def bodyIK(self):
        # Give cm_w and body orientation (variables for optimization), feedback toe XYZ position
        # Toepos_chosen is the toe position in world frame
        # Rbw*(X0_b+staticToePos+dynamicToePos-cm_w(or body shift)-Rwb*X0_b)
        for whichLeg in self.legNo:
            R_bw = xyz_RotationMatrix(self.bodyRotX, self.bodyRotY, self.bodyRotZ + self.dynamicToe[whichLeg].r)
            R_wb = np.transpose(R_bw)
            cm_w = pointXYZ(self.bodyShiftX, self.bodyShiftY, self.bodyShiftZ)
            toepos_chosen = Body.bodyVertices[whichLeg] + self.staticToe[whichLeg] + self.dynamicToe[whichLeg]
            pArr = np.dot(R_bw, (
            pointXYZ.arr(toepos_chosen - cm_w) - np.dot(R_wb, pointXYZ.arr(Body.bodyVertices[whichLeg]))))
            self.toeKeyFrame[whichLeg] = pointXYZ(pArr[0], pArr[1], pArr[2])

    def legIK(self):
        #From self.toeKeyFrame calculate self.motorKeyFrame
        for whichLeg in self.legNo:
            self.motorKeyFrame[whichLeg] = leg_IK(self.toeKeyFrame[whichLeg], self.legInstallation_0[whichLeg])

    def runToeKeyFrame(self, tranTime):

        if tranTime == 0:
            #run until next key frame, takes tranTime long
            #To figure out speed, need current position, initialization doesn't have read -> no current position info
            #So just set speed to be very slow
            if Run == 1:
                for i in self.legNo:
                    dcmc.set_moving_speed((3 * i + 1, 8),
                                          (3 * i + 2, 8),
                                          (3 * i + 3, 8))
                for i in self.legNo:
                    dcmc.set_command_position((3 * i + 1, self.motorKeyFrame[i].coxa),
                                              (3 * i + 2, self.motorKeyFrame[i].femur),
                                              (3 * i + 3, self.motorKeyFrame[i].tibia))

                for i in self.legNo:
                    self.motorCurrentPos[i].coxa = self.motorKeyFrame[i].coxa
                    self.motorCurrentPos[i].femur = self.motorKeyFrame[i].femur
                    self.motorCurrentPos[i].tibia = self.motorKeyFrame[i].tibia

            else:
                pass

        else:

            #Acceleration is 40
            a = 45*8.583/180*np.pi
            v_resol = 0.0119

            for i in self.legNo:
                dcoxa = abs(self.motorKeyFrame[i].coxa - self.motorCurrentPos[i].coxa)
                dfemur = abs(self.motorKeyFrame[i].femur - self.motorCurrentPos[i].femur)
                dtibia = abs(self.motorKeyFrame[i].tibia - self.motorCurrentPos[i].tibia)

                vcoxa_raw = (a * (tranTime) - np.sqrt((a ** 2) * ((tranTime) ** 2) - 4 * dcoxa * a)) / (2 * v_resol)
                if math.isnan(vcoxa_raw):
                    vcoxa = 1024
                else:
                    vcoxa = int(vcoxa_raw)

                if vcoxa > 1023:
                    vcoxa = 1023
                    print("leg", i, "coxa motor reach max speed!")
                elif vcoxa < 1:
                    vcoxa = 1
                    #print("leg", i, "coxa motor reach zero speed!")

                vfemur_raw = (a * (tranTime) - np.sqrt((a ** 2) * ((tranTime) ** 2) - 4 * dfemur * a)) / (2 * v_resol)
                if math.isnan(vfemur_raw):
                    vfemur = 1024
                else:
                    vfemur = int(vfemur_raw)

                if vfemur > 1023:
                    vfemur = 1023
                    print("leg", i, "femur motor reach max speed!")
                elif vfemur < 1:
                    vfemur = 1
                    #print("leg", i, "femur motor reach zero speed!")

                vtibia_raw = (a * (tranTime) - np.sqrt((a ** 2) * ((tranTime) ** 2) - 4 * dtibia * a)) / (2 * v_resol)
                if math.isnan(vtibia_raw):
                    vtibia = 1024
                else:
                    vtibia = int(vtibia_raw)

                if vtibia > 1023:
                    vtibia = 1023
                    print("leg", i, "tibia motor reach max speed!")
                elif vtibia < 1:
                    vtibia = 1
                    #print("leg", i, "tibia motor reach zero speed!")

                self.motorSpeed[i] = motor123(vcoxa, vfemur, vtibia)
                if Run == 1:
                    dcmc.set_moving_speed((3*i+1, self.motorSpeed[i].coxa),
                                          (3*i+2, self.motorSpeed[i].femur),
                                          (3*i+3, self.motorSpeed[i].tibia))
                else:
                    pass

            for i in self.legNo:
                if Run == 1:
                    dcmc.set_command_position((3*i+1, self.motorKeyFrame[i].coxa),
                                              (3*i+2, self.motorKeyFrame[i].femur),
                                              (3*i+3, self.motorKeyFrame[i].tibia))
                else:
                    pass



        for i in self.legNo:
            m = motor123(self.motorKeyFrame[i].coxa, self.motorKeyFrame[i].femur, self.motorKeyFrame[i].tibia)
            self.motorCurrentPos[i] = m


    ###################### The following functions accepts outside variables ###########################################
    def readMotorPos(self):
        if Run == 1:
            q = dcmc.get_current_position(range(1,19))
            for i in range(6):
                self.motorCurrentPos = motor123(q[i*3+1], q[i*3+2], q[i*3+3])

    def walkOneStep(self, tranTime):
        self.generateGait()        #Set dynamicToe
        self.bodyIK()              #Set toeKeyFrame
        self.legIK()               #Set motorKeyFrame
        self.runToeKeyFrame(tranTime)   #Run at transition time 100ms
        self.step+=1
        if self.step >= self.stepsInCycle:
            self.step = 0

    def setAngles(self, angle_sprawl, angle_femur, angle_tibia):
        self.angleSprawl = angle_sprawl
        self.angleFemur = angle_femur
        self.angleTibia = angle_tibia
        self.setStaticToePos()

    def setLiftHeight(self, liftHeight):
        self.liftHeight = liftHeight

    def liftHeightadd5(self):
        self.liftHeight = self.liftHeight + 5

    def liftHeightsub5(self):
        self.liftHeight = self.liftHeight - 5

    def setToeKeyFrame(self, tranTime):
        #This is the user interface function which directly set toe keyframe (skipping gererateGait and bodyIK)
        #self.toeKeyFrame = [] #Accept keyframe parameter from outside
        self.legIK()
        self.runToeKeyFrame(tranTime)

    def setMotorKeyFrame(self, motorKeyFrame):
        #This is the user interface which directly accept self.motorKeyFrame from outside
        #The motor keyframe should be a 2-dimensional array, the first dimension is leg number (0-5), the second dimension is motor angles (0-2)
        for whichLeg in self.legNo:
            self.motorKeyFrame[whichLeg] = motor123(motorKeyFrame[whichLeg][0], motorKeyFrame[whichLeg][1], motorKeyFrame[whichLeg][2])

    def setSpeed(self, xSpeed, ySpeed, rSpeed):
        self.xSpeed = xSpeed
        self.ySpeed = ySpeed
        self.rSpeed = rSpeed

    def setOrientation(self, bodyShiftX, bodyShiftY, bodyShiftZ, bodyRotX, bodyRotY, bodyRotZ):
        self.bodyShiftX = bodyShiftX
        self.bodyShiftY = bodyShiftY
        self.bodyShiftZ = bodyShiftZ
        self.bodyRotX = bodyRotX
        self.bodyRotY = bodyRotY
        self.bodyRotZ = bodyRotZ

    def changeOrientationRun(self, tranTime):
        #This function is only for changing orientation, dynamic toe pos should be zero
        self.bodyIK()  # Set toeKeyFrame
        self.legIK()  # Set motorKeyFrame
        self.runToeKeyFrame(tranTime)  # Run at transition time 100ms

    def setGait(self, gait):
        self.gait = gait
        self.gaitParameterGenerator()

    def changeGait(self):
        self.gait = self.gait+1
        if self.gait >= 3:
            self.gait = 0
        self.gaitParameterGenerator()

    def angleSprawladd10(self):
        self.angleSprawl = self.angleSprawl + 10

    def angleSprawlsub10(self):
        self.angleSprawl = self.angleSprawl - 10

    def increaseWalkTime(self):
        #increase the walking loop time by 0.1 second
        self.walkTime = self.walkTime + 0.1

    def decreaseWalkTime(self):
        #increase the walking loop time by 0.1 second
        self.walkTime = self.walkTime - 0.1
        if self.walkTime < 0.1:
            self.walkTime = 0.1
