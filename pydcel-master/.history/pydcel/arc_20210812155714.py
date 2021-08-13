import math

class circle(object):

    def __init__(self, center_x, center_y, radius, arcLength):
        self.center_x = center_x
        self.center_y = center_y
        self.radius = radius   
        self.arcLength = arcLength 


class Arc(object):
    
    def __init__(self, deg3Start, deg3End):
        self.apollonisCircle = None
        self.deg3Start = deg3Start
        self.deg3End = deg3End
        # self.radian = self.calRadian()
        self.startCircleTangency = None  
        self.endCircleTangency = None
        self.deg3StartFlag = None  # 0-pre, 1-post
        self.deg3EndFlag = None
        self.tangencyPointStart = []
        self.tangencyPointEnd = []


    # calculate the arc, that:
    # the ratio of the points on the circle to the centroid = the ratio of the radius
    def calApollonisCircle(self, face_1_centroid, face_2_centroid, radius_1, radius_2):
        ratio = radius_1 / radius_2
        centroidDis = math.sqrt((face_1_centroid.x - face_2_centroid.x) ** 2 + (face_1_centroid.y - face_2_centroid.y) ** 2)
        
        radius = abs(ratio / (1 - ratio ** 2)) * centroidDis
        center_x = (face_1_centroid.x - ratio ** 2 * face_2_centroid.x) / (1 - ratio ** 2)
        center_y = (face_1_centroid.y - ratio ** 2 * face_2_centroid.y) / (1 - ratio ** 2)
        # TODO
        self.apollonisCircle = circle(center_x, center_y, radius, 0)
        # return radius, centerPoint.x, centerPoint.y

               
     # flag of the two deg3 vertices, one rotate clockwise, the other ccw
    def calFlagForDeg3(self):
        dx1 = self.deg3Start.x - self.apollonisCircle.center_x
        dy1 = self.deg3Start.y - self.apollonisCircle.center_y
        dx2 = self.deg3End.x - self.apollonisCircle.center_x
        dy2 = self.deg3End.y - self.apollonisCircle.center_y
        angle1 = math.atan2(dy1, dx1)
        angle1 = int(angle1 * 180/math.pi)
        angle2 = math.atan2(dy2, dx2)
        angle2 = int(angle2 * 180/math.pi)
        included_angle = 0

        self.deg3StartFlag = 0 if angle1 > angle2 else 1
        self.deg3EndFlag = 0 if angle1 <= angle2 else 1
        
        if angle1*angle2 < 0:
            included_angle = abs(angle1) + abs(angle2)
            if included_angle > 180:
                self.deg3StartFlag = 0 if angle1 < angle2 else 1
                self.deg3EndFlag = 0 if angle1 >= angle2 else 1
       

    # Deg3Circle's radius
    def calTangencyCircle(self):
        self.calFlagForDeg3()

        dis1 = math.sqrt( (self.apollonisCircle.center_x - self.deg3Start.x)**2 + (self.apollonisCircle.center_y - self.deg3Start.y)**2)
        deg3StartCircleRadius = dis1 - self.apollonisCircle.radius
        dis2 = math.sqrt( (self.apollonisCircle.center_x - self.deg3End.x)**2 + (self.apollonisCircle.center_y - self.deg3End.y)**2)
        deg3EndCircleRadius = dis2 - self.apollonisCircle.radius

        # angle want to rotate
        alpha = (2 * dis1**2 - deg3StartCircleRadius**2) / 2 * dis1**2
        deg3StartRadian = math.acos(alpha)
        beta = (2 * dis2**2 - deg3EndCircleRadius**2) / 2 * dis2**2
        deg3EndRadian = math.acos(beta)

        if self.deg3StartFlag == 0:
            # clockwise
            tangencyStartCircleCenter_x = self.apollonisCircle.center_x + math.cos(alpha) * (self.deg3Start.x - self.apollonisCircle.center_x) + math.sin(alpha) * (self.deg3Start.y - self.apollonisCircle.center_y)
            tangencyStartCircleCenter_y = self.apollonisCircle.center_y - math.sin(alpha) * (self.deg3Start.x - self.apollonisCircle.center_x) + math.cos(alpha) * (self.deg3Start.y - self.apollonisCircle.center_y)
            # ccw
            tangencyEndCircleCenter_x = self.apollonisCircle.center_x + math.cos(beta) * (self.deg3End.x - self.apollonisCircle.center_x) - math.sin(beta) * (self.deg3End.y - self.apollonisCircle.center_y)
            tangencyEndCircleCenter_y = self.apollonisCircle.center_y + math.sin(beta) * (self.deg3End.x - self.apollonisCircle.center_x) + math.cos(beta) * (self.deg3End.y - self.apollonisCircle.center_y)
        else:
            # ccw
            tangencyStartCircleCenter_x = self.apollonisCircle.center_x + math.cos(alpha) * (self.deg3Start.x - self.apollonisCircle.center_x) - math.sin(alpha) * (self.deg3Start.y - self.apollonisCircle.center_y)
            tangencyStartCircleCenter_y = self.apollonisCircle.center_y + math.sin(alpha) * (self.deg3Start.x - self.apollonisCircle.center_x) + math.cos(alpha) * (self.deg3Start.y - self.apollonisCircle.center_y)
            # clockwise
            tangencyEndCircleCenter_x = self.apollonisCircle.center_x + math.cos(beta) * (self.deg3End.x - self.apollonisCircle.center_x) + math.sin(beta) * (self.deg3End.y - self.apollonisCircle.center_y)
            tangencyEndCircleCenter_y = self.apollonisCircle.center_y - math.sin(beta) * (self.deg3End.x - self.apollonisCircle.center_x) + math.cos(beta) * (self.deg3End.y - self.apollonisCircle.center_y)

        return deg3StartCircleRadius, tangencyStartCircleCenter_x, tangencyStartCircleCenter_y, deg3EndCircleRadius, tangencyEndCircleCenter_x, tangencyEndCircleCenter_y


    # 先求出这两点间的弦长（设为d）：d＝根号下[(x2－x1)²＋(y2－y1)²]
    # 圆心角θ＝2arcsin(d/2r)
    # 弧长L＝rθ＝2r·arcsin(d/2r)
    def calArcLength(self, x1, x2, y1, y2, radius):
        #  二度点个数，弧度，半径，起点，终点
        chordLength = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        chordRadian = 2 * math.asin(chordLength/(2 * radius)) 
        arcLength = chordLength * chordRadian
        return arcLength
        


    def distributeDeg2Vertices(self, chain):
        # Pull the second point past the arc according to the arc and the center of the circle
        deg3StartCircleRadius, tangencyStartCircleCenter_x, tangencyStartCircleCenter_y, deg3EndCircleRadius, tangencyEndCircleCenter_x, tangencyEndCircleCenter_y = self.calTangencyCircle()
        
        apollonisArcLength = self.calArcLength(tangencyStartCircleCenter_x, tangencyEndCircleCenter_x, tangencyStartCircleCenter_y, tangencyEndCircleCenter_y, self.apollonisCircle.radius)
        startArcLength = self.calArcLength(self.deg3Start.x, tangencyStartCircleCenter_x, self.deg3Start.y, tangencyStartCircleCenter_y, deg3StartCircleRadius)
        endArcLength = self.calArcLength(self.deg3End.x, tangencyEndCircleCenter_x, self.deg3End.y, tangencyEndCircleCenter_y, deg3EndCircleRadius)

        totalArcLength = apollonisArcLength + startArcLength + endArcLength

        vertexNumberOnChain = len(chain) - 1
        interval = totalArcLength / vertexNumberOnChain

        # 然后求单位旋转角度：interval/2*Pi

        #########################
        for deg2 in chain:
            leftInterval = 0
            leftLength = startArcLength
            while leftLength >= interval:
                # 计算移动的位置,顺逆时针
                angle = interval / math.pi * 2 * startRadius
                if startFlag == 0:
                    pass
                else:
                    pass
                leftLength = startArcLength - interval
            leftInterval = interval - leftLength

            if leftInterval > 0:
                # 计算移动的位置,顺逆时针
                if startFlag == 0:
                    pass
                else:
                    pass
                leftLength = apollonisArcLength - leftInterval

            while apollonisArcLength >= interval:
                # 计算移动的位置,顺逆时针
                angle = interval / math.pi * 2 * apolloRadius
                if startFlag == 0:
                    pass
                else:
                    pass
                leftLength = apollonisArcLength - interval
            leftInterval = interval - leftLength

            if leftInterval > 0:
                # 计算移动的位置,顺逆时针
                if startFlag == 0:
                    pass
                else:
                    pass
                leftLength = endArcLength - leftInterval

            while endArcLength >= interval:
                # 计算移动的位置,顺逆时针
                angle = interval / math.pi * 2 * endRadius
                if startFlag == 0:
                    pass
                else:
                    pass
                leftLength = apollonisArcLength - interval
            leftInterval = interval - leftLength
        ############################





    # # We have already know the position for the two deg3+ vertices (use force-directed)
    # # move the ApollonisCircle that make the circle pass the two deg3+ vertices
    # # we know the radius and the two points of the arc, want to know the center and the radian
    # def calCenterOfArc(self, Deg3A, Deg3B, r, centroidOfSmallerCircle):
    #     a = Deg3A.x
    #     b = Deg3A.y
    #     c = Deg3B.x
    #     d = Deg3B.y
    #     c1 = math.sqrt(4*r*r - ((a-c) ** 2) + ((b-d) ** 2))
    #     c2 = 2 * math.sqrt((a - c) ** 2 + (b - d) ** 2)
    #     # the two centers calculated by the two deg3 vertices, we only choose one
    #     x0 = (a + c) / 2 + (b - d) * c1 / 2 / c2
    #     y0 = (b + d) / 2 + (c - a) * c1 / 2 / c2
    #     x1 = (a + c) / 2 - (b - d) * c1 / 2 / c2
    #     y1 = (b + d) / 2 - (c - a) * c1 / 2 / c2

    #     # The way we choose the center is, Draw a straight line over Deg3A and Deg3B, 
    #     # and select the center where the arc bends to that side
    #     aLinear, bLinear, cLinear = getLinearEquation(Deg3A, Deg3B)
    #     # check which side of the line has the smaller radius circle
    #     comparatorcentroidOfSmallerCircle = calComparator(aLinear, bLinear, cLinear, centroidOfSmallerCircle.x, centroidOfSmallerCircle.y)
    #     # Determine which center is on the same side as centredOfSmallerCircle
    #     comparatorPoint0 = calComparator(aLinear, bLinear, cLinear, x0, y0)
    #     if comparatorcentroidOfSmallerCircle == comparatorPoint0:
    #         return x0, y0
    #     else:
    #         return x1, y1


    # # Given the coordinates of two points, the equation of the line pass the two points： a*x+b*y+c = 0  (a >= 0)
    # # Deg3A and Deg3B are used to draw a straight line to determine which side of the line centerPoint is on          
    # def getLinearEquation(self, Deg3A, Deg3B):
    #     p1x = Deg3A.x
    #     p1y = Deg3A.y
    #     p2x = Deg3B.x
    #     p2y = Deg3B.y
    #     sign = 1
    #     a = p2y - p1y
    #     if a < 0:
    #         sign = -1
    #         a = sign * a
    #     b = sign * (p1x - p2x)
    #     c = sign * (p1y * p2x - p1x * p2y)
    #     return a, b, c

    # def calComparator(self, a, b, c, x, y):
    #     comparator = 0
    #     calResult = a * x + b * y + c
    #     if calResult != 0:
    #         comparator = 1 if calResult > 0 else -1
    #     return comparator
  
        




    