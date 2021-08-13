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
        self.calFlagForDeg3ByDis(dx1, dy1, dx2, dy2)

    def calFlagForDeg3ByDis(self, dx1, dy1, dx2, dy2):
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
        alpha = (dis1**2 + dis1**2 - deg3StartCircleRadius**2) / (2 * dis1 * dis1) # cos C
        deg3StartRadian = math.acos(alpha)
        beta = (dis2**2 + dis2**2 - deg3EndCircleRadius**2) / (2 * dis2 * dis2)
        deg3EndRadian = math.acos(beta)

        if self.deg3StartFlag == 0:
            # clockwise
            tangencyStartCircleCenter_x = self.apollonisCircle.center_x + math.cos(deg3StartRadian) * (self.deg3Start.x - self.apollonisCircle.center_x) + math.sin(deg3StartRadian) * (self.deg3Start.y - self.apollonisCircle.center_y)
            tangencyStartCircleCenter_y = self.apollonisCircle.center_y - math.sin(deg3StartRadian) * (self.deg3Start.x - self.apollonisCircle.center_x) + math.cos(deg3StartRadian) * (self.deg3Start.y - self.apollonisCircle.center_y)
            # ccw
            tangencyEndCircleCenter_x = self.apollonisCircle.center_x + math.cos(deg3EndRadian) * (self.deg3End.x - self.apollonisCircle.center_x) - math.sin(deg3EndRadian) * (self.deg3End.y - self.apollonisCircle.center_y)
            tangencyEndCircleCenter_y = self.apollonisCircle.center_y + math.sin(deg3EndRadian) * (self.deg3End.x - self.apollonisCircle.center_x) + math.cos(deg3EndRadian) * (self.deg3End.y - self.apollonisCircle.center_y)
        else:
            # ccw
            tangencyStartCircleCenter_x = self.apollonisCircle.center_x + math.cos(deg3StartRadian) * (self.deg3Start.x - self.apollonisCircle.center_x) - math.sin(deg3StartRadian) * (self.deg3Start.y - self.apollonisCircle.center_y)
            tangencyStartCircleCenter_y = self.apollonisCircle.center_y + math.sin(deg3StartRadian) * (self.deg3Start.x - self.apollonisCircle.center_x) + math.cos(deg3StartRadian) * (self.deg3Start.y - self.apollonisCircle.center_y)
            # clockwise
            tangencyEndCircleCenter_x = self.apollonisCircle.center_x + math.cos(deg3EndRadian) * (self.deg3End.x - self.apollonisCircle.center_x) + math.sin(deg3EndRadian) * (self.deg3End.y - self.apollonisCircle.center_y)
            tangencyEndCircleCenter_y = self.apollonisCircle.center_y - math.sin(deg3EndRadian) * (self.deg3End.x - self.apollonisCircle.center_x) + math.cos(deg3EndRadian) * (self.deg3End.y - self.apollonisCircle.center_y)
    
        return deg3StartCircleRadius, tangencyStartCircleCenter_x, tangencyStartCircleCenter_y, deg3EndCircleRadius, tangencyEndCircleCenter_x, tangencyEndCircleCenter_y, dis1, dis2


    # First calculate the chord length between these two points (d): d = sqrt[(x2－x1)²＋(y2－y1)²]
    # angle: θ = 2arcsin(d/2r)
    # arc length: L = rθ = 2r·arcsin(d/2r)
    def calArcLength(self, x1, x2, y1, y2, radius):
        chordLength = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        chordRadian = 2 * math.asin(chordLength/(2 * radius)) 
        arcLength = radius * chordRadian
        return arcLength
        

    def calDividePointCCW(self, ox, oy, px, py, angle):
        qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
        qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
        return qx, qy
        
    def calDividePointCW(self, ox, oy, px, py, angle):  
        qx = ox + math.cos(angle) * (px - ox) + math.sin(angle) * (py - oy)
        qy = oy - math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
        return qx, qy


    def distributeDeg2Vertices(self, chain):
        deg3StartCircleRadius, tangencyStartCircleCenter_x, tangencyStartCircleCenter_y, deg3EndCircleRadius, tangencyEndCircleCenter_x, tangencyEndCircleCenter_y, dis1, dis2 = self.calTangencyCircle()
        
        # calculate the two tangencyPoint
        tangencyPointStart_x = (self.apollonisCircle.radius * (tangencyStartCircleCenter_x - self.apollonisCircle.center_x))/dis1 + self.apollonisCircle.center_x
        tangencyPointStart_y = (self.apollonisCircle.radius * (tangencyStartCircleCenter_y - self.apollonisCircle.center_y))/dis1 + self.apollonisCircle.center_y
        tangencyPointEnd_x = (self.apollonisCircle.radius * (tangencyEndCircleCenter_x - self.apollonisCircle.center_x))/dis2 + self.apollonisCircle.center_x
        tangencyPointEnd_y = (self.apollonisCircle.radius * (tangencyEndCircleCenter_y - self.apollonisCircle.center_y))/dis2 + self.apollonisCircle.center_y
        print("-------------------两个切点-----------------------")
        print(tangencyPointStart_x,tangencyPointStart_y,tangencyPointEnd_x,tangencyPointEnd_y)
        print("------------------------------两个三度点-----------------")
        print(self.deg3Start.x,self.deg3Start.y,self.deg3End.x,self.deg3End.y)
       
        # calculate the total three arc length
        apollonisArcLength = self.calArcLength(tangencyPointStart_x, tangencyPointEnd_x, tangencyPointStart_y, tangencyPointEnd_y, self.apollonisCircle.radius)
        startArcLength = self.calArcLength(self.deg3Start.x, tangencyPointStart_x, self.deg3Start.y, tangencyPointStart_y, deg3StartCircleRadius)
        endArcLength = self.calArcLength(self.deg3End.x, tangencyPointEnd_x, self.deg3End.y, tangencyPointEnd_y, deg3EndCircleRadius)

        totalArcLength = apollonisArcLength + startArcLength + endArcLength
        print("-----------------------三段弧长-----------------------")
        print(startArcLength, apollonisArcLength, endArcLength)
        vertexNumberOnChain = len(chain) - 1
        # the arc length between every two vertices on the chain
        interval = totalArcLength / vertexNumberOnChain

        #########################
        index = 1
        leftInterval = 0
        leftLength = startArcLength
        chain[index-1].x = self.deg3Start.x
        chain[index-1].y = self.deg3Start.y
        while leftLength >= interval: # still on first arc
            # Calculate the position of the movement, cw or ccw
            angle = interval / math.pi * 2
            if self.deg3StartFlag == 0:   # ccw
                chain[index].x, chain[index].y = self.calDividePointCCW(tangencyStartCircleCenter_x, tangencyStartCircleCenter_y, chain[index-1].x, chain[index-1].y, angle)
                index += 1
            else:   # clockwise
                chain[index].x, chain[index].y = self.calDividePointCW(tangencyStartCircleCenter_x, tangencyStartCircleCenter_y, chain[index-1].x, chain[index-1].y, angle)
                index += 1
            leftLength = leftLength - interval

        leftInterval = interval - leftLength  # second arc start

        if leftInterval > 0 and index < len(chain) - 1: # only one step, rotate arc length is leftInterval
            # Calculate the position of the movement, cw or ccw
            angle = leftInterval / math.pi * 2
            if self.deg3StartFlag == 0:  # cw
                chain[index].x, chain[index].y = self.calDividePointCW(self.apollonisCircle.center_x, self.apollonisCircle.center_y, tangencyPointStart_x, tangencyPointStart_y, angle)
                index += 1
            else:   # ccw
                chain[index].x, chain[index].y = self.calDividePointCCW(self.apollonisCircle.center_x, self.apollonisCircle.center_y, tangencyPointStart_x, tangencyPointStart_y, angle)
                index += 1
        leftLength = apollonisArcLength - leftInterval

        while leftLength > interval and index < len(chain) - 1: # on second arc
            # Calculate the position of the movement, cw or ccw
            angle = interval / math.pi * 2
            if self.deg3StartFlag == 0: # cw
                chain[index].x, chain[index].y = self.calDividePointCW(self.apollonisCircle.center_x, self.apollonisCircle.center_y, chain[index].x, chain[index].y, angle)
                index += 1
            else:  # ccw
                chain[index].x, chain[index].y = self.calDividePointCCW(self.apollonisCircle.center_x, self.apollonisCircle.center_y, chain[index].x, chain[index].y, angle)
                index += 1
            leftLength = leftLength - interval

        leftInterval = interval - leftLength # start third arc

        if leftInterval > 0 and index < len(chain) - 1: # only one step, rotate arc length is leftInterval
            # Calculate the position of the movement, cw or ccw
            angle = leftInterval / math.pi * 2
            if self.deg3StartFlag == 0: # ccw
                chain[index].x, chain[index].y = self.calDividePointCCW(tangencyEndCircleCenter_x, tangencyEndCircleCenter_y, tangencyPointEnd_x, tangencyPointEnd_x, angle)
                index += 1
            else:    # cw
                chain[index].x, chain[index].y = self.calDividePointCW(tangencyEndCircleCenter_x, tangencyEndCircleCenter_y, tangencyPointEnd_x, tangencyPointEnd_x, angle)
                index += 1
        leftLength = endArcLength - leftInterval

        while leftLength > interval and index < len(chain) - 1:
            # Calculate the position of the movement, cw or ccw
            angle = interval / math.pi * 2
            if self.deg3StartFlag == 0: # ccw
                chain[index].x, chain[index].y = self.calDividePointCCW(tangencyEndCircleCenter_x, tangencyEndCircleCenter_y, chain[index-1].x, chain[index-1].y, angle)
                # index += 1
            else: # cw
                chain[index].x, chain[index].y = self.calDividePointCW(tangencyEndCircleCenter_x, tangencyEndCircleCenter_y, chain[index-1].x, chain[index-1].y, angle)
                # index += 1
            leftLength = leftLength - interval
        # leftInterval = interval - leftLength # should be 0
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
  
        




    