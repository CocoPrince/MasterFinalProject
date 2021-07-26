import pydcel

# pydcel.io.xml2ply("sampledata/washington.xml", "sampledata/washington.ply")

d = pydcel.io.ply2dcel('sampledata/sydney.ply')
# d.handleFaces('off')
d.handleFaces('on')


# sydney
# montreal
# wien                
# washington
# karlsruhe
# london1


# 求分割后的单位角度
def calUnitAngle(self, deg3Start, deg3End, centerPoint, count2Deg):
    dx1 = deg3Start.x - centerPoint.x
    dy1 = deg3Start.y - centerPoint.y
    dx2 = deg3End.x - centerPoint.x
    dy2 = deg3End.y - centerPoint.y
    angle1 = math.atan2(dy1, dx1)
    angle2 = math.atan2(dy2, dx2)
    if angle1*angle2 >= 0:
        included_angle = abs(angle1-angle2)
    else:
        included_angle = abs(angle1) + abs(angle2)
        if included_angle > 180:
            included_angle = 360 - included_angle

    deg3StartANgleFromXAxis = angle1 if angle1 >= 0 else 
    return included_angle / count2Deg, 

def calDeg3AngleFromPositiveXAxis(self, deg3Start, deg3End):





