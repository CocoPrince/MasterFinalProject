import math

class arc(object):
    
    def __init__(self, deg3_start, deg3_end):
        self.radius = None
        self.centerPoint = None
        self.deg3_start = deg3_start
        self.deg3_end = deg3_end
        self.radian = self.calRadian()

    # calculate the arc, that:
    # the ratio of the points on the circle to the centroid = the ratio of the radius
    def calApollonisCircle(self, circle_1_centroid, circle_2_centroid, radius_1, radius_2):
        ratio = radius_1 / radius_2
        centroidDis = math.sqrt((circle_1_centroid.x - circle_2_centroid.x) ** 2 + (circle_1_centroid.y - circle_2_centroid.y) ** 2)
        self.radius = abs(ratio / (1 - ratio ** 2)) * centroidDis
        self.centerPoint.x = (circle_1_centroid.x - ratio ** 2 * circle_2_centroid.x) / (1 - ratio ** 2)
        self.centerPoint.y = (circle_1_centroid.y - ratio ** 2 * circle_2_centroid.y) / (1 - ratio ** 2)


    # We have already know the position for the two deg3+ vertices (use force-directed)
    # move the ApollonisCircle that make the circle pass the two deg3+ vertices
    # we know the radius and the two points of the arc, want to know the center and the radian

    def calArc(self, centroid3DegA, centroid3DegB, r, centroidSmallerCircle):

        a = centroid3DegA.x
        b = centroid3DegA.y
        c = centroid3DegB.x
        d = centroid3DegB.y
        c1 = math.sqrt(4*r*r - ((a-c) ** 2) + ((b-d) ** 2))
        c2 = 2 * math.sqrt((a - c) ** 2 + (b - d) ** 2)
        x0 = (a + c) / 2 + (b - d) * c1 / 2 / c2
        y0 = (b + d) / 2 + (c - a) * c1 / 2 / c2
        x1 = (a + c) / 2 - (b - d) * c1 / 2 / c2
        y1 = (b + d) / 2 - (c - a) * c1 / 2 / c2

        # 在两个三度点centroid3DegA, centroid3DegB之间画一条线（直线方程）
        aLinear, bLinear, cLinear = getLinearEquation(centroid3DegA, centroid3DegB)
        # 判断半径更小的圆centroidSmallerCircle的圆心在线的哪一侧
        comparatorSmallerCircle = calComparator(aLinear, bLinear, cLinear, centroidSmallerCircle.x, centroidSmallerCircle.y)

        # 判断哪个点和centroidSmallerCircle在同一侧
        comparatorPoint0 = calComparator(aLinear, bLinear, cLinear, x0, y0)

        if comparatorSmallerCircle == comparatorPoint0:
            return x0, y0
        else:
            return x1, y1


    # 根据已知两点坐标，求过这两点的直线解析方程： a*x+b*y+c = 0  (a >= 0)
    # centroid3DegA与centroid3DegB用于画直线，判断centerPoint在直线的哪一测

           
    def getLinearEquation(self, centroid3DegA, centroid3DegB):
        p1x = centroid3DegA.x
        p1y = centroid3DegA.y
        p2x = centroid3DegB.x
        p2y = centroid3DegB.y
        sign = 1
        a = p2y - p1y
        if a < 0:
            sign = -1
            a = sign * a
        b = sign * (p1x - p2x)
        c = sign * (p1y * p2x - p1x * p2y)

        return a, b, c

    def calComparator(self, a, b, c, x, y):
        comparator = 0
        calResult = a * x + b * y + c

        if calResult != 0:
            comparator = 1 if calResult > 0 else -1
        return comparator
        pass
        



    def splitArc(self):
        # , 二度点个数，弧度，半径，起点，终点
        pass


    def distributeDeg2Vertices(self):
        # 根据弧度和圆心均分弧，把二度点拉过去
        pass
    