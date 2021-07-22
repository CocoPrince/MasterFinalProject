class arc(object):
    
    def __init__(self, deg3_start, deg3_end):
        self.radius = None
        self.centerPoint = None
        self.deg3_start = deg3_start
        self.deg3_end = deg3_end
        self.radian = self.calRadian()

    def calApollonisCircle(self, circle_1_centroid, circle_2_centroid, radius_1, radius_2):
        # TODO
        # 改北京题中的公式
        self.radius = None
        self.centerPoint = None 



    # move the ApollonisCircle that make the circle pass the two deg3+ vertices
    def calArc(self):
        # we know the radius and the two points of the arc, want to know the center and the radian
        # 现在准备用画垂线的方式，经过两个三度点的中点画垂线
        pass
        



    def splitArc(self):
        # , 二度点个数，弧度，半径，起点，终点
        pass


    def distributeDeg2Vertices(self):
        # 根据弧度和圆心均分弧，把二度点拉过去
        pass
    