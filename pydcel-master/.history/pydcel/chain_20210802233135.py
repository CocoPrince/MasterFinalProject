import math

class WrapperChain(object):
    
    def __init__(self, chain):
        self.chain = chain
        self.chainType = self.calChainType() # two types: 1-inside, 2-outside
        self.optimalLength = self.calOptimalLength()
        self.threeSectionArc = None # record the correspond relationship between the three arcs and the chain
        self.locateCircleRadius = None # record the correspond relationship between the locatecircle and the chain
        self.locateCircleCenter = []
        
        self.incidentFacesOfChain = []
        # self.deg3Type = self.calDeg3Type()

    def calOptimalLength(self):
        optimalLength = len(self.chain) - 1
        return optimalLength

    def calLocateCircle(self, circle_1_centroid, circle_2_centroid, radius_1, radius_2):
        locateCircleCenter_x = abs((circle_2_centroid.x - circle_1_centroid.x) * radius_1 / (radius_1 + radius_2) + circle_1_centroid.x)
        locateCircleCenter_y = abs((circle_2_centroid.y - circle_1_centroid.y) * radius_1 / (radius_1 + radius_2) + circle_1_centroid.y)        
        self.locateCircleRadius = 1/2 * self.optimalLength
        self.locateCircleCenter.append(locateCircleCenter_x)
        self.locateCircleCenter.append(locateCircleCenter_y)

    def setThreeSectionArc(self, arc):
        self.threeSectionArc = arc


    def calChainType(self):
        iteratorFlag = True
        for hedgeOfFirstPoint in self.chain[0].incidentEdges:
            if not iteratorFlag:
                break
            for hedgeOfSecondPoint in self.chain[1].incidentEdges:
                iteratorFlag = True
                if hedgeOfFirstPoint.twins.identifier == hedgeOfSecondPoint.identifier:
                    # calculate outside chain type
                    if hedgeOfFirstPoint.incidentFace.identifier == 'i' or hedgeOfSecondPoint.incidentFace.identifier == 'i':
                        self.chainType = 2 
                    else:
                        self.chainType = 1
                        self.incidentFacesOfChain.append(hedgeOfFirstPoint.incidentFace)
                        self.incidentFacesOfChain.append(hedgeOfSecondPoint.incidentFace)
                    iteratorFlag = False
                    break


    def calDeg3Type(self):
        pass


    # calculate the intersection of locateCircle and appoloCircle
    def calLocateIntersection(self, deg3vertex):
        x = self.locateCircleCenter[0]   # center of locateCircle
        y = self.locateCircleCenter[1]   
        R = self.locateCircleRadius      # radius of locateCircle
        a = self.arc.appoloCircle.center_x   # center of appoloCircle
        b = self.arc.appoloCircle.center_y
        S = self.arc.appoloCircle.radius      # radius of appoloCircle
        d = math.sqrt((abs(a - x)) ** 2 + (abs(b - y)) ** 2)    
        A = (R ** 2 - S ** 2 + d ** 2) / (2 * d)
        h = math.sqrt(R ** 2 - A ** 2)
        x2 = x + A * (a - x) / d
        y2 = y + A * (b - y) / d
        x3 = round(x2 - h * (b - y) / d, 2)
        y3 = round(y2 + h * (a - x) / d, 2)
        x4 = round(x2 + h * (b - y) / d, 2)
        y4 = round(y2 - h * (a - x) / d, 2)
        c1 = [x3, y3]
        c2 = [x4, y4]

        return c1, c2
