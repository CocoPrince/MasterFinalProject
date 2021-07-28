class WrapperChain(object):
    
    def __init__(self, chain):
        self.chain = chain
        self.chainType = self.calChainType() # three types: 0-none, 1-inside, 2-outside
        self.optimalLength = self.calOptimalLength()
        self.threeSectionArc = None # record the correspond relationship between the three arcs and the chain
        self.locateCircleRadius = None # record the correspond relationship between the locatecircle and the chain
        self.locateCircleCenter = []
        self.deg3Type

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
        if len(chain) == 2:
            chainType = 0
        elif len(chain[1].incidentFaces) > 1:
            chainType = 1
        else:
            chainType = 2

       
    def calLocateIntersection(self, deg3vertex):
        
