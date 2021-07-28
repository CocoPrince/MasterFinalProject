class chain(object):
    
    def __init__(self, chain):
        self.chain = chain
        self.chainType = chainType # three types: 0-none, 1-inside, 2-outside
        self.optimalLength = self.calOptimalLength()
        self.threeSectionArc = None # record the correspond relationship between the three arcs and the chain
        self.locateCircle = None # record the correspond relationship between the locatecircle and the chain


    def calOptimalLength(self):



    def calLocateCircle(self, circle_1_centroid, circle_2_centroid, radius_1, radius_2):


    def setThreeSectionArc(self, arc):
        self.threeSectionArc = arc


    def calChainType(self)
