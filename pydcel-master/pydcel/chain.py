import math
from pydcel.arc import Arc
from .arc import Arc

class WrapperChain(object):
    
    def __init__(self, chain, edgeDict):
        self.chainId = 0
        self.chain = chain
        self.incidentFacesOfChain = []
        self.chainType = self.calChainType(edgeDict) # two types: 1-inside, 2-outside
        self.optimalLength = self.calOptimalLength()
        self.threeSectionArc = None # record the correspond relationship between the three arcs and the chain
        self.locateCircleRadius = None # record the correspond relationship between the locatecircle and the chain
        self.locateCircleCenter = []
        self.deg3Type = self.calDeg3Type() # two types: 1-inside(except cases of outside type), 2-outside(as long as one of the deg3 vertex is adjacent to the infinate face)

    def geneChainId(self, id):
        self.chainId = id

    def calOptimalLength(self):
        optimalLength = (len(self.chain) - 1) * 10
        return optimalLength

    def calLocateCircleOfInsideChain(self, circle_1_centroid, circle_2_centroid, radius_1, radius_2):
        locateCircleCenter_x = abs((circle_2_centroid.x - circle_1_centroid.x) * radius_1 / (radius_1 + radius_2) + circle_1_centroid.x)
        locateCircleCenter_y = abs((circle_2_centroid.y - circle_1_centroid.y) * radius_1 / (radius_1 + radius_2) + circle_1_centroid.y)        
        self.locateCircleRadius = 1/2 * self.optimalLength
        self.locateCircleCenter.append(locateCircleCenter_x)
        self.locateCircleCenter.append(locateCircleCenter_y)

    def setThreeSectionArc(self, arc):
        self.threeSectionArc = arc


    def calChainType(self, edgeDict):
        for hedgeIdOfFirstPoint in self.chain[0].incidentEdges:
            hedgeOfFirstPoint = edgeDict[hedgeIdOfFirstPoint]
            for hedgeIdOfSecondPoint in self.chain[1].incidentEdges:
                hedgeOfSecondPoint = edgeDict[hedgeIdOfSecondPoint]
                if hedgeOfFirstPoint.twin.identifier == hedgeOfSecondPoint.identifier:
                    # calculate outside chain type
                    if hedgeOfFirstPoint.incidentFace.identifier == 'i' or hedgeOfSecondPoint.incidentFace.identifier == 'i':
                        return 2 
                    else:
                        self.incidentFacesOfChain.append(hedgeOfFirstPoint.incidentFace)
                        self.incidentFacesOfChain.append(hedgeOfSecondPoint.incidentFace)
                        return 1


    def calDeg3Type(self):
        deg3Start = self.chain[0]
        deg3End = self.chain[-1]
        if len(deg3Start.incidentEdges) > len(deg3Start.incidentFaces) or len(deg3End.incidentEdges) > len(deg3End.incidentFaces):
            return 2
        else:
            return 1


    # calculate the intersection of locateCircle and apollonisCircle
    def calLocateIntersection(self, centroid):
        x = self.locateCircleCenter[0]   # center of locateCircle
        y = self.locateCircleCenter[1]   
        R = self.locateCircleRadius      # radius of locateCircle
        if self.threeSectionArc is not None and self.threeSectionArc.apollonisCircle is not None:
            a = self.threeSectionArc.apollonisCircle.center_x   # center of apollonisCircle
            b = self.threeSectionArc.apollonisCircle.center_y
            S = self.threeSectionArc.apollonisCircle.radius      # radius of apollonisCircle
            d = math.sqrt((abs(a - x)) ** 2 + (abs(b - y)) ** 2)    
            A = (R ** 2 - S ** 2 + d ** 2) / (2 * d)
            h = math.sqrt(abs(R ** 2 - A ** 2))
            x2 = x + A * (a - x) / d
            y2 = y + A * (b - y) / d
            x3 = round(x2 - h * (b - y) / d, 2)
            y3 = round(y2 + h * (a - x) / d, 2)
            x4 = round(x2 + h * (b - y) / d, 2)
            y4 = round(y2 - h * (a - x) / d, 2)
            c1 = [x3, y3]
            c2 = [x4, y4]

            
        else:
            x3 = x + math.cos(math.pi / 2) * (centroid.x - x) - math.sin(math.pi / 2) * (centroid.y - y)
            y3 = y + math.sin(math.pi / 2) * (centroid.x - x) + math.cos(math.pi / 2) * (centroid.y - y)
            
            x4 = x + math.cos(math.pi / 2) * (centroid.x - x) + math.sin(math.pi / 2) * (centroid.y - y)
            y4 = y - math.sin(math.pi / 2) * (centroid.x - x) + math.cos(math.pi / 2) * (centroid.y - y)

            ratio = R / math.sqrt((centroid.x - x)**2 + (centroid.y - y)**2)
            x3 = (x3 - x) * ratio + x
            y3 = (y3 - y) * ratio + y
            x4 = (x4 - x) * ratio + x
            y4 = (y4 - y) * ratio + y

        c1 = [x3, y3]
        c2 = [x4, y4]
        return c1, c2
    def getSmallerFace(self, faceCentroidDict, centroidRadiusDict):
        centroid0 = faceCentroidDict[self.incidentFacesOfChain[0].identifier]
        radius0 = centroidRadiusDict[centroid0.identifier]
        centroid1 = faceCentroidDict[self.incidentFacesOfChain[1].identifier]
        radius1 = centroidRadiusDict[centroid1.identifier]
        samllerFace = self.incidentFacesOfChain[0] if radius0 < radius1 else self.incidentFacesOfChain[1]
        # biggerFace = self.incidentFacesOfChain[0] if samllerFace.identifier != self.incidentFacesOfChain[0].identifier else self.incidentFacesOfChain[1]
        biggerFace = self.incidentFacesOfChain[0] if radius0 >= radius1 else self.incidentFacesOfChain[1]
        return samllerFace, biggerFace


    def distributeInsideDeg2Chain(self):
        self.threeSectionArc.distributeDeg2Vertices(self.chain)

    def distributeOutsideDeg2Chain(self, centroid, radius):
        self.threeSectionArc = Arc(self.chain[0], self.chain[-1])
        self.threeSectionArc.calTangencyCircleOutside(self.chain, centroid, radius)
        self.threeSectionArc.distributeDeg2VerticesOutside(self.chain, centroid)


