import math

class edge(object):
    
    def __init__(self, start, end):
        self.start = start
        self.end = end


class force_directed(object):
    
    def __init__(self, centroidList, centroidRadiusDict, edges):
        self.centroidList = centroidList
        self.centroidRadiusDict = centroidRadiusDict
        self.edges = edges
        self.nodeDict = self.geneNodeDict(centroidList)
        self.xDisDit = {}
        self.yDisDit = {}

    def geneNodeDict(self, centroidList):
        nodeDict = {}
        for vertex in centroidList:
            nodeDict[vertex.identifier] = vertex
    
    def handler(self):
        for i in range(10):
            # self.calRepulsiveForce()
            self.calAttractiveForce()
            self.updateCoordinates()
            print(self.checkTotalEnergy())



    def checkTotalEnergy(self):
        total = 0
        for centroid1 in self.centroidList:
            for centroid2 in self.centroidList:
                if centroid1.identifier == centroid2.identifier:
                    continue
                distX = centroid1.x - centroid2.x
                distY = centroid1.y - centroid2.y
                dist = math.sqrt(distX**2 + distY**2)
                idealDis = self.centroidRadiusDict.get(centroid1.identifier) + self.centroidRadiusDict.get(centroid2.identifier)
                total = total + (abs(dist) - idealDis) * (abs(dist) - idealDis)
        return total


    def calRepulsiveForce(self):
            for edge in self.edges:
                centroid1 = edge.start
                centroid2 = edge.end
                distX = distY = dist = 0.0
                self.xDisDit[centroid1.identifier] = 0.0
                self.yDisDit[centroid1.identifier] = 0.0
        # distX = distY = dist = 0.0
        # for centroid1 in self.centroidList:
        #     self.xDisDit[centroid1.identifier] = 0.0
        #     self.yDisDit[centroid1.identifier] = 0.0
        #     for centroid2 in self.centroidList:
        #         if centroid1.identifier == centroid2.identifier:
        #             continue
                distX = centroid1.x - centroid2.x
                distY = centroid1.y - centroid2.y
                dist = math.sqrt(distX**2 + distY**2)
                idealDis = self.centroidRadiusDict.get(centroid1.identifier) + self.centroidRadiusDict.get(centroid2.identifier)
                self.xDisDit[centroid1.identifier] = self.xDisDit[centroid1.identifier] + distX / dist * self.repulsiveForce(idealDis, dist) * dist / abs(dist)
                self.yDisDit[centroid1.identifier] = self.yDisDit[centroid1.identifier] + distY / dist * self.repulsiveForce(idealDis, dist) * dist / abs(dist)
            

    def repulsiveForce(self, idealDis, dist):
        print("real dis:" + str(dist) + "   ideal dis：" + str(idealDis))
        dist = abs(dist)
        return - idealDis / dist

    def attractiveForce(self, idealDis, dist):
        print("real dis:" + str(dist) + "   ideal dis：" + str(idealDis))
        dist = abs(dist)
        return dist / idealDis
        

    def calAttractiveForce(self):
        for edge in self.edges:
            startcentroid = edge.start
            endcentroid = edge.end
            distX = distY = dist = 0.0
            distX = startcentroid.x - endcentroid.x
            distY = startcentroid.y - endcentroid.y
            dist = math.sqrt(distX**2 + distY**2)

            idealDis = self.centroidRadiusDict.get(startcentroid.identifier) + self.centroidRadiusDict.get(endcentroid.identifier)

            self.xDisDit[startcentroid.identifier] = self.xDisDit[startcentroid.identifier] - distX / dist * self.repulsiveForce(idealDis, dist) * dist / abs(dist)
            self.yDisDit[startcentroid.identifier] = self.yDisDit[startcentroid.identifier] - distY / dist * self.repulsiveForce(idealDis, dist) * dist / abs(dist)
            self.xDisDit[endcentroid.identifier] = self.xDisDit[endcentroid.identifier] + distX / dist * self.repulsiveForce(idealDis, dist) * dist / abs(dist)
            self.yDisDit[endcentroid.identifier] = self.yDisDit[endcentroid.identifier] + distX / dist * self.repulsiveForce(idealDis, dist) * dist / abs(dist)


    def updateCoordinates(self):
        for centroid in self.centroidList:
            disp_x = self.xDisDit[centroid.identifier]
            disp_y = self.yDisDit[centroid.identifier]
            centroid.x = centroid.x + disp_x
            centroid.y = centroid.y + disp_y
