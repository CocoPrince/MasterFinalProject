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
        for i in range(5):
            self.calRepulsiveForce()
            self.calAttractiveForce()
            self.updateCoordinates()




    def calRepulsiveForce(self):
        distX = distY = dist = 0.0
        for centroid1 in self.centroidList:
            self.xDisDit[centroid1.identifier] = 0.0
            self.yDisDit[centroid1.identifier] = 0.0
            for centroid2 in self.centroidList:
                if centroid1.identifier == centroid2.identifier:
                    continue
                distX = centroid1.x - centroid2.x
                distY = centroid1.y - centroid2.y
                dist = math.sqrt(distX**2 + distY**2)
                idealDis = self.centroidRadiusDict.get(centroid1.identifier) + self.centroidRadiusDict.get(centroid2.identifier)
                self.xDisDit[centroid1.identifier] = self.xDisDit[centroid1.identifier] + distX / dist * idealDis * idealDis / dist
                self.yDisDit[centroid1.identifier] = self.yDisDit[centroid1.identifier] + distY / dist * idealDis * idealDis / dist


    def calAttractiveForce(self):
        for edge in self.edges:
            startcentroid = edge.start
            endcentroid = edge.end
            distX = distY = dist = 0.0
            distX = startcentroid.x - endcentroid.x
            distY = startcentroid.y - endcentroid.y
            dist = math.sqrt(distX**2 + distY**2)

            idealDis = self.centroidRadiusDict.get(startcentroid.identifier) + self.centroidRadiusDict.get(endcentroid.identifier)

            self.xDisDit[startcentroid.identifier] = self.xDisDit[startcentroid.identifier] - distX * dist / idealDis
            self.yDisDit[startcentroid.identifier] = self.yDisDit[startcentroid.identifier] - distY * dist / idealDis
            self.xDisDit[endcentroid.identifier] = self.xDisDit[endcentroid.identifier] + distX * dist / idealDis
            self.yDisDit[endcentroid.identifier] = self.yDisDit[endcentroid.identifier] + distY * dist / idealDis


    def updateCoordinates(self):
        for centroid in self.centroidList:
            disp_x = self.xDisDit[centroid.identifier]
            disp_y = self.yDisDit[centroid.identifier]
            centroid.x = centroid.x + disp_x
            centroid.y = centroid.y + disp_y









			node.setXPosition((node.getX() + dx) >= CANVAS_WIDTH || (node.getX() + dx) <= 0 ? node.getX() - dx : node.getX() + dx);
			node.setYPosition((node.getY() + dy) >= CANVAS_HEIGHT || (node.getY() + dy <= 0) ? node.getY() - dy : node.getY() + dy);
		}



