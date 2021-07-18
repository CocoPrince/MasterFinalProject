import math
import networkx as nx
import numpy as np

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
        self.networkGraph = self.generateNetwork()

    # construct the graph for dijkstra
    def generateNetwork(self):
        G = nx.Graph()
        
        for centroid in self.centroidList:
            G.add_node(centroid.identifier)
        for edge in self.edges:
            start = edge.start.identifier
            end = edge.end.identifier
            weight = self.centroidRadiusDict[start] + self.centroidRadiusDict[end]
            G.add_weighted_edges_from([(start, end, weight)])
        return G

    def geneNodeDict(self, centroidList):
        nodeDict = {}
        for vertex in centroidList:
            nodeDict[vertex.identifier] = vertex
    

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

    def calIdealDis(self, startcentroid, endcentroid):
        # idealDis = self.centroidRadiusDict.get(startcentroid.identifier) + self.centroidRadiusDict.get(endcentroid.identifier)
        idealDis = nx.dijkstra_path_length(self.networkGraph, source=startcentroid.identifier, target=endcentroid.identifier)
        return idealDis

    def calRepulsiveForce(self):
        distX = distY = dist = 0.0
        for centroid1 in self.centroidList:
            self.xDisDit[centroid1.identifier] = 0.0
            self.yDisDit[centroid1.identifier] = 0.0
            for centroid2 in self.centroidList:
                if centroid1.identifier == centroid2.identifier:
                    continue
                distX = centroid2.x - centroid1.x
                distY = centroid2.y - centroid1.y
                dist = math.sqrt(distX**2 + distY**2)
                idealDis = self.calIdealDis(centroid1, centroid2)
                self.xDisDit[centroid1.identifier] = self.xDisDit[centroid1.identifier] + distX / dist * self.repulsiveForce(idealDis, dist) * distX / abs(distX)
                self.yDisDit[centroid1.identifier] = self.yDisDit[centroid1.identifier] + distY / dist * self.repulsiveForce(idealDis, dist) * distY / abs(distY)

    def repulsiveForce(self, idealDis, dist):
        # print("real dis:" + str(dist) + "   ideal dis：" + str(idealDis))
        dist = abs(dist)
        return - idealDis / dist 

    def attractiveForce(self, idealDis, dist):
        # print("real dis:" + str(dist) + "   ideal dis：" + str(idealDis))
        dist = abs(dist)
        # print("坸引力："+str(dist / idealDis))
        return dist / idealDis
        

    # def calRepulsiveForce(self):
    #     for edge in self.edges:
    #         startcentroid = edge.start
    #         endcentroid = edge.end
    #         self.xDisDit[startcentroid.identifier] = 0.0
    #         self.yDisDit[startcentroid.identifier] = 0.0
    #         self.xDisDit[endcentroid.identifier] = 0.0
    #         self.yDisDit[endcentroid.identifier] = 0.0
    #         distX = distY = dist = 0.0
    #         distX = startcentroid.x - endcentroid.x
    #         distY = startcentroid.y - endcentroid.y
    #         dist = math.sqrt(distX**2 + distY**2)

    #         idealDis = self.centroidRadiusDict.get(startcentroid.identifier) + self.centroidRadiusDict.get(endcentroid.identifier)

    #         self.xDisDit[startcentroid.identifier] = self.xDisDit[startcentroid.identifier] - distX / dist * self.repulsiveForce(idealDis, dist) * dist / abs(dist)
    #         self.yDisDit[startcentroid.identifier] = self.yDisDit[startcentroid.identifier] - distY / dist * self.repulsiveForce(idealDis, dist) * dist / abs(dist)
    #         self.xDisDit[endcentroid.identifier] = self.xDisDit[endcentroid.identifier] + distX / dist * self.repulsiveForce(idealDis, dist) * dist / abs(dist)
    #         self.yDisDit[endcentroid.identifier] = self.yDisDit[endcentroid.identifier] + distX / dist * self.attractiveForce(idealDis, dist) * dist / abs(dist)

        
    def calAttractiveForce(self):
        for edge in self.edges:
            startcentroid = edge.start
            endcentroid = edge.end
            # self.xDisDit[startcentroid.identifier] = 0.0
            # self.yDisDit[startcentroid.identifier] = 0.0
            # self.xDisDit[endcentroid.identifier] = 0.0
            # self.yDisDit[endcentroid.identifier] = 0.0
            distX = distY = dist = 0.0
            distX = startcentroid.x - endcentroid.x
            distY = startcentroid.y - endcentroid.y
            dist = math.sqrt(distX**2 + distY**2)

            idealDis = self.calIdealDis(startcentroid, endcentroid)

            self.xDisDit[startcentroid.identifier] = self.xDisDit[startcentroid.identifier] - distX / dist * self.attractiveForce(idealDis, dist) * distX / abs(distX)
            self.yDisDit[startcentroid.identifier] = self.yDisDit[startcentroid.identifier] - distY / dist * self.attractiveForce(idealDis, dist) * distY / abs(distY)
            self.xDisDit[endcentroid.identifier] = self.xDisDit[endcentroid.identifier] + distX / dist * self.attractiveForce(idealDis, dist) * distX / abs(distX)
            self.yDisDit[endcentroid.identifier] = self.yDisDit[endcentroid.identifier] + distY / dist * self.attractiveForce(idealDis, dist) * distY / abs(distY)


    def updateCoordinates(self):
        for centroid in self.centroidList:
            disp_x = self.xDisDit[centroid.identifier]
            disp_y = self.yDisDit[centroid.identifier]
            centroid.x = centroid.x + disp_x
            centroid.y = centroid.y + disp_y

    # Main function
    def handler(self):
        self.calRepulsiveForce()
        self.calAttractiveForce()
        self.updateCoordinates()
        print("total energy: " + str(self.checkTotalEnergy()))