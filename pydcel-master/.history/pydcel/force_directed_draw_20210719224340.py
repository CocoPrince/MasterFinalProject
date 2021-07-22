'''

'''

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
        self.lastTimeEnergy = 0

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
        idealDis = self.centroidRadiusDict.get(startcentroid.identifier) + self.centroidRadiusDict.get(endcentroid.identifier)
        # idealDis = nx.dijkstra_path_length(self.networkGraph, source=startcentroid.identifier, target=endcentroid.identifier)
        return idealDis

    # ********For every two vertices***********
    # def calRepulsiveForce(self):
    #     distX = distY = dist = 0.0
    #     for centroid1 in self.centroidList:
    #         self.xDisDit[centroid1.identifier] = 0.0
    #         self.yDisDit[centroid1.identifier] = 0.0
    #         for centroid2 in self.centroidList:
    #             if centroid1.identifier == centroid2.identifier:
    #                 continue
    #             # D := v.pos - u.pos;
    #             distX = centroid2.x - centroid1.x
    #             distY = centroid2.y - centroid1.y
    #             dist = math.sqrt(distX**2 + distY**2)
    #             idealDis = self.calIdealDis(centroid1, centroid2)
    #             # v.disp := v.disp + ( D /| D |) * fr (| D |)
    #             self.xDisDit[centroid1.identifier] = self.xDisDit[centroid1.identifier] + distX / abs(distX) * self.repulsiveForce(idealDis, dist) * abs(distX) / abs(dist)
    #             self.yDisDit[centroid1.identifier] = self.yDisDit[centroid1.identifier] + distY / abs(distY) * self.repulsiveForce(idealDis, dist) * abs(distY) / abs(dist)


    # *********For only adjacent vertices*********
    def calRepulsiveForce(self):
        for edge in self.edges:
            startcentroid = edge.start
            endcentroid = edge.end
            self.xDisDit[startcentroid.identifier] = 0.0
            self.yDisDit[startcentroid.identifier] = 0.0
            self.xDisDit[endcentroid.identifier] = 0.0
            self.yDisDit[endcentroid.identifier] = 0.0
            
        for edge in self.edges:
            startcentroid = edge.start
            endcentroid = edge.end
            distX = distY = dist = 0.0
            distX = startcentroid.x - endcentroid.x
            distY = startcentroid.y - endcentroid.y
            dist = math.sqrt(distX**2 + distY**2)

            idealDis = self.centroidRadiusDict.get(startcentroid.identifier) + self.centroidRadiusDict.get(endcentroid.identifier)

            self.xDisDit[startcentroid.identifier] = self.xDisDit[startcentroid.identifier] - distX / abs(distX) * self.repulsiveForce(idealDis, dist) * abs(distX) / abs(dist)
            self.yDisDit[startcentroid.identifier] = self.yDisDit[startcentroid.identifier] - distY / abs(distY) * self.repulsiveForce(idealDis, dist) * abs(distY) / abs(dist)
            self.xDisDit[endcentroid.identifier] = self.xDisDit[endcentroid.identifier] + distX / abs(distX) * self.repulsiveForce(idealDis, dist) * abs(distX) / abs(dist)
            self.yDisDit[endcentroid.identifier] = self.yDisDit[endcentroid.identifier] + distY / abs(distY) * self.repulsiveForce(idealDis, dist) * abs(distY) / abs(dist)

    
    # For only adjacent vertices
    def calAttractiveForce(self):
        for edge in self.edges:
            startcentroid = edge.start
            endcentroid = edge.end
            distX = distY = dist = 0.0
            distX = startcentroid.x - endcentroid.x
            distY = startcentroid.y - endcentroid.y
            dist = math.sqrt(distX**2 + distY**2)

            idealDis = self.calIdealDis(startcentroid, endcentroid)

            self.xDisDit[startcentroid.identifier] = self.xDisDit[startcentroid.identifier] - distX / abs(distX) * self.attractiveForce(idealDis, dist) * abs(distX) / abs(dist)
            self.yDisDit[startcentroid.identifier] = self.yDisDit[startcentroid.identifier] - distY / abs(distY) * self.attractiveForce(idealDis, dist) * abs(distY) / abs(dist)
            self.xDisDit[endcentroid.identifier] = self.xDisDit[endcentroid.identifier] + distX / abs(distX) * self.attractiveForce(idealDis, dist) * abs(distX) / abs(dist)
            self.yDisDit[endcentroid.identifier] = self.yDisDit[endcentroid.identifier] + distY / abs(distY) * self.attractiveForce(idealDis, dist) * abs(distY) / abs(dist)


    def repulsiveForce(self, idealDis, dist):
        dist = abs(dist)
        return - idealDis / dist 

    def attractiveForce(self, idealDis, dist):
        dist = abs(dist)
        return dist / idealDis

    def updateCoordinates(self):
        for centroid in self.centroidList:
            disp_x = self.xDisDit[centroid.identifier]
            disp_y = self.yDisDit[centroid.identifier]
            centroid.x = centroid.x + disp_x
            centroid.y = centroid.y + disp_y

    # Main function
    def handler(self):
        currentEnergy = self.checkTotalEnergy()
        if currentEnergy > self.lastTimeEnergy and 0 != self.lastTimeEnergy:
            return False
        self.calRepulsiveForce()
        self.calAttractiveForce()
        self.updateCoordinates()
        print("total energy: ", currentEnergy)
        self.lastTimeEnergy = currentEnergy
        return True


    # ********For every two vertices***********
    def cal3DimRepulsiveForce(self, threeDimVertex, centroidList):
        distX = distY = dist = 0.0
        xDisDit = 0.0
        yDisDit = 0.0
        for centroid in centroidList:
            if threeDimVertex.identifier == centroid.identifier:
                continue
            # D := v.pos - u.pos;
            distX = centroid.x - threeDimVertex.x
            distY = centroid.y - threeDimVertex.y
            dist = math.sqrt(distX**2 + distY**2)
            idealDis = self.centroidRadiusDict.get(centroid.identifier)
            # v.disp := v.disp + ( D /| D |) * fr (| D |)
            xDisDit = xDisDit + distX / abs(distX) * self.repulsiveForce(idealDis, dist) * abs(distX) / abs(dist)
            yDisDit = yDisDit + distY / abs(distY) * self.repulsiveForce(idealDis, dist) * abs(distY) / abs(dist)      
        
        return xDisDit, yDisDit

    
    def cal3DimAttractiveForce(self, threeDimVertex, centroidList):
        distX = distY = dist = 0.0
        xDisDit = 0.0
        yDisDit = 0.0
        for centroid in centroidList:
            if threeDimVertex.identifier == centroid.identifier:
                continue
            # D := v.pos - u.pos;
            distX = centroid.x - threeDimVertex.x
            distY = centroid.y - threeDimVertex.y
            dist = math.sqrt(distX**2 + distY**2)
            idealDis = self.centroidRadiusDict.get(centroid.identifier)
            # v.disp := v.disp + ( D /| D |) * fr (| D |)
            xDisDit = xDisDit + distX / abs(distX) * self.attractiveForce(idealDis, dist) * abs(distX) / abs(dist)
            yDisDit = yDisDit + distY / abs(distY) * self.attractiveForce(idealDis, dist) * abs(distY) / abs(dist)

        return xDisDit, yDisDit

    def handle3DimVertex(self, threeDimVertex, centroidList):
        flag = True
        total = 0.0
        while flag:
            xDisRepulsive, yDisRepulsive = self.cal3DimRepulsiveForce(self, threeDimVertex, centroidList)
            xDisAttractive, yDisAttractive = self.cal3DimAttractiveForce(self, threeDimVertex, centroidList)
            for centroid in centroidList:
                threeDimVertex.x = threeDimVertex.x + xDisAttractive + xDisRepulsive
                threeDimVertex.y = threeDimVertex.y + yDisAttractive + yDisRepulsive
                distX = centroid.x - threeDimVertex.x
                distY = centroid.y - threeDimVertex.y
                dist = math.sqrt(distX**2 + distY**2)
                idealDis = self.centroidRadiusDict.get(centroid.identifier)
                tmp_total = total
                total = total + (abs(dist) - idealDis) * (abs(dist) - idealDis)
                print(total)
                    
                if tmp_total != 0.0 and tmp_total <= total:
                    print("\n")
                    flag = False
                        

