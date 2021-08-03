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
    
    def __init__(self, centroidList, centroidRadiusDict, edges, centroidFaceDict):
        self.centroidList = centroidList
        self.centroidRadiusDict = centroidRadiusDict
        self.edges = edges
        self.nodeDict = self.geneNodeDict(centroidList)
        self.xDisDit = {}
        self.yDisDit = {}
        self.networkGraph = self.generateNetwork()
        self.lastTimeEnergy = 0
        self.centroidEdgeDict = self.buildCentroidEdgeDict() # Add relationship between the centroid and the edge
        self.centroidFaceDict = centroidFaceDict # add relationship between the centroid and the face, from dcel


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
    
    '''--------------The whole system stops calculating when the energy is minimal-----------'''
    # calculate the degree of unbalance of the whole map
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
                # total = total + (abs(dist) - idealDis)
        return total


    # ideal distance between two centroids: just add the two radius up
    def calIdealDis(self, startcentroid, endcentroid):
        idealDis = self.centroidRadiusDict.get(startcentroid.identifier) + self.centroidRadiusDict.get(endcentroid.identifier)
        # idealDis = nx.dijkstra_path_length(self.networkGraph, source=startcentroid.identifier, target=endcentroid.identifier)
        return idealDis


    '''----------------The formula for calculating forces-------------'''
    # repulsive force between two centroids: linear
    def repulsiveForce(self, idealDis, dist):
        dist = abs(dist)
        return - idealDis / dist 

    # attractive force between two centroids: linear
    def attractiveForce(self, idealDis, dist):
        dist = abs(dist)
        return dist / idealDis





    '''--------------calculate RepulsiveForce-----------------'''
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

    
    '''--------------------calculate AttractiveForce------------------'''
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

    '''--------------move the '''
    # Add the repulsive and attractive forces calculated above and calculate the coordinates of the new position
    def updateCoordinates(self):
        for centroid in self.centroidList:
            disp_x = self.xDisDit[centroid.identifier]
            disp_y = self.yDisDit[centroid.identifier]
            centroid.x = centroid.x + disp_x
            centroid.y = centroid.y + disp_y



    def buildCentroidEdgeDict(self):
        dict = {}
        for edge in self.edges:
            if edge.start not in dict:
                dict[edge.start] = []
            dict[edge.start].append(edge)
            if edge.end not in dict:
                dict[edge.end] = []
            dict[edge.end].append(edge)
        return dict

    # 为了方便，这里直接移动了
    def handleRotateRepusive(self):
        for kCentroid, vEdgeList in self.centroidEdgeDict.items():
            index = 0
            # TODO 暂时两两之间都计算切向排斥力力 
            while index < len(vEdgeList):
                startEdge = vEdgeList[index]
                endEdge = vEdgeList[(index+1) % len(vEdgeList)] # 最后一组是下标n-1 到 0，完成循环

                # 找到要旋转的两个质心
                centroid1 = startEdge.start if startEdge.start.identifier != kCentroid.identifier else startEdge.end
                centroid2 = endEdge.start if endEdge.start.identifier != kCentroid.identifier else endEdge.end
                # TODO 这里要算出一个需要旋转的角度
                vertexListOfKCentroid = self.centroidFaceDict[kCentroid.identifier].vertexList
                radian = self.rotateRepulsive(kCentroid, centroid1, centroid2, vertexListOfKCentroid)
                angle = math.radians(radian)
                

                # TODO 移动centroid1， centroid2，一个顺时针，一个逆时针,具体怎么判断还没写

                centroid1.x = (centroid1.x-kCentroid.x)*math.cos(angle) - (centroid1.y-kCentroid.y)*math.sin(angle)+kCentroid.x
                centroid1.y = (centroid1.y-kCentroid.y)*math.cos(angle) + (centroid1.x-kCentroid.x)*math.sin(angle)+kCentroid.y

                centroid2.x = (centroid2.x-kCentroid.x)*math.cos(angle) - (centroid2.y-kCentroid.y)*math.sin(angle)+kCentroid.x
                centroid2.y = (centroid2.y-kCentroid.y)*math.cos(angle) + (centroid2.x-kCentroid.x)*math.sin(angle)+kCentroid.y
            #                 


    

    def rotateRepulsive(self, kCentroid, centroid1, centroid2, vertexListOfKCentroid):
        # 遍历KCentroid所属m面的所有vertex 找到位于夹角内的点的数量
        count = 0
        for vertex in vertexListOfKCentroid:
            if 1==1:
                count += 1

        dealRadian = count / len(vertexListOfKCentroid) * 360

        # TODO 求实际夹角 kCentroid——centroid1与kCentroid——centroid2
        realRadian = 0

        # 要移动的角度应该是实际夹角与理想夹角之差的一半，因为两个centroid都要移动，或者只移动一个，那就不用除2了
        return realRadian - dealRadian # TODO 这里不知道需不需要带正负号





    # Main function
    def handler(self):
        currentEnergy = self.checkTotalEnergy()
        if currentEnergy > self.lastTimeEnergy and 0 != self.lastTimeEnergy:
            return False
        self.calRepulsiveForce()
        self.calAttractiveForce()
        self.updateCoordinates()
        self.handleRotateRepusive()
        print("total energy: ", currentEnergy)
        self.lastTimeEnergy = currentEnergy
        return True



    '''--------------------------------------------------------------------------------------------------
                part 1
    ################################## Deal with deg3+ vertices #####################################
    For every inside deg3+ vertex, first move it to the centroid of the polygon of adajcent centroids,
    then use force-directed method to determine its position.
    --------------------------------------------------------------------------------------------------------'''

    # all the centroids of deg3+'s incidentFace will have repulsive force to it
    def cal3DegRepulsiveForce(self, threeDegVertex, centroidsOfIncidentFace):
        xDisDit = 0.0
        yDisDit = 0.0
        for centroid in centroidsOfIncidentFace:
            # D := v.pos - u.pos;
            distX = centroid.x - threeDegVertex.x
            distY = centroid.y - threeDegVertex.y
            dist = math.sqrt(distX**2 + distY**2) # real distance
            idealDis = self.centroidRadiusDict.get(centroid.identifier) # ideal distance
            # v.disp := v.disp + ( D /| D |) * fr (| D |)
            xDisDit = xDisDit + distX / abs(distX) * self.repulsiveForce(idealDis, dist) * abs(distX) / abs(dist)
            yDisDit = yDisDit + distY / abs(distY) * self.repulsiveForce(idealDis, dist) * abs(distY) / abs(dist)      
        
        return xDisDit, yDisDit


    # all the centroids of deg3+'s incidentFace will have attractive force to it
    def cal3DegAttractiveForce(self, threeDegVertex, centroidsOfIncidentFace):
        xDisDit = 0.0
        yDisDit = 0.0
        for centroid in centroidsOfIncidentFace:
            # D := v.pos - u.pos;
            distX = centroid.x - threeDegVertex.x
            distY = centroid.y - threeDegVertex.y
            dist = math.sqrt(distX**2 + distY**2)
            idealDis = self.centroidRadiusDict.get(centroid.identifier)
            # v.disp := v.disp + ( D /| D |) * fr (| D |)
            xDisDit = xDisDit + distX / abs(distX) * self.attractiveForce(idealDis, dist) * abs(distX) / abs(dist)
            yDisDit = yDisDit + distY / abs(distY) * self.attractiveForce(idealDis, dist) * abs(distY) / abs(dist)

        return xDisDit, yDisDit


    # Calculate the final force and move the deg3+ points
    def handle3DegVertex_inside(self, threeDegVertex, centroidsOfIncidentFace):
        flag = True
        total = 0.0
        while flag:
            xDisRepulsive, yDisRepulsive = self.cal3DegRepulsiveForce(threeDegVertex, centroidsOfIncidentFace)
            xDisAttractive, yDisAttractive = self.cal3DegAttractiveForce(threeDegVertex, centroidsOfIncidentFace)
            # Using the repulsive and attractive force calculated above, move of deg3+ points
            threeDegVertex.x = threeDegVertex.x + xDisAttractive + xDisRepulsive
            threeDegVertex.y = threeDegVertex.y + yDisAttractive + yDisRepulsive
            # store last time total energy
            last_time_total_energy = total
            total = 0
            # calculate current total energy
            for centroid in centroidsOfIncidentFace:
                distX = centroid.x - threeDegVertex.x
                distY = centroid.y - threeDegVertex.y
                dist = math.sqrt(distX**2 + distY**2)
                idealDis = self.centroidRadiusDict.get(centroid.identifier)
                total = total + (abs(dist) - idealDis) * (abs(dist) - idealDis)
                # total = total + (abs(dist) - idealDis)
                print(total)                               
            # check if the total energy is the minimum
            if last_time_total_energy != 0.0 and last_time_total_energy <= total:
                print("\n")
                flag = False
                        