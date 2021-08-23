'''
'''

import math
# import networkx as nx
# import numpy as np

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
        self.centroidChainDict = {}
        self.centroidEdgeDict = self.buildCentroidEdgeDict() # Add relationship between the centroid and the edges between centroids
        self.centroidFaceDict = centroidFaceDict # add relationship between the centroid and the face, from dcel
        self.xRotateDict = {}
        self.yRotateDict = {}

    # If we know two centroids, then we can know the chain between them
    def setCentroidChainDict(self, centroidChainDict):
        self.centroidChainDict = centroidChainDict


    # construct the graph for dijkstra
    # def generateNetwork(self):
    #     G = nx.Graph()

    #     for centroid in self.centroidList:
    #         G.add_node(centroid.identifier)
    #     for edge in self.edges:
    #         start = edge.start.identifier
    #         end = edge.end.identifier
    #         weight = self.centroidRadiusDict[start] + self.centroidRadiusDict[end]
    #         G.add_weighted_edges_from([(start, end, weight)])
    #     return G


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


    
    '''--------------------------------------------------------------------------------------------------
                    part 1:   Rearrange the centroids to the ideal position
    ------------------------------------------------------------------------------------------------------
    Rearrange the centroids, that adjacent circles are tangent to each other. The formula for calculating 
    force is taken from the paper. Forces only exist between adajcent centroids. The ideal distance between 
    two adajcent centroids is adding the two radius up. 
    --------------------------------------------------------------------------------------------------------'''


    '''--------------calculate RepulsiveForce-----------------'''
    # ********For every two centroids***********
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


    # *********For only adjacent centroids*********
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
    # For only adjacent centroids
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


    '''--------------move the centroids------------'''
    # Add the repulsive and attractive forces calculated above and calculate the coordinates of the new position
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
        self.handleRotateRepusive()
        self.updateCoordinates()
                
        print("total energy: ", currentEnergy)
        self.lastTimeEnergy = currentEnergy
        return True

    '''--------------------------------------------------------------------------------------------------
                    part 2:   Rotate the centroids to improve short-chain problem
    ------------------------------------------------------------------------------------------------------
    After rearranging the centroids, there sometimes exists very short chains with many vertices. Thus,
    rotate the centroid according to the number of vertices along an arc. 
    --------------------------------------------------------------------------------------------------------'''

    '''----------------The formula for calculating forces-------------'''
    # rotate force between two centroid lines
    def rotateRepulsiveForce(self, idealRadian, realRadian):
        # realRadian = abs(realRadian)
        return - idealRadian / 1 * realRadian 


    # A dictionary, find all the adjacent edges(orange lines between centroids) for every centroid
    # key: a centroid, value: all incident edges of this centroid
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


    '''------main function---'''
    # main function for ratation, with movement for centroids
    def handleRotateRepusive(self):
        for kCentroid, vEdgeList in self.centroidEdgeDict.items():
            index = 0
            # TODO each pair of neighbour orange edge 
            # while index < 1:
            while index < len(vEdgeList):
                startEdge = vEdgeList[index]
                endEdge = vEdgeList[(index+1) % len(vEdgeList)] # 最后一组是下标n-1 到 0，完成循环
                index += 1

                # Find the two centroids to rotate
                centroid1 = startEdge.start if startEdge.start.identifier != kCentroid.identifier else startEdge.end
                centroid2 = endEdge.start if endEdge.start.identifier != kCentroid.identifier else endEdge.end
                
                # calculate the real radian
                # kCentroid——centroid1 and kCentroid——centroid2
                # And determine which side is more forward counterclockwise and which side is more back
                realRadian, preCentroid, postCentroid = self.calRandianBetween2Vector(kCentroid, centroid1, centroid2)

                # kCentroid is the center of rotation
                chainBetweenKCentroidAndPre = self.getChainBetweenCentroids(kCentroid, preCentroid)
                chainBetweenKCentroidAndPost = self.getChainBetweenCentroids(kCentroid, postCentroid)

                # TODO               
                # only for adajcent faces with common chains, not with only one vertex
                if centroid1.identifier not in self.centroidChainDict or centroid2.identifier not in self.centroidChainDict:
                    continue
                if len(chainBetweenKCentroidAndPre) < 2 or len(chainBetweenKCentroidAndPost) < 2:
                    continue
                
                # calculate the rotation angle
                # We expect that half of the vertices on the chain belong to this arc               
                vertexListOfKCentroid = [v for v in self.centroidFaceDict[kCentroid.identifier].loopOuterVertices()]
                rotateRadian = 0.5*self.rotateRepulsive(realRadian, chainBetweenKCentroidAndPre, chainBetweenKCentroidAndPost, vertexListOfKCentroid)
                angle = math.degrees(rotateRadian)
                #print("angle:" + str(angle))

                # TODO done move centroid1 and centroid2
                # prex = (preCentroid.x-kCentroid.x)*math.cos(angle) - (preCentroid.y-kCentroid.y)*math.sin(angle)+kCentroid.x
                # prey = (preCentroid.y-kCentroid.y)*math.sin(angle) + (preCentroid.x-kCentroid.x)*math.cos(angle)+kCentroid.y
                # postx = (postCentroid.x-kCentroid.x)*math.cos(angle) + (postCentroid.y-kCentroid.y)*math.sin(angle)+kCentroid.x
                # posty = (postCentroid.y-kCentroid.y)*math.cos(angle) - (postCentroid.x-kCentroid.x)*math.sin(angle)+kCentroid.y
                
                # self.xDisDit[preCentroid.identifier] = self.xDisDit[preCentroid.identifier] + (prex - preCentroid.x)
                # self.yDisDit[preCentroid.identifier] = self.yDisDit[preCentroid.identifier] + (prey - preCentroid.y)
                # self.xDisDit[postCentroid.identifier] = self.xDisDit[postCentroid.identifier] + (postx - postCentroid.x)
                # self.yDisDit[postCentroid.identifier] = self.yDisDit[postCentroid.identifier] + (posty - postCentroid.y)
                
                # small step
                factor = 0.001
                self.xDisDit[postCentroid.identifier] = self.xDisDit[postCentroid.identifier] + ((postCentroid.x-kCentroid.x)*math.cos(factor*rotateRadian) + (postCentroid.y-kCentroid.y)*math.sin(factor*rotateRadian)+kCentroid.x - postCentroid.x)
                self.yDisDit[postCentroid.identifier] = self.yDisDit[postCentroid.identifier] + ((postCentroid.y-kCentroid.y)*math.cos(factor*rotateRadian) - (postCentroid.x-kCentroid.x)*math.sin(factor*rotateRadian)+kCentroid.y - postCentroid.y)
                self.xDisDit[preCentroid.identifier] = self.xDisDit[preCentroid.identifier] + ((preCentroid.x-kCentroid.x)*math.cos(-factor*rotateRadian) + (preCentroid.y-kCentroid.y)*math.sin(-factor*rotateRadian)+kCentroid.x - preCentroid.x)
                self.yDisDit[preCentroid.identifier] = self.yDisDit[preCentroid.identifier] + ((preCentroid.y-kCentroid.y)*math.cos(-factor*rotateRadian) - (preCentroid.x-kCentroid.x)*math.sin(-factor*rotateRadian)+kCentroid.y - preCentroid.y)

                


    '''--------------calculate the rotation angle-------------------'''
    def rotateRepulsive(self, realRadian, chainBetweenKCentroidAndPre, chainBetweenKCentroidAndPost, vertexListOfKCentroid):       
        countPre = len(chainBetweenKCentroidAndPre) if chainBetweenKCentroidAndPre is not None else 0
        countPost = len(chainBetweenKCentroidAndPost) if chainBetweenKCentroidAndPost is not None else 0
        count = (countPre + countPost) / 2 - 1
        #print("countPre: " + str(countPre))
        #print("countPost: " + str(countPost))
        #print("count: " + str(count))
        # ideal radian for this arc with this number of vertices
        idealRadian = count / len(vertexListOfKCentroid) * 2 * math.pi
        #print("idealRadian: " + str(idealRadian))
        # The angle to be moved should be half the difference between the actual angle and the ideal angle
        return abs(idealRadian - realRadian) # TODO done. positive or negative sign
        # return self.rotateRepulsiveForce(idealRadian, realRadian)
         

    def getChainBetweenCentroids(self, kCentroid, vCentroid):
        for chainOfKCentroid in self.centroidChainDict[kCentroid.identifier]:
            for chainOfVCentroid in self.centroidChainDict[vCentroid.identifier]:
                if chainOfKCentroid.chainId == chainOfVCentroid.chainId:
                    # if two faces have a common chain
                    return chainOfKCentroid.chain
        # if two faces do not have a common chain (with only one common vertex)
        return []


    # calculate the real radian
    def calRandianBetween2Vector(self, kCentroid, centroid1, centroid2):
        dx1 = centroid1.x - kCentroid.x
        dy1 = centroid1.y - kCentroid.y
        dx2 = centroid2.x - kCentroid.x
        dy2 = centroid2.y - kCentroid.y
        angle1 = math.atan2(dy1, dx1)
        angle1 = int(angle1 * 180/math.pi)
        # print(angle1)
        angle2 = math.atan2(dy2, dx2)
        angle2 = int(angle2 * 180/math.pi)
        # print(angle2)
        included_angle = 0
        # From the x-positive half axis, the centroid(line) at the front 
        # Check which should rotate clockwise and which should counterclockwise
        preCentroid = centroid1 if angle1 > angle2 else centroid2
        postCentroid = centroid1 if angle1 <= angle2 else centroid2
        if angle1*angle2 >= 0:
            included_angle = abs(angle1-angle2)
        else:
            included_angle = abs(angle1) + abs(angle2)
            if included_angle > 180:
                preCentroid = centroid1 if angle1 < angle2 else centroid2
                postCentroid = centroid1 if angle1 >= angle2 else centroid2
                included_angle = 360 - included_angle
        included_angle = math.radians(included_angle)
        return included_angle, preCentroid, postCentroid

    # calculating the angle of the vector to the positive half axis of x
    def calRandianBetweenVectorAndX(self, kCentroid, vertex):
        dx1 = vertex.x - kCentroid.x
        dy1 = vertex.y - kCentroid.y
        angle1 = math.atan2(dy1, dx1)
        angle1 = int(angle1 * 180/math.pi)
        return angle1





    '''--------------------------------------------------------------------------------------------------
                            part 3:   Deal with deg3+ vertices
    ------------------------------------------------------------------------------------------------------
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