'''
Construct DCEL and some related calculations
'''

import math
import pydcel
from .force_directed_draw import *
import queue
from .chain import *
from .arc import *

'''---------------------------------------------------------------------
data structure: DCEL
---------------------------------------------------------------------'''
class vertex(object):   
    def __init__(self, px, py, pz, identifier):
        self.identifier = identifier  # id ofthe vertex
        self.x = px    # coordinates of the vertex
        self.y = py
        self.z = pz
        self.incidentEdge = None    # one arbitrary edge that is incident to the vertex
        self.incidentFaces = set()   # all faces that are incident to the vertex
        self.incidentEdges = set()   # all edges that are incident to the vertex
        
    def setTopology(self, newIncedentEdge):
        self.incidentEdge = newIncedentEdge
        
    def p(self):
        return (self.x,self.y,self.z)

    def __repr__(self):
        return "v{} ({}, {}, {})".format(self.identifier, self.x, self.y, self.z)


class hedge(object):    
    def __init__(self, identifier):
        self.identifier = identifier   # id ofthe hedge
        self.origin = None    # The origin vertex of the hedge
        self.twin = None    # twin of the hedge
        self.incidentFace = None    # the face to its left
        self.next = None
        self.previous = None

    def setTopology(self, newOrigin, newTwin, newIncindentFace, newNext, newPrevious):
        self.origin = newOrigin
        self.twin = newTwin
        self.incidentFace = newIncindentFace
        self.next = newNext
        self.previous = newPrevious
        
    def loop(self):
        """ Tranvers face f:
        Loop from this hedge to the next ones. Stops when we are at the current one again."""
        yield self
        e = self.next
        while e is not self:
            yield e
            e = e.next
            
    def wind(self):
        """ Tranvers the hedges of a face:
        iterate over hedges emerging from vertex at origin in ccw order"""
        yield self
        e = self.previous.twin
        while e is not self:
            yield e
            e = e.previous.twin

    def __repr__(self):
        return "he{}".format(self.identifier)


class face(object):    
    def __init__(self, identifier):
        self.identifier = identifier   # id ofthe face
        self.outerHedges = None      # all outer hedges of the face 
        self.innerHedges = None      # all inner hedges of the face

    def setTopology(self, newOuterHedges, newInnerHedges=None):
        self.outerHedges = newOuterHedges
        self.innerHedges = newInnerHedges
        
    def loopOuterVertices(self):
        '''Tranvers the vertices of a face'''
        for e in self.outerHedges.loop():
            yield e.origin

    def loopOuterEdges(self):
        '''Tranvers the hedges of a face'''
        for e in self.outerHedges.loop():
            yield e

    def __repr__(self):
        return "f{}".format(self.identifier)


class DCEL(object):    
    def __init__(self):
        self.vertexList = []    # all the vertices of the graph
        self.hedgeList = []     # all the hedges of the graph
        self.faceList = []      # all the hedges of the graph
        self.centroidList = []  # all the centroids of the graph
        self.vertexOnCircleList = []   # equal division points on the circle
        self.faceVerticesDict = {}  # Dictionary: the correspondence between the face and its vertices 
        self.faceCentroidDict = {}    # Dictionary: the correspondence between the face and its centroid 
        self.centroidRadiusDict = {}  # Dictionary: the correspondence between the centroid and the radius
        self.centroidFaceDict = {}
        self.radiusList = []    # all the radius of the circle
        self.infiniteFace = None
        self.centroidEdges = [] # all the edges of the centroid graph
        self.apolloCenterList = []
        self.apolloRadiusList = []
        self.face_dict = {} # dictionary that finds faces by id
        self.edge_dict = {} # dictionary that finds hedges by id

    def getNewId(self, L):
        if len(L) == 0:
            return 0
        else:
            return L[-1].identifier + 1
        
    def createVertex(self, px, py, pz):
        '''Used to create all the vertices of the graph (vertexList)'''
        identifier = self.getNewId(self.vertexList)
        v = vertex(px,py,pz, identifier)
        self.vertexList.append(v)
        return v
    
    def createCentroid(self, px, py, pz):
        '''Used to create centroid for each face (centroidList)'''
        identifier = self.getNewId(self.centroidList)
        v = vertex(px,py,pz, identifier)
        self.centroidList.append(v)
        return v

    def createVertexPure(self, px, py, pz):
        '''Used to create auxiliary points, such as equal division points on the arc'''
        identifier = self.getNewId(self.vertexList)
        v = vertex(px,py,pz, identifier)
        return v
        
    def createHedge(self):
        '''Used to create Hedges of the graph (hedgeList)'''
        identifier = self.getNewId(self.hedgeList)
        e = hedge(identifier)
        self.hedgeList.append(e)
        return e
        
    def createFace(self):
        '''Used to create Faces of the graph (faceList)'''
        identifier = self.getNewId(self.faceList)
        f = face(identifier)
        self.faceList.append(f)
        return f

    def createInfFace(self):
        '''Used to create the infinite face'''
        identifier = "i"
        f = face(identifier)
        self.infiniteFace = f
        return f
        
    def remove(self, element):
        # TODO: verify topology? i.e. make sure no reference to element exists...
        if type(element) is vertex:
            self.vertexList.remove(element)
            del element
        elif type(element) is hedge:
            self.hedgeList.remove(element)
            del element
        elif type(element) is face:
            self.faceList.remove(element)
            del element
        else:
            print("illegal element type")

    def __repr__(self):
        s = ""
        for v in self.vertexList:
            s += "{} : \t{}\n".format(v, v.incidentEdge)
        for e in self.hedgeList:
            s += "{} : \t{} \t{} \t{} \t{} \t{}\n".format(e, e.origin, e.twin, e.incidentFace, e.next, e.previous)
        for f in self.faceList + [self.infiniteFace]:
            s += "{} : \t{} \t{}\n".format(f, f.outerHedges, f.innerHedges)
        return s

    def remove_vertex(self, vertex):
        # remove this vertex from the DCEL and keep the topo right. Assumes no dangling edges.
        e_0 = vertex.incidentEdge       
        # we don't want to come near the infiniteFace
        for e in e_0.wind():
            print((e, e.twin, e.incidentFace))
            if e.incidentFace == self.infiniteFace:
                print("refusing to remove vertex incident to infiniteFace...")
                return                
        # we also don't want to create any dangling edges
        for e in e_0.wind():
            if e.previous == e.twin.next.twin:
                print("refusing to remove this vertex because it will create dangling edge(s)")
                return
            for e_neighbor in e.next.wind():
                if e_neighbor.previous == e_neighbor.twin.next.twin:
                    print("refusing to remove this vertex because it might cause dangling edge(s) in future")
                    return       
        #This face we like so we keep it.
        nice_face = e_0.incidentFace       
        toRemove = [vertex]
        current_edge = e_0.twin
        if current_edge.incidentFace != nice_face:
            toRemove.append(current_edge.incidentFace)
        # update all face references to nice face
        while current_edge != e_0.previous:
            # loop backwards over face that must be removed and set incidentface fields to nice face
            while current_edge.origin != vertex:
                current_edge = current_edge.previous
                # if current_edge == None: return
                current_edge.incidentFace = nice_face
            current_edge = current_edge.twin
            # this face must be gone
            if current_edge.incidentFace != nice_face:
                toRemove.append(current_edge.incidentFace)
        # update prev and next fields
        edges = [e for e in e_0.wind()]
        for e in edges:
            e.next.previous = e.twin.previous
            e.twin.previous.next = e.__next__
            e.twin.origin.incidentEdge = e.__next__
            toRemove.append(e)
            toRemove.append(e.twin)
        # now we can finally get rid of this stuff
        nice_face.outerHedges = e_0.__next__
        for element in toRemove:
            self.remove(element)



    '''============================Now, let's start!!!====================================='''

    '''------------------------------------------------------------------------------------
    calculate the centroid of a face
    ---------------------------------------------------------------------------------------'''
    # calculate the area of each triangle
    def calArea(self, p1, p2):
        # p0 is the origin, so in order to simplify the calculation, we use only the other two vertices
        return (p1.x*p2.y - p2.x*p1.y) / 2.0
          
    # calculate the centroid of a face (polygon)
    def calCentroid(self, face_vertices):
        x = 0
        y = 0
        totalArea = 0
        eachArea = 0
        n = len(face_vertices)
        for i in range(n):
            eachArea = self.calArea(face_vertices[i], face_vertices[(i + 1) % n]) # A triangle formed by two vertices and the origin
            x += (face_vertices[i].x + face_vertices[(i + 1) % n].x)*eachArea # The centroid divided by 3 for each step is dealt with at the end
            y += (face_vertices[i].y + face_vertices[(i + 1) % n].y)*eachArea
            totalArea += eachArea # total area of the face
        x = x / (totalArea * 3) #the coordinate for the centroid of the face
        y = y / (totalArea * 3)
        centroid = self.createCentroid(round(x, 3), round(y, 3), 0)
        return centroid, abs(totalArea)


    '''------------------------------------------------------------------------------------
    calculate the distance
    ---------------------------------------------------------------------------------------'''
    # calculate the edge length (distance between two vertices)
    def calDistance(self, point1, point2):
        return math.sqrt((point1.x - point2.x)**2 + (point1.y - point2.y)**2)
    
    # calculate the perimeter of the face
    def calPerimeter(self, face_vertices):
        perimeter = 0
        for i in range(len(face_vertices)):
            perimeter += self.calDistance(face_vertices[i], face_vertices[(i+1) % (len(face_vertices))])
        return perimeter


    '''------------------------------------------------------------------------------------
    calculate the radius of the circle
    radius = sqrt(face_area/equilateral_area) * stand_radius.
    stand_radius = The radius of a circle that is as large as the area of the equilateral polygon.
    ---------------------------------------------------------------------------------------'''
    # Calculate the area of the standard equilateral polygon with the same vertex number as the face(edge length = 1)
    def calEqualiteral(self, face_vertices):
        vertex_num = len(face_vertices) # The number of vertices in a face, determine the optimal area of the face
        edge_length = 1 # Default length of each edge is 1
        area = edge_length**2 * vertex_num / (4 * math.tan(math.pi/vertex_num)) #formula
        print('Equaliteral:', area)
        return area
    
    # Calculate the radius of the circle
    def calCentroidAndRadius(self, face_vertices):
        centroid, real_area = self.calCentroid(face_vertices) # centroid of the circle is the same as the face 
        optimal_area = self.calEqualiteral(face_vertices) # real_area: current area of the face, optimal_area: area of the equilateral polygon
        stand_radius = math.pow(optimal_area/math.pi, 1/2) * 1 # The radius of a circle that is as large as the area of the equilateral polygon
        radius = math.sqrt(real_area/optimal_area) * stand_radius # radius of this face's circle
        # radius = 1.2*stand_radius
        self.radiusList.append(radius)
        return centroid, radius


    # Record all the incident edges of a vertex
    # NOTE: Cannot use .twin, .previous or .next!!!
    def getAllIncidentEdge(self, vertex, edge_dict): 
        edge_list = []   # All edges that incident to a vertex
        edge_identifier_set = vertex.incidentEdges  # identifier set of all incidentEdges
        for edge_identifier in edge_identifier_set:               
            edge = edge_dict[edge_identifier]
            edge_list.append(edge)
        return edge_list, edge_identifier_set
    

    '''------------------------------------------------------------------------------------
    Pull the points on the face directly to the corresponding position on the circle for each face
    --------------------------------------------------------------------------------------'''
   # Divide the circle equally
    def splitCircle(self, centroid, radius, section_count):
        vertices_in_circle = []  # section_count = number of vertices on the face
        vertices_in_circle.append(self.createVertexPure(centroid.x + radius, centroid.y, 0))
        radian = 2 * math.pi / section_count  # The radian between two equal division points
        count = 1
        while count < section_count:
            x = radius * math.cos(radian * count) + centroid.x # coordinates of equal division points
            y = radius * math.sin(radian * count) + centroid.y
            vertex = self.createVertexPure(x, y, 0)
            vertices_in_circle.append(vertex)
            count = count + 1
        return vertices_in_circle

    # pull the vertices on the face directly onto the circle
    def getDestinationVerticesOnCircle(self, face_dict, vertex):
        destination_vertices = [] # the corresponding position for a vertex on the circle
        distinct_set = set()
        for incident_face_id in vertex.faces: # for all incidientFace to a vertex 
            incident_face = face_dict[incident_face_id]
            if incident_face.identifier == "i":  # skip the infinite face
                continue
            vertices_set_dict = self.faceVerticesDict.get(incident_face.identifier)
            index_circle_vertices = vertices_set_dict.get("cur_vertices").get(vertex.identifier)
            destination_vertex = vertices_set_dict.get("circle_vertices")[index_circle_vertices]
            distinct_key = str(destination_vertex.x) + str(destination_vertex.y)
            if distinct_key in distinct_set:
                continue
            distinct_set.add(distinct_key)    
            destination_vertices.append(destination_vertex)
        return destination_vertices

    # the vertex with largest x-coordinate on the face corresponds 
    # to the point with largest x-coordinates on the circle
    def verticesSort(self, face_vertices):
        max = face_vertices[0].x
        max_index = 0
        index = 0
        for v in face_vertices:
            if max < v.x:
                max = v.x
                max_index = index
            index += 1
        slice_prefix = face_vertices[0:max_index]
        result_list = face_vertices[max_index:len(face_vertices)]
        result_list = result_list + slice_prefix
        identifer_dict = {}
        index = 0
        for v in result_list:
            identifer_dict[v.identifier] = index
            index += 1
        return identifer_dict

    
    # The target position of each adjacent surface is regarded as a vector,
    #  and the final displacement (not target position) is obtained by vector addition
    def calNewPosition(self, vertex, destination_vertices):
        vector_list = []       
        x_distance, y_distance = 0, 0
        for destination in destination_vertices:
            vector_x = destination.x - vertex.x
            vector_y = destination.y - vertex.y
            x_distance = vector_x + x_distance
            y_distance = vector_y + y_distance 
        # x_distance = x_distance / len(destination_vertices) + vertex.x
        # y_distance = y_distance / len(destination_vertices) + vertex.y  
        # used for crossover check  
        x_distance = x_distance / len(destination_vertices)
        y_distance = y_distance / len(destination_vertices)
        return x_distance, y_distance



    def findChain(self, vertexList, edge_dict):
        chainList = []    # All chains that start from a deg3+ vertex and end at a deg3+ vertex
        
        q = queue.Queue() # store deg3+ vertices
        
        # find a start deg3+ vertex
        for vertex in vertexList:
            if len(vertex.incidentEdges) > 2:
                # 
                q.put(vertex)
                
        # find chain
        TraveledVertexSet = set() # record if a vertex has been traveled
        while not q.empty():
            vertex = q.get()           

            edge_list, edge_identifier_set = self.getAllIncidentEdge(vertex, edge_dict) 

            for edge in edge_list:
                chain = []     # one of the chain in chainList
                chain.append(vertex) 
                edge = edge.next    
                nextVertex = edge.origin  
                isTraveledChain = False      # # record if a chain has been traveled  
                while len(nextVertex.incidentEdges) < 3:
                    # TODO
                    if nextVertex.identifier in TraveledVertexSet:
                        isTraveledChain = True
                        break
                    chain.append(nextVertex)
                    TraveledVertexSet.add(nextVertex.identifier)
                    edge = edge.next
                    nextVertex = edge.origin               
                chain.append(nextVertex)
                # TODO
                if isTraveledChain is True:
                    continue
                wrapperChain = WrapperChain(chain, self.edge_dict) # class Wrapperchain
                chainList.append(wrapperChain)
                print(chain) 
        return chainList

   

    # build the dictionary of the chains, used to calculate the optimal distance between two deg3+ vertices
    # the optimal distance is proportional to the number of deg2 vertices on the chain
    # the optimal distance is used to determine the radius of positioning-circle for deg3+ vertices with infinit face
    # key: deg3 vertex, value: all chains incident with this vertex
    def buildDeg3ChainsDict(self, chainList):
        deg3ChainDict = {} # a dict for deg3+ vertices
        for wrapperchain in chainList:
            keyDeg3Start = wrapperchain.chain[0] # start deg3+ vertex
            keyDeg3End = wrapperchain.chain[-1]  # end deg3+ vertex

            # calculate location circle of inside chain for outside deg3+
            if 1 == wrapperchain.chainType:
                incidentFace1 = wrapperchain.incidentFacesOfChain[0]
                incidentFace2 = wrapperchain.incidentFacesOfChain[1]
                centroid1 = self.faceCentroidDict[incidentFace1.identifier]
                centroid2 = self.faceCentroidDict[incidentFace2.identifier]
                radius1 = self.centroidRadiusDict[centroid1.identifier]
                radius2 = self.centroidRadiusDict[centroid2.identifier]
                arc = Arc(keyDeg3Start, keyDeg3End)
                arc.calApollonisCircle(centroid1, centroid2, radius1, radius2)
                
                # TODO other properties of three section arc
                wrapperchain.setThreeSectionArc(arc)
                # calculate location circle of inside chain for outside deg3+
                if wrapperchain.deg3Type == 2:
                    wrapperchain.calLocateCircleOfInsideChain(centroid1, centroid2, radius1, radius2)

        
            if keyDeg3Start not in deg3ChainDict:
                deg3ChainDict[keyDeg3Start] = [] # value(a list of all incident chains)
            deg3ChainDict[keyDeg3Start].append(wrapperchain)

            if keyDeg3End not in deg3ChainDict:
                deg3ChainDict[keyDeg3End] = []
            deg3ChainDict[keyDeg3End].append(wrapperchain)
        return deg3ChainDict
                


    
    
    # 
    def relateDeg3WithIntersection(self, kDeg3, intersection1, intersection2):
        distance1 = math.sqrt((kDeg3.x - intersection1[0]) ** 2 + (kDeg3.y - intersection1[1]) ** 2)
        distance2 = math.sqrt((kDeg3.x - intersection2[0]) ** 2 + (kDeg3.y - intersection2[1]) ** 2)
        return intersection1 if distance1 < distance2 else intersection2

    
    def handleDeg3Vertex_outside(self, deg3ChainDict):
        for kDeg3, vChainList in deg3ChainDict.items():
            # only handle the outside deg3+ vertex
            if len(kDeg3.incidentFaces) == len(kDeg3.incidentEdges):
                continue

            # 2 or more inside chains
            xCoordSum = 0
            yCoordSum = 0
            count = 0
            for wrapperChain in vChainList:
                # TODO whether empty chain should be deal with
                if wrapperChain.chainType == 2:
                    continue
                if wrapperChain.deg3Type == 1:
                    continue
                intersection1, intersection2 = wrapperChain.calLocateIntersection()
                intersection = self.relateDeg3WithIntersection(kDeg3, intersection1, intersection2)
                xCoordSum += intersection[0]
                yCoordSum += intersection[1]
                count += 1
            
            kDeg3.x = xCoordSum / count
            kDeg3.y = yCoordSum / count

        



    '''------------------------------------------------------------------------------------
    Preserve the topology: not cross
    ------------------------------------------------------------------------------------'''
    # Check if the target location will break the topology
    def checkNewPosition(self, vertex, x_distance, y_distance, incident_edge_list, edge_identifier_set):
        # Temporarily modify the point coordinates, used to detect whether the topology is broken
        
        # if(vertex.identifier == 58):
        #     print("烦死了。")
        
        reset_x = vertex.x
        reset_y = vertex.y 
        vertex.x = x_distance + vertex.x
        vertex.y = y_distance + vertex.y

        # Check whether there are edges cross each other
        isPassIntersecCheck = False
        for incident_edge in incident_edge_list:
            for edge in self.hedgeList:
                if edge.identifier not in edge_identifier_set:
                    count_set = set()
                    p1 = edge.origin
                    p2 = edge.next.origin
                    p3 = incident_edge.origin
                    p4 = incident_edge.next.origin
                    # if(p1.identifier == 56 or p1.identifier == 57):
                    #     print("????????????????")
                
                    count_set.add(p1.identifier)
                    count_set.add(p2.identifier)
                    count_set.add(p3.identifier)
                    count_set.add(p4.identifier)
                    if len(count_set) == 4:
                        isPassIntersecCheck = self.isIntersec(p1,p2,p3,p4)
                        
                        if isPassIntersecCheck is True:
                            # print("交叉四点：",p1,p2,p3,p4)
                            break
                    # TODO: if len(count_set) == 3, the overlap of the edges, need to be identified                
            if isPassIntersecCheck is True:
                break     
        # Coordinate reset
        vertex.x = reset_x
        vertex.y = reset_y       
        return not isPassIntersecCheck
    

    # Determines whether the two segments intersect
    def isIntersec(self, p1, p2, p3, p4): 
    #Fast rejection, l1, l2 as diagonal rectangles must intersect, otherwise the two segments do not intersect
        if(max(p1.x,p2.x)>=min(p3.x,p4.x)    #Rectangle 1 Far right Largest Rectangle 2 Leftmost Edge
        and max(p3.x,p4.x)>=min(p1.x,p2.x)   #The far right end of rectangle 2 is larger than the far left end of the rectangle
        and max(p1.y,p2.y)>=min(p3.y,p4.y)   #The highest end of rectangle 1 is greater than the lowest end of the rectangle
        and max(p3.y,p4.y)>=min(p1.y,p2.y)): #Rectangular 2 top end is greater than the lowest end of the rectangle
        #If fast rejection is passed, the Cross-over experiment is carried out
            if(self.cross(p1,p2,p3) * self.cross(p1,p2,p4)<=0
               and self.cross(p3,p4,p1) * self.cross(p3,p4,p2)<=0):
                return True
            else:
                return False
        else:
            return False

    # Cross-over experiment
    def cross(self, p1,p2,p3):
        x1=p2.x-p1.x
        y1=p2.y-p1.y
        x2=p3.x-p1.x
        y2=p3.y-p1.y
        return x1*y2-x2*y1 

    def getCentroidAndRadius(self, face, face_vertices):
        if face.identifier in self.faceCentroidDict.keys():
            centroid = self.faceCentroidDict.get(face.identifier)
            radius = self.centroidRadiusDict.get(centroid.identifier)
        else:
            centroid, radius = self.calCentroidAndRadius(face_vertices)
            self.faceCentroidDict[face.identifier] = centroid
            self.centroidFaceDict[centroid.identifier] = face
            self.centroidRadiusDict[centroid.identifier] = radius
        return centroid, radius


    def calEdgesOfCentroid(self, face_dict, distinct_edge, edges):
        for face in self.faceList:
            face_vertices = [v for v in face.loopOuterVertices()] 
            centroid, radius = self.getCentroidAndRadius(face, face_vertices)
        
            # Find adjacent faces and construct the graph for centroids
            # !!! We think faces with only one common vertex are adajcent
            for vertex in face_vertices:
                for neighbour_face_id in vertex.incidentFaces:
                    neighbour_face = face_dict[neighbour_face_id]
                    # The following two filters
                    if neighbour_face.identifier == face.identifier: # remove the face itself
                        continue
                    distinct_key = "%d-%d" % (face.identifier, neighbour_face.identifier)
                    distinct_key_reverse = "%d-%d" % (neighbour_face.identifier, face.identifier)
                    if distinct_key in distinct_edge or distinct_key_reverse in distinct_edge: # remove ba for ab
                        continue

                    neighbour_face_vertices = [v for v in neighbour_face.loopOuterVertices()]
                    neighbour_centroid, neighbour_radius = self.getCentroidAndRadius(neighbour_face, neighbour_face_vertices)
                    self.faceCentroidDict[neighbour_face.identifier] = neighbour_centroid
                    self.centroidRadiusDict[neighbour_centroid.identifier] = neighbour_radius
                    centroid_edge = edge(centroid, neighbour_centroid)
                    edges.append(centroid_edge)
                    distinct_edge.add(distinct_key)
            face_vertices = [v for v in face.loopOuterVertices()] 
            face_hedges = [e for e in face.loopOuterEdges()] 
            centroid, radius = self.getCentroidAndRadius(face, face_vertices)
         

    '''------------------------------------------------------------------------------------
    Main function: handle face one by one
    ---------------------------------------------------------------------------------------'''
    def handleFaces(self, switch):
        i = 1

        for repeat in range(1):
            edges = []
            distinct_edge = set()

            
            for face in self.faceList:
                self.face_dict[face.identifier] = face
            for hedge in self.hedgeList:
                self.edge_dict[hedge.identifier] = hedge

            #-------------------------------------------------------------------    
            # Traversal faces 1st time: calculate the centroids and radius of all faces
            # and calculate all the adjacent faces of this face(centroid graph)
            self.calEdgesOfCentroid(self.face_dict, distinct_edge, edges)

            # draw original map  
            self.centroidEdges = edges          
            # gui = pydcel.dcelVis(self)  
            
            # use force_directed approach to rearrange the centroid.  
            # note:Initialize the centroid to arrange the processor
            force_directed_draw = force_directed(self.centroidList, self.centroidRadiusDict, edges, self.centroidFaceDict)
            isHandle = True
            for i in range(1000):
                if switch == 'on' and isHandle:
                    # Rearrange the centroid
                    isHandle = force_directed_draw.handler()
                    # gui = pydcel.dcelVis(self)


            #----------------------------------------------------------------------------
            # Traversal faces 2nd time: dictionary to store face and vertices on it, and vertices and target
            # locations on the circle, then we can use it to drag the vertices to the circle
            for face in self.faceList:
                face_vertices = [v for v in face.loopOuterVertices()]
                centroid = self.faceCentroidDict.get(face.identifier)
                radius = self.centroidRadiusDict.get(centroid.identifier)

                face_vertices_dict = self.verticesSort(face_vertices) 
                vertices_in_circle = self.splitCircle(centroid, radius, len(face_vertices))
                self.faceVerticesDict[face.identifier] = {"cur_vertices":face_vertices_dict, "circle_vertices":vertices_in_circle}
               
            chainList = self.findChain(self.vertexList, self.edge_dict)
            print(chainList)

            # Move point to target position
            for vertex in self.vertexList:
                if len(vertex.incidentEdges) < 3:
                    continue
                # edge_list, edge_identifier_set = self.getAllIncidentEdge(vertex, edge_dict)
                # destination_vertices = self.getDestinationVerticesOnCircle(face_dict, vertex)
                centroidsOfIncidentFace = []
                # Find the centroids of the circle around the current deg3+ point
                for faceIdentifier in vertex.incidentFaces:
                    centroid = self.faceCentroidDict.get(faceIdentifier)
                    centroidsOfIncidentFace.append(centroid)

                
                # ----------deg3+ inside
                # move the deg3+ vertex into the ploygon formed by centroids of the around circle
                if len(centroidsOfIncidentFace) < 3:
                    continue


                centroid_of_centroids, area_of_centroids = self.calCentroid(centroidsOfIncidentFace)
                vertex.x = centroid_of_centroids.x
                vertex.y = centroid_of_centroids.y
                # Use the centroids to calculate the attraction and repulsive force and move the current deg3+ point               
                force_directed_draw.handle3DegVertex_inside(vertex, centroidsOfIncidentFace)
                
                chainList = self.findChain(self.vertexList, self.edge_dict)
                print(chainList)
                
            # before rotate           
            # just use straight line to show the results
            for wrapperchain in chainList:
                first_deg2 = wrapperchain.chain[1]
                if len(first_deg2.incidentEdges) > 2:
                    continue
                face = first_deg2.incidentFaces
                # if len(face) < 2:
                #     continue

                xdis = wrapperchain.chain[-1].x - wrapperchain.chain[0].x
                ydis = wrapperchain.chain[-1].y - wrapperchain.chain[0].y
                xUnitDis = xdis / (len(wrapperchain.chain) - 1)
                yUnitDis = ydis / (len(wrapperchain.chain) - 1)

                multi = 1
                for deg2 in wrapperchain.chain[1:-1]:                    
                    deg2.x = wrapperchain.chain[0].x + multi * xUnitDis
                    deg2.y = wrapperchain.chain[0].y + multi * yUnitDis
                    multi += 1
            gui = pydcel.dcelVis(self)

            ####################-----------------------------##########################
            # after the force-directed method, if there existes very short chains, do the next rotate part
            ########## rotate the circle to make  short chain longer
            # for vertex in self.vertexList:
            #     if len(vertex.incidentEdges) < 3:
            #         continue
            #     centroidsOfIncidentFace = []
            #     # Find the centroids of the circle around the current deg3+ point
            #     for faceIdentifier in vertex.incidentFaces:
            #         centroid = self.faceCentroidDict.get(faceIdentifier)
            #         centroidsOfIncidentFace.append(centroid)
            #     if len(centroidsOfIncidentFace) < 3:
            #         continue

            #     # rotate
            #     for wrapperChain in chainList:
            #         chain = wrapperChain.chain
            #         start3Deg = chain[0]
            #         end3Deg = chain[-1]
            #         if start3Deg.identifier == vertex.identifier or end3Deg.identifier == vertex.identifier:
            #             disBetween3Deg = math.sqrt((start3Deg.x - end3Deg.x) ** 2 + (start3Deg.y - end3Deg.y) ** 2) # real length of the chain

            #             # If chain length less than threshold (fractional times optimal distance), rotate
            #             if (len(chain) - 1) * 5 > disBetween3Deg:
            #                 faceId1 = list(chain[1].incidentFaces)[0] # two incident faces of the chain
            #                 faceId2 = list(chain[1].incidentFaces)[1]
            #                 centroid1 = self.faceCentroidDict.get(faceId1) # two centroids of the two incident faces
            #                 centroid2 = self.faceCentroidDict.get(faceId2)
            #                 # rotation angle, can be adjusted
            #                 angle = math.radians(30)
            #                 # formula of rotation
            #                 # take theh cetroid1 as center, rotate the centroid2
            #                 centroid2.x = (centroid2.x-centroid1.x)*math.cos(angle) - (centroid2.y-centroid1.y)*math.sin(angle)+centroid1.x
            #                 centroid2.y = (centroid2.y-centroid1.y)*math.cos(angle) + (centroid2.x-centroid1.x)*math.sin(angle)+centroid1.y
            #                 # put the deg3+ vertex at the centroid of the centroids of the rotated circles(same as before) 
            #                 centroid_of_centroids, area_of_centroids = self.calCentroid(centroidsOfIncidentFace)
            #                 vertex.x = centroid_of_centroids.x
            #                 vertex.y = centroid_of_centroids.y
            #                 # Use the centroids to calculate the attraction and repulsive force and move the current deg3+ point
            #                 force_directed_draw.handle3DegVertex_inside(vertex, centroidsOfIncidentFace)
                           

            # ----------deg3+ outside

            deg3ChainDict = self.buildDeg3ChainsDict(chainList)
            self.handleDeg3Vertex_outside(deg3ChainDict)


                


                # TODO: Check for crossovers
                # x_distance, y_distance = self.calNewPosition(vertex, destination_vertices)
                # proportion = 0
                # for i in range(1, 10):
                #     if not self.checkNewPosition(vertex, x_distance * (i / 10), y_distance * (i / 10), edge_list, edge_identifier_set):
                #         break
                #     proportion = i / 10

                # vertex.x = x_distance * proportion + vertex.x
                # vertex.y = y_distance * proportion + vertex.y
                # # -------------------------
                # count = 0
                # while not self.checkNewPosition(vertex, x_distance, y_distance, edge_list, edge_identifier_set):
                #     x_distance = x_distance / 2
                #     y_distance = y_distance / 2
                #     count +=1
                # vertex.x = x_distance + vertex.x
                # vertex.y = y_distance + vertex.y


            

            # just use straight line to show the results
            for wrapperchain in chainList:
                first_deg2 = wrapperchain.chain[1]
                if len(first_deg2.incidentEdges) > 2:
                    continue
                face = first_deg2.incidentFaces
                # if len(face) < 2:
                #     continue

                xdis = wrapperchain.chain[-1].x - wrapperchain.chain[0].x
                ydis = wrapperchain.chain[-1].y - wrapperchain.chain[0].y
                xUnitDis = xdis / (len(wrapperchain.chain) - 1)
                yUnitDis = ydis / (len(wrapperchain.chain) - 1)

                multi = 1
                for deg2 in wrapperchain.chain[1:-1]:                    
                    deg2.x = wrapperchain.chain[0].x + multi * xUnitDis
                    deg2.y = wrapperchain.chain[0].y + multi * yUnitDis
                    multi += 1

            # draw the adjusted graph   
            gui = pydcel.dcelVis(self)  
             


        gui.mainloop()

        
        



                # another method to calculate the target position, not use vector addition
                # if len(destination_vertices) == 1:
                #     vertex.x = destination_vertices[0].x
                #     vertex.y = destination_vertices[0].y
                # elif len(destination_vertices) == 2:
                #     vertex.x = (destination_vertices[0].x + destination_vertices[1].x) / 2
                #     vertex.y = (destination_vertices[0].y + destination_vertices[1].y) / 2
                # else:
                #     destination_centroid, totalArea = self.calCentroid(destination_vertices)
                #     vertex.x = destination_centroid.x
                #     vertex.y = destination_centroid.y


            # Find adjacent faces and construct the graph for centroids
            # !!! We think faces with only one common vertex are NOT adajcent
            # for hedge in face_hedges:
            #     neighbour_face = hedge.twin.incidentFace
            #     if "i" == neighbour_face.identifier:
            #         continue
                
            #     # The following two filters
            #     if neighbour_face.identifier == face.identifier: # remove the face itself
            #         continue
                
            #     distinct_key = "%d-%d" % (face.identifier, neighbour_face.identifier)
            #     distinct_key_reverse = "%d-%d" % (neighbour_face.identifier, face.identifier)
            #     if distinct_key in distinct_edge or distinct_key_reverse in distinct_edge: # remove ba for ab
            #         continue

            #     neighbour_face_vertices = [v for v in neighbour_face.loopOuterVertices()]
            #     neighbour_centroid, neighbour_radius = self.getCentroidAndRadius(neighbour_face, neighbour_face_vertices)
            #     self.faceCentroidDict[neighbour_face.identifier] = neighbour_centroid
            #     self.centroidRadiusDict[neighbour_centroid.identifier] = neighbour_radius
            #     centroid_edge = edge(centroid, neighbour_centroid)
            #     edges.append(centroid_edge)
            #     distinct_edge.add(distinct_key)