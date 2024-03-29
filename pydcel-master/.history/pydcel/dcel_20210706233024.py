'''
Construct DCEL and some related calculations
'''

import math
import pydcel
from .force_directed_draw import *


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
        self.radiusList = []    # all the radius of the circle
        self.infiniteFace = None

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
        '''Used to create auxiliary points, such as equal division points on the circle'''
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
            if incident_face.identifier is "i":  # skip the infinite face
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
        x_distance = x_distance / len(destination_vertices)
        y_distance = y_distance / len(destination_vertices)
        return x_distance, y_distance



    '''------------------------------------------------------------------------------------
    Preserve the topology: not cross
    ------------------------------------------------------------------------------------'''
    # Check if the target location will break the topology
    def checkNewPosition(self, vertex, x_distance, y_distance, incident_edge_list, edge_identifier_set):
        # Temporarily modify the point coordinates, used to detect whether the topology is broken
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
                    count_set.add(p1.identifier)
                    count_set.add(p2.identifier)
                    count_set.add(p3.identifier)
                    count_set.add(p4.identifier)
                    if len(count_set) == 4:
                        isPassIntersecCheck = self.isIntersec(p1,p2,p3,p4)
                        if isPassIntersecCheck is True:
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
            self.centroidRadiusDict[centroid.identifier] = radius
        return centroid, radius




    '''------------------------------------------------------------------------------------
    Main function: handle face one by one
    ---------------------------------------------------------------------------------------'''
    def handleFaces(self, isHandle):
        i = 1

        for repeat in range(1):
            edges = []
            distinct_edge = set()

            face_dict = {} # dictionary that finds faces by id
            for face in self.faceList:
                face_dict[face.identifier] = face
            edge_dict = {} # dictionary that finds hedges by id
            for hedge in self.hedgeList:
                edge_dict[hedge.identifier] = hedge

            #-------------------------------------------------------------------    
            # Traversal faces 1st time: calculate the centroids and radius of all faces
            # and calculate all the adjacent faces of this face(centroid graph)
            for face in self.faceList:
                face_vertices = [v for v in face.loopOuterVertices()] 
                centroid, radius = self.getCentroidAndRadius(face, face_vertices)
         
                # Find adjacent faces and construct the graph for centroids
                for vertex in face_vertices:
                    for neighbour_face_id in vertex.faces:
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

            # draw original map            
            gui = pydcel.dcelVis(self)  
            
            # use force_directed approach to rearrange the centroid
            force_directed_draw = force_directed(self.centroidList, self.centroidRadiusDict, edges)
            for i in range(1000):
                if isHandle:
                    force_directed_draw.handler()
                    # gui = pydcel.dcelVis(self)


            #----------------------------------------------------------------------------
            # Traversal faces 2nd time: dictionary to store face and vertices on it
            for face in self.faceList:
                face_vertices = [v for v in face.loopOuterVertices()]
                # print("************************************************************")
                # print("Face: %d" %(i))
                # print("Vertices:", end = "")
                # print(face_vertices)
                # centroid, area = self.calCentroid(face_vertices)
                centroid = self.faceCentroidDict.get(face.identifier)
                radius = self.centroidRadiusDict.get(centroid.identifier)

                face_vertices_dict = self.verticesSort(face_vertices) 
                vertices_in_circle = self.splitCircle(centroid, radius, len(face_vertices))
                self.faceVerticesDict[face.identifier] = {"cur_vertices":face_vertices_dict, "circle_vertices":vertices_in_circle}
               
            # Move point to target position
            for vertex in self.vertexList:
                edge_list, edge_identifier_set = self.getAllIncidentEdge(vertex, edge_dict)
                destination_vertices = self.getDestinationVerticesOnCircle(face_dict, vertex)
                
                x_distance, y_distance = self.calNewPosition(vertex, destination_vertices)
                # TODO: Check for crossovers
                # while not self.checkNewPosition(vertex, x_distance, y_distance, edge_list, edge_identifier_set):
                #     x_distance = x_distance / 2
                #     y_distance = y_distance / 2
                vertex.x = x_distance + vertex.x
                vertex.y = y_distance + vertex.y

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

