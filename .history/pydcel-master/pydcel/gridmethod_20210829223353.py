import math
import pydcel

class vertex(object):
    
    def __init__(self, px, py, pz, identifier):
        self.identifier = identifier
        self.x = px
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
        self.identifier = identifier
        self.origin = None
        self.twin = None
        self.incidentFace = None
        self.next = None
        self.previous = None

    def setTopology(self, newOrigin, newTwin, newIncindentFace, newNext, newPrevious):
        self.origin = newOrigin
        self.twin = newTwin
        self.incidentFace = newIncindentFace
        self.next = newNext
        self.previous = newPrevious
        
    def loop(self):
        """Loop from this hedge to the next ones. Stops when we are at the current one again."""
        yield self
        e = self.next
        while e is not self:
            yield e
            e = e.next
            
    def wind(self):
        """iterate over hedges emerging from vertex at origin in ccw order"""
        yield self
        e = self.previous.twin
        while e is not self:
            yield e
            e = e.previous.twin

    def __repr__(self):
        return "he{}".format(self.identifier)


class face(object):
    
    def __init__(self, identifier):
        self.identifier = identifier
        self.outerComponent = None
        self.innerComponent = None

    def setTopology(self, newOuterComponent, newInnerComponent=None):
        self.outerComponent = newOuterComponent
        self.innerComponent = newInnerComponent
        
    def loopOuterVertices(self):
        for e in self.outerComponent.loop():
            yield e.origin

    def __repr__(self):
        # return "face( innerComponent-{}, outerComponent-{} )".format(self.outerComponent, self.innerComponent)
        return "f{}".format(self.identifier)


class DCEL(object):
    
    def __init__(self):
        self.vertexList = []
        self.vertexSet = set()
        self.hedgeList = []
        self.faceList = []
        self.centriodList = []
        self.radiusList = []
        self.infiniteFace = None
        self.incidentEdgeAll = {}               
        self.vertexOnCircleList = []   # equal division points on the circle
        self.faceVerticesDict = {}  # Dictionary: the correspondence between the face and its vertices 
        self.faceCentroidDict = {}    # Dictionary: the correspondence between the face and its centroid 
        self.centroidRadiusDict = {}  # Dictionary: the correspondence between the centroid and the radius
        self.centroidFaceDict = {}        
        self.centroidEdges = [] # all the edges of the centroid graph
        self.apollonisCenterList = []
        self.apollonisRadiusList = []
        self.face_dict = {} # dictionary that finds faces by id
        self.edge_dict = {} # dictionary that finds hedges by id

    def getNewId(self, L):
        if len(L) == 0:
            return 0
        else:
            return L[-1].identifier + 1
        
    def createVertex(self, px, py, pz):
        identifier = self.getNewId(self.vertexList)
        v = vertex(px,py,pz, identifier)
        self.vertexList.append(v)
        self.vertexSet.add("%d-%d" % (v.x, v.y))
        return v
    
    def createVertexPure(self, px, py, pz):
        '''Used to create centroid for each face and five possible positions for each vertex, 
        and not append centroid into vertexList and vertexSet'''
        identifier = self.getNewId(self.vertexList)
        v = vertex(px,py,pz, identifier)
        return v
        
    def createHedge(self):
        identifier = self.getNewId(self.hedgeList)
        e = hedge(identifier)
        self.hedgeList.append(e)
        return e
        
    def createFace(self):
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
            s += "{} : \t{} \t{}\n".format(f, f.outerComponent, f.innerComponent)
        return s

    def checkEdgeTwins(self):
        for e in self.hedgeList:
            if not e == e.twin.twin:
                print(("this edge has a problem with its twin:"), end=' ')
                print(e)

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
        nice_face.outerComponent = e_0.__next__
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
    def calAllCentroid(self, face_vertices):
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
        centroid = self.createVertexPure(round(x, 3), round(y, 3), 0)
        return centroid, abs(totalArea)


    '''------------------------------------------------------------------------------------
    calculate the roundness of a face
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

    # calculate the roundness of the face
    def calRoundness(self, face_vertices):
        centroid, area = self.calAllCentroid(face_vertices)
        perimeter = self.calPerimeter(face_vertices)
        roundness = 4 * math.pi * area / perimeter**2
        return roundness


    '''------------------------------------------------------------------------------------
    calculate the radius of the circle
    radius = sqrt(face_area/equilateral_area) * stand_radius.
    stand_radius = The radius of a circle that is as large as the area of the equilateral polygon.
    ---------------------------------------------------------------------------------------'''
    # Calculate the area and roundness of the standard equilateral polygon with the same vertex number as the face(edge length = 1)
    def calEqualiteral(self, face_vertices):
        vertex_num = len(face_vertices)
        edge_length = 1
        area = edge_length**2 * vertex_num / (4 * math.tan(math.pi/vertex_num))
        roundness = 4 * math.pi * area / (vertex_num * edge_length)**2
        print(math.tan(math.pi/vertex_num))
        print(len(face_vertices))
        print(area)
        return area, roundness
    
    # Calculate the radius of the circle
    def calCircleRadius(self, face_vertices):
        centroid, real_area = self.calAllCentroid(face_vertices)
        optimal_area, max_roundness = self.calEqualiteral(face_vertices)
        stand_radius = math.sqrt(optimal_area/math.pi)
        radius = math.sqrt(real_area/optimal_area) * stand_radius
        return radius


    '''------------------------------------------------------------------------------------
    List the vertices inside the circle both on the left chain and the right chain 
    ---------------------------------------------------------------------------------------'''
    # Determine which points are in the circle
    def checkCloseVertices(self, face_vertices, centroid, radius):
        close_set = set()
        max_vertex = face_vertices[0]
        min_vertex = face_vertices[0]
        for vertex in face_vertices:
            distance = self.calDistance(vertex, centroid)
            if distance < radius:
                close_set.add(vertex.identifier)
            if(vertex.y > max_vertex.y):
                max_vertex = vertex
            if(vertex.y < min_vertex.y):
                min_vertex = vertex
        return close_set, max_vertex, min_vertex

    # divide the vertices of a face into left chain and right chain
    def cal2WayChain(self, face_vertices, centroid, radius):
        close_set, max_vertex, min_vertex = self.checkCloseVertices(face_vertices, centroid, radius)
        print("--------------------------")
        print(close_set)
        print(max_vertex)
        print(min_vertex)
        print("--------------------------")
        prefix_vertices = []
        left_chain_vertices = []
        right_chain_vertices = []
        # The maximum point and the minimun point divide the array into up to three segments, with check_point record the current position, 0, 1, 2
        check_point = 0
        # Record whether part 2 and 3 belong to the left chain or the right chain
        is_left = 0
        index = 0
        incidentEdge = face_vertices[0].incidentEdge
        
        for vertex in face_vertices:
            # The first part belongs to the right chain
            if vertex.identifier == max_vertex.identifier and check_point == 0:
                left_chain_vertices = prefix_vertices
                check_point = 1
                is_left = 0
                continue
            # The first paragraph belongs to the left chain
            if vertex.identifier == min_vertex.identifier and check_point == 0:
                right_chain_vertices = prefix_vertices
                check_point = 1
                is_left = 1
                continue
            if check_point == 0 and vertex.identifier in close_set:
                prefix_vertices.append(vertex)
            
            if vertex.identifier == max_vertex.identifier and check_point == 1:
                check_point = 2
                is_left = 0
                continue
            if vertex.identifier == min_vertex.identifier and check_point == 1:
                check_point = 2
                is_left = 1
                continue
            if check_point == 1 and vertex.identifier in close_set:
                if is_left == 1:
                    left_chain_vertices.append(vertex)
                else:
                    right_chain_vertices.append(vertex)

            if check_point == 2 and vertex.identifier in close_set:
                if is_left == 1:
                    left_chain_vertices.append(vertex)
                else:
                    right_chain_vertices.append(vertex)
        print("Left Chain: ", end='')
        print(left_chain_vertices)
        print("Right Chain: ", end='')
        print(right_chain_vertices)
        return left_chain_vertices, right_chain_vertices

    # 5 possible positions for vertices on left chain inside the circle
    def listLeftChainGridPoint(self, vertex):
        vertex_l = self.createVertexPure(vertex.x-1, vertex.y, vertex.z)
        vertex_ul = self.createVertexPure(vertex.x-1, vertex.y+1, vertex.z)
        vertex_bl = self.createVertexPure(vertex.x-1, vertex.y-1, vertex.z)
        vertex_u = self.createVertexPure(vertex.x, vertex.y+1, vertex.z)
        vertex_b = self.createVertexPure(vertex.x, vertex.y-1, vertex.z)
        left_chain_point_list = [vertex_u, vertex_ul, vertex_l, vertex_bl, vertex_b]
        return left_chain_point_list

    # 5 possible positions for vertices on right chain inside the circle
    def listRightChainGridPoint(self, vertex):
        vertex_u = self.createVertexPure(vertex.x, vertex.y+1, vertex.z)
        vertex_ur = self.createVertexPure(vertex.x+1, vertex.y+1, vertex.z)
        vertex_br = self.createVertexPure(vertex.x+1, vertex.y-1, vertex.z)
        vertex_r = self.createVertexPure(vertex.x+1, vertex.y, vertex.z)
        vertex_b = self.createVertexPure(vertex.x, vertex.y-1, vertex.z)
        right_chain_point_list = [vertex_u, vertex_ur, vertex_r, vertex_br, vertex_b]
        return right_chain_point_list


    '''------------------------------------------------------------------------------------
    Determine which position can a vertex move to
    ---------------------------------------------------------------------------------------'''
    # Octilinear layout, the slope of all edges have to be an integer multiple of 45 degrees
    def checkNeighborSlope(self, vertex, destination_point):
        neighbor_point_list = []
        incidentEdge = vertex.incidentEdge
        incident_identifier = incidentEdge.identifier
        identifier = 0
        while identifier != incident_identifier:
            neighbor_vertex = incidentEdge.next.origin
            if neighbor_vertex.x != destination_point.x:
                slope = (neighbor_vertex.y - destination_point.y) / (neighbor_vertex.x - destination_point.x)
                if 0 != slope and 1 != abs(slope):
                    return False
            twin = incidentEdge.twin
            incidentEdge = twin.next
            identifier = incidentEdge.identifier       
        return True
    
    # Two line segments can not intersect
    def cross(self, p1,p2,p3): # Cross-over experiment
        x1=p2.x-p1.x
        y1=p2.y-p1.y
        x2=p3.x-p1.x
        y2=p3.y-p1.y
        return x1*y2-x2*y1     

    def isIntersec(self, p1, p2, p3, p4): # Repulsion experiment
        # If the rectangle does not intersect, then the two line segments must not intersect
        if(max(p1.x,p2.x)>=min(p3.x,p4.x)    # The rightmost end of Rectangle 1 is greater than the leftmost end of Rectangle 2
        and max(p3.x,p4.x)>=min(p1.x,p2.x)   # The rightmost end of Rectangle 2 is greater than the leftmost end of Rectangle 1
        and max(p1.y,p2.y)>=min(p3.y,p4.y)   # Rectangle 1 at the top is greater than Rectangle 1 at the bottom
        and max(p3.y,p4.y)>=min(p1.y,p2.y)): # Rectangle 2 at the top is greater than Rectangle 1 at the bottom
        # If repulsion experiment is passed, the cross-over experiment is carried out
            if(self.cross(p1,p2,p3) * self.cross(p1,p2,p4)<=0
               and self.cross(p3,p4,p1) * self.cross(p3,p4,p2)<=0):
                return True
            else:
                return False
        else:
            return False
    
    # record the order of the incident edges for each vertex    
    def vertexOrdering(self, vertexList):   
        edge_order_dict = {}
        for vertex in vertexList:
            edge_order_list = []
            incidentEdge = vertex.incidentEdge
            incident_identifier = incidentEdge.identifier
            edge_order_list.append(incident_identifier)
            identifier = 0
            while identifier != incident_identifier:               
                twin = incidentEdge.twin
                incidentEdge = twin.next
                identifier = incidentEdge.identifier 
                edge_order_list.append(identifier)      
            edge_order_dict[vertex.identifier] = edge_order_list
        return edge_order_dict

    # Judge whether the order is the same before and after moving vertex
    def compareOrdering(self, sourceDict, targetDict):
        for key,sourcelist in sourceDict.items():
            targetlist = targetDict[key]
            for i,j in zip(sourcelist, targetlist):
                if i != j:
                    return False
        return True

    # Record all the incident edges of a vertex
    def getAllIncidentEdge(self, vertex): 
        edge_list = []
        edge_ientifier_set = set()
        incidentEdge = vertex.incidentEdge
        incident_identifier = incidentEdge.identifier
        edge_list.append(incidentEdge)
        edge_ientifier_set.add(incident_identifier)
        identifier = 0
        while identifier != incident_identifier:               
            twin = incidentEdge.twin
            incidentEdge = twin.next
            identifier = incidentEdge.identifier 
            edge_list.append(incidentEdge)
            edge_ientifier_set.add(identifier)
        return edge_list, edge_ientifier_set
    
   
    # Determine if a grid can be used
    def getChainGrid(self, vertex, side):
        chain_point_list = []
        if side == "left":
            chain_point_list = self.listLeftChainGridPoint(vertex)
        else:
            chain_point_list = self.listRightChainGridPoint(vertex)
        
        for moving_point in chain_point_list:
            # not overlap with other vertex
            key = "%d-%d" % (moving_point.x, moving_point.y)
            if key in self.vertexSet:
                continue
            # octilinear
            if self.checkNeighborSlope(vertex, moving_point) is False:
                continue
            # not cross, preserve the topology structure            
            reset_x = vertex.x
            reset_y = vertex.y
            vertex.x = moving_point.x
            vertex.y = moving_point.y
            incident_edge_list, edge_ientifier_set = self.getAllIncidentEdge(vertex)
            isPassIntersecCheck = False
            for incident_edge in incident_edge_list:
                if incident_edge.identifier == 3:
                    print()
                for edge in self.hedgeList:
                    if edge.identifier not in edge_ientifier_set:
                        count_set = set()
                        p1 = edge.origin
                        p2 = edge.previous.origin
                        p3 = incident_edge.origin
                        p4 = incident_edge.previous.origin
                        count_set.add(p1.identifier)
                        count_set.add(p2.identifier)
                        count_set.add(p3.identifier)
                        count_set.add(p4.identifier)
                        if len(count_set) == 4:
                            isPassIntersecCheck = self.isIntersec(p1,p2,p3,p4)
                            if isPassIntersecCheck is True:
                                break
                        # todo 3个点时，边重叠的情况，需要识别                        
                if isPassIntersecCheck is True:
                    break                               
            vertex.x = reset_x
            vertex.y = reset_y
            # Update Set
            if isPassIntersecCheck is True:
                continue
            oldkey = "%d-%d" % (vertex.x, vertex.y)
            self.vertexSet.remove(oldkey)
            self.vertexSet.add(key)
            return moving_point        
        return None


    '''------------------------------------------------------------------------------------
    If there a possible position can be used, move the vertex
    ---------------------------------------------------------------------------------------'''   
    # move the vertex
    def moveVertices(self, left_chain_vertices, right_chain_vertices):
        for vertex in left_chain_vertices:
            left_chain_point = self.getChainGrid(vertex, "left")
            if left_chain_point is None:
                continue
            vertex.x = left_chain_point.x
            vertex.y = left_chain_point.y
        for vertex in right_chain_vertices:
            right_chain_point = self.getChainGrid(vertex, "right")
            if right_chain_point is None:
                continue
            vertex.x = right_chain_point.x
            vertex.y = right_chain_point.y
    

    '''------------------------------------------------------------------------------------
    handle face one by one
    ---------------------------------------------------------------------------------------'''
    def handleFaces(self,switch):
        # self.incidentEdgeAll = self.vertexOrdering(self.vertexList)
        # vertex_list = []
        i = 1
        gui = pydcel.dcelVis(self)
        for face in self.faceList:
            face_vertices = [v for v in face.loopOuterVertices()]
            print("************************************************************")
            print("Face: %d" %(i))
            print("Vertices:", end = "")
            print(face_vertices)
            centriod, area = self.calAllCentroid(face_vertices)
            radius = self.calCircleRadius(face_vertices)
            self.radiusList.append(radius)
            self.centriodList.append(centriod)
            roundness = self.calRoundness(face_vertices)
            
            optimal_area, max_roundness = self.calEqualiteral(face_vertices)
            
            #while (area > 1.5 * optimal_area) or (roundness < 0.435 * max_roundness): 
            while (roundness < 0.3 * max_roundness):
                left_chain_vertices, right_chain_vertices = self.cal2WayChain(face_vertices, centriod, radius)
                print("before")
                print(left_chain_vertices)
                optimal_area, max_roundness = self.calEqualiteral(face_vertices)
                self.moveVertices(left_chain_vertices, right_chain_vertices)
                print("after")
                print(left_chain_vertices)
                optimal_area, max_roundness = self.calEqualiteral(face_vertices)
                
                gui = pydcel.dcelVis(self)  
                centriod, area = self.calAllCentroid(face_vertices)
                roundness = self.calRoundness(face_vertices)
                self.centriodList.pop()
                self.radiusList.pop()
                self.radiusList.append(radius)
                self.centriodList.append(centriod)
            i = i + 1
        gui = pydcel.dcelVis(self)  
        gui.mainloop()
