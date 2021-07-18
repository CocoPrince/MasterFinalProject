import math
import pydcel
from .force_directed_draw import *

class vertex(object):
    
    def __init__(self, px, py, pz, identifier):
        self.identifier = identifier
        self.x = px
        self.y = py
        self.z = pz
        self.incidentEdge = None
        
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
        self.centroidList = []
        self.vertexOnCircleList = []
        self.faceVerticesDict = {}
        #记录面与质心的对应关系
        self.faceCentroidDict = {}
        # 记录质心与半径的对应关系
        self.centroidRadiusDict = {}
        self.radiusList = []
        self.infiniteFace = None
        self.incidentEdgeAll = {}

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
    
    def createCentroid(self, px, py, pz):
        '''Used to create centroid for each face and five possible positions for each vertex, 
        and not append centroid into vertexList and vertexSet'''
        identifier = self.getNewId(self.centroidList)
        v = vertex(px,py,pz, identifier)
        self.centroidList.append(v)
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
    def calCentroidAndRadius(self, face_vertices):
        centroid, real_area = self.calCentroid(face_vertices)
        optimal_area, max_roundness = self.calEqualiteral(face_vertices)
        stand_radius = math.pow(optimal_area/math.pi, 1/2) * 1
        radius = math.sqrt(real_area/optimal_area) * stand_radius
        return centroid, radius



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
    
   
    

    def splitCircle(self, centroid, radius, section_count):
        vertices_in_circle = []
        vertices_in_circle.append(self.createVertexPure(centroid.x + radius, centroid.y, 0))
        radian = 2 * math.pi / section_count
        count = 1
        while count < section_count:
            x = radius * math.cos(radian * count) + centroid.x
            y = radius * math.sin(radian * count) + centroid.y
            vertex = self.createVertexPure(x, y, 0)
            vertices_in_circle.append(vertex)
            count = count + 1
        
        print("+++++++++++++++++++++++++++++++分割圆")
        # print(vertices_in_circle)
        return vertices_in_circle


    def getDestinationVerticesOnCircle(self, edge_list, vertex):
        destination_vertices = []
        distinct_set = set()
        for edge in edge_list:
            incident_face = edge.incidentFace
            if incident_face.identifier is "i":
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

        # {identifer_vertex:index_vertex}
        identifer_dict = {}
        index = 0
        for v in result_list:
            identifer_dict[v.identifier] = index
            index += 1
        return identifer_dict

    
    # 向量加法,得到移动的距离，但非最终位置
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

    def checkNewPosition(self, vertex, x_distance, y_distance, incident_edge_list, edge_ientifier_set):
        # 不能交叉

        # 临时修改点坐标
        reset_x = vertex.x
        reset_y = vertex.y
        vertex.x = x_distance + vertex.x
        vertex.y = y_distance + vertex.y

        # 进行交叉检验
        isPassIntersecCheck = False
        for incident_edge in incident_edge_list:
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
        

        # 坐标复位
        vertex.x = reset_x
        vertex.y = reset_y
        
        return not isPassIntersecCheck
    
    def isIntersec(self, p1, p2, p3, p4): #判断两线段是否相交
    #快速排斥，以l1、l2为对角线的矩形必相交，否则两线段不相交
        if(max(p1.x,p2.x)>=min(p3.x,p4.x)    #矩形1最右端大于矩形2最左端
        and max(p3.x,p4.x)>=min(p1.x,p2.x)   #矩形2最右端大于矩形最左端
        and max(p1.y,p2.y)>=min(p3.y,p4.y)   #矩形1最高端大于矩形最低端
        and max(p3.y,p4.y)>=min(p1.y,p2.y)): #矩形2最高端大于矩形最低端

        #若通过快速排斥则进行跨立实验
            if(self.cross(p1,p2,p3) * self.cross(p1,p2,p4)<=0
               and self.cross(p3,p4,p1) * self.cross(p3,p4,p2)<=0):
                return True
            else:
                return False
        else:
            return False

    def cross(self, p1,p2,p3):#跨立实验
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
    handle face one by one
    ---------------------------------------------------------------------------------------'''
    def handleFaces(self, isHandleCentroid):
        # self.incidentEdgeAll = self.vertexOrdering(self.vertexList)
        # vertex_list = []
        i = 1
        gui = pydcel.dcelVis(self)

        for repeat in range(1):
            # face 遍历1——第一次得到计算质心与半径，并计算相邻面的质心与半径
            edges = []
            distinct_edge = set()
            for face in self.faceList:
                face_vertices = [v for v in face.loopOuterVertices()]
                centroid, radius = self.getCentroidAndRadius(face, face_vertices)
                self.radiusList.append(radius)

                # 求相邻面
                for vertex in face_vertices:
                    incident_hedge = vertex.incidentEdge
                    twin = incident_hedge.twin
                    neighbour_face = twin.incidentFace
                    # 以下三个过滤条件
                    if neighbour_face.identifier is 'i':
                        continue
                    # todo 为何存在twin的incidentFace与当前Face相同
                    if neighbour_face.identifier == face.identifier:
                        continue
                    distinct_key = "%d-%d" % (face.identifier, neighbour_face.identifier)
                    distinct_key_reverse = "%d-%d" % (neighbour_face.identifier, face.identifier)
                    if distinct_key in distinct_edge or distinct_key_reverse in distinct_edge:
                        continue

                    neighbour_face_vertices = [v for v in neighbour_face.loopOuterVertices()]
                    neighbour_centroid, neighbour_radius = self.getCentroidAndRadius(neighbour_face, neighbour_face_vertices)
                    self.faceCentroidDict[neighbour_face.identifier] = neighbour_centroid
                    self.centroidRadiusDict[neighbour_centroid.identifier] = neighbour_radius
                    centroid_edge = edge(centroid, neighbour_centroid)
                    edges.append(centroid_edge)
                    distinct_edge.add(distinct_key)

              
            # force_directed
            force_directed_draw = force_directed(self.centroidList, self.centroidRadiusDict, edges)
            print("Original centroid positions：")
            print(self.centroidList)
            if isHandleCentroid:
                force_directed_draw.handler()
            print("Improved centroid positions：")
            print(self.centroidList)
            

            # face 遍历2——构造字典存储面及其下的点集合
            for face in self.faceList:
                face_vertices = [v for v in face.loopOuterVertices()]
                print("************************************************************")
                # print("Face: %d" %(i))
                # print("Vertices:", end = "")
                # print(face_vertices)
                # centroid, area = self.calCentroid(face_vertices)
                centroid = self.faceCentroidDict.get(face.identifier)
                radius = self.centroidRadiusDict.get(centroid.identifier)

                face_vertices_dict = self.verticesSort(face_vertices) 
                # 标准圆分割点集合
                vertices_in_circle = self.splitCircle(centroid, radius, len(face_vertices))

                # self.centroidList.append(centroid)

                #
                # self.vertexOnCircleList = vertices_in_circle

                self.faceVerticesDict[face.identifier] = {"cur_vertices":face_vertices_dict, "circle_vertices":vertices_in_circle}
                # gui = pydcel.dcelVis(self)  

            for vertex in self.vertexList:
                edge_list, edge_ientifier_set = self.getAllIncidentEdge(vertex)
                destination_vertices = self.getDestinationVerticesOnCircle(edge_list, vertex)
                
                x_distance, y_distance = self.calNewPosition(vertex, destination_vertices)
                # while not self.checkNewPosition(vertex, x_distance, y_distance, edge_list, edge_ientifier_set):
                #     x_distance = x_distance / 2
                #     y_distance = y_distance / 2
                vertex.x = x_distance + vertex.x
                vertex.y = y_distance + vertex.y


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
                
            
            
            gui = pydcel.dcelVis(self)  
        gui.mainloop()




