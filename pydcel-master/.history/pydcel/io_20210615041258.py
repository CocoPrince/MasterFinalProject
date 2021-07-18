from .dcel import vertex, hedge, face, DCEL
#from sets import Set
from .writegrid2ply import writegrid2ply
from xml.dom.minidom import parse
import xml.etree.ElementTree as tree
import re
from copy import deepcopy as dc
from functools import reduce
import operator
import math
import numpy as np


def ply2datadict(infile):
    """collect vertex coordinates and normals from input file"""
    datadict = {}
    with open(infile) as f:
        vertexcount = facecount = None
        while True:
            line = f.readline()
            if line.startswith("element vertex"):
                vertexcount = line.split()[-1]
            if line.startswith("element face"):
                facecount = line.split()[-1]
            if line.startswith("end_header"):
                break

        datadict['coords'] = []
        datadict['normals'] = []

        for i in range(int(vertexcount)):
            line = f.readline()
            x, y, z = line.split()
            datadict['coords'].append([float(x), float(y), float(z)])

        if facecount is not None:
            datadict['faces'] = []
            for i in range(int(facecount)):
                line = f.readline().split()
                vertex_ids = [int(x) for x in line[1:]]
                datadict['faces'].append(vertex_ids)

    return datadict


def datadict2dcel(datadict):
    # assume ccw vertex order
    hedges = {}  # he_id: (v_origin, v_end), f, nextedge, prevedge
    vertices = {}  # v_id: (e0,...,en) i.e. the edges originating from this v

    m = len(datadict['coords'])
    for i in range(m):
        vertices[i] = []

    # find all halfedges, keep track of their vertices and faces
    j = 0
    for i, face in enumerate(datadict['faces']):
        # face.reverse()
        n_vertices = len(face)

        for v_i in range(n_vertices):
            # store reference to this hedge in vertex list
            vertices[face[v_i]].append(j)

            if v_i == 0:
                hedges[j] = (face[v_i], face[v_i+1]), i, j+1, j+(n_vertices-1)
                vertices[face[v_i+1]].append(j)
            elif v_i < n_vertices-1:
                hedges[j] = (face[v_i], face[v_i+1]), i, j+1, j-1
                vertices[face[v_i+1]].append(j)
            else:
                hedges[j] = (face[v_i], face[0]), i, j-(n_vertices-1), j-1
                vertices[face[0]].append(j)
            vertices[face[v_i]].append(j)
            j += 1

    D = DCEL()

    # create vertices for all points
    for v in datadict['coords']:
        dcel_v = D.createVertex(v[0], v[1], v[2])

    # create faces
    for f in range(len(datadict['faces'])):
        D.createFace()
    # the last face in the DCEL will be the infinite face:
    infinite_face = D.createInfFace()

    # create all edges except for the ones incident to the infinite face
    for e in range(len(hedges)):
        D.createHedge()

    inf_edge = None
    for this_edge, value in hedges.items():
        v, face, nextedge, prevedge = value
        v_origin, v_end = v

        v_origin_edges = set(vertices[v_origin])
        v_end_edges = set(vertices[v_end])

        # print v_origin_edges, v_end_edges
        twin_edge = v_origin_edges.intersection(v_end_edges)
        twin_edge.discard(this_edge)

        e = D.hedgeList[this_edge]

        if len(twin_edge) == 0:  # oh that must be incident to infinite face...
            # face = infinite_face
            e_twin = D.createHedge()
            # oops, forgetting to set something here...
            e_twin.setTopology(
                D.vertexList[v_end], e, infinite_face, None, None)
            inf_edge = e_twin
        else:
            e_twin = D.hedgeList[twin_edge.pop()]
        D.faceList[face].setTopology(e)

        e.setTopology(D.vertexList[v_origin], e_twin, D.faceList[face],
                      D.hedgeList[nextedge], D.hedgeList[prevedge])
        e.origin.setTopology(e)

    # now fix prev/next refs for all edges incident to inf face
    infinite_face.innerComponent = inf_edge
    current_edge = last_correct_edge = inf_edge

    while inf_edge.previous == None:

        current_edge = last_correct_edge
        while current_edge.twin.incidentFace != infinite_face:
            current_edge = current_edge.twin.previous
        current_edge = current_edge.twin

        last_correct_edge.next = current_edge
        current_edge.previous = last_correct_edge
        last_correct_edge = current_edge

    return D


def ply2dcel(infile):
    datadict = ply2datadict(infile)
    return datadict2dcel(datadict)


def dcel2ply(dcel, outfile):
    vertexcount = len(dcel.vertexList)
    facecount = len(dcel.faceList)
    with open(outfile, 'w') as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("comment python generated\n")
        f.write("element vertex {}\n".format(vertexcount))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("element face {}\n".format(facecount))
        f.write("property list uchar int vertex_indices\n")
        f.write("end_header\n")

        for v in dcel.vertexList:
            f.write("{} {} {}\n".format(v.x, v.y, v.z))

        for face in dcel.faceList:  # don't consider inf face
            vertex_list = [dcel.vertexList.index(
                v) for v in face.loopOuterVertices()]
            f.write(str(len(vertex_list)))
            for v_id in vertex_list:
                f.write(" {}".format(v_id))
            f.write("\n")


def xml2ply(xml_path, ply_path):
    graph, start = xml2graph(xml_path)
    faceSetOrigin = set()
    vertexDict = {}
    dfs(graph, [], start, faceSetOrigin, vertexDict)
    vertexcount = len(vertexDict)
    faceSet = set()
    loop_distinct_set_list = []
    for face in faceSetOrigin:  # don't consider inf face
        distinct_set = set()
        vertex_tag_list = face.split(' ')
        flag = True
        for vertex_tag in vertex_tag_list:
            distinct_set.add(vertexDict[vertex_tag])
        for loop_distinct_set in loop_distinct_set_list:
            if distinct_set  == loop_distinct_set:
                flag = False
        # 去重
        if flag :
            faceSet.add(face)
            loop_distinct_set_list.append(distinct_set)
    facecount = len(faceSet)
    with open(ply_path, 'w') as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("comment python generated\n")
        f.write("element vertex {}\n".format(vertexcount))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("element face {}\n".format(facecount))
        f.write("property list uchar int vertex_indices\n")
        f.write("end_header\n")

        for v in vertexDict:
            x_y = v.split("-")
            x = float(x_y[0])
            y = float(x_y[1])
            f.write("{} {} {}\n".format(x, y, 0))

        
        for face in faceSet:  # don't consider inf face
            
            vertex_tag_list = face.split(' ')
            f.write(str(len(vertex_tag_list)))
            
            for vertex_tag in vertex_tag_list:
                f.write(" {}".format(vertexDict[vertex_tag]))
                
            f.write("\n")


def xml2graph(path):
    text = open(path).read()
    text = re.sub(u"[\x00-\x08\x0b-\x0c\x0e-\x1f]+", u"", text)
    domTree = tree.fromstring(text)
    coord_list_all = []
    graph = {}
    coord_index_dict = {}
    # 文档根元素

    # 所有顾客
    page = domTree.find('page')
    ite = page.iter()
    count = 0
    for path in ite:
        if path.tag == 'path':
            vertex_list = []
            coord_list = path.text.split('\n')
            for coord in coord_list:
                if coord != '':
                    x_y = coord.split(" ")
                    vertex_tag = (x_y[0] + '-' + x_y[1])
                    coord_list_all.append(vertex_tag)
                    if len(vertex_list) != 0:
                        key = vertex_list[len(vertex_list)-1]
                        if key in graph:
                            graph[key].append(vertex_tag)
                        else:
                            graph[key] = [vertex_tag]

                        if vertex_tag in graph:
                            graph[vertex_tag].append(key)
                        else:
                            graph[vertex_tag] = [key]

                    vertex_list.append(vertex_tag)
                    coord_index_dict[vertex_tag] = count
                    count += 1
    return graph, coord_list_all[0]


# 用集合去除重复路径
vertex_index = 0
def dfs(graph, trace, start, ans, vertexDict):
    trace = dc(trace)  # 深拷贝，对不同起点，走过的路径不同
    global vertex_index
    # 如果下一个点在trace中，则返回环
    if start in trace:
        index = trace.index(start)
        tmp = [str(i) for i in trace[index:]]
        tmp_set = {str(i) for i in trace[index:]}
        if len(tmp) == 2:
            return
        if check_big_loop(tmp, tmp_set, graph) is False:
            return
        
        tmp_sorted = sort(tmp)

        ans.add(str(' '.join(tmp_sorted)))
        for tag in tmp_sorted:
            if tag not in vertexDict:
                vertexDict[tag] = vertex_index
                vertex_index += 1
        return

    trace.append(start)

    # 深度优先递归递归
    for i in graph[start]:
        dfs(graph, trace, i, ans, vertexDict)


def check_big_loop(loop_list, loop_set, graph):
    for vertex_loop in loop_list:
        nerbor_count = 0
        inside_count = 0
        for neibor_vertex in graph[vertex_loop]:
            x_y = neibor_vertex.split('-')
            if True == ray_tracing_method(float(x_y[0]), float(x_y[1]), loop_list) and neibor_vertex not in loop_set:
                inside_count += 1
            if neibor_vertex in loop_set:
                nerbor_count += 1
        if inside_count > 0 or nerbor_count > 2:
            return False
    return True


def ray_tracing_method(x,y,poly):
    n = len(poly)
    inside = False
    p1x_str,p1y_str = poly[0].split('-')
    p1x,p1y = float(p1x_str), float(p1y_str)
    for i in range(n+1):
        p2x_str,p2y_str = poly[i % n].split('-')
        p2x,p2y = float(p2x_str), float(p2y_str)
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def sort(loop_list):
    coords = []
    for vertex_loop in loop_list:
        x_y = vertex_loop.split("-")
        coords.append((float(x_y[0]), float(x_y[1])))
    
    # 输入的多边形为二维数组(x, y)，多行两列，可以是二维列表，也可以是npumpy二维数组
    # 方法为先确定最高点，则该点一定是凸边的，再按x增大的方向来开始排序
    
    # 1、转换成numpy的array
    polygon_point = np.array(coords)  
    pp = np.array(polygon_point)  
    
    # 2、如果多边形的头尾相连，去掉重复的一个点
    
    if (pp[0] == pp[-1]).all():
        pp = np.delete(pp, -1, axis=0)
    x = pp[:, 0]
    y = pp[:, 1]
    
    # 3、获取最高点位置, np.argmax是默认取到第一个数，所以两边不可能相等
    max_y_index = np.argmax(y)    
    
    # 4、判断前后相邻点在左右的位置，如果最高点在首个，后一个为最后一个，如果最高点在最后一个，则前一个为首个
    pre_index = max_y_index -1
    next_index = 0 if max_y_index == len(pp) - 1 else max_y_index + 1        
 
    # 5、比较前后两个的X大小来确定顺逆时针
    res = polygon_point
    if x[pre_index] > x[next_index]:
        # X是朝大的方向移动，为顺时针
        res = polygon_point
    else:
        # X是朝小的方向移动，为逆时针，重新排序
        res = polygon_point[::-1]

    result_list = []
    x_p = res[:, 0]
    y_p = res[:, 1]
    for x, y in zip(x_p, y_p):
        result_list.append(str(x) + "-" + str(y))
    return result_list
