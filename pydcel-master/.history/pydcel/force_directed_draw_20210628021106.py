
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

    def geneNodeDict(self, centroidList):
        nodeDict = {}
        for vertex in centroidList:
            nodeDict[vertex.identifier] = vertex
    
    def handler(self):
        pass

    def calK(self, start, end):
        return self.centroidRadiusDict.get(start.identifier) + self.centroidRadiusDict.get(end.identifier)
     

