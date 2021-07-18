import pydcel

# pydcel.io.xml2ply("sampledata/washington.xml", "sampledata/washington.ply")

d = pydcel.io.ply2dcel('sampledata/sydney.ply')
#d.handleFaces(False)
d.handleFaces(True)


# sydney
# montreal
# wien                

# washington
# karlsruhe
# london