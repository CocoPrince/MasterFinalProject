import pydcel

# pydcel.io.xml2ply("sampledata/washington.xml", "sampledata/washington.ply")

d = pydcel.io.ply2dcel('sampledata/wien.ply')
# d.handleFaces('off')
d.handleFaces('on')


# sydney
# montreal
# wien                
# washington
# karlsruhe
# london1