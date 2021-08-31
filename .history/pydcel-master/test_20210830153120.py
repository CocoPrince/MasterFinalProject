import pydcel

# pydcel.io.xml2ply("sampledata/karlsruhe_new.xml", "sampledata/karlsruhe_new.ply")

d = pydcel.io.ply2dcel('sampledata/wien.ply')
# d.handleFaces('off')
d.handleFaces('on')



# sydney
# montreal
# wien                
# washington
# karlsruhe 
# london
# karlsruhe_new