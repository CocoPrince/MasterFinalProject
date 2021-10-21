import pydcel

# pydcel.io.xml2ply("sampledata/karlsruheBox.xml", "sampledata/karlsruheBox.ply")

d = pydcel.io.ply2dcel('sampledata/montrealBox.ply')
# d.handleFaces('off')
d.handleFaces('on')



# sydney
# montreal
# wien                
# washington
# karlsruhe 

# sydneyBox
# montrealBox
# wienBox                
# washingtonBox
# karlsruheBox 