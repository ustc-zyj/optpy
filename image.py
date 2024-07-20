# -*- coding: utf-8 -*-
from PIL import Image
import numpy as np
import gdspy as gp

def convert_image(
        name,  
        center = (0,0),
        size = (256, 256), 
        pixelsize = 1, 
        threshold = 0.5,
        invert = False,  
        layer = 0,
        datatype = 0):
    ld = {"layer": layer, "datatype": datatype}
    p = pixelsize
    threshold = int(threshold * 256)
    
    im = Image.open(name)
    gray = im.convert("L")
    gray.thumbnail((size[0], size[1]), Image.LANCZOS)
    bw = gray.point(lambda x: 0 if x < threshold else 255, "1")
    pix = np.array(bw.getdata(), bool).reshape(bw.size[1], bw.size[0])
    width, height = bw.size

    polygons = []

    for line in range(height):
        x1 = x0 = 0
        lb = lw = 0
        y0 = (height - line) * p
        for pixel in pix[line]:
            if pixel == invert:
                lb += 1
                if lw > 0:
                    x0 += lw * p
                    lw = 0
            else:
                lw += 1
                if lb > 0:
                    x1 = x0 + lb * p
                    xy = [(x0, y0), (x1, y0), (x1, y0 - p), (x0, y0 - p)]
                    polygons.append(xy)
                    x0 = x1
                    lb = 0
        if lb > 0:
            x1 = x0 + lb * p
            xy = [(x0, y0), (x1, y0), (x1, y0 - p), (x0, y0 - p)]
            polygons.append(gp.Polygon(xy))
    
    image = gp.boolean(polygons, None, "or", **ld)
    box = image.get_bounding_box()
    box_center = ((box[0][0] + box[1][0])/2, (box[0][1] + box[1][1])/2)
    image.translate(center[0]-box_center[0], center[1]-box_center[1])
    return image