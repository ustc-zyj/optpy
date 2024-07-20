# -*- coding: utf-8 -*-

from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
import gdspy as gp
from optpy.function import get_ld

def convert_text(
        text, 
        text_pos = (0, 0),
        text_size = 64, 
        font_family = None,
        font_style = None, 
        tolerance = 0.1,
        layer = 0,
        datatype = 0):
    font_prop = FontProperties(family=font_family, style=font_style)
    ld = {"layer": layer, "datatype": datatype}
    path = TextPath((0, 0), text, size=text_size, prop=font_prop)
    polys = []
    xmax = text_pos[0]
    for points, code in path.iter_segments():
        if code == path.MOVETO:
            c = gp.Curve(*points, tolerance=tolerance)
        elif code == path.LINETO:
            c.L(*points)
        elif code == path.CURVE3:
            c.Q(*points)
        elif code == path.CURVE4:
            c.C(*points)
        elif code == path.CLOSEPOLY:
            poly = c.get_points()
            if poly.size > 0:
                if poly[:, 0].min() < xmax:
                    i = len(polys) - 1
                    while i >= 0:
                        if gp.inside(
                            poly[:1], [polys[i]], precision=0.1 * tolerance
                        )[0]:
                            p = polys.pop(i)
                            poly = gp.boolean(
                                [p],
                                [poly],
                                "xor",
                                precision=0.1 * tolerance,
                                max_points=0,
                            ).polygons[0]
                            break
                        elif gp.inside(
                            polys[i][:1], [poly], precision=0.1 * tolerance
                        )[0]:
                            p = polys.pop(i)
                            poly = gp.boolean(
                                [p],
                                [poly],
                                "xor",
                                precision=0.1 * tolerance,
                                max_points=0,
                            ).polygons[0]
                        i -= 1
                xmax = max(xmax, poly[:, 0].max())
                polys.append(poly)
    text = gp.PolygonSet(polys, **ld)
    box = text.get_bounding_box()
    box_pos = (box[0][0], (box[0][1] + box[1][1])/2)
    text.translate(text_pos[0]-box_pos[0], text_pos[1]-box_pos[1])
    return text

def create_name(
        name_pos,
        name_text,
        name_size = 64,
        ):
    ld_text = get_ld('ld_text')
    names = []
    name = convert_text(
        text = name_text,
        text_pos = name_pos,
        text_size = name_size,
        font_family = 'serif',
        font_style = 'normal',
        **ld_text)
    names.append(name)
    return names