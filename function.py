# -*- coding: utf-8 -*-

import numpy as np
import gdspy as gp
import collections

def rotate(point, ref_point, angle_deg):
    angle = angle_deg / 180 * np.pi
    ca = np.cos(angle)
    sa = np.sin(angle)
    dx = point[0] - ref_point[0]
    dy = point[1] - ref_point[1]
    rx = ref_point[0]
    ry = ref_point[1]
    xx = rx + dx * ca - dy * sa
    yy = ry + dx * sa + dy * ca
    return (xx, yy)

def get_ld(ld_key):
    from optpy import layer_datatype_dict
    if ld_key in layer_datatype_dict:
        layer_datatype = layer_datatype_dict[ld_key]
    else:
        raise ValueError('No ' + ld_key + ' in layer-datatype dictionary!')
    return layer_datatype

def create_box(center, size, layer=0, datatype=0):
    ld = {'layer': layer, 'datatype': datatype}
    box = gp.Rectangle((center[0]-size[0]/2, center[1]-size[1]/2), 
                        (center[0]+size[0]/2, center[1]+size[1]/2),
                        **ld)
    return box

def create_posmark(
        center,
        length=50,
        width=2
        ):
    ld_posmark = get_ld('ld_posmark')
    sq1 = create_box(center, (length, width), **ld_posmark)
    sq2 = create_box(center, (width, length), **ld_posmark)
    mark = gp.boolean(sq1, sq2, 'or', **ld_posmark)
    return mark

def create_cutmark(
        center, 
        size=(50, 20)
        ):
    ld_cutmark = get_ld('ld_cutmark')
    mark = create_box(center, size, **ld_cutmark)
    return mark

def create_PhCR(
        center, 
        outer_radius,
        inner_radius,
        iang_ampl,
        iang_freq,
        tolerance,
        layer=0,
        datatype=0
        ):
    ld_ring = {"layer": layer, "datatype": datatype}
    inner_num = int(
        (iang_freq+1)
        * np.pi 
        / np.arccos(1 - tolerance / iang_ampl))
    outer_num = 1 + int(np.pi / np.arccos(1 - tolerance / outer_radius) + 0.5)
    point_num = inner_num + outer_num
    points = np.zeros((point_num, 2))
    theta = np.arange(outer_num) * 2 * np.pi / (outer_num - 1)
    points[:outer_num, 0] = np.cos(theta) * outer_radius + center[0]
    points[:outer_num, 1] = np.sin(theta) * outer_radius + center[0]
    theta = np.arange(inner_num) * -2 * np.pi / (inner_num - 1)
    points[outer_num:, 0] = (
        np.cos(theta) 
        * (inner_radius 
           + iang_ampl * np.cos(iang_freq * theta)) 
        + center[0]
        )
    points[outer_num:, 1] = (
        np.sin(theta) 
        * (inner_radius 
           + iang_ampl * np.cos(iang_freq * theta)) 
        + center[1]
        )
    ring = gp.Polygon(points, **ld_ring)
    return ring

def create_taper(
        posl,
        posr,
        waveguide_width,
        taper_length, 
        tip_length, 
        tip_width_nm, 
        anchor_base,
        anchor_length
        ): 
    ld_taper = get_ld('ld_taper')
    ld_tip = get_ld('ld_tip')
    ld_anchor = get_ld('ld_anchor')
    
    couplers = []
    tip_width = tip_width_nm / 1000
    
    if posl is not None:
        tipl_pos = (posl[0] - taper_length - tip_length, posl[1])
        al_points = [
            (-anchor_length/2, anchor_base/2), (-anchor_length/2, -anchor_base/2),
            (anchor_length/2, -tip_width/2), (anchor_length/2, tip_width/2)]
        al = gp.Polygon(al_points, **ld_anchor)
        al.translate(tipl_pos[0], tipl_pos[1])
        couplers.append(al)
        
        tl = gp.Path(
            width = tip_width,
            initial_point = tipl_pos)
        tl.segment(
            length = tip_length,
            direction = "+x",
            **ld_tip)
        tl.segment(
            length = taper_length,
            direction = "+x",
            final_width = waveguide_width,
            **ld_taper)
        couplers.append(tl)
    
    if posr is not None:   
        tipr_pos = (posr[0] + taper_length + tip_length, posr[1])
        ar_points = [
            (anchor_length/2, anchor_base/2), (anchor_length/2, -anchor_base/2),
            (-anchor_length/2, -tip_width/2), (-anchor_length/2, tip_width/2)]
        ar = gp.Polygon(ar_points, **ld_anchor)
        ar.translate(tipr_pos[0], tipr_pos[1])
        couplers.append(ar)
        
        tr = gp.Path(
            width = tip_width,
            initial_point = tipr_pos)
        tr.segment(
            length = tip_length,
            direction = "-x",
            **ld_tip)
        tr.segment(
            length = taper_length,
            direction = "-x",
            final_width = waveguide_width,
            **ld_taper)
        couplers.append(tr)
    
    return couplers

def create_thermode(
        ring_center,
        ring_radius,
        thermode_type,
        thermode_angle_deg,
        thermode_width,
        thermode_gap,
        hole_radius,
        hole_distance,
        hole_num,
        rotate_angle_deg
        ):    
    def surface_thermode(
            ring_center,
            ring_radius,
            thermode_width,
            thermode_gap,
            hole_radius,
            hole_distance,
            hole_num,
            initial_angle_deg,
            final_angle_deg
            ):
        ld_thermode = get_ld('ld_pad')
        
        initial_angle = initial_angle_deg / 180 * np.pi
        final_angle = final_angle_deg / 180 * np.pi
        hole_pos = (ring_center[0] + hole_distance, ring_center[1])
        thermode = gp.Round(
            center = ring_center, 
            radius = ring_radius + thermode_width / 2 + thermode_gap,
            inner_radius = ring_radius - thermode_width / 2 + thermode_gap,
            initial_angle = initial_angle,
            final_angle = final_angle,
            tolerance = 1e-3,
            max_points = 0,
            **ld_thermode)
        connection = gp.Rectangle(
            point1 = (ring_center[0] + hole_distance, ring_center[1] + hole_radius/2),
            point2 = (ring_center[0] + ring_radius + thermode_gap - thermode_width/2, 
                      ring_center[1] - hole_radius/2),
            **ld_thermode)
        holebase = gp.Round(
            center = hole_pos, 
            radius = hole_radius + 0.5,
            tolerance = 1e-3,
            **ld_thermode)
        arm1 = gp.boolean(connection, holebase, "or", **ld_thermode)
        arm1 = arm1.rotate(initial_angle, ring_center)
        thermode = gp.boolean(thermode, arm1, "or", **ld_thermode)
        arm2 = gp.boolean(connection, holebase, "or", **ld_thermode)
        arm2 = arm2.rotate(final_angle, ring_center)
        thermode = gp.boolean(thermode, arm2, "or", **ld_thermode)
        thermodes = [thermode]
        
        holepos1 = rotate(hole_pos, ring_center, initial_angle_deg)
        holepos2 = rotate(hole_pos, ring_center, final_angle_deg)
        holepos = [holepos2, holepos1]
        
        return thermodes, holepos
    
    def single_thermode(
            ring_center,
            thermode_width,
            thermode_gap,
            ring_radius,
            hole_radius,
            hole_distance,
            hole_num,
            initial_angle_deg,
            final_angle_deg
            ):
        ld_thermode = get_ld('ld_thermode')
        ld_pad = get_ld('ld_pad')
        
        initial_angle = initial_angle_deg / 180 * np.pi
        final_angle = final_angle_deg / 180 * np.pi
        hole_pos = (ring_center[0] + hole_distance, ring_center[1])
        thermode = gp.Round(
            center = ring_center, 
            radius = ring_radius + thermode_gap + thermode_width,
            inner_radius = ring_radius + thermode_gap,
            initial_angle = initial_angle,
            final_angle = final_angle,
            tolerance = 1e-3,
            max_points = 0,
            **ld_thermode)
        connection = gp.Rectangle(
            point1 = (ring_center[0] + hole_distance, ring_center[1] + hole_radius/2),
            point2 = (ring_center[0] + ring_radius + thermode_gap, 
                      ring_center[1] - hole_radius/2),
            **ld_thermode)
        holebase = gp.Round(
            center = hole_pos, 
            radius = hole_radius + 0.3,
            tolerance = 1e-3,
            **ld_thermode)
        arm1 = gp.boolean(connection, holebase, "or", **ld_thermode)
        arm1 = arm1.rotate(initial_angle, ring_center)
        thermode = gp.boolean(thermode, arm1, "or", **ld_thermode)
        arm2 = gp.boolean(connection, holebase, "or", **ld_thermode)
        arm2 = arm2.rotate(final_angle, ring_center)
        thermode = gp.boolean(thermode, arm2, "or", **ld_thermode)
        thermodes = [thermode]
             
        for kk in range(hole_num):
            hole1 = gp.Round(
                center = hole_pos,
                radius = hole_radius + kk * 0.5,
                tolerance = 1e-3,
                **get_ld('ld_hole%1d' %(kk+1)))
            hole1.rotate(initial_angle, ring_center)
            hole2 = gp.Round(
                center = hole_pos,
                radius = hole_radius + kk * 0.5,
                tolerance = 1e-3,
                **get_ld('ld_hole%1d' %(kk+1)))
            hole2.rotate(final_angle, ring_center)
            thermodes.extend([hole1, hole2])
        
        holecap1 = gp.Round(
            center = hole_pos,
            radius = hole_radius + hole_num * 0.5,
            tolerance = 1e-3,
            **ld_pad)
        holecap1.rotate(initial_angle, ring_center)
        holecap2 = gp.Round(
            center = hole_pos,
            radius = hole_radius + hole_num * 0.5,
            tolerance = 1e-3,
            **ld_pad)
        holecap2.rotate(final_angle, ring_center)
        thermodes.extend([holecap1, holecap2])
        
        holepos1 = rotate(hole_pos, ring_center, initial_angle_deg)
        holepos2 = rotate(hole_pos, ring_center, final_angle_deg)
        holepos = [holepos1, holepos2]
        
        return thermodes, holepos
    
    thermodes = []    
    rotate_angle = rotate_angle_deg / 180 * np.pi
    if thermode_type == 'none':
        holepos = []
    elif thermode_type == 'surface':
        thermode, holepos = surface_thermode(
            ring_center = ring_center, 
            thermode_width = thermode_width, 
            thermode_gap = thermode_gap, 
            ring_radius = ring_radius, 
            hole_radius = hole_radius, 
            hole_distance = hole_distance, 
            hole_num = hole_num,
            initial_angle_deg = 90 - thermode_angle_deg / 2, 
            final_angle_deg = 90 + thermode_angle_deg / 2)
        thermodes.extend(thermode)
    elif thermode_type == 'single':
        thermode, holepos = single_thermode(
            ring_center = ring_center, 
            thermode_width = thermode_width, 
            thermode_gap = thermode_gap, 
            ring_radius = ring_radius, 
            hole_radius = hole_radius, 
            hole_distance = hole_distance, 
            hole_num = hole_num,
            initial_angle_deg = 90 - thermode_angle_deg / 2, 
            final_angle_deg = 90 + thermode_angle_deg / 2)
        thermodes.extend(thermode)
    elif thermode_type == 'double':
        thermode1, holepos1 = single_thermode(
            ring_center = ring_center, 
            thermode_width = thermode_width, 
            thermode_gap = thermode_gap, 
            ring_radius = ring_radius, 
            hole_radius = hole_radius, 
            hole_distance = hole_distance, 
            hole_num = hole_num,
            initial_angle_deg = 0 - thermode_angle_deg / 2, 
            final_angle_deg = 0 + thermode_angle_deg / 2)
        thermodes.extend(thermode1)
        thermode2, holepos2 = single_thermode(
            ring_center = ring_center, 
            thermode_width = thermode_width, 
            thermode_gap = thermode_gap, 
            ring_radius = ring_radius, 
            hole_radius = hole_radius, 
            hole_distance = hole_distance, 
            hole_num = hole_num,
            initial_angle_deg = 180 - thermode_angle_deg / 2, 
            final_angle_deg = 180 + thermode_angle_deg / 2,)
        thermodes.extend(thermode2)
        wires = create_wire(
            initial_points = [
                holepos1[0], (holepos1[0][0], holepos1[0][1] - ring_radius)], 
            final_points = [
                (holepos2[1][0], holepos2[1][1] - ring_radius), holepos2[1]], 
            initial_width = hole_radius * 2, 
            final_width = hole_radius * 2, 
            taper_length = 0)
        thermodes.extend(wires)
        if (rotate_angle_deg % 360) <= 90 or (rotate_angle_deg % 360) >= 270:
            holepos = [holepos2[0], holepos1[1]]
        else:
            holepos = [holepos2[1], holepos1[0]]
        
    for thermode in thermodes:
        thermode.rotate(rotate_angle, ring_center)
    for kk, pos in enumerate(holepos):
        holepos[kk] = rotate(pos, ring_center, rotate_angle_deg)
    return thermodes, holepos
    
def create_electrode(
        ring_center,
        ring_radius,
        ring_width,
        electrode_angle_deg,
        electrode_width,
        electrode_gap,
        hole_angle_deg,
        hole_radius,
        hole_num,
        rotate_angle_deg,
        tolerance = 1e-3
        ):
    electrodes = []
    holepos = [
        (ring_center[0], 
         ring_center[1] - ring_radius + ring_width +\
             electrode_gap + electrode_width + hole_radius + 0.5),
        (ring_center[0], 
         ring_center[1] - ring_radius -\
             electrode_gap - electrode_width - hole_radius - 0.5)]
    ld_electrode = get_ld('ld_electrode')
    # ld_holebase = get_ld('ld_holebase')
    ld_pad = get_ld('ld_pad')
    rotate_angle = rotate_angle_deg / 180 * np.pi
    initial_angle = (-90 - electrode_angle_deg / 2) * np.pi / 180
    final_angle = (-90 + electrode_angle_deg / 2) * np.pi / 180
    inner_radius1 = ring_radius - ring_width -\
        electrode_gap - electrode_width - hole_radius*2 - 1
    inner_radius2 = ring_radius + electrode_gap + electrode_width
    angle1 = np.arctan(hole_radius / inner_radius1) * 2
    angle2 = np.arctan(hole_radius / inner_radius2) * 2
    
    electrode1 = gp.Round(
        center = ring_center, 
        radius = ring_radius - ring_width - electrode_gap,
        inner_radius = ring_radius - ring_width - electrode_gap - electrode_width,
        initial_angle = initial_angle,
        final_angle = final_angle,
        tolerance = tolerance,
        **ld_electrode)
    electrodes.append(electrode1)
    
    holebase1 = gp.Round(
        center = ring_center, 
        radius = ring_radius - ring_width - electrode_gap,
        inner_radius = inner_radius1,
        initial_angle = (-np.pi - angle1) / 2,
        final_angle = (-np.pi + angle1) / 2,
        tolerance = tolerance,
        **ld_electrode)
    electrodes.append(holebase1)
    
    electrode2 = gp.Round(
        center = ring_center, 
        radius = ring_radius + electrode_gap + electrode_width,
        inner_radius = ring_radius + electrode_gap,
        initial_angle = initial_angle,
        final_angle = final_angle,
        tolerance = tolerance,
        **ld_electrode)
    electrodes.append(electrode2)
    
    holebase2 = gp.Round(
        center = ring_center, 
        radius = inner_radius2 + hole_radius*2 + 1,
        inner_radius = inner_radius2 - electrode_width,
        initial_angle = (-np.pi - angle2) / 2,
        final_angle = (-np.pi + angle2) / 2,
        tolerance = tolerance,
        **ld_electrode)
    electrodes.append(holebase2)
    
    for hp in holepos:
        for kk in range(hole_num):
            hole = gp.Round(
                center = hp,
                radius = hole_radius + kk * 0.5,
                tolerance = tolerance,
                **get_ld('ld_hole%1d' %(kk+1)))
            electrodes.append(hole)
        holecap = gp.Round(
            center = hp,
            radius = hole_radius + hole_num * 0.5,
            tolerance = tolerance,
            **ld_pad)
        electrodes.append(holecap)
        
    for electrode in electrodes:
        electrode.rotate(rotate_angle, ring_center)
    for kk, pos in enumerate(holepos):
        holepos[kk] = rotate(pos, ring_center, rotate_angle_deg)
    
    return electrodes, holepos

def create_pad(
        pad_center,
        pad_shape = 'round',
        pad_gap = 50,
        pad_size = (200,200),
        port_distance = 50,
        wire_width = 20,
        wire_bias = 70,
        **kwargs
        ):
    ld_pad = get_ld('ld_pad')
    
    def single_pad(pad_center, pad_shape, pad_size, hole_center, wire_width):
        pads = []
        if pad_shape == 'rectangle':
            pad = gp.Rectangle(
                point1 = (pad_center[0] - pad_size[0]/2, pad_center[1] - pad_size[1]/2), 
                point2 = (pad_center[0] + pad_size[0]/2, pad_center[1] + pad_size[1]/2),
                **ld_pad)
        elif pad_shape == 'round':
            pad = gp.Round(pad_center, pad_size[0]/2, **ld_pad)
        else:
            raise ValueError('No pad shape ' + pad_shape)
        pads.append(pad)
        wire = gp.FlexPath(
            points = [
                (hole_center[0], pad_center[1]),
                (hole_center[0], hole_center[1])], 
            width = wire_width,
            corners = 'bevel',
            **ld_pad)
        pads.append(wire.to_polygonset())
        return pads   
        
    if pad_gap > 0:
        pad_center1 = (pad_center[0] - pad_gap/2 - pad_size[0]/2, pad_center[1])
    elif pad_gap == 0:
        pad_center1 = pad_center
    else:
        raise ValueError("Gaps between pads should be positive")
        
    if port_distance > 0:
        port_ypos = pad_center[1] + pad_size[1]/2 + port_distance
    elif port_distance < 0:
        port_ypos = pad_center[1] - pad_size[1]/2 + port_distance
    else:
        port_ypos = pad_center[1]
    
    port_center1 = (pad_center1[0] + wire_bias, port_ypos)   
    portpos = [port_center1]
    
    pads = []
    pad1 = single_pad(pad_center1, pad_shape, pad_size, port_center1, wire_width)
    pads.extend(pad1)
    if pad_gap > 0:
        pad_center2 = (pad_center[0] + pad_gap/2 + pad_size[0]/2, pad_center[1])
        port_center2 = (pad_center2[0] - wire_bias, port_ypos)  
        pad2 = single_pad(pad_center2, pad_shape, pad_size, port_center2, wire_width)
        pads.extend(pad2)
        portpos.append(port_center2)
    return pads, portpos

def create_wire(
        initial_points,
        final_points,
        initial_width,
        final_width,
        taper_length = 100
        ):
    ld_pad = get_ld('ld_pad')
    
    point1 = initial_points[-1]
    point2 = final_points[0]
    dx = point2[0] - point1[0]
    dy = point2[1] - point1[1]
    dl = np.sqrt(dx**2 + dy**2)
    
    wire = gp.FlexPath(
        points = initial_points, 
        width = initial_width,
        corners = 'bevel',
        **ld_pad)   
    if initial_width == final_width:
        wire.segment(point2) 
    elif dl <= taper_length + initial_width + final_width:
        wire.segment(point2, width=final_width)
    elif dl > taper_length*3 + initial_width + final_width:
        ty = dy / dl * taper_length
        tx = dx / dl * taper_length
        point3 = (
            point1[0] + tx, 
            point1[1] + ty)
        point4 = (
            point1[0] + tx*2, 
            point1[1] + ty*2)
        wire.segment(point3)
        wire.segment(point4, width=final_width)
        wire.segment(point2)
    else:
        ty = dy / dl * taper_length / 2
        tx = dx / dl * taper_length / 2
        point3 = (
            (point1[0]+point2[0])/2 - tx, 
            (point1[1]+point2[1])/2 - ty)
        point4 = (
            (point1[0]+point2[0])/2 + tx, 
            (point1[1]+point2[1])/2 + ty)
        wire.segment(point3)
        wire.segment(point4, width=final_width)
        wire.segment(point2)
    for point in final_points[1:]:
        wire.segment(point)
    wires = [wire.to_polygonset()]
    return wires

def create_electest(
        position,
        hole_radius,
        hole_slope,
        hole_num,
        elec_length,
        elec_width,
        wire_length,
        ld_elec,
        tolerance = 1e-3
        ):
    structures = []
    ld_pad = get_ld('ld_pad')
    holepos = [
        (position[0] - elec_length/2, position[1]),
        (position[0] + elec_length/2, position[1])]
    padpos = [
        (position[0] - elec_length/2, position[1] - wire_length),
        (position[0] + elec_length/2, position[1] + wire_length)]
    
    elec = gp.Rectangle(
        (position[0] - elec_length/2, position[1] - elec_width/2), 
        (position[0] + elec_length/2, position[1] + elec_width/2),
        **ld_elec)
    structures.append(elec)
    
    for hp in holepos:
        holebase = gp.Round(
            center = hp,
            radius = hole_radius + hole_slope,
            tolerance = tolerance,
            **ld_elec)
        structures.append(holebase)
        for kk in range(hole_num):
            hole = gp.Round(
                center = hp,
                radius = hole_radius + kk * hole_slope,
                tolerance = tolerance,
                **get_ld('ld_hole%1d' %(kk+1)))
            structures.append(hole)
        holecap = gp.Round(
            center = hp,
            radius = hole_radius + hole_num * hole_slope,
            tolerance = tolerance,
            **ld_pad)
        structures.append(holecap)
    
    for pp in padpos:
        pm = (position[1] - pp[1]) / abs(position[1] - pp[1])
        pad, portpos = create_pad(
            pad_center = pp, 
            pad_shape = 'rectangle', 
            pad_gap = 0, 
            pad_size = (200, 200), 
            port_distance = position[1] - pp[1] - pm * 100, 
            wire_width = (hole_radius + hole_num * hole_slope) * 2, 
            wire_bias = 0)
        structures.extend(pad)
    
    return structures

def bezier_points(bend_radius, direction='l', bezier_parameters=(1, 0.4)):
    t1, t2 = bezier_parameters
    cp4 = (bend_radius, bend_radius)
    if direction == 'l':
        cp1 = (bend_radius * t1 * t2, 0)
        cp2 = (bend_radius * (1 + t1) / 2, bend_radius * (1 - t1) / 2)
        cp3 = (bend_radius, bend_radius * (1 - t1 * t2))
    elif direction == 'r':
        cp1 = (0, bend_radius * t1 * t2)
        cp2 = (bend_radius * (1 - t1) / 2, bend_radius * (1 + t1) / 2)
        cp3 = (bend_radius * (1 - t1 * t2), bend_radius)
    else:
        raise ValueError("Invalid directon input:" + direction)
    control_points = [cp1, cp2, cp3, cp4]
    return control_points

def generate_field(
        xy_range,
        mf_width,
        sf_width,
        mf_size = 500, 
        sf_size = 20, 
        mf_bias = (0, 0), 
        sf_bias = (0, 0),
        mf_ld = {"layer": 0, "datatype": 0},
        sf_ld = {"layer": 0, "datatype": 0}
        ):
    mf_xfield = []
    mf_yfield = []
    sf_xfield = []
    sf_yfield = []
    x_min = min(xy_range[0])
    x_max = max(xy_range[0])
    y_min = min(xy_range[1])
    y_max = max(xy_range[1])
    
    mf_nxmax = (x_max - mf_bias[0]) // mf_size
    mf_nxmin = (x_min - mf_bias[0]) // mf_size
    mf_nymax = (y_max - mf_bias[1]) // mf_size
    mf_nymin = (y_min - mf_bias[1]) // mf_size
    for nn in range(mf_nxmin, mf_nxmax + 1):
        fx = nn * mf_size + mf_bias[0]
        if fx > x_max or fx < x_min:
            continue
        field = gp.Rectangle(
            (fx - mf_width/2, y_min), 
            (fx + mf_width/2, y_max),
            **mf_ld)
        mf_xfield.append(field)
    for nn in range(mf_nymin, mf_nymax + 1):
        fy = nn * mf_size + mf_bias[1]
        if fy > y_max or fy < y_min:
            continue
        field = gp.Rectangle(
            (x_min, fy - mf_width/2), 
            (x_max, fy + mf_width/2),
            **mf_ld)
        mf_yfield.append(field)
    # mf_field = gp.boolean(mf_fields, None, 'or', **mf_ld)
    
    sf_nxmax = (x_max - sf_bias[0]) // sf_size
    sf_nxmin = (x_min - sf_bias[0]) // sf_size
    sf_nymax = (y_max - sf_bias[1]) // sf_size
    sf_nymin = (y_min - sf_bias[1]) // sf_size
    for nn in range(sf_nxmin, sf_nxmax + 1):
        fx = nn * sf_size + sf_bias[0]
        if fx > x_max or fx < x_min:
            continue
        field = gp.Rectangle(
            (fx - sf_width/2, y_min), 
            (fx + sf_width/2, y_max),
            **sf_ld)
        sf_xfield.append(field)
    for nn in range(sf_nymin, sf_nymax + 1):
        fy = nn * sf_size + sf_bias[1]
        if fy > y_max or fy < y_min:
            continue
        field = gp.Rectangle(
            (x_min, fy - sf_width/2), 
            (x_max, fy + sf_width/2),
            **sf_ld)
        sf_yfield.append(field)
    # sf_field = gp.boolean(sf_fields, None, 'or', **sf_ld)
    
    return mf_xfield, mf_yfield, sf_xfield, sf_yfield

def merge_all(cell):
    polygonSet = collections.defaultdict(list)
    for p in cell.polygons:
        for l, d, points in zip(p.layers, p.datatypes, p.polygons):
            polygonSet[(l, d)].append(points)

    cell.polygons = []
    for (l, d), plist in polygonSet.items():
        if len(plist) > 1:
            result = gp.boolean(plist, None, "or", max_points=0, layer=l, datatype=d)
            cell.add(result)
        elif len(plist) == 1:
            cell.add(gp.Polygon(plist[0], layer=l, datatype=d))
    return cell