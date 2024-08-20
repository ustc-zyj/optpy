# -*- coding: utf-8 -*-

import numpy as np
import gdspy as gp
import warnings

from optpy.text import create_name

from optpy.function import(
    get_ld,
    create_box,
    create_taper,
    )

class MMI_MZI(object):
    
    __slots__ = (
        "ww",
        "sw",
        "sl",
        "br",
        "cx",
        "cy",
        "lx",
        "rx",
        "bl",
        "bw",
        "pw",
        "pd",
        "tl",
        "cv",
        "nv",
        "ev",
        )
    
    def __init__(
            self,
            mzi_center,
            waveguide_width_nm,
            port_width_nm,
            saparate_length,
            saparate_width,
            endl_xpos,
            endr_xpos,
            bend_radius,
            box_length,
            box_width,
            port_distance,
            taper_length,
            cpara = None,
            npara = None,
            epara = None,
            ):
        
        # Translate waveguide parameters
        if waveguide_width_nm > 0:
            self.ww = waveguide_width_nm / 1000
        else:
            raise ValueError("Waveguide width should be positive")
            
        if port_width_nm > 0:
            self.pw = port_width_nm / 1000
        else:
            raise ValueError("Port width should be positive")
            
        if bend_radius > 0:
            self.br = bend_radius
        else:
            raise ValueError("Bend radius should be positive")
            
        if taper_length > 0:
            self.tl = taper_length
        else:
            raise ValueError("Taper length should be positive")
            
        # Translate box parameters
        if saparate_length > 0:
            self.bl = box_length
        else:
            raise ValueError("Saparate length should be positive")
            
        if saparate_width > 0:
            self.bw = box_width
        else:
            raise ValueError("Saparate width should be positive")
            
        if port_distance - self.pw >= 0:
            self.pd = port_distance
        elif port_distance > 0:
            self.pd = port_distance
            warnings.warn("Ports of the splitter have overlap")
        else:
            raise ValueError("Port distance should be positive")
            
        # Translate saparation parameters
        if saparate_length > 0:
            self.sl = saparate_length
        else:
            raise ValueError("Saparate length should be positive")
            
        if saparate_width > 0:
            self.sw = saparate_width
        else:
            raise ValueError("Saparate width should be positive")
            
        # Translate positions
        self.cx, self.cy = mzi_center
        
        if endl_xpos < self.cx:
            self.lx = endl_xpos
        else:
            raise ValueError("Center of ring should between two ends of waveguide")
        
        if endr_xpos > self.cx:
            self.rx = endr_xpos
        else:
            raise ValueError("Center of ring should between two ends of waveguide")
        
        if self.rx - self.lx < self.sl + self.bl * 2 + self.tl * 4:
            raise ValueError("Not enough length for the splitter")
            
        # Translate coupler variables
        self.cv = cpara
        
        # Translate name variables
        self.nv = npara
        
        # Translate electrode variables
        self.ev = epara
            
    def create_waveguide(self, tolerance=1e-3):
        ld_waveguide = get_ld('ld_waveguide')
        lx2 = self.cx - self.sl/2
        lx1 = lx2 - self.tl*2 - self.bl
        rx2 = self.cx + self.sl/2
        rx1 = rx2 + self.tl*2 + self.bl
        sd = (self.sw - self.pw) / 2
        if sd < self.br * 2:
            theta = np.arccos(1 - sd / (2 * self.br))
            vl = 0
        else:
            theta = np.pi / 2
            vl = sd - self.br*2
        hl = self.sl - self.br * 4 * np.sin(theta)
            
        waveguides = []
        
        waveguide = gp.Path(
            width = self.ww,
            initial_point = (self.lx, self.cy))
        waveguide.segment(
            length = lx1 - self.lx,
            direction = "+x",
            **ld_waveguide)
        waveguide.segment(
            length = self.tl,
            final_width = self.pw,
            **ld_waveguide)
        waveguides.append(waveguide)
        
        box = gp.Rectangle(
            (lx1 + self.tl, self.cy + self.bw/2), 
            (lx2 - self.tl, self.cy - self.bw/2),
            **ld_waveguide)
        waveguides.append(box)
        
        waveguide = gp.Path(
            width = self.pw,
            initial_point = (lx2 - self.tl, self.cy + self.pd/2))
        waveguide.segment(
            length = self.tl,
            direction = "+x",
            final_width = self.ww,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = vl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = -theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = hl,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = -theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = vl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.tl,
            direction = "+x",
            final_width = self.pw,
            **ld_waveguide)
        waveguides.append(waveguide)
        
        waveguide = gp.Path(
            width = self.pw,
            initial_point = (lx2 - self.tl, self.cy - self.pd/2))
        waveguide.segment(
            length = self.tl,
            direction = "+x",
            final_width = self.ww,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = -theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = vl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = hl,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = vl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = -theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.tl,
            direction = "+x",
            final_width = self.pw,
            **ld_waveguide)
        waveguides.append(waveguide)
        
        box = gp.Rectangle(
            (rx1 - self.tl, self.cy + self.bw/2), 
            (rx2 + self.tl, self.cy - self.bw/2),
            **ld_waveguide)
        waveguides.append(box)
        
        waveguide = gp.Path(
            width = self.pw,
            initial_point = (rx1 - self.tl, self.cy))
        waveguide.segment(
            length = self.tl,
            direction = "+x",
            final_width = self.ww,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - rx1,
            **ld_waveguide)
        waveguides.append(waveguide)
        return waveguides
            
    def create_taper(
            self, 
            taper_length = 100, 
            tip_length = 100, 
            tip_width_nm = 100, 
            anchor_base = 4,
            anchor_length = 4):
        couplers = create_taper(
            posl = (self.lx, self.cy), 
            posr = (self.rx, self.cy), 
            waveguide_width = self.ww, 
            taper_length = taper_length, 
            tip_length = tip_length, 
            tip_width_nm = tip_width_nm, 
            anchor_base = anchor_base, 
            anchor_length = anchor_length)
        return couplers
            
    def create_coupler(self, tolerance=1e-3):
        couplers = []
        if self.cv is not None:
            cpara = self.cv.copy()
            couple_type = cpara.pop("couple_type", None)
            if couple_type is None:
                pass
            elif couple_type == "taper":
                couplers = self.create_taper(**cpara)    
            else:
                raise ValueError('Invalid couple type: ' + couple_type)
        return couplers
    
    def create_name(self):
        names = []
        if self.nv is not None:
            npara = self.nv.copy()
            names = create_name(**npara)
        return names
    
    def create_electrode(self):
        electrodes = []
        hpos_list = []
        if self.ev is None:
            return electrodes, hpos_list
        epara = self.ev.copy()
        ld_electrode = get_ld('ld_electrode')
        ld_holebase = get_ld('ld_holebase')
        ld_pad = get_ld('ld_pad')
        
        eg = epara['electrode_gap_nm'] / 1000
        el = epara['electrode_length']
        ew = epara['electrode_width']
        # ew = self.sw - self.ww - eg * 2
        hn = epara['hole_num']
        hr = epara['hole_radius']
        hx0 = epara['center_hole_xpos']
        hx1 = epara['up_hole_xpos']
        hx2 = epara['down_hole_xpos']
        
        ex = self.cx
        hy0 = self.cy
        hy1 = self.cy + self.sw/2 + self.ww/2 + eg + ew + hr
        hy2 = self.cy - self.sw/2 - self.ww/2 - eg - ew - hr
        ey0a = self.cy + self.sw/2 - self.ww/2 - eg - ew/2
        ey0b = self.cy - self.sw/2 + self.ww/2 + eg + ew/2
        ey1 = self.cy + self.sw/2 + self.ww/2 + eg + ew/2
        ey2 = self.cy - self.sw/2 - self.ww/2 - eg - ew/2
        
        hpos_list = [(hx1, hy1), (hx0, hy0), (hx2, hy2)]
        
        electrode = create_box(
            center = (ex, ey1),
            size = (el, ew),
            **ld_electrode)
        electrodes.append(electrode)
        
        electrode = create_box(
            center = (ex, ey0a),
            size = (el, ew),
            **ld_electrode)
        electrodes.append(electrode)
        
        electrode = create_box(
            center = (ex, ey0b),
            size = (el, ew),
            **ld_electrode)
        electrodes.append(electrode)
        
        electrode = create_box(
            center = (ex, ey2),
            size = (el, ew),
            **ld_electrode)
        electrodes.append(electrode)
        
        holebase = create_box(
            center = (hx1, hy1+0.5), 
            size = (hr*2+2, hr*2+1),
            **ld_holebase)
        electrodes.append(holebase)
        
        holebase = create_box(
            center = (hx0, hy0), 
            size = (hr*2+2, self.sw - self.ww - eg*2 - ew*2),
            **ld_holebase)
        electrodes.append(holebase)
        
        holebase = create_box(
            center = (hx2, hy2-0.5), 
            size = (hr*2+2, hr*2+1),
            **ld_holebase)
        electrodes.append(holebase)
        
        for hpos in hpos_list:
            for kk in range(hn):
                hole = gp.Round(
                    center = hpos,
                    radius = hr + kk * 1,
                    tolerance = 1e-3,
                    **get_ld('ld_hole%1d' %(kk+1)))
                electrodes.append(hole)
            holecap = gp.Round(
                center = hpos,
                radius = hr + hn * 1,
                tolerance = 1e-3,
                **ld_pad)
            electrodes.append(holecap)
        return electrodes, hpos_list
    
    def create(self, bias=(0,0), flip=False, tolerance=1e-3):
        structures = []
        holepos = []
        structures.extend(self.create_waveguide(tolerance))
        structures.extend(self.create_coupler(tolerance))
        structures.extend(self.create_name())
        electrodes, holepos = self.create_electrode()
        structures.extend(electrodes)
        if flip:
            for structure in structures:
                structure.mirror((self.cx, self.cy), (self.cx + 1, self.cy))
            for kk, pos in enumerate(holepos):
                holepos[kk] = (2 * self.cx - pos[0], pos[1])
        for structure in structures:
            structure.translate(bias[0], bias[1])
        for kk in range(len(holepos)):
            pos = holepos[kk]
            holepos[kk] = (pos[0] + bias[0], pos[1] + bias[1])
        return structures, holepos
    
    def __call__(self, cell, bias=(0,0), flip=False, tolerance=1e-3):
        structures, holepos = self.create(bias, flip, tolerance)
        cell.add(structures)
        return holepos


class DC_MZI(object):
    
    __slots__ = (
        "ww",
        "sw",
        "sl",
        "br",
        "cx",
        "cy",
        "lx",
        "rx",
        "dl",
        "dg",
        "dw",
        "cv",
        "nv",
        "ev",
        )
    
    def __init__(
            self,
            mzi_center,
            waveguide_width_nm,
            saparate_length,
            saparate_width,
            endl_xpos,
            endr_xpos,
            bend_radius,
            dc_length,
            dc_gap_nm,
            dc_width_nm,
            cpara = None,
            npara = None,
            epara = None,
            ):
        # Translate waveguide parameters
        if waveguide_width_nm > 0:
            self.ww = waveguide_width_nm / 1000
        else:
            raise ValueError("Waveguide width should be positive")
            
        
            
        if bend_radius > 0:
            self.br = bend_radius
        else:
            raise ValueError("Bend radius should be positive")
            
        # Translate DC parameters
        if dc_length > 0:
            self.dl = dc_length
        else:
            raise ValueError("DC length should be positive")
            
        if dc_gap_nm > 0:
            self.dg = dc_gap_nm / 1000
        else:
            raise ValueError("DC gap should be positive")
            
        if dc_width_nm > 0:
            self.dw = dc_width_nm / 1000
        else:
            raise ValueError("DC width should be positive")    
            
        # Translate saparation parameters
        if saparate_length > 0:
            self.sl = saparate_length
        else:
            raise ValueError("Saparate length should be positive")
            
        if saparate_width > 0:
            self.sw = saparate_width
        else:
            raise ValueError("Saparate width should be positive")
            
        # Translate positions
        self.cx, self.cy = mzi_center
        
        if endl_xpos < self.cx:
            self.lx = endl_xpos
        else:
            raise ValueError("Center of ring should between two ends of waveguide")
        
        if endr_xpos > self.cx:
            self.rx = endr_xpos
        else:
            raise ValueError("Center of ring should between two ends of waveguide")
        
        if self.rx - self.lx < self.sl + self.dl * 2:
            raise ValueError("Not enough length for the splitter")
            
        # Translate coupler variables
        self.cv = cpara
        
        # Translate name variables
        self.nv = npara
        
        # Translate electrode variables
        self.ev = epara
        
    def create_waveguide(self, tolerance=1e-3):
        ld_waveguide = get_ld('ld_waveguide')
        lx2 = self.cx - self.sl/2
        rx2 = self.cx + self.sl/2
        sd = (self.sw - self.dg - self.dw) / 2
        if sd < self.br * 2:
            theta = np.arccos(1 - sd / (2 * self.br))
            vl = 0
        else:
            theta = np.pi / 2
            vl = sd - self.br*2
        dcl = self.br * 4 * np.sin(theta) + self.dl
        lx1 = lx2 - dcl
        rx1 = rx2 + dcl
        
        waveguides = []
        
        waveguide = gp.Path(
            width = self.ww,
            initial_point = (self.lx, self.cy + self.sw / 2))
        waveguide.segment(
            length = lx1 - self.lx,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = -theta,
            tolerance = tolerance,
            final_width = (self.ww + self.dw) / 2,
            **ld_waveguide)
        waveguide.segment(
            length = vl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = theta,
            tolerance = tolerance,
            final_width = self.dw,
            **ld_waveguide)
        waveguide.segment(
            length = self.dl,
            direction = '+x',
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = theta,
            tolerance = tolerance,
            final_width = (self.ww + self.dw) / 2,
            **ld_waveguide)
        waveguide.segment(
            length = vl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = -theta,
            tolerance = tolerance,
            final_width = self.ww,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl,
            direction = '+x',
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = -theta,
            tolerance = tolerance,
            final_width = (self.ww + self.dw) / 2,
            **ld_waveguide)
        waveguide.segment(
            length = vl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = theta,
            tolerance = tolerance,
            final_width = self.dw,
            **ld_waveguide)
        waveguide.segment(
            length = self.dl,
            direction = '+x',
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = theta,
            tolerance = tolerance,
            final_width = (self.ww + self.dw) / 2,
            **ld_waveguide)
        waveguide.segment(
            length = vl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = -theta,
            tolerance = tolerance,
            final_width = self.ww,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - rx1,
            direction = "+x",
            **ld_waveguide)
        waveguides.append(waveguide)
        
        waveguide = gp.Path(
            width = self.ww,
            initial_point = (self.lx, self.cy - self.sw / 2))
        waveguide.segment(
            length = lx1 - self.lx,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = theta,
            tolerance = tolerance,
            final_width = (self.ww + self.dw) / 2,
            **ld_waveguide)
        waveguide.segment(
            length = vl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = -theta,
            tolerance = tolerance,
            final_width = self.dw,
            **ld_waveguide)
        waveguide.segment(
            length = self.dl,
            direction = '+x',
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = -theta,
            tolerance = tolerance,
            final_width = (self.ww + self.dw) / 2,
            **ld_waveguide)
        waveguide.segment(
            length = vl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = theta,
            tolerance = tolerance,
            final_width = self.ww,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl,
            direction = '+x',
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = theta,
            tolerance = tolerance,
            final_width = (self.ww + self.dw) / 2,
            **ld_waveguide)
        waveguide.segment(
            length = vl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = -theta,
            tolerance = tolerance,
            final_width = self.dw,
            **ld_waveguide)
        waveguide.segment(
            length = self.dl,
            direction = '+x',
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = -theta,
            tolerance = tolerance,
            final_width = (self.ww + self.dw) / 2,
            **ld_waveguide)
        waveguide.segment(
            length = vl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = theta,
            tolerance = tolerance,
            final_width = self.ww,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - rx1,
            direction = "+x",
            **ld_waveguide)
        waveguides.append(waveguide)
        
        return waveguides
    
    def create_taper(
            self, 
            taper_length = 100, 
            tip_length = 100, 
            tip_width_nm = 100, 
            anchor_base = 4,
            anchor_length = 4):
        couplers = []
        coupler = create_taper(
            posl = (self.lx, self.cy + self.sw/2), 
            posr = (self.rx, self.cy + self.sw/2), 
            waveguide_width = self.ww, 
            taper_length = taper_length, 
            tip_length = tip_length, 
            tip_width_nm = tip_width_nm, 
            anchor_base = anchor_base, 
            anchor_length = anchor_length)
        couplers.extend(coupler)
        
        coupler = create_taper(
            posl = (self.lx, self.cy - self.sw/2), 
            posr = (self.rx, self.cy - self.sw/2), 
            waveguide_width = self.ww, 
            taper_length = taper_length, 
            tip_length = tip_length, 
            tip_width_nm = tip_width_nm, 
            anchor_base = anchor_base, 
            anchor_length = anchor_length)
        couplers.extend(coupler)
        return couplers
    
    def create_coupler(self, tolerance=1e-3):
        couplers = []
        if self.cv is not None:
            cpara = self.cv.copy()
            couple_type = cpara.pop("couple_type", None)
            if couple_type is None:
                pass
            elif couple_type == "taper":
                couplers = self.create_taper(**cpara)    
            else:
                raise ValueError('Invalid couple type: ' + couple_type)
        return couplers
    
    def create_name(self):
        names = []
        if self.nv is not None:
            npara = self.nv.copy()
            names = create_name(**npara)
        return names
    
    def create_electrode(self):
        electrodes = []
        hpos_list = []
        if self.ev is None:
            return electrodes, hpos_list
        epara = self.ev.copy()
        ld_electrode = get_ld('ld_electrode')
        ld_pad = get_ld('ld_pad')
        
        eg = epara['electrode_gap_nm'] / 1000
        el = epara['electrode_length']
        ew = self.sw - self.ww - eg * 2
        ex = self.cx
        ey0 = self.cy
        ey1 = self.cy + self.sw
        ey2 = self.cy - self.sw
        
        hn = epara['hole_num']
        hr = epara['hole_radius']
        hx0 = epara['center_hole_xpos']
        hx1 = epara['up_hole_xpos']
        hx2 = epara['down_hole_xpos']
        hpos_list = [(hx1, ey1), (hx0, ey0), (hx2, ey2)]
        
        for ey in [ey1, ey0, ey2]:
            electrode = create_box((ex, ey), (el, ew), **ld_electrode)
            electrodes.append(electrode)
        for hpos in hpos_list:
            for kk in range(hn):
                hole = gp.Round(
                    center = hpos,
                    radius = hr + kk * 1,
                    tolerance = 1e-3,
                    **get_ld('ld_hole%1d' %(kk+1)))
                electrodes.append(hole)
            holecap = gp.Round(
                center = hpos,
                radius = hr + hn * 1,
                tolerance = 1e-3,
                **ld_pad)
            electrodes.append(holecap)
        return electrodes, hpos_list
    
    def create(self, bias=(0,0), flip=False, tolerance=1e-3):
        structures = []
        structures.extend(self.create_waveguide(tolerance))
        structures.extend(self.create_coupler(tolerance))
        structures.extend(self.create_name())
        electrodes, holepos = self.create_electrode()
        structures.extend(electrodes)
        if flip:
            for structure in structures:
                structure.mirror((self.cx, self.cy), (self.cx + 1, self.cy))
            for kk, pos in enumerate(holepos):
                holepos[kk] = (2 * self.cx - pos[0], pos[1])
        for structure in structures:
            structure.translate(bias[0], bias[1])
        for kk in range(len(holepos)):
            pos = holepos[kk]
            holepos[kk] = (pos[0] + bias[0], pos[1] + bias[1])
        return structures, holepos
    
    def __call__(self, cell, bias=(0,0), flip=False, tolerance=1e-3):
        structures, holepos = self.create(bias, flip, tolerance)
        cell.add(structures)
        return holepos