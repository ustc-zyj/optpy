# -*- coding: utf-8 -*-

import numpy as np
import gdspy as gp

from optpy.function import(
    get_ld,
    create_thermode,
    )
from optpy.waveguides.basic_waveguide import Basic_Waveguide

class AdjCoupling_Waveguide(Basic_Waveguide):
    
    __slots__ = (
        "nn",
        "rr", 
        "rw", 
        "ww", 
        "cg",
        "ca",
        "al",
        "wy",
        "cx", 
        "cy", 
        "lx",
        "ly",
        "rx",
        "ry",
        "cv",
        "nv",
        "tv"
        )
    
    def __init__(
            self,
            ring_center,
            ring_radius,
            ring_width_nm,
            couple_gap_nm,
            couple_angle_deg,
            endl_xpos,
            endr_xpos,
            waveguide_ypos,
            waveguide_width_nm,
            adjust_length,
            cpara = None,
            npara = None,
            tpara = None):
        
        super().__init__(
            ring_center = ring_center, 
            ring_radius = ring_radius, 
            ring_width_nm = ring_width_nm, 
            couple_gap_nm = couple_gap_nm, 
            couple_angle_deg = couple_angle_deg, 
            endl_xpos = endl_xpos, 
            endr_xpos = endr_xpos, 
            waveguide_width_nm = waveguide_width_nm,
            cpara = cpara,
            npara = npara,
            tpara = tpara)
        
        if adjust_length > 0:
            self.al = adjust_length
        else:
            raise ValueError("Adjust length should be positive")
        
        if waveguide_ypos > ring_center[1] + ring_radius:
            self.wy = waveguide_ypos
        else:
            raise ValueError("Waveguide y-position unsuitable")
            
    def create_waveguide(self, tolerance=1e-3):
        ld_waveguide = get_ld('ld_waveguide')
        waveguides = []
        cr = self.rr + self.cg + self.ww/2
        wx1 = self.cx - cr * (2*np.cos(self.ca/2) - 1)
        wx2 = self.cx + cr * (2*np.cos(self.ca/2) - 1)
        br = (wx2 - wx1) / 2
        waveguide = gp.Path(
            width = self.ww,
            initial_point = (self.lx, self.wy))
        waveguide.segment(
            length = wx1 - self.lx - br,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = br, 
            angle = "r",
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.wy - self.cy - br - cr * 2*np.sin(self.ca/2),
            direction = "-y",
            **ld_waveguide)
        waveguide.turn(
            radius = cr, 
            angle = -self.ca/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr, 
            angle = self.ca,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr, 
            angle = -self.ca/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.al,
            direction = "-y",
            **ld_waveguide)
        waveguide.turn(
            radius = br, 
            angle = "ll",
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.al,
            direction = "+y",
            **ld_waveguide)
        waveguide.turn(
            radius = cr, 
            angle = -self.ca/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr, 
            angle = self.ca,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr, 
            angle = -self.ca/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.wy - self.cy - br - cr * 2*np.sin(self.ca/2),
            direction = "+y",
            **ld_waveguide)
        waveguide.turn(
            radius = br, 
            angle = "r",
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - wx2 - br,
            direction = "+x",
            **ld_waveguide)
        waveguides.append(waveguide)
        return waveguides
    
    def create_thermode(self):
        ld_thermode = get_ld('ld_thermode')
        ld_pad = get_ld('ld_pad')
        thermodes = []
        holepos = []
        cr = self.rr + self.cg + self.ww/2
        wx1 = self.cx - cr * (2*np.cos(self.ca/2) - 1)
        wx2 = self.cx + cr * (2*np.cos(self.ca/2) - 1)
        ty1 = self.cy - self.rr * 2 - cr * 2*np.sin(self.ca/2)
        ty2 = ty1 - self.al + self.rr * 2
        br = (wx2 - wx1) / 2
        if self.tv is not None:
            tw = self.tv["thermode_width"]
            tg = self.tv["thermode_gap"]
            hr = self.tv["hole_radius"]
            hn = self.tv["hole_num"]
            hd = self.tv['hole_distance']
            if self.cg != 0:
                tpara = self.tv.copy()
                thermodes, holepos = create_thermode(
                    ring_center = (self.cx, self.cy),
                    ring_radius = self.rr,
                    **tpara)
            thermode = gp.Rectangle(
                (wx1 - self.ww/2 - tg - tw, ty1), 
                (wx1 - self.ww/2 - tg, ty2),
                **ld_thermode)
            thermodes.append(thermode)
            thermode = gp.Rectangle(
                (wx2 + self.ww/2 + tg + tw, ty1), 
                (wx2 + self.ww/2 + tg, ty2),
                **ld_thermode)
            thermodes.append(thermode)
            thermode = gp.Round(
                center = (self.cx, ty2), 
                radius = br + self.ww/2 + tg + tw,
                inner_radius = br + self.ww/2 + tg,
                initial_angle = np.pi,
                final_angle = np.pi * 2,
                **ld_thermode)
            thermodes.append(thermode)
            arm = gp.Rectangle(
                (wx1 - hd, ty1), 
                (wx1 - self.ww/2 - tg, ty1 - hr*2),
                **ld_thermode)
            thermodes.append(arm)
            arm = gp.Rectangle(
                (wx2 + hd, ty1), 
                (wx2 + self.ww/2 + tg, ty1 - hr*2),
                **ld_thermode)
            thermodes.append(arm)
            hole_pos = [
                (wx1 - hd, ty1 - hr), (wx2 + hd, ty1 - hr)]
            for hp in hole_pos:
                base = gp.Round(
                    center = hp,
                    radius = hr + 0.5,
                    tolerance = 1e-3,
                    **ld_thermode)
                thermodes.append(base)
                for kk in range(hn):
                    hole = gp.Round(
                        center = hp,
                        radius = hr + kk * 0.5,
                        tolerance = 1e-3,
                        **get_ld('ld_hole%1d' %(kk+1)))
                    thermodes.append(hole)
                cap = gp.Round(
                    center = hp,
                    radius = hr + hn * 0.5,
                    tolerance = 1e-3,
                    **ld_pad)
                thermodes.append(cap)
            holepos.extend(hole_pos)
        return thermodes, holepos