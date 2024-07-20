# -*- coding: utf-8 -*-

import numpy as np
import gdspy as gp
import warnings

from optpy.function import(
    get_ld,
    create_taper
    )

from optpy.waveguides.dual_waveguides import BasicDual_Waveguide

#%% Roller Waveguide
class Roller_Waveguide(BasicDual_Waveguide):
    
    __slots__ = (
        "nn",
        "rr", 
        "rw", 
        "ww", 
        "ww2",
        "br",
        "br2",
        "cg",
        "cg2",
        "ca",
        "ca2",
        "rl",
        "sw",
        "sg",
        "sl",
        "cx", 
        "cy", 
        "lx",
        "rx",
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
            couple_gap2_nm,
            couple_angle_deg,
            couple_angle2_deg, 
            endl_xpos, 
            endr_xpos, 
            waveguide_width_nm,
            waveguide_width2_nm,
            reverse_length,
            bend_radius,
            bend_radius2,
            splitter_width_nm,
            splitter_gap_nm,
            splitter_length,
            cpara = None,
            npara = None,
            tpara = None
            ):
        
        super().__init__(
            ring_center = ring_center, 
            ring_radius = ring_radius, 
            ring_width_nm = ring_width_nm, 
            couple_gap_nm = couple_gap_nm, 
            couple_gap2_nm = couple_gap2_nm, 
            couple_angle_deg = couple_angle_deg,
            couple_angle2_deg = couple_angle2_deg,
            endl_xpos = endl_xpos, 
            endr_xpos = endr_xpos, 
            waveguide_width_nm = waveguide_width_nm,
            waveguide_width2_nm = waveguide_width2_nm,
            cpara = cpara,
            npara = npara,
            tpara = tpara)
        
        # Translate bend radius
        if bend_radius > ring_radius:
            self.br = bend_radius
        else:
            raise ValueError("Bend radius should be larger")
        
        if bend_radius2 > ring_radius:
            self.br2 = bend_radius2
        else:
            raise ValueError("Bend radius 2 should be larger")
        
        # Translate reverse length
        if reverse_length > self.rr * 2:
            self.rl = reverse_length
        elif reverse_length <=0:
            raise ValueError("Reverse length should be positive")
        else:
            raise ValueError("Reverse length should be larger")
            
        # Translate splitter width
        if splitter_width_nm > 0:
            self.sw = splitter_width_nm / 1000
        else:
            raise ValueError("Splitter width should be positive")
        
        # Translate splitter gap
        if splitter_gap_nm > 0:
            self.sg = splitter_gap_nm / 1000
        else:
            raise ValueError("Splitter gap should be positive")
        
        # Translate splitter length
        if splitter_length > 0:
            self.sl = splitter_length
        else:
            raise ValueError("Splitter length should be positive") 
            
    def create_waveguide(self, tolerance=1e-3):
        ld_waveguide = get_ld('ld_waveguide')
        cr1 = self.rr + self.cg + self.ww/2
        cr2 = self.rr + self.cg2 + self.ww2/2
        wy1 = self.cy + cr1 * (2*np.cos(self.ca/2) -1)
        wy2 = self.cy - self.br2 * 2 - self.rl - cr2 * (2*np.cos(self.ca2/2) -1)
        
        # main waveguide
        waveguides = []
        ly1 = wy2 + self.sw + self.sg
        waveguide = gp.Path(
            width = self.ww,
            initial_point = (self.lx, ly1))
        waveguide.segment(
            length = self.br2,
            direction = "+x",
            **ld_waveguide)
        waveguide.segment(
            length = self.br2 * 2,
            final_width = self.sw,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'l',
            tolerance = tolerance,
            final_width = (self.sw + self.ww)/2,
            **ld_waveguide)
        waveguide.segment(
            length = wy1 - ly1 - self.br * 2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'r',
            tolerance = tolerance,
            final_width = self.ww,
            **ld_waveguide)
        waveguide.segment(
            length = self.cx - self.lx - self.br2 * 3 - self.br * 2 -\
                self.sl - 2 * cr1 * np.sin(self.ca/2),
            **ld_waveguide)
        waveguide.turn(
            radius = self.rr + self.cg + self.ww/2, 
            angle = self.ca/2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.turn(
            radius = self.rr + self.cg + self.ww/2, 
            angle = -self.ca,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.turn(
            radius = self.rr + self.cg + self.ww/2, 
            angle = self.ca/2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - self.cx - self.br * 3 - self.br2 * 2 -\
                self.sl - 2 * cr1 * np.sin(self.ca/2),
            direction = "+x",
            **ld_waveguide)
        waveguide.segment(
            length = self.br2 * 2,
            final_width = self.sw,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'l',
            tolerance = tolerance,
            final_width = (self.sw + self.ww)/2,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'r',
            tolerance = tolerance,
            final_width = self.ww,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            direction = "+x",
            **ld_waveguide)
        waveguides.append(waveguide)
        
        # auxillary waveguide
        ly2 = wy2 - self.br2 * 3
        waveguide = gp.Path(
            width = self.ww2,
            initial_point = (self.lx, ly2))
        waveguide.segment(
            length = self.br2,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br2,
            angle = "l",
            tolerance = tolerance,
            final_width = (self.sw + self.ww2)/2,
            **ld_waveguide)
        waveguide.segment(
            length = self.br2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br2,
            angle = "r",
            tolerance = tolerance,
            final_width = self.sw,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl,
            **ld_waveguide)
        waveguide.segment(
            length = self.br * 2,
            final_width = self.ww2,
            **ld_waveguide)
        waveguide.segment(
            length = self.cx - self.lx - self.br2 * 3 - self.br*2 \
                - self.sl - self.rl/2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br2,
            angle = np.pi / 4,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl * np.sqrt(2),
            **ld_waveguide)
        waveguide.turn(
            radius = self.br2,
            angle = np.pi * 3 / 4,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl/2 - cr2 * 2*np.sin(self.ca2/2),
            direction = "-x",
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = -self.ca2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl/2 - cr2 * 2*np.sin(self.ca2/2),
            direction = "-x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br2,
            angle = np.pi * 3 / 4,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl * np.sqrt(2),
            **ld_waveguide)
        waveguide.turn(
            radius = self.br2,
            angle = np.pi / 4,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - self.cx - self.br * 3 - self.br2 * 2 \
                - self.sl - self.rl/2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br2,
            angle = "l",
            tolerance = tolerance,
            final_width = (self.sw + self.ww2)/2,
            **ld_waveguide)
        waveguide.segment(
            length = wy1 - wy2 - self.br2 * 2 - self.sw - self.sg,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br2,
            angle = "r",
            tolerance = tolerance,
            final_width = self.sw,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl,
            **ld_waveguide)
        waveguide.segment(
            length = self.br * 2,
            final_width = self.ww2,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            direction = "+x",
            **ld_waveguide)
        waveguides.append(waveguide)
        
        return waveguides
    
    def create_taper(
            self, 
            taper_length = 100, 
            tip_length = 100, 
            tip_width_nm = 100, 
            tip_width2_nm = 100,
            anchor_base = 4,
            anchor_length = 4):         
        cr1 = self.rr + self.cg + self.ww/2
        cr2 = self.rr + self.cg2 + self.ww2/2
        wy1 = self.cy + cr1 * (2*np.cos(self.ca/2) -1)
        wy2 = self.cy - self.br2 * 2 - self.rl - cr2 * (2*np.cos(self.ca2/2) -1)
        ly1 = wy2 + self.sw + self.sg
        ly2 = wy2 - self.br2 * 3
        ry1 = wy1 + self.br * 3
        ry2 = wy1 - self.ww/2 - self.ww2/2 - self.sg
            
        couplers = create_taper(
            posl = (self.lx, ly1), 
            posr = (self.rx, ry1), 
            waveguide_width = self.ww, 
            taper_length = taper_length, 
            tip_length = tip_length, 
            tip_width_nm = tip_width_nm, 
            anchor_base = anchor_base, 
            anchor_length = anchor_length)
        
        couplers.extend(
            create_taper(
                posl = (self.lx, ly2), 
                posr = (self.rx, ry2), 
                waveguide_width = self.ww2, 
                taper_length = taper_length, 
                tip_length = tip_length, 
                tip_width_nm = tip_width2_nm, 
                anchor_base = anchor_base, 
                anchor_length = anchor_length))
        
        return couplers

#%% Roller MultiRing Waveguide
class RollerMR_Waveguide(BasicDual_Waveguide):
    
    __slots__ = (
        "nn",
        "rr", 
        "rw", 
        "ww", 
        "ww2",
        "br",
        "br2",
        "cg",
        "cg2",
        "ca",
        "ca2",
        "cx", 
        "cy", 
        "lx",
        "rx",
        "cv",
        "nv",
        "tv"
        )
    
    def __init__(
            self, 
            ring_xpos,
            ring_ypos,
            ring_radius, 
            ring_width_nm, 
            couple_gap_nm,
            couple_gap2_nm,
            couple_angle_deg,
            couple_angle2_deg, 
            endl_xpos, 
            endr_xpos, 
            waveguide_width_nm,
            waveguide_width2_nm,
            bend_radius,
            bend_radius2,
            cpara = None,
            npara = None,
            ):
        
        for para in [ring_xpos, ring_width_nm]:
            if not isinstance(para, list):
                raise ValueError('Lists should be used for ring parameters')
                
        super().__init__(
            ring_center = (ring_xpos[0], ring_ypos), 
            ring_radius = ring_radius, 
            ring_width_nm = ring_width_nm[0], 
            couple_gap_nm = couple_gap_nm, 
            couple_gap2_nm = couple_gap2_nm, 
            couple_angle_deg = couple_angle_deg,
            couple_angle2_deg = couple_angle2_deg,
            endl_xpos = endl_xpos, 
            endr_xpos = endr_xpos, 
            waveguide_width_nm = waveguide_width_nm,
            waveguide_width2_nm = waveguide_width2_nm,
            cpara = cpara,
            npara = npara,
            tpara = None)
        
        if len(ring_width_nm) != len(ring_xpos):
            raise ValueError('Length of parameter list does not match')
        else:
            self.nn = len(ring_xpos)
            
        if any([ring_xpos[k] >= ring_xpos[k+1] for k in range(self.nn-1)]):
            raise ValueError("X-position of rings are not sorted")
        elif any([ring_xpos[k] >= ring_xpos[k+1] - self.rr * 2 \
                  for k in range(self.nn-1)]):
            raise ValueError("Some rings have overlap")
        else:
            self.cx = ring_xpos
        
        if all([rw > 0 and rw < self.rr * 1000 for rw in ring_width_nm]):
            self.rw = [rw / 1000 for rw in ring_width_nm]
        elif any([rw == self.rr * 1000 for rw in ring_width_nm]):
            warnings.warn("Creating a disk rather than a ring")
            self.rw = [rw / 1000 for rw in ring_width_nm]
        elif any([rw > self.rr * 1000 for rw in ring_width_nm]):
            raise ValueError("Width of rings should be smaller than radius")
        else:
            raise ValueError("Width of rings should be positive")
            
        if bend_radius >= ring_radius:
            self.br = bend_radius
        else:
            raise ValueError("Bending radius should be larger")
        if bend_radius2 >= ring_radius:
            self.br2 = bend_radius2
        else:
            raise ValueError("Bending radius 2 should be larger")
            
    def create_ring(self, tolerance=1e-4):
        ld_ring = get_ld('ld_ring')
        rings = []
        for kk in range(self.nn):
            ring = gp.Round(
                center = (self.cx[kk], self.cy), 
                radius = self.rr,
                inner_radius = self.rr - self.rw[kk],
                tolerance = tolerance,
                max_points = 0,
                **ld_ring)
            rings.append(ring)
        return rings
    
    def create_waveguide(self, tolerance=1e-3):
        ld_waveguide = get_ld('ld_waveguide')
        waveguides = []
        cr1 = self.rr + self.cg + self.ww/2
        cr2 = self.rr + self.cg2 + self.ww2/2
        wy1 = self.cy + cr1 * (2*np.cos(self.ca/2) -1)
        wy2 = self.cy - cr2 * (2*np.cos(self.ca2/2) -1)
        ly1 = wy1 - self.br * 3
        ly2 = wy2 - self.br2 * 4 - self.br * 2
        
        # main waveguide
        waveguide = gp.Path(
            width = self.ww,
            initial_point = (self.lx, ly1))
        waveguide.segment(
            length = self.br,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "l",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "r",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.cx[0] - self.lx - self.br * 3 -\
                cr1 * 2 * np.sin(self.ca/2),
            **ld_waveguide)
        waveguide.turn(
            radius = cr1, 
            angle = self.ca/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr1, 
            angle = -self.ca,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr1, 
            angle = self.ca/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        for kk in range(1, self.nn):
            waveguide.segment(
                length = self.cx[kk] - self.cx[kk-1] - cr1 * 4*np.sin(self.ca/2),
                direction = "+x",
                **ld_waveguide)
            waveguide.turn(
                radius = cr1, 
                angle = self.ca/2,
                tolerance = tolerance * 0.1,
                **ld_waveguide)
            waveguide.turn(
                radius = cr1, 
                angle = -self.ca,
                tolerance = tolerance * 0.1,
                **ld_waveguide)
            waveguide.turn(
                radius = cr1, 
                angle = self.ca/2,
                tolerance = tolerance * 0.1,
                **ld_waveguide)
        waveguide.segment(
            length = self.rx - self.cx[self.nn-1] - self.br * 4 -\
                cr1 * 2 * np.sin(self.ca/2),
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "l",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "r",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br * 2,
            **ld_waveguide)
        waveguides.append(waveguide)
        
        # Auxilary waveguide
        waveguide = gp.Path(
            width = self.ww2,
            initial_point = (self.lx, ly2))
        waveguide.segment(
            length = self.br * 2,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "l",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "r",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.cx[self.nn-1] - self.lx - self.br * 4 + self.br2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br2, 
            angle = "l",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br2 * 2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br2, 
            angle = "l",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br2 - cr2 * 2*np.sin(self.ca2/2),
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = -self.ca2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        for kk in range(1, self.nn):
            waveguide.segment(
                length = self.cx[self.nn - kk] - self.cx[self.nn - kk-1] -\
                    cr2 * 4*np.sin(self.ca2/2),
                direction = "-x",
                **ld_waveguide)
            waveguide.turn(
                radius = cr2, 
                angle = self.ca2/2,
                tolerance = tolerance * 0.1,
                **ld_waveguide)
            waveguide.turn(
                radius = cr2, 
                angle = -self.ca2,
                tolerance = tolerance * 0.1,
                **ld_waveguide)
            waveguide.turn(
                radius = cr2, 
                angle = self.ca2/2,
                tolerance = tolerance * 0.1,
                **ld_waveguide)
        waveguide.segment(
            length = self.br2 - cr2 * 2*np.sin(self.ca2/2),
            direction = "-x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br2, 
            angle = "ll",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - self.cx[0] + self.br2 - self.br * 3,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "l",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "r",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            **ld_waveguide)
        waveguides.append(waveguide)
        return waveguides
    
    def create_taper(
            self, 
            taper_length = 100, 
            tip_length = 100, 
            tip_width_nm = 100,
            tip_width2_nm = 100,
            anchor_base = 4,
            anchor_length = 4):
        cr1 = self.rr + self.cg + self.ww/2
        cr2 = self.rr + self.cg2 + self.ww2/2
        wy1 = self.cy + cr1 * (2*np.cos(self.ca/2) -1)
        wy2 = self.cy - cr2 * (2*np.cos(self.ca2/2) -1)
        ly1 = wy1 - self.br * 3
        ly2 = wy2 - self.br2 * 4 - self.br * 2
        ry1 = wy1 + self.br * 3
        ry2 = wy2 - self.br2 * 2 + self.br * 3
        
        couplers = create_taper(
            posl = (self.lx, ly1), 
            posr = (self.rx, ry1), 
            waveguide_width = self.ww, 
            taper_length = taper_length, 
            tip_length = tip_length, 
            tip_width_nm = tip_width_nm, 
            anchor_base = anchor_base, 
            anchor_length = anchor_length)
        
        couplers.extend(
            create_taper(
                posl = (self.lx, ly2), 
                posr = (self.rx, ry2), 
                waveguide_width = self.ww2, 
                taper_length = taper_length, 
                tip_length = tip_length, 
                tip_width_nm = tip_width2_nm, 
                anchor_base = anchor_base, 
                anchor_length = anchor_length))
        
        return couplers
#%% Balanced MZ-Roller Waveguide
class MZRoller_Waveguide(BasicDual_Waveguide):
    
    __slots__ = (
        "nn",
        "rr", 
        "rw", 
        "ww", 
        "ww2",
        "br",
        "cg",
        "cg2",
        "ca",
        "ca2",
        "rl",
        "pb",
        "sw",
        "sg",
        "sl",
        "cx", 
        "cy", 
        "lx",
        "rx",
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
            couple_gap2_nm,
            couple_angle_deg,
            couple_angle2_deg, 
            endl_xpos, 
            endr_xpos, 
            waveguide_width_nm,
            waveguide_width2_nm,
            reverse_length,
            bend_radius,
            port_bias,
            splitter_width_nm,
            splitter_gap_nm,
            splitter_length,
            cpara = None,
            npara = None,
            tpara = None
            ):
        
        super().__init__(
            ring_center = ring_center, 
            ring_radius = ring_radius, 
            ring_width_nm = ring_width_nm, 
            couple_gap_nm = couple_gap_nm, 
            couple_gap2_nm = couple_gap2_nm, 
            couple_angle_deg = couple_angle_deg,
            couple_angle2_deg = couple_angle2_deg,
            endl_xpos = endl_xpos, 
            endr_xpos = endr_xpos, 
            waveguide_width_nm = waveguide_width_nm,
            waveguide_width2_nm = waveguide_width2_nm,
            cpara = cpara,
            npara = npara,
            tpara = tpara)
        
        # Translate bend radius
        if bend_radius >= ring_radius:
            self.br = bend_radius
        else:
            raise ValueError("Bend radius should be larger")
        
        # Translate reverse length
        if reverse_length > self.rr * 2:
            self.rl = reverse_length
        elif reverse_length <=0:
            raise ValueError("Reverse length should be positive")
        else:
            raise ValueError("Reverse length should be larger")
            
        # Translate splitter width
        if splitter_width_nm > 0:
            self.sw = splitter_width_nm / 1000
        else:
            raise ValueError("Splitter width should be positive")
        
        # Translate splitter gap
        if splitter_gap_nm > 0:
            self.sg = splitter_gap_nm / 1000
        else:
            raise ValueError("Splitter gap should be positive")
            
        # Translate splitter gap
        if port_bias > 0:
            self.pb = port_bias
        else:
            raise ValueError("Port bias should be positive")
        
        # Translate splitter length
        if splitter_length > 0:
            self.sl = splitter_length
        else:
            raise ValueError("Splitter length should be positive") 
            
    def create_waveguide(self, tolerance=1e-3):
        ld_waveguide = get_ld('ld_waveguide')
        cr1 = self.rr + self.cg + self.ww/2
        cr2 = self.rr + self.cg2 + self.ww2/2
        wly1 = self.cy + cr1 * (2*np.cos(self.ca/2) -1)
        wly2 = self.cy - self.br * 2 - cr2 * (2*np.cos(self.ca2/2) -1)
        wry1 = wly1 + self.br * 6
        wry2 = wly2 - self.br * 2
        sly = (wly1 + wly2) / 2
        sry = (wry1 + wry2) / 2
        dy = (wly1 - wly2 - self.sw - self.sg) / 2
        dx = np.sqrt(4 * self.br * dy - dy**2)
        theta = np.arccos((self.br * 2 - dy) / (self.br * 2))
        
        # main waveguide
        waveguides = []
        ly1 = sly + self.sw / 2 + self.sg / 2
        ry1 = sry + self.sw / 2 + self.sg / 2
        waveguide = gp.Path(
            width = self.sw,
            initial_point = (self.lx, ly1))
        waveguide.segment(
            length = self.br * 3 + self.sl,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = -theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.cx - self.lx - self.br*3 - dx - self.sl - self.rl/2,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl / 2 - 2 * cr1 * np.sin(self.ca/2),
            final_width = self.ww,
            **ld_waveguide)
        waveguide.turn(
            radius = cr1, 
            angle = self.ca/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr1, 
            angle = -self.ca,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr1, 
            angle = self.ca/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl / 2 - 2 * cr1 * np.sin(self.ca/2),
            direction = "+x",
            final_width = self.sw,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'll',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl,
            direction = "-x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'r',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br * 2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'r',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl + self.br*2,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'r',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br * 2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'l',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - self.cx - dx - self.br*7 - self.rl/2 - self.sl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = -theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'l',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.pb,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'r',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            **ld_waveguide)
        waveguides.append(waveguide)
        
        # auxillary waveguide
        ly2 = sly - self.sw / 2 - self.sg / 2
        ry2 = sry - self.sw / 2 - self.sg / 2
        waveguide = gp.Path(
            width = self.sw,
            initial_point = (self.lx, ly2 - self.br*2 - self.pb))
        waveguide.segment(
            length = self.br,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "l",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.pb,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "r",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = -theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.cx - self.lx - self.br * 3 - dx - self.sl + self.rl/2,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "ll",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl / 2 - 2 * cr2 * np.sin(self.ca2/2),
            final_width = self.ww2,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = -self.ca2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl / 2 - 2 * cr2 * np.sin(self.ca2/2),
            direction = "-x",
            final_width = self.sw,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "l",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br * 2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "l",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl + self.br*2,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'l',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br * 2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'r',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - self.cx - self.rl / 2 - dx - self.br*7 - self.sl,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = -theta,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl + self.br * 3,
            **ld_waveguide)
        waveguides.append(waveguide)
        
        return waveguides
    
    def create_taper(
            self, 
            taper_length = 100, 
            tip_length = 100, 
            tip_width_nm = 100, 
            tip_width2_nm = 100,
            anchor_base = 4,
            anchor_length = 4):         
        cr1 = self.rr + self.cg + self.ww/2
        cr2 = self.rr + self.cg2 + self.ww2/2
        wly1 = self.cy + cr1 * (2*np.cos(self.ca/2) -1)
        wly2 = self.cy - self.br * 2 - cr2 * (2*np.cos(self.ca2/2) -1)
        wry1 = wly1 + self.br * 6
        wry2 = wly2 - self.br * 2
        sly = (wly1 + wly2) / 2
        sry = (wry1 + wry2) / 2
        ly1 = sly + self.sw / 2 + self.sg / 2
        ry1 = sry + self.sw / 2 + self.sg / 2 + self.br*2 + self.pb
        ly2 = sly - self.sw / 2 - self.sg / 2 - self.br*2 - self.pb
        ry2 = sry - self.sw / 2 - self.sg / 2
            
        couplers = create_taper(
            posl = (self.lx, ly1), 
            posr = (self.rx, ry1), 
            waveguide_width = self.sw, 
            taper_length = taper_length, 
            tip_length = tip_length, 
            tip_width_nm = tip_width_nm, 
            anchor_base = anchor_base, 
            anchor_length = anchor_length)
        
        couplers.extend(
            create_taper(
                posl = (self.lx, ly2), 
                posr = (self.rx, ry2), 
                waveguide_width = self.sw, 
                taper_length = taper_length, 
                tip_length = tip_length, 
                tip_width_nm = tip_width2_nm, 
                anchor_base = anchor_base, 
                anchor_length = anchor_length))
        
        return couplers