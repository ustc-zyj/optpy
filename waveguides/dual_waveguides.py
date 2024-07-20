# -*- coding: utf-8 -*-

import numpy as np
import gdspy as gp
import warnings

from optpy.function import(
    get_ld,
    create_taper
    )

from optpy.waveguides.basic_waveguide import Basic_Waveguide

#%% Dual Waveguide
class BasicDual_Waveguide(Basic_Waveguide):
    
    __slots__ = (
        "nn",
        "rr", 
        "rw", 
        "ww", 
        "ww2",
        "cw",
        "cw2",
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
            ring_center, 
            ring_radius, 
            ring_width_nm, 
            couple_width_nm,
            couple_width2_nm,
            couple_gap_nm,
            couple_gap2_nm,
            couple_angle_deg,
            couple_angle2_deg, 
            endl_xpos, 
            endr_xpos, 
            waveguide_width_nm,
            waveguide_width2_nm,
            cpara = None,
            npara = None,
            tpara = None
            ):
        
        super().__init__(
            ring_center = ring_center, 
            ring_radius = ring_radius, 
            ring_width_nm = ring_width_nm, 
            couple_width_nm = couple_width_nm,
            couple_gap_nm = couple_gap_nm, 
            couple_angle_deg = couple_angle_deg, 
            endl_xpos = endl_xpos, 
            endr_xpos = endr_xpos, 
            waveguide_width_nm = waveguide_width_nm,
            cpara = cpara,
            npara = npara,
            tpara = tpara)
        
        # Translate couple width2
        if couple_width2_nm >= 0:
            self.cw2 = couple_width2_nm / 1000
        else:
            raise ValueError("Couple width should be positive")
        
        # Translate couple gap2
        if couple_gap2_nm >= 0:
            self.cg2 = couple_gap2_nm / 1000
        else:
            raise ValueError("Couple gap should be positive or zero")
            
        # Translate couple angle2
        if couple_angle_deg >= 0:
            self.ca2 = couple_angle2_deg / 180 * np.pi
        else:
            raise ValueError("Couple angle should be positive or zero")
        
        # Translate waveguide width
        if waveguide_width2_nm > 0:
            self.ww2 = waveguide_width2_nm / 1000
        else:
            raise ValueError("Waveguide width should be positive")
        
        
    def create_waveguide(self, tolerance = 1e-3):
        ld_waveguide = get_ld('ld_waveguide')
        
        cr2 = self.rr + self.cg2 + self.cw2/2
        wy2 = self.cy - cr2 * (2*np.cos(self.ca2/2) -1)
        
        waveguides = super().create_waveguide(tolerance)
        
        waveguide = gp.Path(
            width = self.ww2,
            initial_point = (self.lx, wy2))
        waveguide.segment(
            length = self.cx - self.lx - self.rr -\
                2 * cr2 * np.sin(self.ca2/2),
            direction = "+x",
            **ld_waveguide)
        waveguide.segment(
            length = self.rr,
            final_width = self.cw2,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = -self.ca2/2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = self.ca2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = -self.ca2/2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.segment(
            length = self.rr,
            final_width = self.ww2,
            direction = "+x",
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - self.cx - self.rr -\
                2 * cr2 * np.sin(self.ca2/2),
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
            anchor_length = 4,
            istp2hidden = False): 
        wy2 = self.cy - \
            (self.rr + self.cg2 + self.ww2/2) * (2*np.cos(self.ca2/2) -1)
            
        couplers = super().create_taper(
            taper_length = taper_length,
            tip_length = tip_length,
            tip_width_nm = tip_width_nm,
            anchor_base = anchor_base,
            anchor_length = anchor_length)
        
        if istp2hidden:
            return couplers
        
        couplers.extend(
            create_taper(
                posl = (self.lx, wy2), 
                posr = (self.rx, wy2), 
                waveguide_width = self.ww2, 
                taper_length = taper_length, 
                tip_length = tip_length, 
                tip_width_nm = tip_width2_nm, 
                anchor_base = anchor_base, 
                anchor_length = anchor_length))
        
        return couplers

#%% Reverse Waveguide
class Reverse_Waveguide(BasicDual_Waveguide):
    
    __slots__ = (
        "nn",
        "rr", 
        "rw", 
        "ww", 
        "ww2",
        "cw",
        "cw2",
        "cg",
        "cg2",
        "ca",
        "ca2",
        "br",
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
            couple_width_nm,
            couple_width2_nm,
            couple_gap_nm,
            couple_gap2_nm,
            couple_angle_deg,
            couple_angle2_deg, 
            endl_xpos, 
            endr_xpos, 
            waveguide_width_nm,
            waveguide_width2_nm,
            bend_radius,
            cpara = None,
            npara = None,
            tpara = None
            ):
        
        super().__init__(
            ring_center = ring_center, 
            ring_radius = ring_radius, 
            ring_width_nm = ring_width_nm, 
            couple_width_nm = couple_width_nm,
            couple_width2_nm = couple_width2_nm,
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
        if bend_radius > 0:
            self.br = bend_radius
        else:
            raise ValueError("Bend radius should be positive")
            
    def create_waveguide(self, tolerance=1e-3):
        ld_waveguide = get_ld('ld_waveguide')
        cr2 = self.rr + self.cg2 + self.cw2/2
        ly2 = self.cy - self.br * 2 - cr2 * (2*np.cos(self.ca2/2) -1)
        waveguides = []
        waveguide = gp.Path(
            width = self.ww2,
            initial_point = (self.lx, ly2))
        waveguide.segment(
            length = self.cx - self.lx + self.br,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "ll",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.br - self.rr - 2 * cr2 * np.sin(self.ca2/2),
            direction = "-x",
            **ld_waveguide)
        waveguide.segment(
            length = self.rr,
            final_width = self.cw2,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = -self.ca2/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)        
        waveguide.turn(
            radius = cr2, 
            angle = -self.ca2/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.rr,
            final_width = self.ww2,
            **ld_waveguide)
        waveguide.segment(
            length = self.cx - self.lx - self.rr - 2 * cr2 * np.sin(self.ca2/2),
            direction = "-x",
            **ld_waveguide)
        waveguides.append(waveguide)
        
        waveguides.extend(
            super().create_waveguide(tolerance, iswg2hidden=True))
        
        return waveguides
        
    def create_taper(
            self, 
            taper_length = 100, 
            tip_length = 100, 
            tip_width_nm = 100, 
            tip_width2_nm = 100,
            anchor_base = 4,
            anchor_length = 4
            ):
        
        ly2 = self.cy - self.br * 2 -\
            (self.rr + self.cg2 + self.ww2/2) * (2*np.cos(self.ca2/2) -1)
        ry2 = (self.cy - self.rr - self.cg2 - self.ww2/2) * 2 - ly2
            
        couplers = []
        
        if self.LT:
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
        else:
            couplers.extend(
                super().create_taper(
                    taper_length = taper_length,
                    tip_length = tip_length,
                    tip_width_nm = tip_width_nm,
                    anchor_base = anchor_base,
                    anchor_length = anchor_length,
                    istp2hidden = True))
            couplers.extend(
                create_taper(
                    posl = (self.lx, ly2), 
                    posr = None, 
                    waveguide_width = self.ww2, 
                    taper_length = taper_length, 
                    tip_length = tip_length, 
                    tip_width_nm = tip_width2_nm, 
                    anchor_base = anchor_base, 
                    anchor_length = anchor_length))
            couplers.extend(
                create_taper(
                    posl = (self.lx, ly2 + self.br*2), 
                    posr = None, 
                    waveguide_width = self.ww2, 
                    taper_length = taper_length, 
                    tip_length = tip_length, 
                    tip_width_nm = tip_width2_nm, 
                    anchor_base = anchor_base, 
                    anchor_length = anchor_length))
            
        return couplers
            
#%% Multiplexing Waveguide
class Multiplexing_Waveguide(BasicDual_Waveguide):
    
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
        if bend_radius > 0:
            self.br = bend_radius
        else:
            raise ValueError("Bend radius should be positive")
        
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
        # waveguides = super().create_waveguide(tolerance, True)
        ly1 = self.cy +\
            (self.rr + self.cg + self.ww/2) * (2*np.cos(self.ca/2) -1)
        ly2 = self.cy - self.br * 2 -\
            (self.rr + self.cg2 + self.ww2/2) * (2*np.cos(self.ca2/2) -1)
        
        # main waveguide
        waveguides = []
        wy = self.cy + \
            (self.rr + self.cg + self.ww/2) * (2*np.cos(self.ca/2) - 1)
        waveguide = gp.Path(
            width = self.ww,
            initial_point = (self.lx, wy))
        waveguide.segment(
            length = self.br,
            direction = "+x",
            **ld_waveguide)
        waveguide.segment(
            length = self.br * 2,
            final_width = self.sw,
            **ld_waveguide)
        waveguide.segment(
            length = self.cx - self.lx - self.br * 5 - self.rl/2,
            **ld_waveguide)
        waveguide.segment(
            length = self.br * 2,
            final_width = self.ww,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl/2 -\
                2 * (self.rr + self.cg + self.ww/2) * np.sin(self.ca/2),
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
            length = self.rx - self.cx - self.br * 5 - self.sl -\
                2 * (self.rr + self.cg + self.ww/2) * np.sin(self.ca/2),
            direction = "+x",
            **ld_waveguide)
        waveguide.segment(
            length = self.br * 2,
            final_width = self.sw,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl,
            **ld_waveguide)
        waveguide.segment(
            length = self.br * 2,
            final_width = self.ww,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            direction = "+x",
            **ld_waveguide)
        waveguides.append(waveguide)
        
        # auxillary waveguide
        waveguide = gp.Path(
            width = self.ww2,
            initial_point = (self.lx, ly2))
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
            length = ly1 - ly2 - \
                self.br * 2 - self.sw - self.sg,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "r",
            tolerance = tolerance,
            final_width = self.sw,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "r",
            tolerance = tolerance,
            final_width = self.ww2,
            **ld_waveguide)
        waveguide.segment(
            length = ly1 - ly2 - self.br * 2 - self.sw - self.sg,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "l",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.cx - self.lx - self.br * 5 - self.sl + self.rl/2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "ll",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl/2 -\
                (self.rr + self.cg2 + self.ww2/2) * 2*np.sin(self.ca2/2),
            direction = "-x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.rr + self.cg2 + self.ww2/2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = self.rr + self.cg2 + self.ww2/2, 
            angle = -self.ca2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = self.rr + self.cg2 + self.ww2/2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.rl/2 -\
                (self.rr + self.cg2 + self.ww2/2) * 2*np.sin(self.ca2/2),
            direction = "-x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = np.pi / 3,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = -np.pi / 3,
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "r",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = (ly1 - ly2) * 2 - self.br * 3,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "r",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - self.cx - self.sl + self.rl/2\
                - self.br * (5 - 2 * np.sin(np.pi/3)),
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "r",
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = ly1 - ly2 - self.br * 2 - self.sw - self.sg,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "l",
            tolerance = tolerance,
            final_width = self.sw,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = "l",
            tolerance = tolerance,
            final_width = self.ww2,
            **ld_waveguide)
        waveguide.segment(
            length = ly1 - ly2 -\
                self.br * 2 - self.sw - self.sg, 
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
        ly1 = self.cy +\
            (self.rr + self.cg + self.ww/2) * (2*np.cos(self.ca/2) -1)
        ly2 = self.cy - self.br * 2 -\
            (self.rr + self.cg2 + self.ww2/2) * (2*np.cos(self.ca2/2) -1)
        ry2 = ly1 * 2 - ly2
            
        couplers = super().create_taper(
            taper_length = taper_length,
            tip_length = tip_length,
            tip_width_nm = tip_width_nm,
            anchor_base = anchor_base,
            anchor_length = anchor_length,
            istp2hidden = True)
        
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


#%% Triport Waveguide
class Triport_Waveguide(BasicDual_Waveguide):
    
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
        "sx",
        "sg",
        "sr",
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
            bend_radius,
            splitter_xpos,
            splitter_gap_nm,
            splitter_radius,
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
        if bend_radius > 0:
            self.br = bend_radius
        else:
            raise ValueError("Bend radius should be positive")
        
        # Translate splitter position
        if splitter_xpos > self.lx + self.br *3 \
            and splitter_xpos < self.cx - self.rr:
            self.sx = splitter_xpos
        else:
            raise ValueError(
                "Splitter should between ring and waveguide left bend")
        
        # Translate splitter gap
        if splitter_gap_nm > 0:
            self.sg = splitter_gap_nm / 1000
        else:
            raise ValueError("Splitter gap should be positive")
            
        # Translate splitter radius
        if splitter_radius > 0:
            self.sr = splitter_radius
        else:
            raise ValueError("Splitter radius should be positive")
        
        # Translate splitter length
        if splitter_length > 0:
            self.sl = splitter_length
        else:
            raise ValueError("Splitter length should be positive")   
        
    def create_waveguide(self, tolerance = 1e-3):
        ld_waveguide = get_ld('ld_waveguide')
        waveguides = []
        ly1 = self.cy +  2 * self.br +\
            (self.rr + self.cg + self.ww/2) * (2*np.cos(self.ca/2) -1)
        ry2 = self.cy - \
            (self.rr + self.cg2 + self.ww2/2) * (2*np.cos(self.ca2/2) -1)
        ly2 = ry2 - self.sr * 2 - self.ww2 - self.sg
        
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
            angle = "r",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "l",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.segment(
            length = self.cx - self.lx - self.br * 3 -\
                2 * (self.rr + self.cg + self.ww/2) * np.sin(self.ca/2),
            direction = "+x",
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
            length = self.rx - self.cx - self.br * 3 -\
                2 * (self.rr + self.cg + self.ww/2) * np.sin(self.ca/2),
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "l",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "r",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            direction = "+x",
            **ld_waveguide)
        waveguides.append(waveguide)
        
        # Auxilary waveguide left
        waveguide = gp.Path(
            width = self.ww2,
            initial_point = (self.lx, ly2))
        waveguide.segment(
            length = self.sx - self.lx + self.sl/2,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.sr, 
            angle = "rr",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.segment(
            length = self.sx - self.lx + self.sl/2 - self.br *3,
            direction = "-x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "l",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "r",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            direction = "-x",
            **ld_waveguide)
        waveguides.append(waveguide)
        
        # Auxillary waveguide right
        waveguide = gp.Path(
            width = self.ww2,
            initial_point = (self.rx, ry2))
        waveguide.segment(
            length = self.rx - self.cx -\
                (self.rr + self.cg2 + self.ww2/2) * 2*np.sin(self.ca2/2),
            direction = "-x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.rr + self.cg2 + self.ww2/2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.turn(
            radius = self.rr + self.cg2 + self.ww2/2, 
            angle = -self.ca2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.turn(
            radius = self.rr + self.cg2 + self.ww2/2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.segment(
            length = self.cx - self.sx + self.sl/2-\
                (self.rr + self.cg2 + self.ww2/2) * 2*np.sin(self.ca2/2),
            direction = "-x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.sr, 
            angle = "ll",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - self.sx + self.sl/2 - self.br * 3,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "r",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br, 
            angle = "l",
            tolerance = tolerance * 0.1,
            max_points = 0,
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
        ly1 = self.cy +  2 * self.br +\
            (self.rr + self.cg + self.ww/2) * (2*np.cos(self.ca/2) -1)
        ry1 = ly1
        ry2 = self.cy - \
            (self.rr + self.cg2 + self.ww2/2) * (2*np.cos(self.ca2/2) -1)
        ry3 = ry2 - self.sr * 2 - self.br * 2
        ly2 = ry2 - self.sr * 2 - self.ww2 - self.sg
        ly3 = ly2 - self.sr * 2 - self.br * 2
        couplers = []
        couplers.extend(
            create_taper(
                posl = (self.lx, ly1), 
                posr = (self.rx, ry1), 
                waveguide_width = self.ww, 
                taper_length = taper_length, 
                tip_length = tip_length, 
                tip_width_nm = tip_width_nm, 
                anchor_base = anchor_base, 
                anchor_length = anchor_length))
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
        couplers.extend(
            create_taper(
                posl = (self.lx, ly3), 
                posr = (self.rx, ry3), 
                waveguide_width = self.ww2, 
                taper_length = taper_length, 
                tip_length = tip_length, 
                tip_width_nm = tip_width2_nm, 
                anchor_base = anchor_base, 
                anchor_length = anchor_length))
        return couplers
    
#%% Triport MultiRing Waveguide
class TriportMR_Waveguide(Triport_Waveguide):
    
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
        "sx",
        "sg",
        "sr",
        "sl",
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
            splitter_xpos, 
            splitter_gap_nm, 
            splitter_radius, 
            splitter_length,
            cpara = None,
            npara = None
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
            bend_radius = bend_radius, 
            splitter_xpos = splitter_xpos, 
            splitter_gap_nm = splitter_gap_nm, 
            splitter_radius = splitter_radius, 
            splitter_length = splitter_length,
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
        ly1 = self.cy +  2 * self.br +\
            (self.rr + self.cg + self.ww/2) * (2*np.cos(self.ca/2) -1)
        ry2 = self.cy - \
            (self.rr + self.cg2 + self.ww2/2) * (2*np.cos(self.ca2/2) -1)
        ly2 = ry2 - self.sr * 2 - self.ww2 - self.sg
        
        # main waveguide
        waveguide0 = gp.Path(
            width = self.ww,
            initial_point = (self.lx, ly1))
        waveguide0.segment(
            length = self.br,
            direction = "+x",
            **ld_waveguide)
        waveguide0.turn(
            radius = self.br, 
            angle = "r",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide0.turn(
            radius = self.br, 
            angle = "l",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide0.segment(
            length = self.cx[0] - self.lx - self.br * 3 -\
                2 * (self.rr + self.cg + self.ww/2) * np.sin(self.ca/2),
            direction = "+x",
            **ld_waveguide)
        waveguide0.turn(
            radius = self.rr + self.cg + self.ww/2, 
            angle = self.ca/2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide0.turn(
            radius = self.rr + self.cg + self.ww/2, 
            angle = -self.ca,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide0.turn(
            radius = self.rr + self.cg + self.ww/2, 
            angle = self.ca/2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        for kk in range(1, self.nn):
            waveguide0.segment(
                length = self.cx[kk] - self.cx[kk-1] -\
                    4 * (self.rr + self.cg + self.ww/2) * np.sin(self.ca/2),
                direction = "+x",
                **ld_waveguide)
            waveguide0.turn(
                radius = self.rr + self.cg + self.ww/2, 
                angle = self.ca/2,
                tolerance = tolerance * 0.1,
                max_points = 0,
                **ld_waveguide)
            waveguide0.turn(
                radius = self.rr + self.cg + self.ww/2, 
                angle = -self.ca,
                tolerance = tolerance * 0.1,
                max_points = 0,
                **ld_waveguide)
            waveguide0.turn(
                radius = self.rr + self.cg + self.ww/2, 
                angle = self.ca/2,
                tolerance = tolerance * 0.1,
                max_points = 0,
                **ld_waveguide)
        waveguide0.segment(
            length = self.rx - self.cx[self.nn-1] - self.br * 3 -\
                2 * (self.rr + self.cg + self.ww/2) * np.sin(self.ca/2),
            direction = "+x",
            **ld_waveguide)
        waveguide0.turn(
            radius = self.br, 
            angle = "l",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide0.turn(
            radius = self.br, 
            angle = "r",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide0.segment(
            length = self.br,
            direction = "+x",
            **ld_waveguide)
        waveguides.append(waveguide0)
        
        # Auxilary waveguide left
        waveguide1 = gp.Path(
            width = self.ww2,
            initial_point = (self.lx, ly2))
        waveguide1.segment(
            length = self.sx - self.lx + self.sl/2,
            direction = "+x",
            **ld_waveguide)
        waveguide1.turn(
            radius = self.sr, 
            angle = "rr",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide1.segment(
            length = self.sx - self.lx + self.sl/2 - self.br *3,
            direction = "-x",
            **ld_waveguide)
        waveguide1.turn(
            radius = self.br, 
            angle = "l",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide1.turn(
            radius = self.br, 
            angle = "r",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide1.segment(
            length = self.br,
            direction = "-x",
            **ld_waveguide)
        waveguides.append(waveguide1)
        
        # Auxillary waveguide right
        waveguide2 = gp.Path(
            width = self.ww2,
            initial_point = (self.rx, ry2))
        waveguide2.segment(
            length = self.rx - self.cx[self.nn-1] -\
                2 * (self.rr + self.cg2 + self.ww2/2) * np.sin(self.ca2/2),
            direction = "-x",
            **ld_waveguide)
        waveguide2.turn(
            radius = self.rr + self.cg2 + self.ww2/2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide2.turn(
            radius = self.rr + self.cg2 + self.ww2/2, 
            angle = -self.ca2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide2.turn(
            radius = self.rr + self.cg2 + self.ww2/2, 
            angle = self.ca2/2,
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        for kk in range(self.nn-1, 0, -1):
            waveguide2.segment(
                length = self.cx[kk] - self.cx[kk-1] -\
                    4 * (self.rr + self.cg2 + self.ww2/2) * np.sin(self.ca2/2),
                direction = "-x",
                **ld_waveguide)
            waveguide2.turn(
                radius = self.rr + self.cg2 + self.ww2/2, 
                angle = self.ca2/2,
                tolerance = tolerance * 0.1,
                max_points = 0,
                **ld_waveguide)
            waveguide2.turn(
                radius = self.rr + self.cg2 + self.ww2/2, 
                angle = -self.ca2,
                tolerance = tolerance * 0.1,
                max_points = 0,
                **ld_waveguide)
            waveguide2.turn(
                radius = self.rr + self.cg2 + self.ww2/2, 
                angle = self.ca2/2,
                tolerance = tolerance * 0.1,
                max_points = 0,
                **ld_waveguide)
        waveguide2.segment(
            length = self.cx[0] - self.sx + self.sl/2-\
                2 *(self.rr + self.cg2 + self.ww2/2) * np.sin(self.ca2/2),
            direction = "-x",
            **ld_waveguide)
        waveguide2.turn(
            radius = self.sr, 
            angle = "ll",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide2.segment(
            length = self.rx - self.sx + self.sl/2 - self.br * 3,
            direction = "+x",
            **ld_waveguide)
        waveguide2.turn(
            radius = self.br, 
            angle = "r",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide2.turn(
            radius = self.br, 
            angle = "l",
            tolerance = tolerance * 0.1,
            max_points = 0,
            **ld_waveguide)
        waveguide2.segment(
            length = self.br,
            direction = "+x",
            **ld_waveguide)
        waveguides.append(waveguide2)
        return waveguides