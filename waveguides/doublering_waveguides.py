# -*- coding: utf-8 -*-

import numpy as np
import gdspy as gp

from optpy.function import(
    get_ld,
    create_taper,
    create_thermode,
    )

from optpy.waveguides.basic_waveguide import Basic_Waveguide

#%% Balanced MZ-Roller Waveguide
class UPB_Waveguide(Basic_Waveguide):
    
    __slots__ = (
        "nn",
        "rr", 
        "rw", 
        "rg",
        "ra",
        "ww", 
        "ww2",
        "br",
        "cw",
        "cw2",
        "cg",
        "cg2",
        "ca",
        "ca2",
        "al",
        "pb",
        "sw",
        "sg",
        "sl",
        "cx", 
        "cy", 
        "lx",
        "rx",
        "db",
        "dw1",
        "dw2",
        "cv",
        "nv",
        "tv"
        )
    
    def __init__(
            self, 
            ring_center, 
            ring_radius, 
            ring_width_nm, 
            ring_gap_nm,
            ring_angle_deg,
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
            adj_length,
            port_bias,
            splitter_width_nm,
            splitter_gap_nm,
            splitter_length,
            draw_waveguide1 = True,
            draw_waveguide2 = True,
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
        
        # Translate ring distance
        if ring_gap_nm >= 0:
            self.rg = ring_gap_nm / 1000
        else:
            raise ValueError("Ring distance should be positive")
            
        # Translate ring angle
        self.ra = ring_angle_deg / 180 * np.pi
            
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
        if couple_angle2_deg >= 0:
            self.ca2 = couple_angle2_deg / 180 * np.pi
        else:
            raise ValueError("Couple angle should be positive or zero")
            
        # Translate waveguide width
        if waveguide_width2_nm > 0:
            self.ww2 = waveguide_width2_nm / 1000
        else:
            raise ValueError("Waveguide width should be positive")
        
        # Translate bend radius
        if bend_radius >= ring_radius:
            self.br = bend_radius
        else:
            raise ValueError("Bend radius should be larger than ring radius")
        
        # Translate adjusting length
        if adj_length > 0:
            self.al = adj_length
        else:
            raise ValueError("Adjusting length should be positive")
            
        # Translate port bias
        if port_bias >= 0:
            self.pb = port_bias
        else:
            raise ValueError("Port bias should be positive")    
            
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
        
        # Translate flags
        self.dw1 = draw_waveguide1
        self.dw2 = draw_waveguide2
    
    def create_ring(self, tolerance=1e-4):
        ld_ring = get_ld('ld_ring')
        rings = super().create_ring(tolerance)
        rdx = (self.rr*2 + self.rg) * np.cos(self.ra)
        rdy = (self.rr*2 + self.rg) * np.sin(self.ra)
        if self.cg > 0 and self.rg > 0:
            ring = gp.Round(
               center = (self.cx + rdx, self.cy + rdy), 
               radius = self.rr,
               inner_radius = self.rr - self.rw,
               tolerance = tolerance,
               **ld_ring)
            rings.append(ring)
        return rings
    
    def create_waveguide(self, tolerance=1e-3):
        ld_waveguide = get_ld('ld_waveguide')
        ld_couple = get_ld('ld_couple')
        cr1 = self.rr + self.cg + self.cw/2
        cr2 = self.rr + self.cg2 + self.cw2/2
        wy1 = self.cy + cr1
        wy2 = self.cy - cr2
        sly = self.cy
        sry = self.cy
        
        # main waveguide
        waveguides = []
        ly1 = sly + self.sw / 2 + self.sg / 2
        ry1 = sry + self.sw / 2 + self.sg / 2
        dy1 = self.br*2 - cr1 * (2 - 2 * (1-np.cos(self.ca/2))) + ly1 - wy1
        dxl = self.cx - self.lx - self.br * 6 - self.sl - cr1*2
        dxr = self.rx - self.cx - self.br * 6 - self.sl - cr1*2
        waveguide = gp.Path(
            width = self.sw,
            initial_point = (self.lx, ly1))
        waveguide.segment(
            length = self.br * 3 + self.sl,
            direction = "+x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'l',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'r',
            tolerance = tolerance,
            final_width = self.ww
            **ld_waveguide)
        waveguide.segment(
            length = dxl,
            **ld_waveguide)
        waveguide.segment(
            length = self.br - 2 * cr1 * np.sin(self.ca/2),
            final_width = self.cw,
            **ld_waveguide)
        waveguide.turn(
            radius = cr1, 
            angle = 'r',
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = dy1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr1, 
            angle = 'l',
            tolerance = tolerance * 0.1,
            **ld_couple)   
        if self.ca != 0:
            waveguide.turn(
                radius = cr1, 
                angle = self.ca/2,
                tolerance = tolerance * 0.1,
                **ld_couple)
            waveguide.turn(
                radius = cr1, 
                angle = -self.ca,
                tolerance = tolerance * 0.1,
                **ld_couple)
            waveguide.turn(
                radius = cr1, 
                angle = self.ca/2,
                tolerance = tolerance * 0.1,
                **ld_couple)
        waveguide.turn(
            radius = cr1,
            angle = 'l',
            tolerance = tolerance * 0.1,
            **ld_couple)
        waveguide.segment(
            length = dy1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr1,
            angle = 'r',
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.br - 2 * cr1 * np.sin(self.ca/2),
            direction = "+x",
            final_width = self.sw,
            **ld_waveguide)
        waveguide.segment(
            length = dxr,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'r',
            tolerance = tolerance,
            final_width = self.sw,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'l',
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
        if self.dw1:
            waveguides.append(waveguide)
        
        # auxillary waveguide
        ly2 = sly - self.sw / 2 - self.sg / 2
        ry2 = sry - self.sw / 2 - self.sg / 2
        dy2 = self.br*2 - cr2 * (2 - 2 * (1-np.cos(self.ca2/2))) - ly2 + wy2
        cl = self.cy - ly2 + self.br + cr2 * (1 + 2 * np.sin(self.ca2/2))
        ar = cr2 * (2 - 2*(1-np.cos(self.ca2/2)) - 2*np.sin(self.ca2/2)) / 2
        waveguide = gp.Path(
            width = self.sw,
            initial_point = (self.lx, ly2 - self.pb - self.br*2))
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
            angle = 'r',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'l',
            tolerance = tolerance,
            final_width = self.ww,
            **ld_waveguide)
        waveguide.segment(
            length = self.cx - self.lx - self.br*6 - self.sl +\
                cr2 * 2 * np.sin(self.ca2/2),
            direction = "+x",
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            final_width = self.cw
            **ld_waveguide)
        waveguide.turn(
            radius = cr2,
            angle = 'l',
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = dy2,
            **ld_waveguide)
        waveguide.turn(
            radius = cr2,
            angle = 'l',
            tolerance = tolerance * 0.1,
            **ld_couple)
        if self.ca2 != 0:
            waveguide.turn(
                radius = cr2, 
                angle = self.ca2/2,
                tolerance = tolerance * 0.1,
                **ld_couple)
            waveguide.turn(
                radius = cr2, 
                angle = -self.ca2,
                tolerance = tolerance * 0.1,
                **ld_couple)
            waveguide.turn(
                radius = cr2, 
                angle = self.ca2/2,
                tolerance = tolerance * 0.1,
                **ld_couple)
        waveguide.turn(
            radius = cr2,
            angle = 'l',
            tolerance = tolerance * 0.1,
            **ld_couple)
        waveguide.turn(
            radius = cr2,
            angle = 'r',
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            final_width = self.ww
            **ld_waveguide)
        waveguide.segment(
            length = self.al/2 - self.br,
            direction = "-x",
            **ld_waveguide)
        waveguide.turn(
            radius = ar,
            angle = 'rr',
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.br,
            final_width = self.cw
            **ld_waveguide)
        waveguide.segment(
            length = self.al/2 - self.br +\
                cr2 * (2 * (1-np.cos(self.ca2/2)) + 2 * np.sin(self.ca2/2)),
            **ld_waveguide)
        waveguide.turn(
            radius = cr2,
            angle = 'l',
            tolerance = tolerance * 0.1,
            **ld_couple)
        if self.ca2 != 0:
            waveguide.turn(
                radius = cr2, 
                angle = self.ca2/2,
                tolerance = tolerance * 0.1,
                **ld_couple)
            waveguide.turn(
                radius = cr2, 
                angle = -self.ca2,
                tolerance = tolerance * 0.1,
                **ld_couple)
            waveguide.turn(
                radius = cr2, 
                angle = self.ca2/2,
                tolerance = tolerance * 0.1,
                **ld_couple)
        waveguide.turn(
            radius = cr2,
            angle = 'l',
            tolerance = tolerance * 0.1,
            **ld_couple)
        waveguide.segment(
            length = self.br,
            final_width = self.ww
            **ld_waveguide)
        waveguide.segment(
            length = self.al/2 + ar - self.br,
            direction = "-x",
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'l',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = cl * 2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'l',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.al/2 + ar + cr2 * (1 - 2*(1 - np.cos(self.ca2/2))),
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'l',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = cl - self.br,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'r',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - self.cx - self.br * 7 - self.sl + cr2,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'l',
            tolerance = tolerance,
            final_width = self.sw,
            **ld_waveguide)
        waveguide.turn(
            radius = self.br,
            angle = 'r',
            tolerance = tolerance,
            **ld_waveguide)
        waveguide.segment(
            length = self.sl + self.br * 3,
            **ld_waveguide)
        if self.dw2:
            waveguides.append(waveguide)
        
        return waveguides
    
    def create_taper(
            self, 
            taper_length, 
            tip_length, 
            tip_width_nm, 
            tip_width2_nm,
            anchor_base = 4,
            anchor_length = 4):         
        sly = self.cy
        sry = self.cy
        ly1 = sly + self.sw / 2 + self.sg / 2
        ry1 = sry + self.sw / 2 + self.sg / 2 + self.br*2 + self.pb
        ly2 = sly - self.sw / 2 - self.sg / 2 - self.br*2 - self.pb
        ry2 = sry - self.sw / 2 - self.sg / 2
        couplers = []
            
        coupler = create_taper(
            posl = (self.lx, ly1), 
            posr = (self.rx, ry1), 
            waveguide_width = self.sw, 
            taper_length = taper_length, 
            tip_length = tip_length, 
            tip_width_nm = tip_width_nm, 
            anchor_base = anchor_base, 
            anchor_length = anchor_length)
        if self.dw1:
            couplers.extend(coupler)
        
        coupler = create_taper(
            posl = (self.lx, ly2), 
            posr = (self.rx, ry2), 
            waveguide_width = self.sw, 
            taper_length = taper_length, 
            tip_length = tip_length, 
            tip_width_nm = tip_width2_nm, 
            anchor_base = anchor_base, 
            anchor_length = anchor_length)
        if self.dw2:
            couplers.extend(coupler)
        
        return couplers
    
    def create_thermode(self, tolerance=1e-3):
        thermodes = []
        holepos = []
        ld_thermode = get_ld('ld_pad')
        if self.tv is None:
            return thermodes, holepos
        cr1 = self.rr + self.cg + self.cw/2
        cr2 = self.rr + self.cg2 + self.cw2/2
        ar = cr2 * (2 - 2*(1-np.cos(self.ca2/2)) - 2*np.sin(self.ca2/2)) / 2
        tpara = self.tv.copy()
        tw = tpara['thermode_width']
        hr = tpara['hole_radius']
        
        dt0 = tpara['draw_thermode0']
        dt1 = tpara['draw_thermode1']
        dt2 = tpara['draw_thermode2']
        inv = tpara['inverse']
        
        # Thermode 0 (Adjust coupling)
        txr = self.cx - cr2 * (2 - 2 * (1 - np.cos(self.ca2/2)))
        ty1 = self.cy - cr2 * (1 + 2 * np.sin(self.ca2/2))
        ty2 = self.cy - cr2 * (3 - 2 * (1 - np.cos(self.ca2/2)))
        if inv:
            hx1 = txr - self.br*2
            hx2 = txr - self.br
        else:
            hx1 = txr - self.br
            hx2 = txr - self.br*2
        T1 = gp.Rectangle(
            (hx1, ty2 - tw/2), 
            (txr - self.al/2, ty2 + tw/2),
            **ld_thermode)
        T2 = gp.Rectangle(
            (hx2, ty1 + tw/2), 
            (txr - self.al/2, ty1 - tw/2),
            **ld_thermode)
        T3 = gp.Round(
            center = (txr - self.al/2, (ty1 + ty2) / 2), 
            radius = ar + tw/2,
            inner_radius = ar - tw/2,
            initial_angle = np.pi/2,
            final_angle = np.pi/2 * 3,
            tolerance = tolerance,
            **ld_thermode)
        H1 = gp.Round(
            center = (hx1, ty2 - self.cw2/2 + tw/2), 
            radius = hr,
            tolerance = tolerance,
            **ld_thermode)
        H2 = gp.Round(
            center = (hx2, ty1 + self.cw2/2 - tw/2), 
            radius = hr,
            tolerance = tolerance,
            **ld_thermode)
        if dt0:
            thermodes.extend([T1, T2, T3, H1, H2])
            holepos.extend([
                (hx1, ty2 - self.cw2/2 + tw/2), 
                (hx2, ty1 + self.cw2/2 - tw/2)])
           
        # Thermode 1 on ring 1
        T, hp = create_thermode(
            ring_center = (self.cx, self.cy),
            thermode_type = 'single',
            thermode_angle_deg = 180,
            thermode_width = tw,
            thermode_gap = -(self.rw/2 + tw/2),
            ring_radius = self.rr,
            hole_radius = hr,
            hole_distance = self.rr * 1.5,
            hole_num = 0,
            rotate_angle_deg = 90
            )
        if dt1:
            thermodes.extend(T)
            holepos.extend(hp)
        
        # Thermode 2 on ring 2
        rdx = (self.rr*2 + self.rg) * np.cos(self.ra)
        rdy = (self.rr*2 + self.rg) * np.sin(self.ra)
        T, hp = create_thermode(
            ring_center = (self.cx + rdx, self.cy + rdy),
            thermode_type = 'single',
            thermode_angle_deg = 180,
            thermode_width = tw,
            thermode_gap = -(self.rw/2 + tw/2),
            ring_radius = self.rr,
            hole_radius = hr,
            hole_distance = self.rr * 1.5,
            hole_num = 0,
            rotate_angle_deg = -90
            )
        if dt2:
            thermodes.extend(T)
            holepos.extend(hp)
        return thermodes, holepos