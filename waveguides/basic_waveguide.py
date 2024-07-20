# -*- coding: utf-8 -*-

import numpy as np
import gdspy as gp
import warnings

from optpy.text import create_name

from optpy.function import(
    get_ld,
    create_taper,
    create_thermode,
    create_electrode,
    )

class Basic_Waveguide(object):
    
    __slots__ = (
        "nn",
        "rr",
        "rw",
        "cw",
        "ww",
        "cg",
        "ca",
        "cx",
        "cy",
        "lx",
        "rx",
        "cv",
        "nv",
        "tv",
        "ev",
        )
    
    def __init__(
            self,
            ring_center,
            ring_radius,
            ring_width_nm,
            couple_width_nm,
            couple_gap_nm,
            couple_angle_deg,
            endl_xpos,
            endr_xpos,
            waveguide_width_nm,
            cpara = None,
            npara = None,
            tpara = None,
            epara = None,
            ):
        
        self.nn = 1
        
        # Translate positions
        self.cx, self.cy = ring_center
        
        if endl_xpos < self.cx:
            self.lx = endl_xpos
        else:
            raise ValueError("Center of ring should between two ends of waveguide")
        
        if endr_xpos > self.cx:
            self.rx = endr_xpos
        else:
            raise ValueError("Center of ring should between two ends of waveguide")
        
        # Translate ring radius
        if ring_radius > 0:
            self.rr = ring_radius
        else:
            raise ValueError("Ring radius should be positive")
            
        # Translate ring width
        if ring_width_nm < ring_radius * 1000 and ring_width_nm > 0:
            self.rw = ring_width_nm / 1000
        elif ring_width_nm == ring_radius * 1000:
            warnings.warn("Creating a disk rather than a ring")
            self.rw = ring_width_nm / 1000
        elif ring_width_nm <= 0:
            raise ValueError("Ring width should be positive")
        else:
            raise ValueError("Ring width should be smaller than radius")
            
        # Translate couple width
        if couple_width_nm > 0:
            self.cw = couple_width_nm / 1000
        else:
            raise ValueError("Couple width should be positive")
        
        # Translate couple gap
        if couple_gap_nm >= 0:
            self.cg = couple_gap_nm / 1000
        else:
            raise ValueError("Couple gap should be positive or zero")
            
        # Translate couple angle
        if couple_angle_deg >= 0:
            self.ca = couple_angle_deg / 180 * np.pi
        else:
            raise ValueError("Couple angle should be positive or zero")
        
        # Translate waveguide width
        if waveguide_width_nm > 0:
            self.ww = waveguide_width_nm / 1000
        else:
            raise ValueError("Waveguide width should be positive")
        
        # Translate coupler variables
        self.cv = cpara
        
        # Translate name variables
        self.nv = npara
        
        # Translate thermode variables
        self.tv = tpara
        
        # Translate electrode variables
        self.ev = epara
    
    def create_ring(self, tolerance=1e-4):
        ld_ring = get_ld('ld_ring')
        rings = []
        if self.cg > 0:
            ring = gp.Round(
               center = (self.cx, self.cy), 
               radius = self.rr,
               inner_radius = self.rr - self.rw,
               tolerance = tolerance,
               **ld_ring)
            rings.append(ring)
        return rings
    
    def create_waveguide(self, tolerance=1e-3):
        ld_waveguide = get_ld('ld_waveguide')
        waveguides = []
        cr = self.rr + self.cg + self.cw/2
        wy = self.cy + cr * (2*np.cos(self.ca/2) - 1)
        waveguide = gp.Path(
            width = self.ww,
            initial_point = (self.lx, wy))
        waveguide.segment(
            length = self.cx - self.lx - 2 * cr * np.sin(self.ca/2) - self.rr,
            direction = "+x",
            **ld_waveguide)
        waveguide.segment(
            length = self.rr,
            final_width = self.cw,
            **ld_waveguide)
        waveguide.turn(
            radius = cr, 
            angle = self.ca/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr, 
            angle = -self.ca,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.turn(
            radius = cr, 
            angle = self.ca/2,
            tolerance = tolerance * 0.1,
            **ld_waveguide)
        waveguide.segment(
            length = self.rr,
            final_width = self.ww,
            **ld_waveguide)
        waveguide.segment(
            length = self.rx - self.cx - 2 * cr * np.sin(self.ca/2) - self.rr,
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
        if hasattr(self, 'ly') and hasattr(self, 'ry'):
            couplers = create_taper(
                posl = (self.lx, self.ly), 
                posr = (self.rx, self.ry), 
                waveguide_width = self.ww, 
                taper_length = taper_length, 
                tip_length = tip_length, 
                tip_width_nm = tip_width_nm, 
                anchor_base = anchor_base, 
                anchor_length = anchor_length)
        elif hasattr(self, 'wy'):
            couplers = create_taper(
                posl = (self.lx, self.wy), 
                posr = (self.rx, self.wy), 
                waveguide_width = self.ww, 
                taper_length = taper_length, 
                tip_length = tip_length, 
                tip_width_nm = tip_width_nm, 
                anchor_base = anchor_base, 
                anchor_length = anchor_length)
        else:
            wy = self.cy + \
                (self.rr + self.cg + self.ww/2) * (2*np.cos(self.ca/2) - 1)
            couplers = create_taper(
                posl = (self.lx, wy), 
                posr = (self.rx, wy), 
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
    
    def create_thermode(self):
        thermodes = []
        holepos = []
        if self.tv is not None and self.cg != 0:
            tpara = self.tv.copy()
            thermodes, holepos = create_thermode(
                ring_center = (self.cx, self.cy),
                ring_radius = self.rr,
                **tpara)
        return thermodes, holepos
    
    def create_electrode(self):
        electrodes = []
        holepos = []
        if self.ev is not None and self.cg != 0:
            epara = self.ev.copy()
            electrodes, holepos = create_electrode(
                ring_center = (self.cx, self.cy),
                ring_radius = self.rr,
                ring_width = self.rw,
                **epara)
        return electrodes, holepos
    
    def create(self, bias=(0,0), flip=False, tolerance=1e-3):
        structures = []
        structures.extend(self.create_ring(tolerance * 0.1))
        structures.extend(self.create_waveguide(tolerance))
        structures.extend(self.create_coupler(tolerance))
        structures.extend(self.create_name())
        thermodes, therm_holepos = self.create_thermode()
        structures.extend(thermodes)
        electrodes, elec_holepos = self.create_electrode()
        structures.extend(electrodes)
        holepos = therm_holepos + elec_holepos
        if flip:
            for structure in structures:
                structure.mirror((self.cx, self.cy), (self.cx + 1, self.cy))
            for kk, pos in enumerate(holepos):
                holepos[kk] = (2 * self.cx - pos[0], pos[1])
        for structure in structures:
            structure.translate(bias[0], bias[1])
        for kk, pos in enumerate(holepos):
            holepos[kk] = (pos[0] + bias[0], pos[1] + bias[1])
        return structures, holepos
    
    def __call__(self, cell, bias=(0,0), flip=False, tolerance=1e-3):
        structures, holepos = self.create(bias, flip, tolerance)
        cell.add(structures)
        return holepos