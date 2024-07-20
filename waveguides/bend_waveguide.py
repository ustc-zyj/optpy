# -*- coding: utf-8 -*-

import numpy as np
import gdspy as gp

from optpy.function import(
    get_ld,
    bezier_points
    )

from optpy.waveguides.basic_waveguide import Basic_Waveguide

class Bend_Waveguide(Basic_Waveguide):
    
    __slots__ = (
        "nn",
        "rr", 
        "rw", 
        "cw",
        "ww", 
        "br", 
        "bt",
        "cg",
        "ca",
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
            couple_width_nm,
            couple_gap_nm,
            couple_angle_deg,
            endl_pos,
            endr_pos,
            waveguide_width_nm,
            bend_radius,
            bend_type,
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
            endl_xpos = endl_pos[0], 
            endr_xpos = endr_pos[0], 
            waveguide_width_nm = waveguide_width_nm,
            cpara = cpara,
            npara = npara,
            tpara = tpara)
        
        # Re-translate positions
        self.ly = endl_pos[1]
        self.ry = endr_pos[1]
        
        # Translate bend radius
        if bend_radius > 0:
            self.br = bend_radius
        else:
            raise ValueError("Bend radius should be positive")
        
        # Translate bend type
        if bend_type == 'round':
            self.bt = 0
        elif bend_type == 'bezier':
            self.bt = 1
        else:
            raise ValueError("Invalid bend type: " + bend_type)
        
        # Check waveguide geometry
        geo_tol = (
            self.cx >= self.lx + self.br + self.rr and \
            self.cx <= self.rx - self.br - self.rr and \
            self.cy >= self.ly + self.br + self.rr and \
            self.cy <= self.ry - self.br - self.rr
            )
        if not geo_tol:
            raise ValueError("Ring position out of range!")
        
    def create_waveguide(self, tolerance=1e-3):
        ld_waveguide = get_ld('ld_waveguide')
        cr = self.rr + self.cg + self.cw/2
        waveguides = []
        waveguide = gp.Path(
            width = self.ww,
            initial_point = (self.lx, self.ly))
        waveguide.segment(
            length = self.cx - self.lx - self.br -\
                cr * (2*np.cos(self.ca/2) - 1),
            direction = '+x',
            **ld_waveguide)
        if self.bt == 0:
            waveguide.turn(
                radius = self.br, 
                angle = 'l',
                tolerance = tolerance,
                max_points = 0,
                **ld_waveguide)
        elif self.bt == 1:
            waveguide.bezier(
                points = bezier_points(self.br, 'l'),
                tolerance = tolerance,
                max_points = 0,
                **ld_waveguide)
        waveguide.segment(
            length = self.cy - self.ly - self.br -\
                2 * cr * np.sin(self.ca/2) - self.rr,
            direction = "+y",
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
            direction = "+y",
            **ld_waveguide)
        waveguide.segment(
            length = self.ry - self.cy - self.br - \
                2 * cr * np.sin(self.ca/2) - self.rr,
            **ld_waveguide)
        if self.bt == 0:
            waveguide.turn(
                radius = self.br, 
                angle = 'r',
                tolerance = tolerance,
                max_points = 0,
                **ld_waveguide)
        elif self.bt == 1:
            waveguide.bezier(
                points = bezier_points(self.br, 'r'),
                tolerance = tolerance,
                max_points = 0,
                **ld_waveguide)
        waveguide.segment(
            length = self.rx - self.cx - self.br +\
                cr * (2*np.cos(self.ca/2) - 1),
            direction = '+x',
            **ld_waveguide)
        waveguides.append(waveguide)
        return waveguides
    
class Bend_Waveguide_v2(Bend_Waveguide):
    
    # def __init__(
    #         self, 
    #         ring_center,
    #         ring_radius,
    #         ring_width_nm,
    #         couple_width_nm,
    #         couple_gap_nm,
    #         couple_angle_deg,
    #         endl_pos,
    #         endr_pos,
    #         waveguide_width_nm,
    #         bend_radius,
    #         bend_type,
    #         cpara = None,
    #         npara = None,
    #         tpara = None
    #         ):
        
    #     super().__init__(
    #         ring_center = ring_center, 
    #         ring_radius = ring_radius, 
    #         ring_width_nm = ring_width_nm, 
    #         couple_width_nm = couple_width_nm,
    #         couple_gap_nm = couple_gap_nm, 
    #         couple_angle_deg = couple_angle_deg, 
    #         endl_pos = endl_pos, 
    #         endr_pos = endr_pos, 
    #         waveguide_width_nm = waveguide_width_nm,
    #         bend_radius = bend_radius,
    #         bend_type = bend_type,
    #         cpara = cpara,
    #         npara = npara,
    #         tpara = tpara)
        
    def create_sector(self, tolerance=1e-4):
        if self.cg != 0 and self.ca > 0:
            sector = gp.Round(
               center = (self.cx, self.cy), 
               radius = self.rr + self.cg + self.cw * 2,
               initial_angle = np.pi - self.ca/2,
               final_angle = np.pi + self.ca/2,
               tolerance = tolerance)
        else:
            sector = None
        return sector
        
    def create(self, bias=(0,0), tolerance=1e-3):
        ld_couple = get_ld('ld_couple')
        ld_ring = get_ld('ld_ring')
        ld_waveguide = get_ld('ld_waveguide')
        structures = []
        ring = self.create_ring(tolerance * 0.1)
        waveguide = self.create_waveguide(tolerance)
        sector = self.create_sector(tolerance * 0.1)
        
        if sector is not None:
            ring_lap = gp.boolean(ring, sector, 'and', **ld_couple)
            if ring_lap is not None:
                structures.append(ring_lap)
                ring_res = gp.boolean(ring, sector, 'not', **ld_ring)
                structures.append(ring_res)
            else:
                structures.extend(ring)
                
            waveguide_lap = gp.boolean(waveguide, sector, 'and', **ld_couple)
            if waveguide_lap is not None:
                structures.append(waveguide_lap)
                waveguide_res = gp.boolean(waveguide, sector, 'not', **ld_waveguide)
                structures.append(waveguide_res)
            else:
                structures.extend(waveguide)
        else:
            structures.extend(ring)
            structures.extend(waveguide)
            
        # print([(st is None) for st in structures])
            
        structures.extend(self.create_coupler(tolerance))
        structures.extend(self.create_name())
        thermodes, therm_holepos = self.create_thermode()
        structures.extend(thermodes)
        electrodes, elec_holepos = self.create_electrode()
        structures.extend(electrodes)
        holepos = therm_holepos + elec_holepos
        for structure in structures:
            structure.translate(bias[0], bias[1])
        for kk in range(len(holepos)):
            pos = holepos[kk]
            holepos[kk] = (pos[0] + bias[0], pos[1] + bias[1])
        return structures, holepos