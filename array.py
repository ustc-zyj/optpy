# -*- coding: utf-8 -*-
import numpy as np
import gdspy as gp

from optpy.text import create_name

from optpy.function import create_pad, create_wire

from optpy.waveguides import (
    Basic_Waveguide,
    Bend_Waveguide,
    Bend_Waveguide_v2,
    BasicDual_Waveguide,
    Reverse_Waveguide,
    Multiplexing_Waveguide,
    Roller_Waveguide,
    MZRoller_Waveguide,
    Triport_Waveguide,
    TriportMR_Waveguide,
    RollerMR_Waveguide,
    AdjCoupling_Waveguide,
    UPB_Waveguide,
    MMI_MZI,
    DC_MZI,
    )

_type_dict = {
    "BASIC": "Basic_Waveguide",
    "BEND": "Bend_Waveguide",
    "BEND2": "Bend_Waveguide_v2",
    "DUAL": "BasicDual_Waveguide",
    "REVRS": "Reverse_Waveguide",
    "MTPLX": "Multiplexing_Waveguide",
    "ROLLR": "Roller_Waveguide",
    "MZRLR": "MZRoller_Waveguide",
    "3PORT": "Triport_Waveguide",
    "3PMR": "TriportMR_Waveguide",
    "RLRMR": "RollerMR_Waveguide",
    "ADJCP": "AdjCoupling_Waveguide",
    "UPBWG": "UPB_Waveguide",
    "MMZI": "MMI_MZI",
    "DMZI": "DC_MZI",
    }

class Array(object):
    
    __slots__ = (
        "nm",
        "np",
        "ns",
        "nn",
        "wt",
        "cx",
        "cy",
        "dy",
        "wv",
        "nv",
        "cv",
        "tv",
        "pv",
        "ev",
        )
    
    def __init__(
            self,
            array_name,
            array_name_pos,
            array_name_size,
            waveguide_type,
            waveguide_num,
            waveguide_period,
            wpara,
            npara = None,
            cpara = None,
            tpara = None,
            ppara = None,
            epara = None,
            ):
        
        self.nm = array_name            
        self.np = array_name_pos        
        self.ns = array_name_size
        
        if waveguide_type in _type_dict:
            self.wt = _type_dict[waveguide_type]
        else:
            raise ValueError(
                "Invalid waveguide type: " + waveguide_type)
        
        if isinstance(waveguide_num, int) and waveguide_num > 0:
            self.nn = waveguide_num
        else:
            raise ValueError("Number of waveguides should be a positive integer")
        
        if waveguide_period > 0:
            self.dy = waveguide_period
        else:
            raise ValueError('Waveguide period should be positive')
        
        if isinstance(wpara, dict):
            for key, value in wpara.items():
                if isinstance(value, list) and len(value) != self.nn:
                    raise ValueError(
                        "Waveguide variable " + key + \
                            " length should match waveguide number")
            self.wv = wpara  
        else:
            raise ValueError("Waveguide variables should in dict form")
        
        if isinstance(npara, dict):
            for key, value in npara.items():
                if isinstance(value, list) and len(value) != self.nn:
                    raise ValueError(
                        "Name variable " + key + \
                            " length should match waveguide number")
            self.nv = npara
        elif npara is None:
            self.nv = None
        else:
            raise ValueError("Name variables should be a dict or none")
        
        if isinstance(cpara, dict):
            for key, value in cpara.items():
                if isinstance(value, list) and len(value) != self.nn:
                    raise ValueError(
                        "Coupling variable " + key + \
                            " length should match waveguide number")
            self.cv = cpara   
        elif cpara is None:
            self.cv = None
        else:
            raise ValueError("Coupling variables should be a dict or none")
            
        if isinstance(tpara, dict):
            for key, value in tpara.items():
                if isinstance(value, list) and len(value) != self.nn:
                    raise ValueError(
                        "Thermal variable " + key + \
                            " length should match waveguide number")
            self.tv = tpara   
        elif tpara is None:
            self.tv = None
        else:
            raise ValueError("Thermal variables should be a dict or none")
            
        if isinstance(epara, dict):
            for key, value in epara.items():
                if isinstance(value, list) and len(value) != self.nn:
                    raise ValueError(
                        "Electrical variable " + key + \
                            " length should match waveguide number")
            self.ev = epara   
        elif epara is None:
            self.ev = None
        else:
            raise ValueError("Thermal variables should be a dict or none")
        
        if isinstance(ppara, dict):
            self.pv = ppara
        elif ppara is None:
            self.pv = None
        else:
            raise ValueError("Pad variables should be a dict or none")
    
    def get_opt(self):
        opts = []
        for kk in range(self.nn):
            wpara = {}
            for key, value in self.wv.items():
                if isinstance(value, list):
                    wpara[key] = value[kk]
                else:
                    wpara[key] = value
                    
            if self.nv is not None:
                npara = {}
                for key, value in self.nv.items():
                    if isinstance(value, list):
                        npara[key] = value[kk]
                    else:
                        npara[key] = value
                wpara["npara"] = npara
                
            if self.cv is not None:
                cpara = {}
                for key, value in self.cv.items():
                    if isinstance(value, list):
                        cpara[key] = value[kk]
                    else:
                        cpara[key] = value
                wpara["cpara"] = cpara
            
            if self.tv is not None:
                tpara = {}
                for key, value in self.tv.items():
                    if isinstance(value, list):
                        tpara[key] = value[kk]
                    else:
                        tpara[key] = value
                wpara["tpara"] = tpara
                
            if self.ev is not None:
                epara = {}
                for key, value in self.ev.items():
                    if isinstance(value, list):
                        epara[key] = value[kk]
                    else:
                        epara[key] = value
                wpara["epara"] = epara
                
            opt = eval(self.wt)(**wpara)
            opts.append(opt)
        return opts
    
    def create_opt(self, tolerance=1e-3):
        opts = self.get_opt()
        structures = []
        holepos = []
        disp_list = [((self.nn-1)/2 - kk) * self.dy for kk in range(self.nn)]
        for kk in range(self.nn):
            structure, holes = opts[kk].create((0, disp_list[kk]), tolerance)            
            structures.extend(structure)
            holepos.extend(holes)
        if self.ns > 0:
            names = create_name(
                name_pos = self.np,
                name_text = self.nm,
                name_size = self.ns)
            structures.extend(names)
        return structures, holepos
    
    def create_elec(self, holepos):
        ppara = self.pv.copy()
        pad_num = len(self.pv['pad_center'])
        portpos = []
        structures = []
        # w1 = (self.tv['hole_radius'] + self.tv['hole_num'] * 0.5) * 2
        w2 = self.pv['wire_width']
        tl = self.pv['wire_taper_length']
        pg = self.pv['pad_gap']
        distance1 = self.pv['distance1']
        distance2 = self.pv['distance2']
        inverse = self.pv['inverse']
        extra_dx1 = self.pv['extra_dx1']
        extra_dx2 = self.pv['extra_dx2']
        for kk in range(pad_num):
            for key, value in self.pv.items():
                if isinstance(value, list):
                    ppara[key] = value[kk]
                else:
                    ppara[key] = value
            pads, ports = create_pad(**ppara)
            structures.extend(pads)
            portpos.extend(ports)
        if len(portpos) != len(holepos):
            raise ValueError(
                "Number of pads (%d) does not match the need (%d)"\
                    %(len(portpos), len(holepos)))
        else:
            wire_num = len(portpos)
            
        if self.ev is None and self.tv is not None:
            w1 = (self.tv['hole_radius'] + self.tv['hole_num'] * 0.5) * 2
        elif self.ev is not None:
            w1 = (self.ev['hole_radius'] + self.ev['hole_num'] * 0.5) * 2
        
        if pg != 0:
            for kk in range(wire_num // 2):
                if isinstance(distance1, list):
                    d1 = distance1[kk]
                else:
                    d1 = distance1
                    
                if isinstance(distance2, list):
                    d2 = distance2[kk]
                else:
                    d2 = distance2
                    
                if isinstance(extra_dx1, list):
                    edx1 = extra_dx1[kk]
                else:
                    edx1 = extra_dx1
                    
                if isinstance(extra_dx2, list):
                    edx2 = extra_dx2[kk]
                else:
                    edx2 = extra_dx2
                
                if isinstance(inverse, list):
                    inv = inverse[kk]
                else:
                    inv = inverse
                    
                pp = portpos[kk*2 : kk*2+2]
                hp = holepos[kk*2 : kk*2+2]
                if pp[0][0] > pp[1][0]:
                    pp.reverse()
                if hp[0][0] > hp[1][0]:
                    hp.reverse()
                if inv:
                    hp.reverse()
                
                dx1 = pp[0][0] - hp[0][0]
                dx2 = pp[1][0] - hp[1][0]
                ypm = np.sign(pp[0][1] - hp[0][1])
                dhy1 = hp[0][1] + ypm * d1
                dhy2 = hp[1][1] + ypm * d2
                if edx1 == 0:
                    ip1 = [hp[0], (hp[0][0], dhy1)]
                else:
                    ip1 = [hp[0], (hp[0][0]+edx1, hp[0][1]), (hp[0][0]+edx1, dhy1)]
                if edx2 == 0:
                    ip2 = [hp[1], (hp[1][0], dhy2)]
                else:
                    ip2 = [hp[1], (hp[1][0]+edx2, hp[1][1]), (hp[1][0]+edx2, dhy2)]
                    
                if abs(dx1) < tl:
                    dhy1 = dhy1 + ypm * np.sqrt(tl**2 - dx1**2)
                fp1 = [(pp[0][0], dhy1), pp[0]]
                if abs(dx2) < tl:
                    dhy2 = dhy2 + ypm * np.sqrt(tl**2 - dx2**2)
                fp2 = [(pp[1][0], dhy2), pp[1]]
                    
                wires1 = create_wire(
                    initial_points = ip1, 
                    final_points = fp1, 
                    initial_width = w1,
                    final_width = w2,
                    taper_length = tl)
                wires2 = create_wire(
                    initial_points = ip2, 
                    final_points = fp2, 
                    initial_width = w1,
                    final_width = w2,
                    taper_length = tl)
                structures.extend(wires1)
                structures.extend(wires2)   
        else:
            for kk in range(wire_num):
                if isinstance(distance1, list):
                    d1 = distance1[kk]
                else:
                    d1 = distance1
                    
                pp = portpos[kk]
                hp = holepos[kk]
                dx = pp[0] - hp[0]
                ypm = np.sign(pp[1] - hp[1])
                dhy = hp[1] + ypm * d1
                ip = [hp, (hp[0], dhy)]
                if abs(dx) < tl:
                    dhy = dhy + ypm * np.sqrt(tl**2 - dx**2)
                fp = [(pp[0], dhy), pp]
                wires = create_wire(
                    initial_points = ip, 
                    final_points = fp, 
                    initial_width = w1,
                    final_width = w2,
                    taper_length = tl)
                structures.extend(wires)
        return structures    
    
    def __call__(self, cell, bias=(0,0), rotate_angle_deg=0, tolerance=1e-3):
        structures, opt_holepos = self.create_opt(tolerance)
        if self.pv is not None:
            structures.extend(self.create_elec(opt_holepos))
        if rotate_angle_deg != 0:
            rotate_angle = rotate_angle_deg * np.pi / 180
            for structure in structures:
                structure.rotate(rotate_angle)
        for structure in structures:
            structure.translate(*bias)
        cell.add(structures)