# -*- coding: utf-8 -*-

__version__ = "0.0.1"

from optpy.function import(
    get_ld,
    create_box,
    create_posmark,
    create_cutmark,
    create_PhCR,
    create_electest,
    generate_field,
    merge_all
    )

from optpy.array import Array

from optpy.image import convert_image

from optpy.text import convert_text

layer_datatype_dict = {
    "ld_chip": {"layer": 0, "datatype": 0},
    "ld_film": {"layer": 1, "datatype": 0},
    "ld_posmark": {"layer": 2, "datatype": 1},
    "ld_logo": {"layer": 3, "datatype": 1},
    "ld_ring": {"layer": 4, "datatype": 2},
    "ld_waveguide": {"layer": 5, "datatype": 2},
    "ld_tip": {"layer": 6, "datatype": 2},
    "ld_taper": {"layer": 6, "datatype": 2},
    "ld_anchor": {"layer": 7, "datatype": 2},
    "ld_text": {"layer": 8, "datatype": 2},
    "ld_cutmark": {"layer": 9, "datatype": 2},
    "ld_thermode": {"layer": 10, "datatype": 3},
    "ld_hole1": {"layer": 11, "datatype": 4},
    "ld_hole2": {"layer": 12, "datatype": 5},
    "ld_pad": {"layer": 13, "datatype": 6} 
}

create_cutting_mark = create_cutmark