#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later

from pycesam.constants import lsun_dict, rsun_dict, msun_dict, ggrav_dict

class CConstants( object ):

    lsun     = 3.846e33
    msun     = 1.98919e33
    rsun     = 6.9599e10
    gmsun    = 1.32712438e26
    ggrav    = 6.67168e-8

    def __init__( self ):
        super( CConstants, self ).__init__( )

    def set( self, params ):
        self.lsun     = lsun_dict[params.nom_ctes]
        self.msun     = msun_dict[params.nom_ctes]
        self.rsun     = rsun_dict[params.nom_ctes]
        self.ggrav    = ggrav_dict[params.nom_ctes]
        self.gmsun    = self.ggrav * self.msun
