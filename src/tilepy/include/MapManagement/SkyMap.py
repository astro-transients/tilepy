#
# Copyright (C) 2016-2024  tilepy developers (Monica Seglar-Arroyo, Halim Ashkar, Fabian Schussler, Mathieu de Bony)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import astropy.units as u
import healpy as hp
import numpy as np
from astropy.coordinates import SkyCoord
from scipy.stats import norm

from tilepy.include.PointingTools import Tools

__all__ = ["SkyMap"]


class SkyMap:

    def __init__(self, obspar, mapReader):
        self.raw_map_prob_density = mapReader.getMap("prob")
        if obspar.algorithm == "2D":
            self.is3D = False
        else:
            self.is3D = self.determine3D(obspar, mapReader)
        self.mode = getattr(obspar, "mode", None)

        if self.is3D:
            self.raw_map_dist_mean = mapReader.getMap("distMean")
            self.raw_map_dist_sigma = mapReader.getMap("distSigma")
            self.raw_map_dist_norm = mapReader.getMap("distNorm")

        self.minimumProbCutForGalaxyCatalog = obspar.minimumProbCutForCatalogue
        self.rasterized_map_cache = {}
        self.pix_id_area_cache = {}

    def determine3D(self, obspar, mapReader):
        is3D = True
        if mapReader.has3D:
            dist_mean, dist_std = mapReader.getDistance()
            if dist_mean + 2 * dist_std > obspar.distCut:
                is3D = False

            pix_maximum = np.argmax(self.raw_map_prob_density.data)
            lon, lat = self.raw_map_prob_density.pix2ang(pix_maximum, lonlat=True)
            coordinate = SkyCoord(ra=lon * u.deg, dec=lat * u.deg)
            if Tools.GalacticPlaneBorder(coordinate):
                is3D = False
        else:
            is3D = False
        return is3D

    def getPixIdArea(self, fraction_localisation):

        if fraction_localisation in self.pix_id_area_cache.keys():
            return self.pix_id_area_cache[fraction_localisation]

        sorted_pixel_id = np.flipud(np.argsort(self.raw_map_prob_density.data))
        prob = (
            self.raw_map_prob_density[sorted_pixel_id]
            * self.raw_map_prob_density.pixarea(sorted_pixel_id)
        ).value
        summed_probability = np.cumsum(prob)

        self.pix_id_area_cache[fraction_localisation] = sorted_pixel_id[
            summed_probability <= fraction_localisation
        ]

        return self.pix_id_area_cache[fraction_localisation]

    def getArea(self, fraction_localisation):
        area_vals = self.raw_map_prob_density.pixarea(
            self.getPixIdArea(fraction_localisation)
        )
        return np.sum(area_vals.to(u.deg * u.deg))

    def getMap(self, mapType, nside, scheme="ring"):

        cache_entry = mapType + "_" + str(nside) + "_" + scheme
        if cache_entry in self.rasterized_map_cache.keys():
            return self.rasterized_map_cache[cache_entry]

        if mapType == "prob_density":
            self.rasterized_map_cache[cache_entry] = (
                self.raw_map_prob_density.rasterize(nside=nside, scheme=scheme).data
            )
        elif mapType == "prob":
            self.rasterized_map_cache[cache_entry] = self.getMap(
                "prob_density", nside, scheme
            ) * hp.nside2pixarea(nside)
        elif mapType == "coordinate":
            npix = hp.nside2npix(nside)
            id_pix = np.arange(npix)
            lon, lat = hp.pix2ang(nside, id_pix, nest=(scheme == "nested"), lonlat=True)
            self.rasterized_map_cache[cache_entry] = SkyCoord(
                ra=lon * u.deg, dec=lat * u.deg
            )
        else:
            raise Exception("Unknown type of map")

        return self.rasterized_map_cache[cache_entry]

    def computeGalaxyProbability(self, galaxyCatalog, mangrove=False):

        # Identify the pixel associated to each galaxy
        ra = galaxyCatalog["RAJ2000"]
        dec = galaxyCatalog["DEJ2000"]
        dist = galaxyCatalog["Dist"]
        theta = 0.5 * np.pi - np.deg2rad(dec)
        phi = np.deg2rad(ra)
        pix_id = self.raw_map_prob_density.ang2pix(theta, phi)

        # Compute the probability associated to each galaxy
        if self.is3D:
            galaxyCatalog["dp_dV"] = (
                self.raw_map_prob_density[pix_id]
                * self.raw_map_dist_norm[pix_id]
                * norm(
                    self.raw_map_dist_mean[pix_id], self.raw_map_dist_sigma[pix_id]
                ).pdf(dist)
                * self.raw_map_prob_density.pixarea(pix_id)
            )
        else:
            galaxyCatalog["dp_dV"] = self.raw_map_prob_density[
                pix_id
            ] * self.raw_map_prob_density.pixarea(pix_id)

        return galaxyCatalog
