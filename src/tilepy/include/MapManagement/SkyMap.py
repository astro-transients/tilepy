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
    """
    Representation and utility methods for gravitational-wave localization sky maps.

    This class handles both 2D (probability density) and 3D (distance) sky maps,
    provides area and pixel selection utilities, manages rasterization and
    caching, and can compute probabilities for galaxies in a catalog.
    """

    def __init__(self, obspar, mapReader):
        """
        Initialize the SkyMap object.

        Parameters
        ----------
        obspar : object
            Observation parameters.
        mapReader : object
            Interface for sky map data.

        """
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
        """
        Decide whether the sky map should be treated as 3D.

        Checks map properties and observation parameters to set 3D usage.

        Parameters
        ----------
        obspar : object
            Observation parameters.
        mapReader : object
            Map reading interface.

        Returns
        -------
        is3D : bool
            True if 3D mode is used, False otherwise.

        """

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

    def getPixIdArea(self, fraction_localisation, nside=None, scheme="ring"):
        """
        Return pixel indices covering a specified localization probability.

        Selects pixels, ordered by probability, until the given
        `fraction_localisation` of the total probability is included.

        Parameters
        ----------
        fraction_localisation : float
            Probability threshold (e.g., 0.9 for 90% localization).
        nside : int or None, optional
            HEALPix nside to use for rasterization (default: None).
        scheme : str, optional
             HEALPix ordering scheme, either 'ring' (default) or 'nested'.

        Returns
        -------
        pix_id : ndarray of int
            Pixel indices containing the requested probability fraction.

        """

        cache_line = (
            f"{fraction_localisation}_raw"
            if nside is None
            else f"{fraction_localisation}_{nside}_{scheme}"
        )
        if cache_line in self.pix_id_area_cache.keys():
            return self.pix_id_area_cache[cache_line]

        if nside is None:
            sorted_pixel_id = np.flipud(np.argsort(self.raw_map_prob_density.data))
            prob_sorted = (
                self.raw_map_prob_density[sorted_pixel_id]
                * self.raw_map_prob_density.pixarea(sorted_pixel_id)
            ).value
        else:
            prob = self.getMap("prob", nside=nside, scheme=scheme)
            sorted_pixel_id = np.flipud(np.argsort(prob))
            prob_sorted = prob[sorted_pixel_id]
        summed_probability = np.cumsum(prob_sorted)

        self.pix_id_area_cache[cache_line] = sorted_pixel_id[
            summed_probability <= fraction_localisation
        ]

        return self.pix_id_area_cache[cache_line]

    def getArea(self, fraction_localisation):
        area_vals = self.raw_map_prob_density.pixarea(
            self.getPixIdArea(fraction_localisation)
        )
        return np.sum(area_vals.to(u.deg * u.deg))

    def getMap(self, mapType, nside, scheme="ring"):
        """
        Get a rasterized map of the specified type and pixelization.

        Parameters
        ----------
        mapType : str
            Type of map ('prob_density', 'prob', 'coordinate').
        nside : int
            Desired HEALPix nside resolution.
        scheme : str, optional
            HEALPix ordering scheme: 'ring' (default) or 'nested'.

        Returns
        -------
        The requested map: a pixel array for 'prob' or 'prob_density', or an array of sky coordinates for 'coordinate'.

        Raises
        ------
        Exception
            If `mapType` is not recognized.

        """

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

    def getMaximumProbabilityCoordinates(self):
        """
        Returns the sky coordinates (RA, Dec) of the highest-probability pixel
        in the raw probability density map.
        """
        ipix_max = np.argmax(self.raw_map_prob_density.data)
        lon, lat = self.raw_map_prob_density.pix2ang(ipix_max, lonlat=True)
        return SkyCoord(ra=lon * u.deg, dec=lat * u.deg)

    def computeGalaxyProbability(self, galaxyCatalog, mangrove=False):
        """
        Compute localization probability for each galaxy in a catalog.

        Adds a column `dp_dV` to the input table with the probability
        density at each galaxy's position (and distance if in 3D).

        Parameters
        ----------
        galaxyCatalog : astropy.table.Table or pandas.DataFrame
            Table of galaxies with 'RAJ2000', 'DEJ2000', and 'Dist' columns.
        mangrove : bool, optional
            Flag to use the mangrove method of weighting by the mass of the host galaxy.

        Returns
        -------
        galaxyCatalog : same type as input
            Input table with additional 'dp_dV' probability column.

        """

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
