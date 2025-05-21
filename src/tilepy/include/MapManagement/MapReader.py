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


import logging
import os
import time
import traceback
from abc import ABC, abstractmethod
from urllib.parse import urlparse
from urllib.request import urlretrieve
from urllib.error import HTTPError

import astropy.units as u
import healpy as hp
import mhealpy as mh
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

##################################################################################################
#                        Read Healpix map from fits file                                         #
##################################################################################################

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)

__all__ = ["MapReader", "create_map_reader"]


def validate_source_params(obspar):
    if obspar.mode == "gaussian":
        for attr in ["raSource", "decSource", "sigmaSource"]:
            if getattr(obspar, attr, None) is None:
                raise ValueError(f"{attr} must be defined in 'gaussian' mode")


class SimpleHealpixMap:
    def __init__(self, data, nside, ordering="nested"):
        self.data = data
        self.nside = nside
        self.unit = u.Unit("sr^-1")
        self.ordering = ordering

    def pixarea(self, pix_ids):
        area_sr = hp.nside2pixarea(self.nside)
        return np.full_like(pix_ids, area_sr, dtype=float) * u.sr

    def __getitem__(self, key):
        return self.data[key]

    def rasterize(self, nside, scheme="NESTED"):
        downgraded_data = hp.ud_grade(
            self.data, nside_out=nside, order_in="NESTED", order_out=scheme.upper()
        )
        return SimpleHealpixMap(downgraded_data, nside, ordering=scheme)


def create_map_reader(obspar):
    """
    Factory function to create the appropriate MapReader subclass based on the observation parameters.

    Priority:
    1. Use obspar.mode if set (explicit is better)
    2. Otherwise, infer the map type from the skymap filename

    Parameters
    ----------
    obspar : ObservationParameters
        The full configuration object with mode, skymap path, etc.

    Returns
    -------
    MapReader subclass instance
        One of: GaussianMapReader, LocProbMapReader, HealpixMapReader
    """
    mode = getattr(obspar, "mode", None)
    if mode is not None:
        mode = mode.lower()
        if mode == "gaussian":
            return GaussianMapReader(obspar)
        elif mode == "locprob":
            return LocProbMapReader(obspar)
        elif mode == "healpix":
            return HealpixMapReader(obspar)
        else:
            raise ValueError(f"Unknown obspar.mode = '{mode}'")

    if hasattr(obspar, "skymap"):
        if obspar.skymap is not None:
            filename = os.path.basename(obspar.skymap)
            if "glg_locprob_all" in filename:
                return LocProbMapReader(obspar)
            else:
                return HealpixMapReader(obspar)

    raise ValueError(
        "Unable to determine appropriate MapReader. "
        "Set `obspar.mode` to one of ['gaussian', 'locprob', 'healpix'], "
        "or ensure `obspar.skymap` is provided and recognizable."
    )


class MapReader(ABC):
    def __init__(self, obspar):
        self.obspar = obspar
        self.name_event = obspar.event_name or "undefined"
        self.has3D = False
        self.prob_density = False

    @abstractmethod
    def getMap(self, mapType):
        pass

    def getDistance(self):
        raise NotImplementedError("This reader does not support 3D distance.")

    def downloadMap(self, download_max_nb_try=1, time_wait_retry=60):
        """
        Download the map from a url

        :param download_max_nb_try: If more than 1, will attempt multiple time to download the file in case of error the first time
        :type download_max_nb_try: int
        :param time_wait_retry: the time to wait between each attempt at downloading the file
        :type time_wait_retry: int

        :return: filename
        """

        if download_max_nb_try < 1:
            download_max_nb_try = 1

        filename = self.url.split("/")[-1]
        logger.info(f"The filename is {filename}")

        file_exist = os.path.isfile(filename)
        if file_exist:
            logger.warning("The file exists, it will be re-downloaded")

        for i in range(download_max_nb_try):
            if "png" in self.url:  # Change Fermi-GBM url(if png) to fit format
                self.url = self.url.replace("png", "fit")
            try:
                urlretrieve(self.url, filename)
                break

            except HTTPError:

                if i == (download_max_nb_try - 1):
                    logger.error("Issue to download map from url")

                logger.info(
                    f"Map not available, waiting for {time_wait_retry}s before a new attempt"
                )
                time.sleep(time_wait_retry)

            except Exception as e:

                if i == (download_max_nb_try - 1):
                    logger.error("Issue to download map from url")

                traceback.print_exc()
                if not file_exist:
                    raise e
                else:
                    logger.warning("The existing file will be used")

        return filename


class GaussianMapReader(MapReader):
    def __init__(self, obspar):
        self.mode = "gaussian"
        obspar.mode = self.mode
        self.name_event = getattr(obspar, "event_name", "gaussian_event")
        self.has3D = False
        self.prob_density = True

        self.validate_source_params(obspar)

        ra_deg = float(obspar.raSource)
        dec_deg = float(obspar.decSource)
        sigma_deg = float(obspar.sigmaSource)
        self.nside = int(getattr(obspar, "nside", 64))

        npix = hp.nside2npix(self.nside)

        theta_c = np.radians(90.0 - dec_deg)
        phi_c = np.radians(ra_deg)
        center_vec = np.array(
            [
                np.sin(theta_c) * np.cos(phi_c),
                np.sin(theta_c) * np.sin(phi_c),
                np.cos(theta_c),
            ]
        )

        pix_vecs = np.array(hp.pix2vec(self.nside, np.arange(npix), nest=True))
        dots = np.dot(center_vec, pix_vecs)
        ang_dist_rad = np.arccos(np.clip(dots, -1.0, 1.0))
        ang_dist_deg = np.degrees(ang_dist_rad)

        prob_unnorm = np.exp(-0.5 * (ang_dist_deg / sigma_deg) ** 2)

        pixarea_sr = hp.nside2pixarea(self.nside)
        norm_factor = np.sum(prob_unnorm * pixarea_sr)
        prob_density = prob_unnorm / norm_factor

        self.simulated_map = SimpleHealpixMap(
            prob_density, self.nside, ordering="nested"
        )

    def getMap(self, mapType):
        if mapType != "prob":
            raise ValueError("Only 'prob' map type is supported for GaussianMapReader.")
        return self.simulated_map

    def validate_source_params(self, obspar):
        if (
            not hasattr(obspar, "raSource")
            or not hasattr(obspar, "decSource")
            or not hasattr(obspar, "sigmaSource")
        ):
            raise ValueError(
                "Gaussian mode requires 'raSource', 'decSource', and 'sigmaSource' in obspar."
            )


class LocProbMapReader(MapReader):
    def __init__(self, obspar):
        self.mode = "locprob"
        obspar.mode = self.mode
        self.name_event = getattr(obspar, "event_name", "locprob_event")
        self.has3D = False
        self.prob_density = True

        self.nside = getattr(obspar, "nside", 128)

        skymap_localisation = obspar.skymap
        parsed_localisation = urlparse(skymap_localisation)
        self.is_remote = parsed_localisation.scheme != ""

        self.url = None
        self.skymap_filename = skymap_localisation
        if self.is_remote:
            self.url = skymap_localisation
            max_retry = obspar.downloadMaxRetry
            time_wait_retry = obspar.downloadWaitPeriodRetry
            self.skymap_filename = self.downloadMap(
                download_max_nb_try=max_retry + 1, time_wait_retry=time_wait_retry
            )

        if not os.path.isfile(self.skymap_filename):
            raise Exception("Map file not found: " + self.skymap_filename)

        self.simulated_map = self._convert_locprob_to_healpix(
            self.skymap_filename, self.nside
        )

    def getMap(self, mapType):
        if mapType != "prob":
            raise ValueError("Only 'prob' map type is supported for LocProbMapReader.")
        return self.simulated_map

    def _convert_locprob_to_healpix(self, filename, nside):
        hdulist = fits.open(filename)
        if hdulist[1].data is not None:
            hdu = hdulist[1]
        else:
            hdu = hdulist[0]

        data = hdu.data
        if data is None:
            raise ValueError(f"No image data found in {filename}")

        data = data.astype(float)
        wcs = WCS(hdu.header)

        data[np.isnan(data)] = 0.0
        data[data < 0] = 0.0
        total_prob = np.sum(data)
        if total_prob > 0:
            data /= total_prob

        npix = hp.nside2npix(nside)
        healpix_map = np.zeros(npix)

        ny, nx = data.shape
        for y in range(ny):
            for x in range(nx):
                p = data[y, x]
                if p == 0:
                    continue
                world = wcs.pixel_to_world(x, y)
                ra = world.ra.deg
                dec = world.dec.deg
                theta = np.radians(90 - dec)
                phi = np.radians(ra)
                pix = hp.ang2pix(nside, theta, phi, nest=True)
                healpix_map[pix] += p

        pixarea = hp.nside2pixarea(nside)
        norm = np.sum(healpix_map) * pixarea
        healpix_map /= norm

        return SimpleHealpixMap(healpix_map, nside, ordering="nested")


class HealpixMapReader(MapReader):
    def __init__(self, obspar):
        self.mode = "healpix"
        obspar.mode = self.mode
        self.name_event = getattr(obspar, "event_name", "undefined")
        self.has3D = False
        self.prob_density = False

        skymap_localisation = obspar.skymap
        parsed_localisation = urlparse(skymap_localisation)
        self.is_remote = parsed_localisation.scheme != ""

        self.url = None
        self.skymap_filename = skymap_localisation
        if self.is_remote:
            self.url = skymap_localisation
            max_retry = obspar.downloadMaxRetry
            time_wait_retry = obspar.downloadWaitPeriodRetry
            self.skymap_filename = self.downloadMap(
                download_max_nb_try=max_retry + 1, time_wait_retry=time_wait_retry
            )

        if not os.path.isfile(self.skymap_filename):
            raise Exception("Map file not found: " + self.skymap_filename)

        self.fits_map = fits.open(self.skymap_filename)
        self.id_hdu_map = self.getMapHDUId()
        self.identifyColumns()

        if obspar.event_name is None:
            obspar.event_name = self.getSourceName()
        self.name_event = obspar.event_name

    def getSourceName(self):
        if "OBJECT" in self.fits_map[self.id_hdu_map].header:
            return self.fits_map[self.id_hdu_map].header["OBJECT"]
        return "undefined"

    def identifyColumns(self):
        header = self.fits_map[self.id_hdu_map].header
        nb_column = header["TFIELDS"]
        self.offset_column = 1 if "TTYPE0" not in header else 0

        self.id_prob = None
        self.unit_prob = None
        self.id_dist_mean = self.id_dist_sigma = self.id_dist_norm = None
        self.unit_dist_mean = self.unit_dist_sigma = self.unit_dist_norm = None

        for i in range(self.offset_column, nb_column + self.offset_column):
            colname = header.get(f"TTYPE{i}")
            unit = header.get(f"TUNIT{i}")

            if colname in ["PROB", "T", "PROBABILITY", "PROBDENSITY"]:
                self.id_prob = i
                self.unit_prob = u.Unit(unit) if unit else u.dimensionless_unscaled
            elif colname == "DISTMU":
                self.id_dist_mean = i
                self.unit_dist_mean = u.Unit(unit) if unit else u.Mpc
            elif colname == "DISTSIGMA":
                self.id_dist_sigma = i
                self.unit_dist_sigma = u.Unit(unit) if unit else u.Mpc
            elif colname == "DISTNORM":
                self.id_dist_norm = i
                self.unit_dist_norm = u.Unit(unit) if unit else u.Unit("Mpc^-2")

        if self.unit_prob and self.unit_prob.is_equivalent(u.Unit("pix^-1")):
            self.unit_prob = u.dimensionless_unscaled

        self.prob_density = not self.unit_prob.is_equivalent(u.dimensionless_unscaled)

        self.has3D = all(
            x is not None
            for x in [self.id_dist_mean, self.id_dist_sigma, self.id_dist_norm]
        )

    def getMap(self, mapType):
        field_map = {
            "prob": (self.id_prob, self.unit_prob, self.prob_density, u.Unit("sr^-1")),
            "distMean": (self.id_dist_mean, self.unit_dist_mean, False, u.Unit("Mpc")),
            "distSigma": (
                self.id_dist_sigma,
                self.unit_dist_sigma,
                False,
                u.Unit("Mpc"),
            ),
            "distNorm": (
                self.id_dist_norm,
                self.unit_dist_norm,
                False,
                u.Unit("Mpc^-2"),
            ),
        }

        if mapType not in field_map:
            raise Exception(f"Unknown or unsupported map type: {mapType}")

        field_id, unit, is_density, target_unit = field_map[mapType]

        if field_id is None:
            raise Exception(f"Map type '{mapType}' not available in this FITS file")

        raw_map = mh.HealpixMap.read_map(
            self.skymap_filename,
            field=field_id - self.offset_column,
            hdu=self.id_hdu_map,
            density=is_density,
        )

        quantity = raw_map.data * unit
        if not unit.is_equivalent(target_unit):
            # If unit is dimensionless but we're expecting density, compute manually
            if unit.is_equivalent(u.dimensionless_unscaled) and target_unit == u.Unit(
                "1/sr"
            ):
                pixarea = raw_map.pixarea()
                quantity = quantity / pixarea  # Convert to density manually
                raw_map._density = True
        raw_map._data = quantity.to_value(target_unit)
        raw_map._unit = target_unit

        return raw_map

    def getDistance(self):
        if not self.has3D:
            raise Exception("No distance information available")

        header = self.fits_map[self.id_hdu_map].header
        return header["DISTMEAN"], header["DISTSTD"]

    def getMapHDUId(self):
        for i, hdu in enumerate(self.fits_map):
            if hdu.header.get("XTENSION") == "BINTABLE":
                return i
        raise Exception("No valid BINTABLE HDU found for HEALPix map")


class MapReaderLegacy:
    def __init__(self, obspar):
        self.mode = getattr(obspar, "mode", "file")
        self.name_event = "undefined"
        self.has3D = False
        self.prob_density = False

        # Validate the the parameters early so to fail if something is not properly configured:
        validate_source_params(obspar)

        # -------------------------
        # 1. GAUSSIAN MODE
        # -------------------------
        if self.mode == "gaussian":
            """
            Required obspar fields:
            - ra, dec (degrees)
            - sigma_deg
            - nside
            """
            self.generate_gaussian_map(obspar)
            self.name_event = obspar.event_name or "gaussian_event"
            return

        # -------------------------
        # 2. locprob FIT MODE
        # -------------------------
        skymap_localisation = obspar.skymap
        parsed_localisation = urlparse(skymap_localisation)
        self.is_remote = parsed_localisation.scheme != ""

        self.url = None
        self.skymap_filename = skymap_localisation
        if self.is_remote:
            self.url = skymap_localisation
            max_retry = obspar.downloadMaxRetry
            time_wait_retry = obspar.downloadWaitPeriodRetry
            self.skymap_filename = self.downloadMap(
                download_max_nb_try=max_retry + 1, time_wait_retry=time_wait_retry
            )

        if not os.path.isfile(self.skymap_filename):
            raise Exception("Map file not found: " + self.skymap_filename)

        # Early GBM format — convert to HEALPix
        # its important to have this before the healpix fits map mode,
        # because it will be bypassing the normal healpix map reading
        if "glg_locprob_all" in os.path.basename(self.skymap_filename):
            self.simulated_map = self.convert_locprob_to_healpix(
                self.skymap_filename, nside=128
            )
            self.name_event = obspar.event_name or "gbm_locprob"
            self.prob_density = True
            self.has3D = False
            return

        # -------------------------
        # 3. HEALPix FITS MAP MODE
        # -------------------------
        self.fits_map = fits.open(self.skymap_filename)
        self.id_hdu_map = self.getMapHDUId()
        self.name_event = self.getSourceName()
        if obspar.event_name is not None:
            self.name_event = obspar.event_name
        if obspar.event_name is None:
            obspar.event_name = self.name_event
        self.identifyColumns()

    def convert_locprob_to_healpix(self, filename, nside=128):
        """
        Convert a GBM glg_locprob_all_*.fit 2D probability grid to a HEALPix map.
        """
        hdulist = fits.open(filename)
        hdu = hdulist[1] if hdulist[1].data is not None else hdulist[0]

        data = hdu.data.astype(float)
        wcs = WCS(hdu.header)

        # Clean invalid entries
        data[np.isnan(data)] = 0.0
        data[data < 0] = 0.0

        # Normalize (assuming flat prior)
        total = np.sum(data)
        if total > 0:
            data /= total

        ny, nx = data.shape

        # Generate all pixel grid coordinates
        y_indices, x_indices = np.mgrid[0:ny, 0:nx]
        flat_x = x_indices.flatten()
        flat_y = y_indices.flatten()
        flat_p = data.flatten()

        # Remove zero entries early
        mask = flat_p > 0
        flat_x = flat_x[mask]
        flat_y = flat_y[mask]
        flat_p = flat_p[mask]

        # Convert all to sky coordinates at once
        sky_coords = wcs.pixel_to_world(flat_x, flat_y)
        ra = sky_coords.ra.deg
        dec = sky_coords.dec.deg

        theta = np.radians(90 - dec)
        phi = np.radians(ra)
        pix = hp.ang2pix(nside, theta, phi, nest=True)

        # Accumulate into HEALPix map
        healpix_map = np.bincount(pix, weights=flat_p, minlength=hp.nside2npix(nside))

        # Normalize as probability density
        pixarea = hp.nside2pixarea(nside)
        norm = np.sum(healpix_map) * pixarea
        healpix_map /= norm

        return SimpleHealpixMap(healpix_map, nside, ordering="nested")

    def generate_gaussian_map(self, obspar):

        ra_deg = float(obspar.raSource)
        dec_deg = float(obspar.decSource)
        sigma_deg = float(obspar.sigmaSource)

        self.nside = int(getattr(obspar, "nside", 64))

        npix = hp.nside2npix(self.nside)

        # Convert center to unit vector
        theta_c = np.radians(90.0 - dec_deg)
        phi_c = np.radians(ra_deg)
        center_vec = np.array(
            [
                np.sin(theta_c) * np.cos(phi_c),
                np.sin(theta_c) * np.sin(phi_c),
                np.cos(theta_c),
            ]
        )

        # Compute angular distance to each pixel center
        pix_vecs = np.array(hp.pix2vec(self.nside, np.arange(npix), nest=True))
        dots = np.dot(center_vec, pix_vecs)
        ang_dist_rad = np.arccos(np.clip(dots, -1.0, 1.0))
        ang_dist_deg = np.degrees(ang_dist_rad)

        # Build Gaussian profile (not normalized yet)
        prob_unnorm = np.exp(-0.5 * (ang_dist_deg / sigma_deg) ** 2)

        # Normalize to make it a **density**: total int(p,domega) = 1
        pixarea_sr = hp.nside2pixarea(self.nside)  # steradians
        norm_factor = np.sum(prob_unnorm * pixarea_sr)
        prob_density = prob_unnorm / norm_factor  # unit: sr⁻¹

        self.simulated_map = SimpleHealpixMap(
            prob_density, self.nside, ordering="nested"
        )
        self.name_event = f"Gaussian_RA{ra_deg}_Dec{dec_deg}"
        self.prob_density = True
        self.has3D = False

    def getMapHDUId(self):
        id_hdu_map = -1
        for i in range(len(self.fits_map)):
            if (
                "XTENSION" in self.fits_map[i].header
                and self.fits_map[i].header["XTENSION"] == "BINTABLE"
            ):
                if id_hdu_map != -1:
                    raise Exception(
                        "Multiple map detected, please provide a file with only one map"
                    )
                else:
                    id_hdu_map = i
        if id_hdu_map == -1:
            raise Exception("No map detected, please provide a file with a map")
        return id_hdu_map

    def identifyColumns(self):
        nb_column = self.fits_map[self.id_hdu_map].header["TFIELDS"]

        # Identify if column numeration start at 0 or 1
        self.offset_column = 1
        if "TTYPE0" in self.fits_map[self.id_hdu_map].header:
            self.offset_column = 0

        if nb_column == 1:
            self.id_prob = 0
            self.unit_prob = u.dimensionless_unscaled
            self.id_dist_mean = None
            self.unit_dist_mean = None
            self.id_dist_sigma = None
            self.unit_dist_sigma = None
            self.id_dist_norm = None
            self.unit_dist_norm = None
        else:
            self.unit_prob = None
            self.id_dist_mean = None
            self.unit_dist_mean = None
            self.id_dist_sigma = None
            self.unit_dist_sigma = None
            self.id_dist_norm = None
            self.unit_dist_norm = None
            for i in range(self.offset_column, nb_column + self.offset_column):
                columns_name = self.fits_map[self.id_hdu_map].header["TTYPE" + str(i)]
                unit_information = ("TUNIT" + str(i)) in self.fits_map[
                    self.id_hdu_map
                ].header
                if columns_name in ["PROB", "T", "PROBABILITY", "PROBDENSITY"]:
                    self.id_prob = i
                    if unit_information:
                        self.unit_prob = u.Unit(
                            self.fits_map[self.id_hdu_map].header["TUNIT" + str(i)]
                        )
                    else:
                        columns_name = self.fits_map[self.id_hdu_map].header[
                            "TTYPE" + str(i)
                        ]
                        if columns_name in ["T", "PROBABILITY"]:
                            self.unit_prob = u.dimensionless_unscaled
                        else:
                            self.unit_prob = u.Unit("sr^-1")
                elif columns_name == "DISTMU":
                    self.id_dist_mean = i
                    if unit_information:
                        self.unit_dist_mean = u.Unit(
                            self.fits_map[self.id_hdu_map].header["TUNIT" + str(i)]
                        )
                    else:
                        self.unit_dist_mean = u.Mpc
                elif columns_name == "DISTSIGMA":
                    self.id_dist_sigma = i
                    if unit_information:
                        self.unit_dist_sigma = u.Unit(
                            self.fits_map[self.id_hdu_map].header["TUNIT" + str(i)]
                        )
                    else:
                        self.unit_dist_sigma = u.Mpc
                elif columns_name == "DISTNORM":
                    self.id_dist_norm = i
                    if unit_information:
                        self.unit_dist_norm = u.Unit(
                            self.fits_map[self.id_hdu_map].header["TUNIT" + str(i)]
                        )
                    else:
                        self.unit_dist_norm = u.Unit("Mpc^-2")

        # Correct unit if needed
        if self.unit_prob.is_equivalent(u.Unit("pix^-1")):
            self.unit_prob = u.dimensionless_unscaled

        # Check if probabilities is probabilities density
        self.prob_density = True
        if self.unit_prob.is_equivalent(u.dimensionless_unscaled):
            self.prob_density = False

        if (
            self.id_dist_mean is None
            or self.id_dist_mean is None
            or self.id_dist_norm is None
        ):
            self.has3D = False
        else:
            self.has3D = True

    def getSourceName(self):
        """
        Get the source name from the contents of the fits file

        :return: name
        """

        name = "undefined"
        if "OBJECT" in self.fits_map[self.id_hdu_map].header:
            name = self.fits_map[self.id_hdu_map].header["OBJECT"]
        if (
            "SENDER" in self.fits_map[self.id_hdu_map].header
            and self.fits_map[self.id_hdu_map].header["SENDER"]
            == "IceCube Collaboration"
        ):
            name = (
                str(self.fits_map[self.id_hdu_map].header["RUNID"])
                + "_"
                + str(self.fits_map[self.id_hdu_map].header["EVENTID"])
            )
        # if the event if from LVK and the URL is from GraceDB, get the superevent name
        if "ORIGIN" in self.fits_map[self.id_hdu_map].header.keys():
            if (
                "LIGO/Virgo/KAGRA" in self.fits_map[self.id_hdu_map].header["ORIGIN"]
                and self.is_remote
            ):
                if "https://gracedb.ligo.org/api/superevents" in self.url:
                    name = self.url.split("/")[5]
        # if the event if from Fermi-GBM, get the GBM name (i.e. replace GRB with bn)
        if "TELESCOP" in self.fits_map["PRIMARY"].header.keys():
            if self.fits_map["PRIMARY"].header["TELESCOP"] == "GLAST":
                name = (
                    self.fits_map[self.id_hdu_map].header["OBJECT"].replace("GRB", "bn")
                )

        return name

    def downloadMap(self, download_max_nb_try=1, time_wait_retry=60):
        """
        Download the map from a url

        :param download_max_nb_try: If more than 1, will attempt multiple time to download the file in case of error the first time
        :type download_max_nb_try: int
        :param time_wait_retry: the time to wait between each attempt at downloading the file
        :type time_wait_retry: int

        :return: filename
        """

        if download_max_nb_try < 1:
            download_max_nb_try = 1

        filename = self.url.split("/")[-1]
        logger.info(f"The filename is {filename}")

        file_exist = os.path.isfile(filename)
        if file_exist:
            logger.warning("The file exists, it will be re-downloaded")

        for i in range(download_max_nb_try):
            if "png" in self.url:  # Change Fermi-GBM url(if png) to fit format
                self.url = self.url.replace("png", "fit")
            try:
                urlretrieve(self.url, filename)
                break

            except Exception as e:

                print(f"Exception {e == HTTPError}")
                if (
                    i == (download_max_nb_try - 1)
                    or str(e) != "HTTP Error 404: Not Found"
                ):
                    logger.error("Issue to download map from url")
                    traceback.print_exc()
                    if not file_exist:
                        raise e
                    else:
                        logger.warning("The existing file will be used")
                else:
                    logger.info(
                        f"Map not available, waiting for {time_wait_retry} before a new attempt"
                    )
                    time.sleep(time_wait_retry)

        return filename

    def getMap(self, mapType):
        if hasattr(self, "simulated_map"):
            if mapType == "prob":
                return self.simulated_map
            else:
                raise Exception(
                    f"Map type '{mapType}' not available in simulated map mode"
                )

        if self.mode == "gaussian":
            if mapType != "prob":
                raise Exception("Only 'prob' map type supported in gaussian mode.")
            return self.simulated_map

        if mapType == "prob":
            raw_map = mh.HealpixMap.read_map(
                self.skymap_filename,
                field=self.id_prob - self.offset_column,
                hdu=self.id_hdu_map,
                density=True,
            )
            if not self.prob_density:
                raw_map._data = (
                    raw_map.data * self.unit_prob / raw_map.pixarea()
                ).to_value(u.Unit("sr^-1"))
                raw_map._unit = u.Unit("sr^-1")
            else:
                raw_map._data = (raw_map.data * self.unit_prob).to_value(
                    u.Unit("sr^-1")
                )
        elif mapType == "distMean" and self.has3D:
            raw_map = mh.HealpixMap.read_map(
                self.skymap_filename,
                field=self.id_dist_mean - self.offset_column,
                hdu=self.id_hdu_map,
                density=False,
            )
            raw_map._data = (raw_map.data * self.unit_dist_mean).to_value(u.Unit("Mpc"))
        elif mapType == "distSigma" and self.has3D:
            raw_map = mh.HealpixMap.read_map(
                self.skymap_filename,
                field=self.id_dist_sigma - self.offset_column,
                hdu=self.id_hdu_map,
                density=False,
            )
            raw_map._data = (raw_map.data * self.unit_dist_sigma).to_value(
                u.Unit("Mpc")
            )
        elif mapType == "distNorm" and self.has3D:
            raw_map = mh.HealpixMap.read_map(
                self.skymap_filename,
                field=self.id_dist_norm - self.offset_column,
                hdu=self.id_hdu_map,
                density=False,
            )
            raw_map._data = (raw_map.data * self.unit_dist_norm).to_value(
                u.Unit("Mpc^-2")
            )
        elif not self.has3D and (
            mapType == "distMean" or mapType == "distSigma" or mapType == "distNorm"
        ):
            raise Exception("No distance information available")
        else:
            raise Exception("Unknown type of map")

        return raw_map

    def getDistance(self):
        if self.has3D:
            return (
                self.fits_map[self.id_hdu_map].header["DISTMEAN"],
                self.fits_map[self.id_hdu_map].header["DISTSTD"],
            )
        else:
            raise Exception("No distance information available")
