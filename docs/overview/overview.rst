.. _overview:

=======
Purpose
=======


.. dropdown::  tilepy concept

   The growing number and improving sensitivity of gravitational-wave detectors allow us to probe deeper into the universe, increasing the observed volume by a factor proportional
   to the cube of the detection horizon (see :doc:`Observing Capabilities <userguide:capabilities>`).
   As a result, the detection rate of compact binary coalescences (BBH, NSBH, and BNS) by the international gravitational-wave network (IGWN) continues to rise.

   However, gravitational-wave localization regions often remain large, ranging from several hundred to a few thousand square degrees :footcite:`2022ApJ...924...54P`, :footcite:`2023ApJ...958..158K`.
   This poses a major challenge for identifying potential electromagnetic counterparts such as short gamma-ray bursts (GRBs) and kilonovae, which are fast-evolving and require rapid response.

   |tilepy| addresses this challenge by generating optimized observation plans, starting from :doc:`IGWN-provided skymaps <userguide:tutorial/skymaps>` and incorporating realistic telescope-specific constraints.
   It currently includes ground-based facilities such as H.E.S.S., LST, and both CTA-South and CTA-North.

   It accounts for constraints and background noise like:

   .. tab-set::

      .. tab-item:: Current

         **Backgrounds:**

         - Milky Way dust extinction
         - Galactic diffuse ultraviolet background

         **Constraints:**

         - Airmass and altitude limits
         - Twilight phase or solar altitude conditions
         - Sun and Moon exclusion angles

      .. tab-item:: Planned

         **Backgrounds:**

         - Sky background from zodiacal light (sunlight scattered by interplanetary dust)

         **Constraints:**

         - Dynamic constraints such as telescope slew time

   .. admonition:: Telescope-specific optimization
      :class: important

      These constraints are customized per telescope, enabling more focused and efficient tiling of the sky.
      This improves the chances of capturing the early light from fast transients and maximizes the scientific return of follow-up observations.

   .. note::

      For a full description of **tilepy**, please see :footcite:`Seglar-Arroyo_2024`.



.. dropdown::  How far can GW counterparts be difficult to follow-up and rare to detect?


   The first phase of the fourth observing run (O4a) of the LIGO-Virgo-KAGRA Collaboration took place from May 24, 2023,
   to January 16, 2024. During this period, only the two LIGO detectors, the LIGO Hanford Observatory (LHO, H1) and the LIGO Livingston Observatory (LLO, L1), were operational.
   These detectors identified 87 gravitational wave (GW) candidate events with false alarm rates (FAR) below 1 :math:`\mathrm{yr}^{-1}`
   for which detailed source property estimates have been provided. Based on the inferred component masses, these candidates are consistent
   with mergers of binary black holes (BBH) and neutron star–black hole (NSBH) systems. Combined, this brings the total number of compact binary coalescences (CBCs) observed since the first observing run (O1) to 159.
   Notably, no binary neutron star (BNS) mergers were detected during O4a.

   .. admonition:: Motivation
      :class: info

      Those kinds of difficulties motivated the implementation of |tilepy| to optimize GRBs, early ultraviolet (UV) and optical follow-up of gravitational wave events.


==========
References
==========

.. footbibliography::
