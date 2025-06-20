.. _overview:

=======
Purpose
=======

The growing number and improving sensitivity of gravitational-wave detectors allow us to probe deeper into the universe, increasing the observed volume by a factor proportional
to the cube of the detection horizon (see :doc:`Observing Capabilities <userguide:capabilities>`).
As a result, the detection rate of compact binary coalescences (BBH, NSBH, and BNS) by the international gravitational-wave network (IGWN) continues to rise.

However, gravitational-wave localization regions often remain large — ranging from several hundred to a few thousand square degrees :footcite:`2022ApJ...924...54P`, :footcite:`2023ApJ...958..158K`.
This poses a major challenge for identifying potential electromagnetic counterparts such as short gamma-ray bursts (GRBs) and kilonovae, which are fast-evolving and require rapid response.

|tilepy| addresses this challenge by generating optimized observation plans, starting from :doc:`IGWN-provided skymaps <userguide:tutorial/skymaps>` and incorporating realistic telescope-specific constraints.
It currently includes ground-based facilities such as H.E.S.S., LST, and both CTA-South and CTA-North.

It accounts for constraints and bacground noise like:

- Sky background from zodiacal light (sunlight scattered by interplanetary dust),
- Milky Way dust extinction,
- Galactic diffuse ultraviolet background (modeled with a piecewise cosecant profile),
- Airmass and altitude limits,
- Twilight phase or solar altitude conditions,
- Sun and Moon exclusion angles,
- Dynamic constraints such as telescope slew time.

These constraints are customized per telescope, enabling more focused and efficient tiling of the sky.
This improves the chances of capturing the early light from fast transients and maximizes the scientific return of follow-up observations.

.. note::

   For a full description of **tilepy**, please see :footcite:`Seglar-Arroyo_2024`.




==========
References
==========

.. footbibliography::
