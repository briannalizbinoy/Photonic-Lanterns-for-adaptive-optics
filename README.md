Ground-based telescopes are limited by atmospheric turbulence, which distorts incoming optical wavefronts
and degrades image resolution. Traditional wavefront sensors (WFSs) used in adaptive optics (AO), such as
the Shack-Hartmann sensor, face challenges including sensitivity limitations, alignment requirements, and non-
common path aberrations. In this work, I investigate the use of mode-selective photonic lanterns (MSPLs) as
focal-plane WFSs, offering a compact alternative with a potential for multiplexing.
I modeled a six-output MSPL using a Maxwell-equations solver to simulate its response to various input
wavefronts, including individual linearly polarized (LP) fiber modes and Kolmogorov-distorted wavefronts.
The output intensity distributions of the single-mode fibers (SMFs) were shown to encode the spatial structure
of the input wavefront. A dataset of 20,000 simulated distorted wavefronts and their corresponding Zernike
coefficients was generated to later train a feedforward neural network (FFNN) to learn the inverse mapping from
lantern output to wavefront shape.
