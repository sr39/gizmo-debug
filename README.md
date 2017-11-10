Welcome!

This is GIZMO: a flexible, multi-method mhd+gravity+radiation code. The code includes a hydro solver using Lagrangian mesh-free finite-volume Godunov methods (or SPH), and self-gravity solved fast with hybrid PM-Tree methods and fully-adaptive gravitational softenings. Other physics include: magnetic fields (ideal and non-ideal), anisotropic conduction and viscosity, sub-grid turbulent diffusion, radiative cooling, cosmological integration, sink particles, dust-gas mixtures, cosmic rays, degenerate equations of state, radiation-hydrodynamics, galaxy/star/black hole formation and feedback, self-interacting and scalar-field dark matter, on-the-fly structure finding, and more. See the User Guide for an up-to-date list.

The code is descended from P-SPH/P-GADGET, itself descended from GADGET-3 (so a huge debt owes to the work of Volker Springel), and many of the GADGET conventions remain (for compatibility with GADGET-compatible codes). See the source code for appropriate attribution of the code elements. 

The BitBucket site hosts code updates and news, so check back regularly. For code issues, feature requests, or bugs, please post to the GIZMO Google Group: https://groups.google.com/d/forum/gizmo-code

Basic Rules: 

1. The main reference for methods, code setup, use policies and citations is the User Guide, available as `gizmo_documentation.html` in the `scripts` folder as part of this repository, or via download [here](http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html)

Read it! The code is extensively documented. Important things you must know before running are there. Most questions I get asked are already answered there.

2. Access to the development code does **not** imply permission to use any modules identified as proprietary either in the User Guide or `Template_Config.sh` file (or elsewhere). The development code can only be used or distributed with explicit permission from the code authors. Many of the non-public modules are proprietary and developed by students for on-going research; it is not acceptable to use or share these routines without first obtaining the explicit permission of both the lead code author and author(s) of the relevant routines. If in doubt, ask. Anyone violating these terms will have code access immediately revoked.

The public version of the code is free software, distributed under the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html). You may freely distribute and copy the public code. You may also modify it as you wish, and distribute these modified versions as long as you indicate prominently any changes you made in the original code, and as long as you leave the copyright notices, and the no-warranty notice intact. Please read the General Public License for more details. Note that the authors retain their copyright on the code. The public code is available [here](http://www.tapir.caltech.edu/~phopkins/public/gizmo_public.tgz).

If you use any version of the code, please reference the code paper at: http://arxiv.org/abs/1409.7395 (Hopkins 2015); you should also reference Volker Springel's GADGET paper (Springel, 2005, MNRAS, 364, 1105) for the domain decomposition and N-body algorithms. Appropriate citations for specific modules are described in the Users Guide.
