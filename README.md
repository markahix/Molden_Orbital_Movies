# Molden_Orbital_Movies
Molden Orbital Movies takes a set of molden files from a QM/MM MD simulation, checks the phases of each orbital compared to the previous step, and produces a VMD input file that can be used to generate a movie of the MD trajectory with a selected orbital shown.

---

#### Usage

To process Molden files and generate VMD input:

`Molden_Orbital_Movies.py -p </directory/with/molden/files> -n <orbital number> [--color1 <color>] [--color2 <color>]`

To execute VMD script:

`vmd -e orbital_trajectories.vmd`

Generates frames named `temp.#####.tga` using the TachyonLOptixInternal renderer with VMD, numbered by frame.  You can use ImageMagick or GIMP to generate a video/animated GIF from these frames.

---
