# EQdecluster

Declustering of earthquake catalogs using window methods

By: Álvaro González
    alvaro@geonaut.eu
    www.geonaut.eu
    
Coded in Fortran 2008.

Declusters an earthquake catalogue using four algorithms:

  - Gardner, J.K. and Knopoff, L. (1974). “Is the sequence of earthquakes in Southern California,
    with aftershocks removed, Poissonian?”. Bulletin of the Seismological Society of America,
    64, 1363–1367.

  - Grünthal, G. (1985), in: van Stiphout, T.; Zhuang, J. & Marsan, D. (2012). "Theme V -Models and Techniques for
    Analysing Seismicity". Community Online Resource for Statistical Seismicity Analysis. http://www.corssa.org

  - Peláez, J.A.; Chourak, M.; Tadili, B.A.; Brahim, L.A.; Hamdache, M.; López Casado, C.
    & Martínez Solares, J.M. (2007). "A catalog of main Moroccan earthquakes from 1045 to 2005".
    Seismological Research Letters, 78, 614-621, 2007

  - IGN-UPM (2013; digital edition 2017) "Actualización de mapas de peligrosidad sísmica de España 2012".
    Centro Nacional de Información Geográfica (Madrid, Spain).
    https://www.ign.es/resources/acercaDe/libDigPub/ActualizacionMapasPeligrosidadSismica2012.pdf


Notes:

  * The coefficients of Gardner & Knopoff (1974) and Peláez et al. (2007)
    have been here interpolated for increased precision.

  * For the method of IGN-UPM (2013), the magnitude ranges used are the recommended ones
    to avoid large gaps in the curves defining the windows.

  * Declustering is done according to the "mainshock" method, detailed by Luen & Stark (2012):
    "Consider the events in chronological order. If the i-th event is in the window
    of a preceding larger shock that has not already been deleted, delete it.
    If a larger shock is in the window of the ith event, delete the ith event.
    Otherwise, retain the ith event." (also detailed by Peláez et al., 2007).
    (Luen, B. & Stark, P.B. (2012). “Poisson tests of declustered catalogues”.
    Geophysical Journal International, 189, 691–700.)

Input file:

  * catalogue.dat, with nine columns and no headers:
    Year, month, day, hour, minute, second, longitude, latitude, magnitude

Output file:

  * results.dat, with four columns, corresponding to each algorithm above,
    with one line per earthquake, in the same order as the original catalogue
    (1=mainshock, 0=aftershock or foreshock).
    Includes one line of headers.  

Compilation:

  * Suggested compilation with gfortran:
    $ gfortran -static-libgfortran -O3 EQdecluster.f08 -o EQdecluster.exe
    (The flag -O3 optimizes for speed.)

Execution:

  * Copy and run the .exe file in the directory containing the file catalogue.dat

__________________________________________________________________________________________________________________________________

     Copyright (c) 2024 by Álvaro González

     This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with this program.  If not, see <http://www.gnu.org/licenses/>.

     Contact info: alvaro@geonaut.eu
__________________________________________________________________________________________________________________________________
