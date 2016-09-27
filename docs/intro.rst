Introduction
============


Eye of Thundera
---------------

.. image:: Figures/eye_thundera.png
   :width: 150pt


"Eye of Thundera - give me sight beyond sight".
 
When Lion-O's friends were in danger, he invoked the power of Thundera through the Sword of Omens. The sword has the mystical Eye of Thundera, that gave to Lion-O the sight beyond sight.

Wavelet analysis is similar to the Eye of Thundera, in the sense that it'll give you the power to localize a pulse in frequency and time domain - sight beyond stationarity.

This guide includes a description of a Python module named waipy that calculates the Continuous Wavelet Transform (CWT) and significance tests based on Torrence and Compo (1998). It provides further a Cross Wavelet Analysis (CWA) based on Maraun and Kurths(2004).


Building the puzzle ...

Why/when should I use the wavelet analysis?
-------------------------------------------

How can anyone turn a 1D to 2D information?

The code will explain to you!

The wavelet analysis is used for detecting and characterizing its possible singularities, and in particular the continuous wavelet transform is well suited for analyzing the local differentiability of a function (Farge, 1992).

"Therefore the wavelet analysis or synthesis can be performed locally on the signal, as opposed to the Fourier transform which is inherently nonlocal due to the space-filling nature of the trigonometric functions. " (Farge,1992).

                +----------+-----------------+----------------+
                |          |    Fourier      |     Wavelet    |
                +==========+=================+================+
                |context   |    sound        |     image      |
                +----------+-----------------+----------------+ 
                |character |   stationary    |  nonstationary |
                +----------+-----------------+----------------+
                |to see    | global features |  singularities |
                +----------+-----------------+----------------+


Choose the right glasses for what you want to see !

.. image:: Figures/lennon_glasses.png
   :width: 150pt



.. note::
    The Morlet wavelet is used as default in this code.

.. image:: Figures/wavelet.png
   :width: 250pt
   :alt: Morlet wavelet with real and imaginary part.

For more Information:
---------------------

.. image:: Figures/info.png
   :width: 350pt

Fonte: Domingues (2012)


Papers:
-------

Farge, M. 1992. Wavelet transforms and their applications to turbulence. Annu. Rev. Mech., 24: 395-457

Domingues, M. O.; Kaibar, M.K. 2012. Wavelet biortogonais. Revista brasileira de Ensino de Física,n.3, 34: 3701
