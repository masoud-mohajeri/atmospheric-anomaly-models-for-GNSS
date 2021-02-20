# atmospheric-anomaly-models-for-GNSS
klobuchar and collins models calculation

#Klobuchar Ionospheric Model

![Klobuchar image](https://gssc.esa.int/navipedia/images/e/ed/NeQuickIonoVTECmap.jpeg)

GPS satellites broadcast the parameters of the Klobuchar ionospheric model for single frequency users.
The Klobuchar model was designed to minimise user computational complexity and user computer storage as
far as to keep a minimum number of coefficients to transmit on satellite-user link.

This broadcast model is based on an empirical approach and is estimated to reduce
about the 50% RMS ionospheric range error worldwide. It is assumed that the electron content is concentrated 
in a thin layer at 350 kilometres in height. Thence, the slant delay is computed from the vertical delay at the 
Ionospheric Pierce Point (IPP) multiplying by a obliquity factor, i.e., the mapping function .

#Tropospheric Delay

Troposphere is the atmospheric layer placed between earth's surface and an altitude of about 60 kilometres.

The effect of the troposphere on the GNSS signals appears as an extra delay in the measurement of the signal 
traveling from the satellite to receiver. This delay depends on the temperature, pressure, humidity as well as
the transmitter and receiver antennas location and, according to :

Δ=∫straight line(n−1)dl(1)

it can be written as:

T=∫(n−1)dl=10−6∫Ndl(2)

where n is the refractive index of air and N=10−6(n−1) is the refractivity.
The refractivity can be divided in hydrostatic, i.e., Dry gases (mainly N2 and O2), and wet, i.e., Water vapour, 
components N=Nhydr+Nwet.

Each of these components has different effects on GNSS signals. The main feature of the troposphere is that it is
a non dispersive media with respect to electromagnetic waves up to 15GHz, i.e., the tropospheric effects are not frequency
dependent for the GNSS signals. Thence, the carrier phase and code measurements are affected by the same delay.
