# This file is part of a variant of version 5.2 of the Zilany, Bruce, and Carney (2014)
# auditory-nerve model. Specifically, this file is part of a Python 3 wrapper for the
# original model code written in the C programming language. Licensing information can be
# found below. See README.md for a list of changes made versus the original version, which
# is available at:
# https://www.urmc.rochester.edu/labs/carney/publications-code/auditory-models.aspx.


# Papers related to this code are cited below:

# Zilany, M. S., & Bruce, I. C. (2006). Modeling auditory-nerve responses for high sound
# pressure levels in the normal and impaired auditory periphery. The Journal of the
# Acoustical Society of America, 120(3), 1446-1466.

# * Zilany, M.S.A., Bruce, I.C., Nelson, P.C., and Carney, L.H. (2009). "A
# Phenomenological model of the synapse between the inner hair cell and auditory
# nerve : Long-term adaptation with power-law dynamics," Journal of the
# Acoustical Society of America 126(5): 2390-2412.

# Ibrahim, R. A., and Bruce, I. C. (2010). "Effects of peripheral tuning
# on the auditory nerve's representation of speech envelope and temporal fine
# structure cues," in The Neurophysiological Bases of Auditory Perception, eds.
# E. A. Lopez-Poveda and A. R. Palmer and R. Meddis, Springer, NY, pp. 429-438.

# Zilany, M.S.A., Bruce, I.C., Ibrahim, R.A., and Carney, L.H. (2013).
# "Improved parameters and expanded simulation options for a model of the
# auditory periphery," in Abstracts of the 36th ARO Midwinter Research Meeting.

# * Zilany, M. S., Bruce, I. C., & Carney, L. H. (2014). Updated parameters and expanded
# simulation options for a model of the auditory periphery. The Journal of the Acoustical
# Society of America, 135(1), 283-286.

# Please cite these papers marked with asterisks (*) if you publish any research results
# obtained with this code or any modified versions of this code.


# Zilany, Bruce, and Carney (2014) auditory-nerve model Python 3 wrapper, Copyright (C) 2024
#   Daniel R. Guest <daniel_guest@urmc.rochester.edu>

# This program is free software: you can redistribute it and/or modify it under the terms of
# the GNU Affero General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License along with this
# program. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import scipy as sp
import ctypes
from numpy.ctypeslib import ndpointer
import warnings


def get_lib_path():
    from importlib.resources import files

    lib_names = []
    for name in files('pyzbc2014._lib').iterdir():
        if name.name.startswith('libzbc2014') and (
                name.name.endswith('.so') or
                name.name.endswith('.pyd') or
                name.name.endswith('.dll')):
            lib_names.append(str(name))

    if not lib_names:
        raise FileNotFoundError("Compiled library not found.")

    # If multiple matches, return the first
    return lib_names[0]



def sim_ihc_zbc2014(
        px,
        cf=1e3,
        nrep=1,
        fs=100e3,
        cohc=1.0,
        cihc=1.0,
        species="human",
):
    """ Simulates inner-hair-cell response to sound-pressure waveform.

    Performs some basic input checking and then passes sound-pressure-waveform and
    pre-allocated output array to the C IHCAN function via the ctypes library.

    Args:
      px:
        Sound pressure waveform (1D vector of floats, units of Pa)
      cf:
        Characteristic frequency (Hz)
      nrep:
        Number of reps to generate (only a single rep is simulated, others are appended copies)
      fs:
        Sampling rate (Hz)
      cohc:
        Status of the outer hair cells, in range of [0, 1] with 1 being fully normal/healthy
      cihc:
        Status of the inner hair cells, in range of [0, 1] with 1 being fully normal/healthy
      species:
        String, either "cat", "human", or "human-glasberg" (corresponding to MATLAB/Mex wrapper species=1, species=2, and species=3 respectively)

    """
    # First, enforce assumptions about inputs
    assert np.ndim(px) == 1  # input is 1D vector
    assert nrep >= 1  # number of reps is geq 1
    assert 0 <= cohc <= 1  # C_OHC is between zero and one
    assert 0 <= cihc <= 1  # C_IHC is between zero and one
    assert species in ["cat", "human", "human-glasberg"]  # species is one of available options
    if species == "cat":
        assert .125e3 <= cf <= 40e3  # for cat, CF must be between 0.125 and 40 kHz
    else:
        assert .125e3 <= cf <= 20e3  # for human, CF must be between 0.125 and 20 kHz

    # Emit warnings
    if fs < 100e3:
        warnings.warn("Note, time-domain resolution is less than recommended (sampling rate > 100 kHz, sampling period < 1e-5 s)")

    # Allocate empty storage for output
    totalstim = len(px)
    ihcout = np.zeros(totalstim*nrep)  # we need a vector of size totalstim*nrep to store results

    # Map from species string to species integer
    match species:
        case "cat":
            species_int = 1
        case "human":
            species_int = 2
        case "human-glasberg":
            species_int = 3

    # Open library, fetch IHCAN function, declare input types for call
    lib_path = get_lib_path()
    lib = ctypes.cdll.LoadLibrary(lib_path)
    fun = lib.IHCAN
    fun.argtypes = [
        ndpointer(ctypes.c_double),
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_int,
        ndpointer(ctypes.c_double),
    ]

    # Place function call and return IHC output waveform
    fun(px, cf, nrep, 1/fs, totalstim, cohc, cihc, species_int, ihcout)
    return ihcout


def sim_anrate_zbc2014(
        ihc,
        cf=1e3,
        nrep=1,
        fs=100e3,
        fibertype="hsr",
        powerlaw="true",
        noisetype="fresh",
):
    """ Simulates AN firing rate response to inner-hair-cell potential

    Performs some basic input checking and then passes IHC waveform and inputs to
    the C Synapse function via the ctypes library.

    Args:
      ihc:
        Inner-hair-cell potential (1D vector of floats, units a.u.)
      cf:
        Characteristic frequency (Hz)
      nrep:
        Number of reps to simulate (must match input to sim_ihc_zbc2014)
      fibertype:
        Whether to simulate high-spont ("hsr"), medium-spont ("msr"), or low-spont ("lsr") fibers
      powerlaw:
        Whether to use true ("true") or approximate ("approx") implementation of powerlaw adaptation
      noisetype:
        Whether to use no fractional Gaussian noise ("none") or fresh fractional Gaussian noise ("fresh")
    """
    # First, enforce assumptions about inputs
    assert np.ndim(ihc) == 1  # input is 1D vector
    assert nrep >= 1  # number of reps is geq 1
    assert fibertype in ["hsr", "msr", "lsr"]
    assert powerlaw in ["true", "approx"]
    assert noisetype in ["none", "fresh"]

    # Second, map from fibertype string to spont value
    match fibertype:
        case "hsr":
            spont = 100.0
        case "msr":
            spont = 4.0
        case "lsr":
            spont = 0.1

    # Third, map from implnt string to integer value
    match powerlaw:
        case "true":
            implnt = 1.0
        case "approx":
            implnt = 0.0

    # Emit warnings
    if fs < 100e3:
        warnings.warn("Note, time-domain resolution is less than recommended (sampling rate > 100 kHz, sampling period < 1e-5 s)")

    # Allocate empty storage for output and for noise input
    # (length for noise input is determined based on magic equation extracted from source code)
    synout = np.zeros(len(ihc))
    len_noise = int(np.ceil((len(ihc) + 2 * np.floor(7500 / (cf / 1e3))) * 1/fs * 10e3))

    # Synthesize fGn based on noisetype param
    match noisetype:
        case "none":
            fGn = np.zeros(len_noise)
        case "fresh":
            fGn = ffGn(len(ihc), 1/fs, 0.9, fibertype)

    # Open library, fetch IHCAN function, declare input types for call
    lib_path = get_lib_path()
    lib = ctypes.cdll.LoadLibrary(lib_path)
    fun = lib.Synapse
    fun.argtypes = [
        ndpointer(ctypes.c_double),  # ihcout
        ndpointer(ctypes.c_double),  # randNums
        ctypes.c_double,             # tdres
        ctypes.c_double,             # cf
        ctypes.c_int,                # totalstim
        ctypes.c_int,                # nrep
        ctypes.c_double,             # spont
        ctypes.c_double,             # implnt
        ctypes.c_double,             # sampFreq
        ndpointer(ctypes.c_double)   # synout
    ]

    # Place function call and return output waveform
    fun(ihc, fGn, 1/fs, cf, int(len(ihc)/nrep), nrep, spont, implnt, 10e3, synout)

    # Return synout passed through the pointwise nonlinearity mapping from pre-refractory
    # to post-refractory rates
    return synout / (1.0 + 0.75e-3 * synout)


def ffGn(N, tdres, Hinput, fibertype):
    from scipy.signal import resample
    """
    Generates fractional Gaussian noise (fGn) based on the specified parameters.

    Adapted from original code provided in the 2014 model release, as well as the Python
    translation in the defunct cochlea Python package. Documentation does not extend much
    beyond that offered by the original code; for a more detailed explanation of the
    function and a better implementation of it, see the function `ffGn_rochester` at
    https://osf.io/6bsnt/. An important note is that the parameter values listed in the 2009
    model paper (https://doi.org/10.1121/1.3238250) are incorrect, and instead the values
    below should be use (and match was was present in the 2009/2014 code releases).

    Parameters:
    N (int): Number of points to generate.
    tdres (float): Time-domain resolution (i.e., 1/fs; s)
    Hinput (float): Hurst parameter, must be in the range [0, 2].
    spont (float): Fiber type/spont group of the AN fiber to be simulated in [hsr, msr, lsr], used to determine noise sigma

    Returns:
    np.ndarray: Array of generated noise values, of size (N, )
    """
    assert (N > 0)
    assert (tdres < 1)
    assert (Hinput >= 0) and (Hinput <= 2)

    # Downsampling No. of points to match with those of Scott Jackson (tau 1e-1)
    resamp = int(np.ceil(1e-1 / tdres))
    nop = N
    N = int(np.ceil(N / resamp) + 1)
    if N < 10:
        N = 10

    # Determine whether fGn or fBn should be produced.
    if Hinput <= 1:
        H = Hinput
        fBn = 0
    else:
        H = Hinput - 1
        fBn = 1

    # Calculate the fGn.
    if H == 0.5:
        # If H=0.5, then fGn is equivalent to white Gaussian noise.
        y = np.random.randn(N)
    else:
        Nfft = int(2 ** np.ceil(np.log2(2*(N-1))))
        NfftHalf = np.round(Nfft / 2)

        k = np.concatenate( (np.arange(0,NfftHalf), np.arange(NfftHalf,0,-1)) )
        Zmag = 0.5 * ( (k+1)**(2*H) -2*k**(2*H) + np.abs(k-1)**(2*H) )

        Zmag = np.real(np.fft.fft(Zmag))
        assert np.all(Zmag >= 0)

        Zmag = np.sqrt(Zmag)

        Z = Zmag * (np.random.randn(Nfft) + 1j*np.random.randn(Nfft))

        y = np.real(np.fft.ifft(Z)) * np.sqrt(Nfft)

        y = y[0:N]

        # Convert the fGn to fBn, if necessary.
        if fBn == 1:
            y = np.cumsum(y)

        # Resampling to match with the AN model
        y = resample(y, resamp*len(y))

        if fibertype == "lsr":
            sigma = 3
        elif fibertype == "msr":
            sigma = 30
        elif fibertype == "hsr":
            sigma = 200
        else:
            ValueError("fibertype not recognized, must be in [hsr, msr, lsr]")

        y = y*sigma

        return y[0:nop]
