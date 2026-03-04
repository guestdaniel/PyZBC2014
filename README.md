# Introduction
PyZBC2014 is a barebones Python 3 wrapper around the Zilany, Bruce, and Carney (2014) auditory-nerve (AN) model, which is written in C.
Currently, PyZBC2014 is set up as a Python package (i.e., it can be imported like `import pyzbc2014 as zbc` in Python scripts).

PyZBC2014 is licensed under the GNU AGPLv3 license (see `LICENSE`).
Please send any questions or issues to `daniel_guest@urmc.rochester.edu`, or raise an issue on the [GitHub page](https://github.com/guestdaniel/PyZBC2014).

# Install
## Requirements
- Python>=3.10
- numpy>=1.21.0
- scipy>=1.7.2

Ensure you have the required dependencies installed before trying to install this package — pip will likely produce an error if you don't.

## Install package
Use pip to install
```
pip install pyzbc2014
```

## Uninstall package
Use pip to uninstall:
```
pip uninstall pyzbc2014
```

# Usage
PyZBC2014 defines two functions, `sim_ihc_zbc2014` and `sim_anrate_zbc2014`, to simulate inner-hair-cell (IHC) potentials from sound-pressure waveforms or AN instantaneous rates from IHC potentials, respectively.
Comparing to the origianl MATLAB/Mex implementation of the 2014 model, these behave similarly to the `IHCAN` and `Synapse` functions, respectively, with some small differences in behavior and some improvements in quality-of-life.

## Simulate IHC response
To simulate an IHC response, you need to provide a sound-pressure waveform in 1D Numpy vector format to the `sim_ihc_zbc2014` function; it will return an IHC potential waveform of exactly the same size and format.
Parameters can be modified by keyword arguments:
- `cf`: Characteristic frequency (Hz)
- `nrep`: Number of reps to generate (only a single rep is simulated, others are appended copies)
- `fs`: Sampling rate (Hz), should generally be at least 100 kHz, and higher for cat/very-high-CF simulations; note that your stimulus must be sampled at this frequency!
- `cohc`: Status of the outer hair cells, in range of [0, 1] with 1 being fully normal/healthy
- `cihc`: Status of the inner hair cells, in range of [0, 1] with 1 being fully normal/healthy
- `species`: String encoding species to use in [`"cat"`, `"human"`, `"human-glasberg"`]. In the MATLAB/Mex wrapper, `species=1` corresponds to `"cat"`, `species=2` corresponds to `"human"` (and indicates human with the sharp peripheral tuning suggested by Shera and Oxenham), and `species=3` corresponds to `"human-glasberg"` (and indicates broader peripheral tuning suggested by Glasberg and Moore)

The full function call using default keyword arguments is thus:
```
ihc_potential = sim_ihc_zbc2014(px, cf=1e3, nrep=1, fs=100e3, cohc=1.0, cihc=1.0, species="human")
```

## Simulate AN response
To simulate an AN response, you need to provide an IHC potential (from the output of `sim_ihc_zbc2014`) in 1D Numpy vector format to the `sim_anrate_zbc2014` function; it will return a prediction of instantaneous AN firing rate rate of exactly the same size and format.
Parameters can be modified by keyword arguments:
- `cf`: Characteristic frequency (Hz)
- `nrep`: Number of reps to simulate, should always match value used in `sim_ihc_zbc2014`
- `fs`: Sampling rate (Hz), should generally be at least 100 kHz, and higher for cat/very-high-CF simulations, should always match value used for `sim_ihc_zbc2014`
- `fibertype`: String encoding which fiber type to simulate in [`"hsr"`, `"msr"`, `"lsr"`]
- `powerlaw`: String encoding whether to use true powerlaw adaptation (`"true"`) or approximate powerlaw adaptation (`"approx"`). In the MATLAB/Mex wrapper, `implnt=1` corresponds to `"true"` and `implnt=0` corresponds to `"approx"`
- `noisetype`: String encoding whether to use no fractional Gaussian noise (`"none"`), or fresh noise (`"fresh"`). In the MATLAB/Mex wrapper, `noiseType=1` corresponds to `"fresh"`

The full function call using default keyword arguments is thus:
```
anrate = sim_anrate_zbc2014(ihc, cf=1e3, nrep=1, fs=100e3, fibertype="hsr", powerlaw="approx", noisetype="fresh")
```

# Citing
If you use this code, please cite the model paper:

Zilany, M. S., Bruce, I. C., & Carney, L. H. (2014). Updated parameters and expanded simulation options for a model of the auditory periphery. The Journal of the Acoustical Society of America, 135(1), 283-286.


# Developers
This package is setup to build wheels whenever the version is incremented, but you have to manually trigger the action to publish to PyPI. To make any updates, follow these steps:
1. Edit the code, fix bugs, update docs, etc.
2. Bump the version in `src/pyzbc2014/__init__.py`
3. Commit the changes and push to GitHub. For example:
```
git add --all
git commit -m "[commit message]"
git push origin main
```
4. Create and push a matching git tag when you want to build the wheels.
```
git tag v0.0.1  # this tag needs to match the version in `pyzbc2014/__init__.py`
git push origin v0.0.1  # again, match the tag
```
GitHub Actions should then check the version and build the wheels automatically. Make sure it completes these steps successfully. While you're debugging, you may not wish to bump the version every time. You can manually trigger the wheel building action workflow "Check version and publish wheels" in GitHub actions.

To publish the new version to PyPI, the version must have been incremented. After testing the code and checking the wheels all built without error, run the GitHub action "Publish to PyPI".