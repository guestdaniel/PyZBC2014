from setuptools import setup, Extension
from setuptools.command.build_py import build_py as build_py_orig
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
import os
import sys
import subprocess

so_name = "libzbc2014"
ext = {
    "linux": ".so",
    "darwin": ".so",
    "win32": ".dll"
}[sys.platform if sys.platform != "darwin" else "darwin"]
# Where the compiled lib will reside after building:
source_lib_path = os.path.join("src", "pyzbc2014", "_lib", so_name + ext)
ext_modules=[Extension('pyzbc2014._lib.dummy_ext_for_platlib', 
                       sources=[os.path.join('src', 'pyzbc2014', '_lib', 'dummy.c')])]

class build_py(build_py_orig):
    def run(self):
        this_dir = os.path.abspath(os.path.dirname(__file__))
        lib_dir = os.path.join(this_dir, "src", "pyzbc2014", "_lib")
        # C sources are still in model/
        model_dir = os.path.join(this_dir, "src", "pyzbc2014", "model")
        sources = [
            os.path.join(model_dir, "complex.c"),
            os.path.join(model_dir, "model_IHC.c"),
            os.path.join(model_dir, "model_Synapse.c"),
        ]
        full_lib_path = os.path.join(lib_dir, so_name + ext)
        if not os.path.exists(full_lib_path):
            print(f"Compiling C library: {' '.join(sources)} -> {full_lib_path}")
            if sys.platform == "win32":
                cmd = [
                    "gcc", "-shared", "-O3", "-o", full_lib_path, *sources
                ]
            else:
                cmd = [
                    "gcc", "-fPIC", "-O3", "-shared", "-o", full_lib_path, *sources
                ]
            subprocess.check_call(cmd)
        super().run()

class bdist_wheel(_bdist_wheel):
    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False

setup(
    package_dir={"": "src"},
    cmdclass={
        "build_py": build_py,
        "bdist_wheel": bdist_wheel,
    },
    ext_modules=ext_modules,
    package_data={"pyzbc2014._lib": ["libzbc2014.*"]},
    zip_safe=False,
)