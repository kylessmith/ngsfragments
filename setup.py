from __future__ import division, print_function, absolute_import
import numpy as np
import ailist
import glob
import platform
import re
import subprocess
import sys
import sysconfig
import os
import collections
from contextlib import contextmanager
from distutils import log
from setuptools import setup, Command
from cy_build import CyExtension as Extension, cy_build_ext as build_ext
try:
    import cython
    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False

IS_PYTHON3 = sys.version_info.major >= 3


@contextmanager
def changedir(path):
    save_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(save_dir)


def run_configure(option):
    sys.stdout.flush()
    try:
        retcode = subprocess.call(
            " ".join(("./configure", option)),
            shell=True)
        if retcode != 0:
            return False
        else:
            return True
    except OSError as e:
        return False


def run_make_print_config():
    stdout = subprocess.check_output(["make", "-s", "print-config"])
    if IS_PYTHON3:
        stdout = stdout.decode("ascii")

    make_print_config = {}
    for line in stdout.splitlines():
        if "=" in line:
            row = line.split("=")
            if len(row) == 2:
                make_print_config.update(
                    {row[0].strip(): row[1].strip()})
    return make_print_config


@contextmanager
def set_compiler_envvars():
    tmp_vars = []
    for var in ['CC', 'CFLAGS', 'LDFLAGS']:
        if var in os.environ:
            print ("# ngsfragments: (env) {}={}".format(var, os.environ[var]))
        elif var in sysconfig.get_config_vars():
            value = sysconfig.get_config_var(var)
            print ("# ngsfragments: (sysconfig) {}={}".format(var, value))
            os.environ[var] = value
            tmp_vars += [var]

    try:
        yield
    finally:
        for var in tmp_vars:
            del os.environ[var]


def configure_library(library_dir, env_options=None, options=[]):

    configure_script = os.path.join(library_dir, "configure")

    on_rtd = os.environ.get("READTHEDOCS") == "True"
    # RTD has no bzip2 development libraries installed:
    if on_rtd:
        env_options = "--disable-bz2"

    if not os.path.exists(configure_script):
        raise ValueError(
            "configure script {} does not exist".format(configure_script))

    with changedir(library_dir), set_compiler_envvars():
        if env_options is not None:
            if run_configure(env_options):
                return env_options

        for option in options:
            if run_configure(option):
                return option

    return None


def distutils_dir_name(dname):
    """Returns the name of a distutils build directory
    see: http://stackoverflow.com/questions/14320220/
               testing-python-c-libraries-get-build-path
    """
    f = "{dirname}.{platform}-{version[0]}.{version[1]}"
    return f.format(dirname=dname,
                    platform=sysconfig.get_platform(),
                    version=sys.version_info)


class clean_ext(Command):
    description = "clean up Cython temporary files"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        objs = glob.glob(os.path.join("ngsfragments", "libc*.c"))
        if objs:
            self.announce("removing 'ngsfragments/libc*.c' ({} Cython objects)".format(len(objs)), level=log.INFO)
        for obj in objs:
            os.remove(obj)


# Python 3
MyFileNotFoundError = FileNotFoundError

# Library name
libname = "ngsfragments"
# Build type
build_type = "optimized"
#build_type="debug"

# Descriptions of package
SHORTDESC = "Python package for Next Generation Sequencing fragment manipulation"
#DESC = """A python package for Next Generation Sequencing fragment manipulation."""
with open("README.md", "r") as fh:
    long_description = fh.read()

# Directories (relative to the top-level directory where setup.py resides) in which to look for data files.
datadirs  = ("tests", "data")
# File extensions to be considered as data files. (Literal, no wildcards.)
dataexts  = (".py", ".pyx", ".pxd", ".c", ".h", ".h5", ".txt", ".bed", ".gz", ".bed.gz")
# Standard documentation to detect (and package if it exists).
standard_docs     = ["README", "LICENSE", "TODO", "CHANGELOG", "AUTHORS"]
standard_doc_exts = [".md", ".rst", ".txt", ""]

#########################################################
# Init
#########################################################

# check for Python 2.7 or later
if sys.version_info < (2,7):
    sys.exit('Sorry, Python < 2.7 is not supported')

# Gather user-defined data files
datafiles = []
getext = lambda filename: os.path.splitext(filename)[1]
for datadir in datadirs:
    datafiles.extend( [(root, [os.path.join(root, f) for f in files if getext(f) in dataexts])
                       for root, dirs, files in os.walk(datadir)] )

# Add standard documentation (README et al.), if any, to data files
detected_docs = []
for docname in standard_docs:
    for ext in standard_doc_exts:
        filename = "".join( (docname, ext) )  # relative to the directory in which setup.py resides
        if os.path.isfile(filename):
            detected_docs.append(filename)
datafiles.append( ('.', detected_docs) )


# Extract __version__ from the package __init__.py
import ast
init_py_path = os.path.join(libname, '__init__.py')
version = '0.0.0'
try:
    with open(init_py_path) as f:
        for line in f:
            if line.startswith('__version__'):
                version = ast.parse(line).body[0].value.s
                break
        else:
            print( "WARNING: Version information not found in '%s', using placeholder '%s'" % (init_py_path, version), file=sys.stderr )
except MyFileNotFoundError:
    print( "WARNING: Could not find file '%s', using placeholder version information '%s'" % (init_py_path, version), file=sys.stderr )


#########################################################
# Set up modules
#########################################################


# How to link against HTSLIB
# shared:   build shared chtslib from builtin htslib code.
# external: use shared libhts.so compiled outside of
#           ngsfragments
# separate: use included htslib and include in each extension
#           module. No dependencies between modules and works with
#           setup.py install, but wasteful in terms of memory and
#           compilation time. Fallback if shared module compilation
#           fails.

HTSLIB_MODE = os.environ.get("HTSLIB_MODE", "shared")
HTSLIB_LIBRARY_DIR = os.environ.get("HTSLIB_LIBRARY_DIR", None)
HTSLIB_INCLUDE_DIR = os.environ.get("HTSLIB_INCLUDE_DIR", None)
HTSLIB_CONFIGURE_OPTIONS = os.environ.get("HTSLIB_CONFIGURE_OPTIONS", None)
HTSLIB_SOURCE = None

package_list = ['ngsfragments',
                'ngsfragments.read_sam',
                'ngsfragments.peak_calling',
                'ngsfragments.segment',
                'ngsfragments.plot']
package_dirs = {'ngsfragments': 'ngsfragments',
                'ngsfragments.read_sam': 'ngsfragments/read_sam',
                'ngsfragments.peak_calling': 'ngsfragments/peak_calling',
                'ngsfragments.segment': 'ngsfragments/segment',
                'ngsfragments.plot': 'ngsfragments/plot'}

# list of config files that will be automatically generated should
# they not already exist or be created by configure scripts in the
# subpackages.
config_headers = []

# If cython is available, the ngsfragments will be built using cython from
# the .pyx files. If no cython is available, the C-files included in the
# distribution will be used.
if HAVE_CYTHON:
    print ("# ngsfragments: cython is available - using cythonize if necessary")
    source_pattern = "ngsfragments/libc%s.pyx"
else:
    print ("# ngsfragments: no cython available - using pre-compiled C")
    source_pattern = "ngsfragments/libc%s.c"

    # Exit if there are no pre-compiled files and no cython available
    fn = source_pattern % "htslib"
    if not os.path.exists(fn):
        raise ValueError(
            "no cython installed, but can not find {}."
            "Make sure that cython is installed when building "
            "from the repository"
            .format(fn))

print ("# ngsfragments: htslib mode is {}".format(HTSLIB_MODE))
print ("# ngsfragments: HTSLIB_CONFIGURE_OPTIONS={}".format(
    HTSLIB_CONFIGURE_OPTIONS))
htslib_configure_options = None

if HTSLIB_MODE in ['shared', 'separate']:
    package_list += ['ngsfragments.htslib',
                     'ngsfragments.htslib.htslib']
    package_dirs.update({'ngsfragments.htslib':'htslib'})

    htslib_configure_options = configure_library(
        "htslib",
        HTSLIB_CONFIGURE_OPTIONS,
        ["--enable-libcurl",
         "--disable-libcurl"])

    HTSLIB_SOURCE = "builtin"
    print ("# ngsfragments: htslib configure options: {}".format(
        str(htslib_configure_options)))

    config_headers += ["htslib/config.h"]
    if htslib_configure_options is None:
        # create empty config.h file
        with open("htslib/config.h", "w") as outf:
            outf.write(
                "/* empty config.h created by ngsfragments */\n")
            outf.write(
                "/* conservative compilation options */\n")

    with changedir("htslib"):
        htslib_make_options = run_make_print_config()

    for key, value in htslib_make_options.items():
        print ("# ngsfragments: htslib_config {}={}".format(key, value))

    external_htslib_libraries = ['z']
    if "LIBS" in htslib_make_options:
        external_htslib_libraries.extend(
            [re.sub("^-l", "", x) for x in htslib_make_options["LIBS"].split(" ") if x.strip()])

    shared_htslib_sources = [re.sub("\.o", ".c", os.path.join("htslib", x))
                             for x in
                             htslib_make_options["LIBHTS_OBJS"].split(" ")]

    htslib_sources = []

if HTSLIB_LIBRARY_DIR:
    # linking against a shared, externally installed htslib version, no
    # sources required for htslib
    htslib_sources = []
    shared_htslib_sources = []
    chtslib_sources = []
    htslib_library_dirs = [HTSLIB_LIBRARY_DIR]
    htslib_include_dirs = [HTSLIB_INCLUDE_DIR]
    external_htslib_libraries = ['z', 'hts']
elif HTSLIB_MODE == 'separate':
    # add to each ngsfragments component a separately compiled
    # htslib
    htslib_sources = shared_htslib_sources
    shared_htslib_sources = htslib_sources
    htslib_library_dirs = []
    htslib_include_dirs = ['htslib']
elif HTSLIB_MODE == 'shared':
    # link each ngsfragments component against the same
    # htslib built from sources included in the ngsfragments
    # package.
    htslib_library_dirs = [
        "ngsfragments",  # when using setup.py develop?
        ".",  # when using setup.py develop?
        os.path.join("build", distutils_dir_name("lib"), "ngsfragments")]

    htslib_include_dirs = ['htslib']
else:
    raise ValueError("unknown HTSLIB value '%s'" % HTSLIB_MODE)

# build config.py
with open(os.path.join("ngsfragments", "config.py"), "w") as outf:
    outf.write('HTSLIB = "{}"\n'.format(HTSLIB_SOURCE))
    config_values = collections.defaultdict(int)

    if HTSLIB_SOURCE == "builtin":
        with open(os.path.join("htslib", "config.h")) as inf:
            for line in inf:
                if line.startswith("#define"):
                    key, value = re.match(
                        "#define (\S+)\s+(\S+)", line).groups()
                    config_values[key] = value
            for key in ["ENABLE_GCS",
                        "ENABLE_PLUGINS",
                        "ENABLE_S3",
                        "HAVE_COMMONCRYPTO",
                        "HAVE_HMAC",
                        "HAVE_LIBBZ2",
                        "HAVE_LIBCURL",
                        "HAVE_LIBDEFLATE",
                        "HAVE_LIBLZMA",
                        "HAVE_MMAP"]:
                outf.write("{} = {}\n".format(key, config_values[key]))
                print ("# ngsfragments: config_option: {}={}".format(key, config_values[key]))

# create empty config.h files if they have not been created automatically
# or created by the user:
for fn in config_headers:
    if not os.path.exists(fn):
        with open(fn, "w") as outf:
            outf.write(
                "/* empty config.h created by ngsfragments */\n")
            outf.write(
                "/* conservative compilation options */\n")

#######################################################
# Windows compatibility - untested
if platform.system() == 'Windows':
    include_os = ['win32']
    os_c_files = ['win32/getopt.c']
    extra_compile_args = []
else:
    include_os = []
    os_c_files = []
    # for python 3.4, see for example
    # http://stackoverflow.com/questions/25587039/
    # error-compiling-rpy2-on-python3-4-due-to-werror-
    # declaration-after-statement
    extra_compile_args = [
        "-Wno-unused",
        "-Wno-strict-prototypes",
        "-Wno-sign-compare",
        "-Wno-error=declaration-after-statement"]

define_macros = []

suffix = sysconfig.get_config_var('EXT_SUFFIX')
if not suffix:
    suffix = sysconfig.get_config_var('SO')


libraries_for_ngs_module = external_htslib_libraries

# The list below uses the union of include_dirs and library_dirs for
# reasons of simplicity.

modules = [
    dict(name="ngsfragments.read_sam.ReadSam",
         sources=["ngsfragments.read_sam.ReadSam".replace(".", os.path.sep)+".pyx"] + shared_htslib_sources + os_c_files,
         libraries=external_htslib_libraries),
    dict(name="ngsfragments.peak_calling.RunningMedian",
         sources=["ngsfragments.peak_calling.RunningMedian".replace(".", os.path.sep)+".pyx"] + htslib_sources + os_c_files,
         libraries=external_htslib_libraries),
    dict(name="ngsfragments.peak_calling.RunningMean",
         sources=["ngsfragments.peak_calling.RunningMean".replace(".", os.path.sep)+".pyx"] + htslib_sources + os_c_files,
         libraries=external_htslib_libraries),
    dict(name="ngsfragments.peak_calling.CallPeaks",
         sources=["ngsfragments.peak_calling.CallPeaks".replace(".", os.path.sep)+".pyx"] + htslib_sources + os_c_files,
         libraries=external_htslib_libraries),
    dict(name="ngsfragments.segment.smooth_cnv.smooth_cnv",
         sources=["ngsfragments.segment.smooth_cnv.smooth_cnv".replace(".", os.path.sep)+".pyx"] + htslib_sources + os_c_files,
         libraries=external_htslib_libraries),
]

common_options = dict(
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros,
    # for out-of-tree compilation, use absolute paths
    library_dirs=[os.path.abspath(x) for x in ["ngsfragments"] + htslib_library_dirs],
    include_dirs=[os.path.abspath(x) for x in htslib_include_dirs + \
                  ["ngsfragments", ".", np.get_include(), ailist.get_include()] + include_os])

# add common options (in python >3.5, could use n = {**a, **b}
for module in modules:
    module.update(**common_options)

metadata = dict(
    name = "ngsfragments",
    version = version,
    author = "Kyle S. Smith",
    author_email = "kyle.smith@stjude.org",
    url = "https://github.com/kylessmtih/ngsfragments",
    description = SHORTDESC,
    long_description = long_description,
    long_description_content_type = "text/markdown",
    # CHANGE THIS
    license = "GPL2",
    # free-form text field
    platforms = ["POSIX", "UNIX", "MacOS"],
    classifiers = [ "Development Status :: 4 - Beta",
                    "Environment :: Console",
                    "Intended Audience :: Developers",
                    "Intended Audience :: Science/Research",
                    "Operating System :: POSIX :: Linux",
                    "Programming Language :: Cython",
                    "Programming Language :: Python",
                    "Programming Language :: Python :: 3",
                    "Programming Language :: Python :: 3.4",
                    "Topic :: Scientific/Engineering",
                    "Topic :: Scientific/Engineering :: Mathematics",
                    "Topic :: Software Development :: Libraries",
                    "Topic :: Software Development :: Libraries :: Python Modules"
                    ],
    setup_requires = ["cython", "numpy","ailist","cbseg","bcpseg", "pysam", "matplotlib","seaborn","bokeh","scipy"],
    install_requires = ["numpy","ailist","cbseg","bcpseg","intervalframe", "statsmodels", "h5py","pysam"],
    provides = ["ngsfragments"],
    keywords = ["next generation sequencing fragment"],
    ext_modules = [Extension(**opts) for opts in modules],
    cmdclass = {'build_ext': build_ext, 'clean_ext': clean_ext},
    packages = package_list,
    package_dir = package_dirs,
    package_data={'ngsfragments': ['*.pxd', '*.pyx', '*.c', '*.h'],
                  '': ['*.pxd', '*.h']},
    # Disable zip_safe
    zip_safe = False,
    # Custom data files not inside a Python package
    data_files = datafiles,
)

if __name__ == '__main__':
    dist = setup(**metadata)







































