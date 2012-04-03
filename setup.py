#!/usr/bin/env python
"""Riemann: Clawpack-style Riemann solvers

Riemann is a central repository for all Riemann solvers compatible with the Clawpack family of computational tools

"""

# much of the functionality of this file was taken from the scipy setup.py script.

DOCLINES = __doc__.split("\n")

import os
import sys
import warnings
import subprocess
import shutil
import re

if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: C
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS

"""

MAJOR               = 0
MINOR               = 1
MICRO               = 0
ISRELEASED          = False
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

def write_version_py(filename='src/python/riemann/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM RIEMANN SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    # Adding the git rev number needs to be done inside
    # write_version_py(), otherwise the import of scipy.version messes
    # up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('src/python/riemann/version.py'):
        # must be a source distribution, use existing version file
        from riemann.version import git_revision as GIT_REVISION
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev-' + GIT_REVISION[:7]

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version' : FULLVERSION,
                       'git_revision' : GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()

 
one_d_riemann =   ['acoustics',
                   'advection',
                   'burgers',
                   'euler_with_efix',
                   'nonlinear_elasticity_fwave',
                   'reactive_euler_with_efix',
                   'shallow_roe_with_efix']

two_d_riemann =   ['acoustics',
                   'advection',
                   'euler_5wave',
                   'psystem',
                   'shallow_roe_with_efix',
                   'shallow_sphere',
                   'vc_acoustics',
                   'vc_advection']

three_d_riemann = ['vc_acoustics']
                   
# special rules for rp2_kpp, rp2_euler_mapgrid, and rp3acv

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    package_path=os.path.join(os.path.dirname(__file__),'src','python','riemann')
    config = Configuration('riemann', parent_package, top_path, package_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.get_version('src/python/riemann/version.py')

    src_dir = os.path.join(os.path.dirname(__file__),'src')

    for rp in one_d_riemann:
        rp_ext = 'rp1_'+rp
        rp_src = [os.path.join(src_dir,rp_ext+'.f90')]
        config.add_extension(rp_ext,rp_src)

    for rp in two_d_riemann:
        rp_ext = 'rp2_'+rp
        rp_src = [os.path.join(src_dir,prefix+rp+'.f90')
                  for prefix in ['rpn2_','rpt2_']]
        config.add_extension(rp_ext,rp_src)

    for rp in three_d_riemann:
        rp_ext = 'rp3_'+rp
        rp_src = [os.path.join(src_dir,prefix+rp+'.f90')
                  for prefix in ['rpn3_','rpt3_','rptt3_']]
        config.add_extension(rp_ext,rp_src)

    # special targets
    special_target_list = \
    [{'ext' :'rp2_kpp',
      'srcs':['rpn2_kpp.f90','rpt2_dummy.f90']},
     {'ext' :'rp2_euler_mapgrid',
      'srcs':['rpn2_euler_mapgrid.f90','rpt2_euler_mapgrid.f90',
              'euler_roe_solver_mapgrid.f90','getquadinfo_mapgrid.f90']},
     {'ext' :'rp3acv',
      'srcs':['rpn3acv.f90','rpt3acv.f90','rptt3acv.f90']}]
    
    for rp_dict in special_target_list:
        rp_ext = rp_dict['ext']
        rp_src = [os.path.join(src_dir,src) for src in rp_dict['srcs']]
        config.add_extension(rp_ext,rp_src)
    return config


def setup_package():
    from numpy.distutils.core import setup

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    src_path = local_path

    os.chdir(local_path)
    sys.path.insert(0, local_path)
    sys.path.insert(0, os.path.join(local_path, 'src/python/riemann'))  # to retrieve version

    # Run build
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    # Rewrite the version file everytime
    write_version_py()

    try:
        setup(
            name = 'riemann',
            maintainer = "Clawpack Developers",
            maintainer_email = "claw-dev@googlegroups.com",
            description = DOCLINES[0],
            long_description = "\n".join(DOCLINES[2:]),
            url = "http://www.clawpack.org",
            download_url = "https://github.com/clawpack/riemann/tarball/master", 
            license = 'BSD',
            classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
            platforms = ["Linux", "Solaris", "Mac OS-X", "Unix"],
            configuration=configuration )
    finally:
        del sys.path[0]
        os.chdir(old_path)
    return


if __name__ == '__main__':
    setup_package()
