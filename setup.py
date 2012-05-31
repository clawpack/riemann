one_d_riemann =   ['acoustics',
                   'advection',
                   'burgers',
                   'euler_with_efix',
                   'nonlinear_elasticity_fwave',
                   'reactive_euler_with_efix',
                   'shallow_roe_with_efix',
                   'layered_shallow_water']

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

    package_path=os.path.join(os.path.dirname(__file__),'src','python')
    config = Configuration('riemann', parent_package, top_path, package_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

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
     # Not working yet!
     # {'ext':'rp2_layered_shallow_water',
     #  'srcs':['rpn2_layered_shallow_water.f90','rpt2_layered_shallow_water.f90','geoclaw_rp.f']}]
    
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
            package_dir={'':os.path.join('src','python')},
            configuration=configuration )
    finally:
        del sys.path[0]
        os.chdir(old_path)
    return


if __name__ == '__main__':
    setup_package()
