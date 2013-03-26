one_d_riemann =   ['acoustics',
                   'advection',
                   'burgers',
                   'traffic',
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
                   
# special rules for rp2_kpp, rp2_euler_mapgrid

import os

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('riemann', parent_package, top_path)

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
     {'ext' :'rp2_euler_4wave',
      'srcs':['rpn2_euler_4wave.f90','rpt2_euler.f90']}]
    
    for rp_dict in special_target_list:
        rp_ext = rp_dict['ext']
        rp_src = [os.path.join(src_dir,src) for src in rp_dict['srcs']]
        config.add_extension(rp_ext,rp_src)
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
