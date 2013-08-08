import os
from numpy.distutils.misc_util import Configuration
from numpy.distutils import __config__ as numpy_config
from collections import defaultdict

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
                   'burgers',
                   'euler_5wave',
                   'psystem',
                   'shallow_roe_with_efix',
                   'shallow_sphere',
                   'vc_acoustics',
                   'vc_advection',
                   'vc_elasticity']

three_d_riemann = ['vc_acoustics',
                   'euler',
                   'burgers',
                   'vc_advection']

# special rules for rp2_kpp, rp2_euler_mapgrid


def configuration(parent_package='',top_path=None):
    
    config = Configuration('riemann', parent_package, top_path)

    # configuration updates for lapack/blas
    config_kwargs = defaultdict(list)

    def update_dict_of_lists(base_dict, update_dict):
        """Appends list items from update_dict into defaultdict base_dict"""

        for key, row in update_dict.iteritems():
            config_kwargs[key] += row
            

    blas_config = numpy_config.blas_opt_info
    lapack_config = numpy_config.lapack_opt_info

    update_dict_of_lists(config_kwargs, blas_config)
    update_dict_of_lists(config_kwargs, lapack_config)
    
    src_dir = os.path.join(os.path.dirname(__file__),'src')

    for rp in one_d_riemann:
        rp_ext = rp+'_1D'
        rp_src = [os.path.join(src_dir,'rp1_'+rp+'.f90')]
        config.add_extension(rp_ext,rp_src, **config_kwargs)

    for rp in two_d_riemann:
        rp_ext = rp+'_2D'
        rp_src = [os.path.join(src_dir,prefix+rp+'.f90')
                  for prefix in ['rpn2_','rpt2_']]
        config.add_extension(rp_ext,rp_src, **config_kwargs)

    for rp in three_d_riemann:
        rp_ext = rp+'_3D'
        rp_src = [os.path.join(src_dir,prefix+rp+'.f90')
                  for prefix in ['rpn3_','rpt3_','rptt3_']]
        config.add_extension(rp_ext,rp_src, **config_kwargs)

    # special targets
    special_target_list = \
    [{'ext' :'kpp_2D',
      'srcs':['rpn2_kpp.f90','rpt2_dummy.f90']},
     {'ext' :'euler_mapgrid_2D',
      'srcs':['rpn2_euler_mapgrid.f90','rpt2_euler_mapgrid.f90',
              'euler_roe_solver_mapgrid.f90','getquadinfo_mapgrid.f90']},
     {'ext' :'euler_4wave_2D',
      'srcs':['rpn2_euler_4wave.f90','rpt2_euler.f90']}]

    for rp_dict in special_target_list:
        rp_ext = rp_dict['ext']
        rp_src = [os.path.join(src_dir,src) for src in rp_dict['srcs']]
        config.add_extension(rp_ext,rp_src, **config_kwargs)
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
