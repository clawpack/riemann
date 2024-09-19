package = 'clawpack.riemann'
pkg_dir = join_paths(package.split('.'))

python_sources = [
  '__init__.py',
  'acoustics_1D_constants.py',
  'acoustics_1D_py.py',
  'acoustics_2D_constants.py',
  'acoustics_mapped_2D_constants.py',
  'acoustics_ptwise_2D_constants.py',
  'acoustics_variable_1D_constants.py',
  'advection_1D_constants.py',
  'advection_1D_py.py',
  'advection_2D_constants.py',
  'advection_color_1D_constants.py',
  'burgers_1D_constants.py',
  'burgers_1D_py.py',
  'burgers_2D_constants.py',
  'burgers_3D_constants.py',
  'cubic_1D_constants.py',
  'euler_1D_py.py',
  'euler_3D_constants.py',
  'euler_4wave_2D_constants.py',
  'euler_5wave_2D_constants.py',
  'euler_mapgrid_2D_constants.py',
  'euler_with_efix_1D_constants.py',
  'kpp_2D_constants.py',
  'layered_shallow_water_1D_constants.py',
  'mhd_1D_constants.py',
  'nonlinear_elasticity_1D_py.py',
  'nonlinear_elasticity_fwave_1D_constants.py',
  'psystem_2D_constants.py',
  'reactive_euler_with_efix_1D_constants.py',
  'riemann_interactive.py',
  'riemann_tools.py',
  'setup.py',
  'shallow_1D_py.py',
  'shallow_roe_tracer_1D_constants.py',
  'shallow_roe_with_efix_1D_constants.py',
  'shallow_roe_with_efix_2D_constants.py',
  'shallow_sphere_2D_constants.py',
  'static.py',
  'traffic_1D_constants.py',
  'traffic_vc_1D_constants.py',
  'traffic_vc_tracer_1D_constants.py',
  'vc_acoustics_2D_constants.py',
  'vc_acoustics_3D_constants.py',
  'vc_advection_1D_py.py',
  'vc_advection_2D_constants.py',
  'vc_advection_3D_constants.py',
  'vc_elasticity_2D_constants.py',
  'write_constants.py',
]

py.install_sources(
  python_sources,
  subdir: pkg_dir,
)

srcdir = '..' / 'src'
sep = '_'

riemann_1D_ptwise = [
  'acoustics',
  'advection',
]

riemann_2D_ptwise = [
  'acoustics',
]

riemann_1D = [
  'acoustics',
  'acoustics_variable',
  'advection',
  'advection_color',
  'burgers',
  'cubic',
  'traffic',
  'traffic_vc',
  'traffic_vc_fwave',
  'traffic_vc_tracer',
  'euler_with_efix',
  'euler_hlle',
  'mhd_roe',
  'nonlinear_elasticity_fwave',
  'psystem_fwave',
  'reactive_euler_with_efix',
  'shallow_hlle',
  'shallow_roe_with_efix',
  'shallow_bathymetry_fwave',
  'shallow_roe_tracer',
]

riemann_2D = [
  'acoustics',
  'acoustics_mapped',
  'advection',
  'burgers',
  'euler_5wave',
  'psystem',
  'shallow_roe_with_efix',
  'shallow_bathymetry_fwave',
  'sw_aug',
  'shallow_sphere',
  'vc_acoustics',
  'vc_advection',
  'vc_elasticity',
]

riemann_3D = [
  'vc_acoustics',
  'euler',
  'euler_mapgrid',
  'burgers',
  'vc_advection',
]

riemann_ptwise = {
  '1D': riemann_1D_ptwise,
  '2D': riemann_2D_ptwise,
}

riemann = {
  '1D': riemann_1D,
  '2D': riemann_2D,
  '3D': riemann_3D,
}

prefixes = {
  '1D': ['rp1'],
  '2D': ['rpn2', 'rpt2'],
  '3D': ['rpn3', 'rpt3', 'rptt3'],
}

foreach suffix, names: riemann_ptwise
  foreach name: names
    ext_name = sep.join(name, suffix, 'ptwise')
    ext_srcs = []
    foreach prefix: prefixes[suffix]
      ext_srcs += srcdir / sep.join(prefix, 'ptwise.f90')
      ext_srcs += srcdir / sep.join(prefix, name, 'ptwise.f90')
    endforeach
    f2py_srcs = custom_target(
      command: [f2py, ext_name],
      input: ext_srcs,
      output: [ext_name + 'module.c', ext_name + '-f2pywrappers.f'],
    )
    py.extension_module(
      ext_name, [ext_srcs, f2py_srcs],
      incdir_f2py / 'fortranobject.c',
      include_directories: inc_np,
      subdir: pkg_dir,
      install : true
    )
  endforeach
endforeach

foreach suffix, names: riemann
  foreach name: names
    ext_name = sep.join(name, suffix)
    ext_srcs = []
    foreach prefix: prefixes[suffix]
      ext_srcs += srcdir / sep.join(prefix, name + '.f90')
    endforeach
    f2py_srcs = custom_target(
      command: [f2py, ext_name],
      input: ext_srcs,
      output: [ext_name + 'module.c', ext_name + '-f2pywrappers.f'],
    )
    py.extension_module(
      ext_name, [ext_srcs, f2py_srcs],
      incdir_f2py / 'fortranobject.c',
      include_directories: inc_np,
      subdir: pkg_dir,
      install : true
    )
  endforeach
endforeach


special_extensions = {
  'kpp_2D': [
    'rpn2_kpp.f90',
    'rpt2_dummy.f90',
  ],
  'shallow_hlle_2D': [
    'rpn2_shallow_hlle.f90',
    'rpt2_dummy.f90',
  ],
  'euler_hlle_2D': [
    'rpn2_euler_hlle.f90',
    'rpt2_dummy.f90',
  ],
  'euler_hlle_with_walls_2D': [
    'rpn2_euler_hlle_with_walls.f90',
    'rpt2_dummy.f90',
  ],
  'euler_mapgrid_2D': [
    'rpn2_euler_mapgrid.f90',
    'rpt2_euler_mapgrid.f90',
    'euler_roe_solver_mapgrid.f90',
    'getquadinfo_mapgrid.f90',
  ],
  'euler_4wave_2D': [
    'rpn2_euler_4wave.f90',
    'rpt2_euler.f90',
  ],
}

foreach name, sources: special_extensions
  ext_name = name
  ext_srcs = []
  foreach src: sources
    ext_srcs += srcdir / src
  endforeach
  f2py_srcs = custom_target(
    command: [f2py, ext_name],
    input: ext_srcs,
    output: [ext_name + 'module.c', ext_name + '-f2pywrappers.f'],
  )
  py.extension_module(
    ext_name, [ext_srcs, f2py_srcs],
    incdir_f2py / 'fortranobject.c',
    include_directories: inc_np,
    subdir: pkg_dir,
    install : true
  )
endforeach

ext_name = 'mappedGrid'
ext_srcs = [
  srcdir / 'euler3d_mappedGrid.f90',
]
f2py_srcs = custom_target(
  'f2py_mappedGrid',
  command: [f2py, ext_name],
  input: ext_srcs,
  output: [ext_name + 'module.c', ext_name + '-f2pywrappers2.f90'],
  build_by_default: true,
)

py.extension_module(
  ext_name, [ext_srcs, f2py_srcs],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  subdir: pkg_dir,
  install : true
)
