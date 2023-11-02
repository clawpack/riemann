"""
This script reads constants from comments in the Fortran Riemann solvers and
writes the values to miniature modules so they can be imported and used in
Python scripts.  The constants that are read and written include:

- num_eqn
- num_waves
- names/indices of q fields
- names/indices of aux fields
"""
import os

for filename in os.listdir('.'):
    # Only process normal Riemann solver files
    if filename.startswith('rp1') or filename.startswith('rpn'):
        f = open(filename,'r')
        q_indices = []; q_names = []
        aux_indices = []; aux_names = []
        q_block   = False
        aux_block = False
        num_aux = None
        num_eqn = None
        # Read variable names and other constants
        for line in f.readlines():
            if len(line.split())<2:
                q_block   = False
                aux_block = False
            elif q_block:
                q_indices.append(int(line.split()[-2])-1)
                q_names.append(line.split()[-1])
            elif aux_block:
                aux_indices.append(int(line.split()[-2])-1)
                aux_names.append(line.split()[-1])
            elif 'waves:' in line.lower():
                num_waves = int(line.split()[-1])
            elif 'equations:' in line.lower():
                num_eqn = int(line.split()[-1])
            elif 'aux fields:' in line.lower():
                num_aux = int(line.split()[-1])
            elif 'conserved quantities:' in line.lower():
                q_block = True
            elif 'auxiliary variables:' in line.lower():
                aux_block = True
        f.close()

        if num_eqn is None: # No information from comments
            continue
        # Write constants and variable names to file
        if filename.startswith('rp1'):
            out_name = filename.split('.')[0][4:]+'_1D_constants.py'
        elif filename.startswith('rpn2'):
            out_name = filename.split('.')[0][5:]+'_2D_constants.py'
        elif filename.startswith('rpn3'):
            out_name = filename.split('.')[0][5:]+'_3D_constants.py'
        out_name = './'+out_name
        print('writing ' + out_name)
        f = open(out_name,'w')
        f.write('num_waves = '+str(num_waves)+'\n')
        f.write('num_eqn   = '+str(num_eqn)+'\n')
        if num_aux is not None:
            f.write('num_aux   = '+str(num_aux)+'\n')
        f.write('\n')
        f.write('# Conserved quantities\n')
        for variable, index in zip(q_names,q_indices):
            f.write(variable + ' = ' + str(index) + '\n')
        f.write('\n')
        if len(aux_names)>0:
            f.write('# Auxiliary variables\n')
        for variable, index in zip(aux_names,aux_indices):
            f.write(variable + ' = ' + str(index) + '\n')
        f.close()
    else:
        pass
