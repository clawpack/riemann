#!/usr/bin/env python
# encoding: utf-8
r"""
Data Module

Contains the general class definition and the subclasses of the Clawpack data
objects.

Rewritten for 5.0:
 Stripped down version of Data object, no longer keeps track of "owners".
 data classes are now subclasses of ClawData, which checks if attributes
   already exist before setting.


:Authors:
    Kyle T. Mandli and Randall J. LeVeque 

"""
# ============================================================================
#  Distributed under the terms of the Berkeley Software Distribution (BSD)
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import os
import logging
        


# ========================
class Data(object):
# ========================
    """
    Unrestricted object, can assign to new attributes.
    """
    
    def __init__(self, attributes=None):
        
        # Initialize from attribute list provided
        if attributes:
            for attr in attributes:
                self.__setattr__(attr,None)


# ========================
class ClawData(object):
# ========================
    r"""
    Class to be subclassed when defining data objects that should have
    a limited set of allowed attributes.  Useful to guard against
    typos or misrembering the names of expected attributes.

    Resetting values of existing attributes is allowed as usual,
    but new attributes can only be added using the method add_attribute.

    Trying to set a nonexistent attribute will raise an AttributeError
    exception, except for those starting with '_'.   
    """


    def __init__(self, attributes=None):
        
        # Attribute to store a list of the allowed attributes, 
        # appended to when add_attribute is used: 
        object.__setattr__(self,'_attributes',[])

        # Initialize from attribute list provided
        if attributes:
            for attr in attributes:
                self.add_attribute(attr,None)

    def __setattr__(self,name,value):
        r"""
        Check that attribute exists before setting it.
        If not, raise an AttributeError.
        Exception: attributes starting with '_' are ok to set.
        """

        if (not name in self._attributes) and (name[0] != '_'):
            print "*** Unrecognized attribute: ",name
            print "*** Perhaps a typo?"
            print "*** Add new attributes using add_attribute method"
            raise AttributeError("Unrecognized attribute: %s" % name)
        
        # attribute exists, ok to set:
        object.__setattr__(self,name,value)

    def add_attribute(self, name, value=None, add_to_list=True):
        r"""
        Adds an attribute called name to the data object

        If an attribute needs to be added to the object, this routine must be
        called or the attribute will not be written out.

        :Input:
         - *name* - (string) Name of the data attribute
         - *value* - (id) Value to set *name* to, defaults to None
        """
        if (name not in self._attributes) and add_to_list:
            self._attributes.append(name)
        object.__setattr__(self,name,value)

    def add_attributes(self, arg_list, value=None):
        r"""
        Add a list of attributes, each initialized to *value*.
        """
        for name in arg_list:
            self.add_attribute(name, value)

    def remove_attributes(self, arg_list):
        r"""
        Remove the listed attributes.
        """

        # Convert to list if args is not already a list
        if not isinstance(arg_list,list):
            arg_list = [arg_list]

        for arg in arg_list:
            self._attributes.remove(arg)
            delattr(self,arg)

    #def attributes():
    #    def fget(self): return self._attributes
    #    return locals()
    #attributes = property(**attributes())

    def has_attribute(self,name):
        r"""
        Check if this data object has the given attributes

        :Input:
         - *name* - (string) Name of attribute

        :Output:
         - (bool) - True if data object contains a data attribute name
        """
        return name in self._attributes
        
    def iteritems(self):
        r"""
        Returns an iterator of attributes and values from this object

        :Output:
         - (Iterator) Iterator over attributes and values
        """
        return [(k,getattr(self,k)) for k in self._attributes]



# ========== Parse Value Utility Function ====================================

def _parse_value(value):
    r"""
    Attempt to make sense of a value string from a config file.  If the
    value is not obviously an integer, float, or boolean, it is returned as
    a string stripped of leading and trailing whitespace.

    :Input:
        - *value* - (string) Value string to be parsed

    :Output:
        - (id) - Appropriate object based on *value*
    """
    value = value.strip()
    if not value:
        return None

    # assume that values containing spaces are lists of values
    if len(value.split()) > 1:
        return [_parse_value(vv) for vv in value.split()]

    try:
        # see if it's an integer
        value = int(value)
    except ValueError:
        try:
            # see if it's a float
            value = float(value)
        except ValueError:
            # see if it's a bool
            if value[0] == 'T':
                value = True
            elif value[0] == 'F':
                value = False

    return value


#-----------------------------------------------------------------------

# New classes and functions for dealing with data in setrun function.

class ClawInputData(ClawData):
    r"""
    Object that will be written out to claw.data.
    """
    def __init__(self, ndim):
        super(ClawInputData,self).__init__()
        self.add_attribute('ndim',ndim)

        # Set default values:
        self.add_attribute('mx',100)
        self.add_attribute('meqn',1)
        self.add_attribute('mwaves',1)
        self.add_attribute('maux',0)
        self.add_attribute('output_style',1)
        self.add_attribute('output_ntimes',1)
        self.add_attribute('output_times',[1.])
        self.add_attribute('output_step_interval',10)
        self.add_attribute('total_steps',5)     
        self.add_attribute('output_time_interval',1.)
        self.add_attribute('output_format',1)
        self.add_attribute('output_q_components','all')
        self.add_attribute('output_aux_components',[])
        self.add_attribute('output_aux_onlyonce',True)
        self.add_attribute('tfinal',1.0)
        
        self.add_attribute('dt_initial',1.e-5)
        self.add_attribute('dt_max',1.e99)
        self.add_attribute('dt_variable',1)
        self.add_attribute('cfl_desired',0.9)
        self.add_attribute('cfl_max',1.0)
        self.add_attribute('steps_max',50000)
        self.add_attribute('order',2)
        self.add_attribute('order_trans',0)
        self.add_attribute('dimensional_split',0)
        self.add_attribute('verbosity',0)
        self.add_attribute('verbosity_regrid',0)
        self.add_attribute('src_split',0)
        self.add_attribute('mcapa',0)
        self.add_attribute('limiter',[4])
        self.add_attribute('t0',0.)
        self.add_attribute('xlower',0.)
        self.add_attribute('xupper',1.)
        self.add_attribute('num_ghost',2)
        self.add_attribute('bc_xlower',1)
        self.add_attribute('bc_xupper',1)
        self.add_attribute('restart',0)
        self.add_attribute('restart_frame',0)
        self.add_attribute('fwave',False)
        self.add_attribute('restart',False)
        self.add_attribute('restart_file','')
        self.add_attribute('regions',[])
        self.add_attribute('gauges',[])


        if ndim >= 2:
            self.add_attribute('my',100)
            self.add_attribute('ylower',0.)
            self.add_attribute('yupper',1.)
            self.add_attribute('bc_ylower',1)
            self.add_attribute('bc_yupper',1)       

        if ndim == 3:
            self.add_attribute('mz',100)
            self.add_attribute('zlower',0.)
            self.add_attribute('zupper',1.)
            self.add_attribute('bc_zlower',1)
            self.add_attribute('bc_zupper',1)       

        if ndim not in [1,2,3]:
            raise ValueError("Only ndim=1, 2, or 3 supported ")

    def write(self):
        print 'Creating data file claw.data for use with xclaw'
        make_clawdatafile(self)


class AmrclawInputData(ClawInputData):
    r"""
    Object that will be written out to amrclaw.data.
    """
    def __init__(self, ndim):

        # Set default values:
        
        # Some defaults are inherited from ClawInputData:
        super(AmrclawInputData,self).__init__(ndim)
        
        self.add_attribute('amrlevels_max',1)
        self.add_attribute('refinement_ratio_x',[1])
        self.add_attribute('refinement_ratio_t',[1])
        self.add_attribute('auxtype',[])

        self.add_attribute('checkpt_style',1)
        self.add_attribute('checkpt_interval',1000)
        self.add_attribute('checkpt_time_interval',1000.)
        self.add_attribute('checkpt_times',[1000.])
        self.add_attribute('checkpt_ntimes',1)
        
        self.add_attribute('flag_richardson',False)
        self.add_attribute('flag_richardson_tol',1.0)
        self.add_attribute('flag_gradient',True)
        self.add_attribute('flag_gradient_tol',0.05)
        self.add_attribute('regrid_interval',2)
        self.add_attribute('regrid_buffer_width',3)
        self.add_attribute('clustering_cutoff',0.7)

        
        # debugging flags:
        self.add_attribute('dprint',False)
        self.add_attribute('eprint',False)
        self.add_attribute('edebug',False)
        self.add_attribute('gprint',False)
        self.add_attribute('nprint',False)
        self.add_attribute('pprint',False)
        self.add_attribute('rprint',False)
        self.add_attribute('sprint',False)
        self.add_attribute('tprint',False)
        self.add_attribute('uprint',False)

        
        
        if ndim == 1:
            
            # attributes needed only because 1d AMR is done using 2d amrclaw:
            self.add_attribute('my',1)
            self.add_attribute('ylower',0.)
            self.add_attribute('yupper',1.)
            self.add_attribute('bc_ylower',1)
            self.add_attribute('bc_yupper',1)
            self.add_attribute('refinement_ratio_y',[1,1,1,1,1,1])

        elif ndim >= 2:

            self.add_attribute('refinement_ratio_y',[1])


        if ndim not in [1,2]:
            print '*** Error: only ndim=1 or 2 supported so far ***'
            raise AttributeError("Only ndim=1 or 2 supported so far")

    def write(self):
        print 'Creating data file amrclaw.data for use with xamr'
        make_amrclawdatafile(self)
        #make_setgauges_datafile(self)



def open_datafile(name, datasource='setrun.py'):
    """
    Open a data file and write a warning header.
    Warning header starts with '#' character.  These lines are skipped if
    data file is opened using the library routine opendatafile.

    :Input:
     - *name* - (string) Name of data file
     - *datasource* - (string) Source for the data

    :Output:
     - (file) - file object
    """

    import string

    source = string.ljust(datasource,25)
    file = open(name, 'w')
    file.write('########################################################\n')
    file.write('### DO NOT EDIT THIS FILE:  GENERATED AUTOMATICALLY ####\n')
    file.write('### To modify data, edit  %s ####\n' % source)
    file.write('###    and then "make .data"                        ####\n')
    file.write('########################################################\n\n')

    return file


def data_write(file, dataobj, name=None, descr=''):
    r"""
    Write out value to data file, in the form ::

       value =: name  descr

    Remove brackets and commas from lists, and replace booleans by T/F.

    :Input:
     - *name* - (string) normally a string defining the variable,
       ``if name==None``, write a blank line.
     - *descr* - (string) A short description to appear on the line
    """

    import string
    if name is None:
        file.write('\n')
    else:
        try:
            value = getattr(dataobj, name)
        except:
            print "Variable missing: ",name
            print "  from dataobj = ", dataobj
            raise
        # Convert value to an appropriate string repr
        if isinstance(value,tuple) | isinstance(value,list):
            # Remove [], (), and ','
            string_value = repr(value)[1:-1]
            string_value = string_value.replace(',','')
        elif isinstance(value,bool):
            if value:
                string_value = 'T'
            else:
                string_value = 'F'
        else:
            string_value = repr(value)
        padded_value = string.ljust(string_value, 25)
        padded_name = string.ljust(name, 12)
        file.write('%s =: %s %s\n' % (padded_value, padded_name, descr))


def make_clawdatafile(clawdata):
    r"""
    Take the data specified in clawdata and write it to claw.data in the
    form required by the Fortran code classic/src/Nd/??.
    """


    # open file and write a warning header:
    file = open_datafile('claw.data')

    ndim = clawdata.ndim
    data_write(file, clawdata, 'ndim', '(number of dimensions)')
    data_write(file, clawdata, 'mx', '(cells in x direction)')
    if ndim > 1:
        data_write(file, clawdata, 'my', '(cells in y direction)')
    if ndim == 3:
        data_write(file, clawdata, 'mz', '(cells in z direction)')
    data_write(file, clawdata, None)  # writes blank line
    data_write(file, clawdata, 'meqn', '(number of equations)')
    data_write(file, clawdata, 'mwaves', '(number of waves)')
    data_write(file, clawdata, 'maux', '(number of aux variables)')
    data_write(file, clawdata, None)  # writes blank line

    data_write(file, clawdata, 'output_style', '(style of specifying output times)')
    if clawdata.output_style == 1:
        data_write(file, clawdata, 'output_ntimes', '(number of output times)')
        data_write(file, clawdata, 'tfinal', '(final time)')
    elif clawdata.output_style == 2:
        clawdata.output_ntimes = len(clawdata.output_times)
        data_write(file, clawdata, 'output_ntimes', '(number of output times)')
        data_write(file, clawdata, 'output_times', '(output times)')
    elif clawdata.output_style == 3:
        data_write(file, clawdata, 'output_step_interval', '(timesteps between output)')
        data_write(file, clawdata, 'total_steps', '(number of output times)')
    elif clawdata.output_style == 4:
        data_write(file, clawdata, 'output_time_interval', '(output interval)')
        data_write(file, clawdata, 'tfinal', '(final time)')
    else:
        print '*** Error: unrecognized output_style'
        raise
        return

    data_write(file, clawdata, None)
    if clawdata.output_format in [1,'ascii']:
        clawdata.output_format = 1
    elif clawdata.output_format in [2,'binary']:
        clawdata.output_format = 2
    else:
        raise ValueError("*** Error in data parameter: " + \
              "output_format unrecognized: ",clawdata.output_format)
        
    data_write(file, clawdata, 'output_format', '(output format)')
    clawdata._iout_q = clawdata.meqn * [1]
    if clawdata.output_q_components != 'all':
        for i in range(clawdata.meqn):
            if i+1 not in clawdata.output_q_components:
                clawdata._iout_q[i] = 0
    data_write(file, clawdata, '_iout_q', '(which components of q)')
    if clawdata.maux > 0:
        clawdata._iout_aux = clawdata.maux * [1]
        if clawdata.output_aux_components != 'all':
            for i in range(clawdata.maux):
                if i+1 not in clawdata.output_aux_components:
                    clawdata._iout_aux[i] = 0
        data_write(file, clawdata, '_iout_aux', '(which components of aux)')
        data_write(file, clawdata, 'output_aux_onlyonce', '(only once?)')
    data_write(file, clawdata, None)
    data_write(file, clawdata, 'dt_initial', '(initial time step dt)')
    data_write(file, clawdata, 'dt_max', '(max allowable dt)')
    data_write(file, clawdata, 'cfl_max', '(max allowable Courant number)')
    data_write(file, clawdata, 'cfl_desired', '(desired Courant number)')
    data_write(file, clawdata, 'steps_max', '(max time steps per call to claw)')
    data_write(file, clawdata, None)
    data_write(file, clawdata, 'dt_variable', '(1 for variable dt, 0 for fixed)')
    data_write(file, clawdata, 'order', '(1 or 2)')
    if ndim == 1:
        #data_write(file, clawdata, 'order_trans', '(not used in 1d)')
        pass
    else:
        data_write(file, clawdata, 'order_trans', '(transverse order)')
        data_write(file, clawdata, 'dimensional_split', '(use dimensional splitting?)')
        
    data_write(file, clawdata, 'verbosity', '(verbosity of output)')
    data_write(file, clawdata, 'src_split', '(source term splitting)')
    data_write(file, clawdata, 'mcapa', '(aux index for capacity fcn)')
    data_write(file, clawdata, None)

    data_write(file, clawdata, 'limiter', '(limiter choice for each wave)')
    data_write(file, clawdata, None)

    data_write(file, clawdata, 't0', '(initial time)')
    data_write(file, clawdata, 'xlower', '(xlower)')
    data_write(file, clawdata, 'xupper', '(xupper)')
    if ndim > 1:
        data_write(file, clawdata, 'ylower', '(ylower)')
        data_write(file, clawdata, 'yupper', '(yupper)')
    if ndim == 3:
        data_write(file, clawdata, 'zlower', '(zlower)')
        data_write(file, clawdata, 'zupper', '(zupper)')
    data_write(file, clawdata, None)

    data_write(file, clawdata, 'num_ghost', '(number of ghost cells)')
    
    for bdry in ['xlower','xupper','ylower','yupper','zlower','zupper']:
        bc = getattr(clawdata, 'bc_'+bdry, None) 
        if bc == 'user':
            setattr(clawdata, 'bc_'+bdry, 0)
        if bc == 'extrap':
            setattr(clawdata, 'bc_'+bdry, 1)
        if bc == 'periodic':
            setattr(clawdata, 'bc_'+bdry, 2)
        if bc == 'wall':
            setattr(clawdata, 'bc_'+bdry, 3)
            
    data_write(file, clawdata, 'bc_xlower', '(type of BC at xlower)')
    data_write(file, clawdata, 'bc_xupper', '(type of BC at xupper)')
    if ndim > 1:
        data_write(file, clawdata, 'bc_ylower', '(type of BC at ylower)')
        data_write(file, clawdata, 'bc_yupper', '(type of BC at yupper)')
    if ndim == 3:
        data_write(file, clawdata, 'bc_zlower', '(type of BC at zlower)')
        data_write(file, clawdata, 'bc_zupper', '(type of BC at zupper)')

    data_write(file, clawdata, None)
    data_write(file, clawdata, 'restart', '(T to restart from a past run)')
    data_write(file, clawdata, 'restart_frame', '(which frame to restart from)')
    data_write(file, clawdata, None)

    file.close()

def make_amrclawdatafile(clawdata):
    r"""
    Take the data specified in clawdata and write it to claw.data in the
    form required by the Fortran code lib/main.f95.
    """


    # open file and write a warning header:
    file = open_datafile('amrclaw.data')

    ndim = clawdata.ndim
    data_write(file, clawdata, 'ndim', '(number of dimensions)')
    data_write(file, clawdata, 'mx', '(cells in x direction)')
    data_write(file, clawdata, 'my', '(cells in y direction)')
    if ndim == 3:
        data_write(file, clawdata, 'mz', '(cells in z direction)')

    data_write(file, clawdata, 'amrlevels_max', '(max number of grid levels)')
    if len(clawdata.refinement_ratio_x) < max(abs(clawdata.amrlevels_max)-1, 1):
        raise ValueError("*** Error in data parameter: " + \
              "require len(refinement_ratio_x) >= %s " % max(abs(clawdata.amrlevels_max) - 1, 1))
    if len(clawdata.refinement_ratio_y) < max(abs(clawdata.amrlevels_max)-1, 1):
        raise ValueError("*** Error in data parameter: " + \
              "require len(refinement_ratio_y) >= %s " % max(abs(clawdata.amrlevels_max) - 1, 1))
    data_write(file, clawdata, 'refinement_ratio_x', '(refinement ratios)')
    data_write(file, clawdata, 'refinement_ratio_y', '(refinement ratios)')
    if ndim == 3:
        if len(clawdata.refinement_ratio_z) < max(abs(clawdata.amrlevels_max)-1, 1):
                raise ValueError("*** Error in data parameter: " + \
                  "require len(refinement_ratio_z) >= %s " % max(abs(clawdata.amrlevels_max) - 1, 1))
        data_write(file, clawdata, 'refinement_ratio_z', '(refinement ratios)')
    if len(clawdata.refinement_ratio_t) < max(abs(clawdata.amrlevels_max)-1, 1):
        raise ValueError("*** Error in data parameter: " + \
              "require len(refinement_ratio_t) >= %s " % max(abs(clawdata.amrlevels_max) - 1, 1))
    data_write(file, clawdata, 'refinement_ratio_t', '(refinement ratios)')

    data_write(file, clawdata, None)  # writes blank line
    data_write(file, clawdata, 't0', '(initial time)')

    data_write(file, clawdata, 'output_style', '(style of specifying output times)')
    if clawdata.output_style == 1:
        data_write(file, clawdata, 'output_ntimes', '(number of output times)')
        data_write(file, clawdata, 'tfinal', '(final time)')
    elif clawdata.output_style == 2:
        clawdata.output_ntimes = len(clawdata.output_times)
        data_write(file, clawdata, 'output_ntimes', '(number of output times)')
        data_write(file, clawdata, 'output_times', '(output times)')
    elif clawdata.output_style == 3:
        data_write(file, clawdata, 'output_step_interval', '(output every output_step_interval steps)')
        data_write(file, clawdata, 'total_steps', '(number of steps to take)')
    elif clawdata.output_style == 4:
        data_write(file, clawdata, 'output_time_interval', '(time between outputs)')
        data_write(file, clawdata, 'tfinal', '(final time)')
    else:
        print '*** Error: unrecognized output_style = ',output_style
        raise
        return

    data_write(file, clawdata, None)
    data_write(file, clawdata, 'dt_initial', '(initial time step dt)')
    data_write(file, clawdata, 'dt_max', '(max allowable dt)')
    data_write(file, clawdata, 'cfl_max', '(max allowable Courant number)')
    data_write(file, clawdata, 'cfl_desired', '(desired Courant number)')
    data_write(file, clawdata, 'steps_max', '(max time steps per call to claw)')
    data_write(file, clawdata, None)
    data_write(file, clawdata, 'dt_variable', '(1 for variable dt, 0 for fixed)')
    data_write(file, clawdata, 'order', '(1 or 2)')
    if ndim == 1:
        data_write(file, clawdata, 'order_trans', '(not used in 1d)')
    else:
        data_write(file, clawdata, 'order_trans', '(transverse order)')
    data_write(file, clawdata, 'dimensional_split', '(use dimensional splitting?)')
        
    data_write(file, clawdata, 'verbosity', '(verbosity of output)')
    data_write(file, clawdata, 'src_split', '(source term splitting)')
    data_write(file, clawdata, 'mcapa', '(aux index for capacity fcn)')
    data_write(file, clawdata, 'maux', '(number of aux variables)')
    if len(clawdata.auxtype) != clawdata.maux:
        file.close()
        print "*** Error: An auxtype array must be specified of length maux"
        raise AttributeError, "require len(clawdata.auxtype) == clawdata.maux"
    for i in range(clawdata.maux):
        file.write("'%s'\n" % clawdata.auxtype[i])
    data_write(file, clawdata, None)

    data_write(file, clawdata, 'meqn', '(number of equations)')
    data_write(file, clawdata, 'mwaves', '(number of waves)')
    data_write(file, clawdata, 'limiter', '(limiter choice for each wave)')
    if clawdata.fwave:
        clawdata.add_attribute('ifwave',1)
    else:
        clawdata.add_attribute('ifwave',0)
    data_write(file, clawdata, 'ifwave', '(use f-wave form?)')
    data_write(file, clawdata, None)

    data_write(file, clawdata, 'xlower', '(xlower)')
    data_write(file, clawdata, 'xupper', '(xupper)')
    data_write(file, clawdata, 'ylower', '(ylower)')
    data_write(file, clawdata, 'yupper', '(yupper)')
    if ndim == 3:
        data_write(file, clawdata, 'zlower', '(zlower)')
        data_write(file, clawdata, 'zupper', '(zupper)')
    data_write(file, clawdata, None)

    for bdry in ['xlower','xupper','ylower','yupper','zlower','zupper']:
        bc = getattr(clawdata, 'bc_'+bdry, None) 
        if bc == 'user':
            setattr(clawdata, 'bc_'+bdry, 0)
        if bc == 'extrap':
            setattr(clawdata, 'bc_'+bdry, 1)
        if bc == 'periodic':
            setattr(clawdata, 'bc_'+bdry, 2)
        if bc == 'wall':
            setattr(clawdata, 'bc_'+bdry, 3)

    data_write(file, clawdata, 'num_ghost', '(number of ghost cells)')
    data_write(file, clawdata, 'bc_xlower', '(type of BC at xlower)')
    data_write(file, clawdata, 'bc_xupper', '(type of BC at xupper)')
    data_write(file, clawdata, 'bc_ylower', '(type of BC at ylower)')
    data_write(file, clawdata, 'bc_yupper', '(type of BC at yupper)')
    if ndim == 3:
        data_write(file, clawdata, 'bc_zlower', '(type of BC at zlower)')
        data_write(file, clawdata, 'bc_zupper', '(type of BC at zupper)')
    data_write(file, clawdata, None)

    data_write(file, clawdata, 'restart', '(1 to restart from a past run)')
    data_write(file, clawdata, 'checkpt_style', '(how checkpoints specified)')
    
    if clawdata.checkpt_style == 2:  
        clawdata.checkpt_ntimes = len(clawdata.checkpt_times)
        data_write(file, clawdata, 'checkpt_ntimes', '(number of checkpoint times)')
        data_write(file, clawdata, 'checkpt_times', '(checkpoint times)')
    elif clawdata.checkpt_style == 3:        
        data_write(file, clawdata, 'checkpt_interval', '(step interval for checkpoint)')
    elif clawdata.checkpt_style == 4:  
        print "*** Error: checkpt_style==4 not yet implemented"
        raise ValueError("Invalid checkpt_style")      
        data_write(file, clawdata, 'checkpt_time_interval', '(time interval for checkpoint)')
    else:
        print "*** Error, unrecognized checkpt_style = ",clawdata.checkpt_style
        print "*** Require: 0,2,3, or 4"
        raise ValueError("Unrecognized checkpt_style")
        return
        
    data_write(file, clawdata, None)

    data_write(file, clawdata, 'flag_richardson', '(use Richardson extrap?)')
    data_write(file, clawdata, 'flag_richardson_tol', '(tolerance for Richardson)')
    data_write(file, clawdata, 'flag_gradient', '(use gradient flagging?)')
    data_write(file, clawdata, 'flag_gradient_tol', '(tolerance used for gradient)')
    data_write(file, clawdata, 'regrid_interval', '(how often to regrid)')
    data_write(file, clawdata, 'regrid_buffer_width', '(buffer zone around flagged pts)')
    data_write(file, clawdata, 'clustering_cutoff', '(efficiency cutoff for clustering)')
    data_write(file, clawdata, 'verbosity_regrid', '(what levels to print grid info)')
    data_write(file, clawdata, None)

    if clawdata.output_format == 'ascii':
        clawdata.output_format = 1
    elif clawdata.output_format == 'binary':
        clawdata.output_format = 2
    else:
        if clawdata.output_format not in [1,2]:
            print "*** Unrecognized output_format: ",clawdata.output_format
            raise
            return
    data_write(file, clawdata, 'output_format', '(format for fort.q files)')
    
    if clawdata.output_q_components=='all':
        clawdata.output_q_components = range(1,clawdata.meqn+1)

    nq_components = len(clawdata.output_q_components)
    clawdata.add_attribute('output_nq_components',nq_components)
    data_write(file, clawdata, 'output_nq_components', '(number of q vals to print)')
    if clawdata.output_nq_components > 0:
        data_write(file, clawdata, 'output_q_components', '(which components of q)')
    
    if clawdata.output_aux_components=='all':
        clawdata.output_aux_components = range(1,clawdata.maux+1)
        
    naux_components = len(clawdata.output_aux_components)
    clawdata.add_attribute('output_naux_components',naux_components)
    data_write(file, clawdata, 'output_naux_components', '(number of aux vals to print)')
    if clawdata.output_naux_components > 0:
        data_write(file, clawdata, 'output_aux_components', '(which components of aux)')
    data_write(file, clawdata, 'output_aux_onlyonce', '(only at t0?)')
        
             
        
    data_write(file, clawdata, None)

    data_write(file, clawdata, 'dprint', '(print domain flags)')
    data_write(file, clawdata, 'eprint', '(print err est flags)')
    data_write(file, clawdata, 'edebug', '(even more err est flags)')
    data_write(file, clawdata, 'gprint', '(grid bisection/clustering)')
    data_write(file, clawdata, 'nprint', '(proper nesting output)')
    data_write(file, clawdata, 'pprint', '(proj. of tagged points)')
    data_write(file, clawdata, 'rprint', '(print regridding summary)')
    data_write(file, clawdata, 'sprint', '(space/memory output)')
    data_write(file, clawdata, 'tprint', '(time step reporting each level)')
    data_write(file, clawdata, 'uprint', '(update/upbnd reporting)')
    data_write(file, clawdata, None)

    clawdata.add_attribute('nregions', len(clawdata.regions))
    data_write(file, clawdata, 'nregions', '(nregions)')
    for regions in clawdata.regions:
        file.write(8*"   %g" % tuple(regions) +"\n")

    clawdata.add_attribute('ngauges', len(clawdata.gauges))
    data_write(file, clawdata, 'ngauges', '(ngauges)')
    gaugeno_used = []
    for gauge in clawdata.gauges:
        gaugeno = gauge[0]
        if gaugeno in gaugeno_used:
            print "*** Gauge number %s used more than once! " % gaugeno
            raise Exception("Repeated gauge number")
        else:
            gaugeno_used.append(gauge[0])
        #file.write("%4i %19.10e  %17.10e  %13.6e  %13.6e\n" % tuple(gauge))
        file.write(5*"   %g" % tuple(gauge) +"\n")

    file.close()

def make_userdatafile(userdata):
    r"""
    Create the data file using the parameters in userdata.
    The parameters will be written to this file in the same order they were
    specified using userdata.add_attribute.
    Presumably the user will read these in using a Fortran routine, such as
    setprob.f95, and the order is important.
    """

    # open file and write a warning header:
    file = open_datafile(userdata.__fname__)

    # write all the parameters:
    for param in userdata._attributes:
        data_write(file, userdata, param, \
                   userdata.__descr__[param])

    file.close()

def make_setgauges_datafile(clawdata):
    """
    Create setgauges.data using gauges attribute of clawdata.
    """
    gauges = getattr(clawdata,'gauges',[])
    ngauges = len(gauges)

    print 'Creating data file setgauges.data'
    # open file and write a warning header:
    file = open_datafile('setgauges.data')
    file.write("%4i   =: ngauges\n" % ngauges)
    gaugeno_used = []
    for gauge in gauges:
        gaugeno = gauge[0]
        if gaugeno in gaugeno_used:
            print "*** Gauge number %s used more than once! " % gaugeno
            raise Exception("Repeated gauge number")
        else:
            gaugeno_used.append(gauge[0])
        file.write("%4i %19.10e  %17.10e  %13.6e  %13.6e\n" % tuple(gauge))
        # or use this variant with =:
        #gauge.append(gaugeno)
        #file.write("%4i %19.10e  %17.10e  %13.6e  %13.6e  =: gauge%s\n" % tuple(gauge))
    file.close()


#-----------------------------------------------------
# New version 6/30/09

class ClawRunData(ClawData):
    r"""
    Object that will be written out to claw.data.
    """
    def __init__(self, pkg, ndim):
        super(ClawRunData,self).__init__()
        self.add_attribute('pkg',pkg)
        self.add_attribute('ndim',ndim)
        self.add_attribute('datalist',[])


        if pkg.lower() in ['classic', 'classicclaw']:
            self.add_attribute('xclawcmd', 'xclaw')

            # Required data set for basic run parameters:
            clawdata = ClawInputData(ndim)
            self.add_attribute('clawdata', clawdata)
            self.datalist.append(clawdata)

        elif pkg.lower() in ['amrclaw', 'amr']:
            self.add_attribute('xclawcmd', 'xamr')

            # Required data set for basic run parameters:
            clawdata = AmrclawInputData(ndim)
            self.add_attribute('clawdata', clawdata)
            self.datalist.append(clawdata)

        elif pkg.lower() in ['geoclaw']:
            self.add_attribute('xclawcmd', 'xgeoclaw')

            # Required data set for basic run parameters:
            clawdata = AmrclawInputData(ndim)
            self.add_attribute('clawdata', clawdata)
            self.datalist.append(clawdata)
            geodata = GeoclawInputData(ndim)
            self.add_attribute('geodata', geodata)
            self.datalist.append(geodata)

        else:
            raise AttributeError("Unrecognized Clawpack pkg = %s" % pkg)

    def new_UserData(self,name,fname):
        r"""
        Create a new attribute called name
        for application specific data to be written
        to the data file fname.
        """
        userdata = UserData(fname)
        self.datalist.append(userdata)
        self.add_attribute(name,userdata)
        exec('self.%s = userdata' % name)
        return userdata

    def add_GaugeData(self):
        r"""
        Create a gaugedata attribute for writing to gauges.data.
        """
        gaugedata = GaugeData(self.ndim)
        self.datalist.append(gaugedata)
        self.gaugedata = gaugedata
        return gaugedata

    def write(self):
        for d in self.datalist:
            d.write()

class UserData(ClawData):
    r"""
    Object that will be written out to user file such as setprob.data, as
    determined by the fname attribute.
    """

    def __init__(self, fname):

        super(UserData,self).__init__()

        # Create attributes without adding to attributes list:

        # file to be read by Fortran for this data:
        object.__setattr__(self,'__fname__',fname)

        # dictionary to hold descriptions:
        object.__setattr__(self,'__descr__',{})

    def add_param(self,name,value,descr=''):
         self.add_attribute(name,value)
         descr_dict = self.__descr__
         descr_dict[name] = descr

    def write(self):
         print 'Creating data file %s' % self.__fname__
         make_userdatafile(self)


class GeoclawInputData(ClawData):
    r"""
    Object that will be written out to the various GeoClaw data files.
    """
    def __init__(self, ndim):
        super(GeoclawInputData,self).__init__()

        # Set default values:
        self.add_attribute('igravity',1)
        self.add_attribute('iqinit',0)
        self.add_attribute('icoriolis',1)
        self.add_attribute('Rearth',6367500.0)
        # NEED TO CONTINUE!

    def write(self):

        print 'Creating data file setgeo.data'
        # open file and write a warning header:
        file = open_datafile('setgeo.data')
        data_write(file, self, 'igravity')
        data_write(file, self, 'gravity')
        data_write(file, self, 'icoordsys')
        data_write(file, self, 'icoriolis')
        data_write(file, self, 'Rearth')
        file.close()

        print 'Creating data file settsunami.data'
        # open file and write a warning header:
        file = open_datafile('settsunami.data')
        data_write(file, self, 'sealevel')
        data_write(file, self, 'drytolerance')
        data_write(file, self, 'wavetolerance')
        data_write(file, self, 'depthdeep')
        data_write(file, self, 'maxleveldeep')
        data_write(file, self, 'ifriction')
        data_write(file, self, 'coeffmanning')
        data_write(file, self, 'frictiondepth')
        file.close()

        print 'Creating data file settopo.data'
        # open file and write a warning header:
        file = open_datafile('settopo.data')
        self.ntopofiles = len(self.topofiles)
        data_write(file, self, 'ntopofiles')
        for tfile in self.topofiles:
            try:
                fname = os.path.abspath(tfile[-1])
            except:
                print "*** Error: file not found: ",tfile[-1]
                raise MissingFile("file not found")
            file.write("\n'%s' \n " % fname)
            file.write("%3i %3i %3i %20.10e %20.10e \n" % tuple(tfile[:-1]))
        file.close()

        print 'Creating data file setdtopo.data'
        # open file and write a warning header:
        file = open_datafile('setdtopo.data')
        self.mdtopofiles = len(self.dtopofiles)
        data_write(file, self, 'mdtopofiles')
        data_write(file, self, None)
        for tfile in self.dtopofiles:
            try:
                fname = "'%s'" % os.path.abspath(tfile[-1])
            except:
                print "*** Error: file not found: ",tfile[-1]
                raise MissingFile("file not found")
            file.write("\n%s \n" % fname)
            file.write("%3i %3i %3i\n" % tuple(tfile[:-1]))
        file.close()

        print 'Creating data file setqinit.data'
        # open file and write a warning header:
        file = open_datafile('setqinit.data')
        # self.iqinit tells which component of q is perturbed!
        data_write(file, self, 'iqinit')
        data_write(file, self, None)
        for tfile in self.qinitfiles:
            try:
                fname = "'%s'" % os.path.abspath(tfile[-1])
            except:
                print "*** Error: file not found: ",tfile[-1]
                raise MissingFile("file not found")
            file.write("\n%s  \n" % fname)
            file.write("%3i %3i \n" % tuple(tfile[:-1]))
        file.close()

        make_setgauges_datafile(self)


        print 'Creating data file setfixedgrids.data'
        # open file and write a warning header:
        file = open_datafile('setfixedgrids.data')
        self.nfixedgrids = len(self.fixedgrids)
        data_write(file, self, 'nfixedgrids')
        data_write(file, self, None)
        for fixedgrid in self.fixedgrids:
            file.write(11*"%g  " % tuple(fixedgrid) +"\n")
        file.close()


        print 'Creating data file setregions.data'
        # open file and write a warning header:
        file = open_datafile('setregions.data')
        self.nregions = len(self.regions)
        data_write(file, self, 'nregions')
        data_write(file, self, None)
        for regions in self.regions:
            file.write(8*"%g  " % tuple(regions) +"\n")
        file.close()

