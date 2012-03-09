#!/usr/bin/env python
# encoding: utf-8
r'''
Generate ClawPack environment variables for bash and csh

Checks to see which Clawpack projects/repositories have been downloaded
by either checking in the CLAW directory (which defaults to the current
directory) or by overridden paths set in the command line arguments.
'''

# ============================================================================
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import sys
import os
import getopt

# ============================================================================
# Parameters for the projects and git repositories

# Available git repositories to be considered
git_repos = ["classic","amrclaw","geoclaw","pyclaw","clawapps","clawutil",
             "clawpack-4.x","riemann","doc","visclaw","sharpclaw"]

# Project dependencies
project_dependencies = {"classic":["riemann","clawutil"],
                        "amrclaw":["riemann","clawutil"],
                        "geoclaw":["riemann","clawutil","amrclaw"],
                        "pyclaw":["riemann","clawutil"],
                        "clawapps":["riemann","clawutil"],
                        "clawutil":[],
                        "clawpack-4.x":[],
                        "riemann":["clawutil"],
                        "doc":[],
                        "visclaw":["clawutil"],
                        "sharpclaw":["riemann","clawutil"]}

# This was intended to be a prefered means for adding and deleting paths for
# the resulting output script but we need a csh version before this will really
# work correctly without adding a lot of complexity to this utility (KTM)
bash_path_modification_functions = """
# These are utility functions for manipulating paths
var_append () {
    # Check to see if the variable exists
    if [ -z "${1}" ]; then
        export ${1}="${2}"
    else
        var_remove $1 $2
        export ${1}="`/usr/bin/printenv $1`:${2}"
    fi
}
var_prepend () {
    # Check to see if variable exists
    if [ -z "${1}" ]; then
        export ${1}="${2}"
    else
        var_remove $1 $2
        export ${1}="${2}:`/usr/bin/printenv $1`"
    fi
}
var_remove () {
    VAR_CONTENTS=`/usr/bin/printenv $1`
    NEW_VAR=`echo -n $VAR_CONTENTS | awk -v RS=: -v ORS=: '$0 != "'$2'"' | sed 's/:$//'`
    export ${1}=${NEW_VAR}
} 

path_append () { var_append PATH $1; }
path_prepend () { var_prepend PATH $1; }
path_remove () { var_remove PATH $1; }
python_append () { var_append PYTHONPATH $1;}
python_prepend () { var_prepend PYTHONPATH $1;}
python_remove () { var_remove PYTHONPATH $1;}  

"""

# ============================================================================
#  Help display
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

repo_options = ''
for repo in git_repos:
    repo_options += "  --%s= (string) Path to %s project.  (default == $CLAW/%s)\n" % (repo,repo,repo)
help_message = '''
Generate ClawPack environment variables for bash and csh

Checks to see which Clawpack projects/repositories have been downloaded
by either checking in the CLAW directory (which defaults to the current
directory) or by overridden paths set in the command line arguments.  The 
currently accepted project/repository names are

  %s

Command line script parameters:
  -v, --verbose - Verbose output (default == True)
  -h, --help - Display help message
  -o, --output= (string) - The base name for the output bash and csh files 
                           (default == "setenv")
  -s, --shell= (string) - Type of shell script to output, valid options include
                         'csh', 'bash', 'sh', or 'both'.  The option 'both'
                         will output both a "csh" and "sh" compatable file.
                         (default == 'both')
                           
Project path options:
  -c, --claw= (string) - Path to base CLAW directory.  If this option is choosen
                         and a path provided then all project paths will be set
                         to as this as it's base, regardless of what has been 
                         set before this flag.  If you want to have a specific
                         path to the CLAW directory, give this flag first and 
                         then specify specific paths to the projects.  
                         (default == cwd)
%s 
''' % ('\n  '.join(git_repos),repo_options)

# ============================================================================
#  Helper functions
def write_environment_variable(csh_handle,bash_handle,var,value):
    if csh_handle is not None:
        csh_handle.write('setenv %s "%s"\n' % (var.upper(),value))
    if bash_handle is not None:
        bash_handle.write('export %s="%s"\n' % (var.upper(),value))

def check_repos_dependencies(project_name,available_projects):
    r"""Checks that required repositories of project_name are present"""
    missing_projects = []
    for dependency in project_dependencies[project_name]:
        if dependency not in available_projects:
            missing_projects.append(dependency)
    if len(missing_projects) > 0:
        return missing_projects
    return None

# ============================================================================
def write_env_files(claw_path,verbose=True,outfile_base="setenv",
                                                shell_type='both',**kargs):
    # Find projects
    available_projects = {}
    print "Found the following Clawpack projects:"
    for repo in git_repos:
        if repo in kargs.keys():
            if os.path.exists(kargs[repo]):
                available_projects[repo] = kargs[repo]
                if verbose:
                    print "  %s -- %s" % (repo,available_projects[repo])
                else:
                    print "  %s" % (repo)
            else:
                raise Exception("Project %s not found at custom path %s" % (repo,kargs[repo]))
        else:
            if os.path.exists(os.path.join(claw_path,repo)):
                available_projects[repo] = os.path.join(claw_path,repo)
                if verbose:
                    print "  %s -- %s" % (repo,available_projects[repo])
                else:
                    print "  %s" % (repo)
            else:
                if verbose:
                    print "  %s -- *** Not found ***" % repo
    
    # Check here for environment dependencies
    for project in available_projects.keys():
        missing_projects = check_repos_dependencies(project,available_projects.keys())
        if missing_projects is not None:
            error_msg = "The project %s depends on the following missing projects:" % project_name
            error_msg += ("\n %s" % name for name in missing_projects)
            
    # =========================================================================
    #  Write out out_file_base.csh and out_file_base.sh
    # Open output files
    boiler_plate = ("# Clawpack environment settings\n")
    if "csh" in shell_type:
        csh_file = open(os.path.join(claw_path,".".join((outfile_base,"csh"))),'w')
        csh_file.write(boiler_plate)
    else:
        csh_file = None
    if "bash" == shell_type or "sh" == shell_type:
        bash_file = open(os.path.join(claw_path,".".join((outfile_base,"bash"))),'w')
        bash_file.write(boiler_plate)
    else:
        bash_file = None
    
    # Write out variables
    python_path = "${PYTHONPATH}"
    matlab_path = "${MATLABPATH}"
    
    print ""
    print "The following variables will be set:"
    print "  CLAW = %s" % claw_path
    write_environment_variable(csh_file,bash_file,"CLAW",claw_path)
    
    if "clawutil" in available_projects:
        print "  CLAWUTIL = %s" % available_projects["clawutil"]
        write_environment_variable(csh_file,bash_file,"CLAWUTIL",available_projects["clawutil"])
        python_path = ":".join((os.path.join(available_projects["clawutil"],"src","python"),python_path))
        
    if "classic" in available_projects:
        print "  CLASSICCLAW = %s" % available_projects["classic"]
        write_environment_variable(csh_file,bash_file,"CLASSICCLAW",available_projects["classic"])

    if "amrclaw" in available_projects:
        print "  AMRCLAW = %s" % available_projects["amrclaw"]
        write_environment_variable(csh_file,bash_file,"AMRCLAW",available_projects["amrclaw"])

    if "geoclaw" in available_projects:
        python_path = ":".join((os.path.join(available_projects["geoclaw"],"src","python"),python_path))
        print "  GEOCLAW = %s" % available_projects["geoclaw"]
        write_environment_variable(csh_file,bash_file,"GEOCLAW",available_projects["geoclaw"])

    if "pyclaw" in available_projects:
        python_path = ":".join((os.path.join(available_projects["pyclaw"],"src"),python_path))
        print "  PYCLAW = %s" % available_projects["pyclaw"]
        write_environment_variable(csh_file,bash_file,"PYCLAW",available_projects["pyclaw"])

    if "clawapps" in available_projects:
        print "  CLAWAPPS = %s" % available_projects["clawapps"]
        write_environment_variable(csh_file,bash_file,"CLAWAPPS",available_projects["clawapps"])
    if "doc" in available_projects:
        pass
        
    if "riemann" in available_projects:
        python_path = ":".join((os.path.join(available_projects["riemann"],"src","python"),python_path))
        print "  RIEMANN = %s" % available_projects["riemann"]
        write_environment_variable(csh_file,bash_file,"RIEMANN",available_projects["riemann"])

    if "visclaw" in available_projects:
        python_path = ":".join((os.path.join(available_projects["visclaw"],"src","python"),python_path))
        matlab_path = ":".join((os.path.join(available_projects["visclaw"],"src","matlab"),matlab_path))
        print "  VISCLAW = %s" % available_projects["visclaw"]
        write_environment_variable(csh_file,bash_file,"VISCLAW",available_projects["visclaw"])

    if "shaprclaw" in available_projects:
        python_path = ":".join((os.path.join(available_projects["sharpclaw"],"src","python"),python_path))
        print "  SHARPCLAW = %s" % available_projects["sharpclaw"]
        write_environment_variable(csh_file,bash_file,"SHARPCLAW",available_projects["sharpclaw"])

    if "clawpack-4.x" in available_projects:
        # python_path = ":".join((os.path.join(available_projects["clawpack-4.x"],"python"),python_path))
        print "  CLAW_4 = %s" % available_projects["clawpack-4.x"]
        write_environment_variable(csh_file,bash_file,"CLAW_4",available_projects["clawpack-4.x"])

    if len(python_path) > 13:
        print "  PYTHONPATH = %s" % python_path
        write_environment_variable(csh_file,bash_file,"PYTHONPATH",python_path)
    if len(matlab_path) > 13:
        print "  MATLABPATH = %s" % matlab_path
        write_environment_variable(csh_file,bash_file,"MATLABPATH",matlab_path)
    print ""
        
    # Close output files
    if csh_file is not None:
        csh_file.close()
    if bash_file is not None:
        bash_file.close()

if __name__ == "__main__":    
    # Parse input arguments
    argv = sys.argv
    project_paths = {}
    try:
        try:
            long_options = ["help","output=","verbose","shell=",
                 "claw="]
            for proj_name in git_repos:
                long_options.append("%s=" % proj_name)
            opts, args = getopt.getopt(argv[1:], "ho:vs:c",long_options)
        except getopt.error, msg:
            raise Usage(msg)
            
        # Default script parameter values
        verbose = False
        out_file_base = "setenv"
        shell_type = 'both'
        
        # Default claw path
        claw_path = os.path.abspath(os.curdir)
    
        # option processing
        for option, value in opts:
            # Script parameters
            if option in ("-v","--verbose"):
                verbose = True
            if option in ("-o","--output"):
                out_file_base = value
            if option in ("-s","--shell"):
                shell_type = value
            if option in ("-h","--help"):
                raise Usage(help_message)
                                
            # Project path overrides
            if option in ("-c","--claw"):
                claw_path = os.path.abspath(value)
            for proj_name in git_repos:
                if option == "--%s" % proj_name:
                    project_paths[proj_name] = os.path.abspath(value)
                            
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        sys.exit(2)

    sys.exit(write_env_files(claw_path,verbose=verbose,
                             outfile_base=out_file_base,shell_type=shell_type,
                             **project_paths))
                
