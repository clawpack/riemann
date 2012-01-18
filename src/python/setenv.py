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
    csh_handle.write('setenv %s "%s"\n' % (var.upper(),value))
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
def write_env_files(claw_path,verbose=True,outfile_base="setenv",**kargs):
    
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
            for project in missing_projects:
                error_msg += "\n  %s" % project
            
    # =========================================================================
    #  Write out out_file_base.csh and out_file_base.sh
    
    # Open output files
    csh_file = open(os.path.join(claw_path,".".join((outfile_base,"csh"))),'w')
    bash_file = open(os.path.join(claw_path,".".join((outfile_base,"bash"))),'w')
    
    # Write out boiler plate
    boiler_plate = ("# Clawpack environment settings\n")
    csh_file.write(boiler_plate)
    bash_file.write(boiler_plate)
    
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
        raise NotImplementedError("Environment settings not implemented for clawapps!")

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
        # python_path
        # = ":".join((os.path.join(available_projects["clawpack-4.x"],"python"),python_path))
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
    csh_file.close()
    bash_file.close()

if __name__ == "__main__":    
    # Parse input arguments
    argv = sys.argv
    project_paths = {}
    try:
        try:
            long_options = ["help","output=","verbose",
                 "claw="]
            for proj_name in git_repos:
                long_options.append("%s=" % proj_name)
            opts, args = getopt.getopt(argv[1:], "ho:vc:",long_options)
        except getopt.error, msg:
            raise Usage(msg)
            
        # Default script parameter values
        verbose = False
        out_file_base = "setenv"
        
        # Default claw path
        claw_path = os.path.abspath(os.curdir)
    
        # option processing
        for option, value in opts:
            # Script parameters
            if option in ("-v","--verbose"):
                 verbose = True
            if option in ("-o","--output"):
                out_file_base = value
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
                outfile_base=out_file_base,**project_paths))
                
