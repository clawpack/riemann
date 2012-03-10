#!/usr/bin/env python
r"""
Clawpack GitHub utility functions module

Contains functions for cloning the repositories at www.github.com/clawpack
and checking the status of the local repositories.

:Available Functions:

    status - Checks the status of the local repositories
    clone - Clones git repositories from GitHub
"""

import sys
import os
import subprocess
import tempfile
import getopt

# These environment variables, if they exist, point to the paths
# of the repositories
env_variables = {"classic":"CLASSICCLAW",
                 "amrclaw":"AMRCLAW",
                 "geoclaw":"GEOCLAW",
                 "pyclaw":"PYCLAW",
                 "clawapps":"CLAWAPPS",
                 "clawutil":"CLAWUTIL",
                 "clawpack-4.x":"CLAW_4",
                 "riemann":"RIEMANN",
                 "doc":"CLAWDOCS",
                 "visclaw":"VISCLAW",
                 "sharpclaw":"SHARPCLAW"}

def status(path=None):
    r"""
    This function checks the status of the repository at path.
    
    """
    
    if path is None:
        path = os.getcwd()
    
    expanded_path = os.path.expandvars(os.path.expanduser(path))
    cmd = "cd %s; git status" % path
    output = tempfile.TemporaryFile("rw")

    # Get repos status
    subprocess.Popen(cmd,shell=True,stdout=output,stderr=output).wait()

    # Read back in the status that was captured from the command
    output.seek(0)
    status = output.read()
    output.close()

    return status


def clone(repo,force=False,verbose=False):
    r"""
    This function clones git repositories located at GitHub.

    It first checks to make sure a directory of each name does not already 
    exist unless force == True.

    :Input:
     - repos_list (list) - List of repositories to fetch, if empty then fetch
                           all.
     - force (bool) - Whether to overwrite existing directories if present,
                      (`default == False`).
    """

    if (os.path.isdir(repo) or os.path.isfile(repo)) and force:
        print >> sys.stderr, "*** %s already exists, not cloning! ***" % repo
    else:
        git_cmd = "git clone git@github.com:clawpack/%s.git" % repo
        if verbose:
            print "Running: %s" % git_cmd
        return subprocess.call(git_cmd,shell=True)
                
                

# Command line interface to the git functions here
# These work as git does such that calling this module directly will
# perform the desired function
class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg
    
if __name__ == "__main__":

    argv = sys.argv[1:]
    try:
        if len(argv) == 0:
            raise Usage("")
        elif len(argv) == 1:
            function = argv[0]
            arguments = []
        else:
            function = argv[0]
            arguments = argv[1:]
            
        # Parse each function's argument list
        if function == "status":
            try:
                opts,args = getopt.getopt(arguments,"hv",["help","verbose"])
            except getopt.error, msg:
                raise Usage(msg)
        elif function == "clone":
            try:
                opts,args = getopt.getopt(arguments,"hvf",["help","verbose",
                                            "force"])
            except getopt.error, msg:
                raise Usage(msg)
        else:
            raise Usage("Unsupported (or unknown) command %s." % function)
            
        # Default values
        verbose = False
        force = False
            
        # Option parsing
        for option,value in opts:
            if option in ('-v','--verbose'):
                verbose = True
            if option in ('-h','--help'):
                raise Usage(__doc__)
            if option in ('-f','--force'):
                force = True
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        sys.exit(2)
            

    if function == "status":
        if len(args) > 0:
            for arg in args:
                if env_variables.has_key(arg):
                    repo_checks.append(arg)
                else:
                    print ("*** WARNING *** %s may not be a Clawpack repos." 
                            % argument)
        else:
            repo_checks = env_variables.keys()

        # Check if the environment variables exist for each requested repos
        available_projects = {}
        for repo in repo_checks:
            if os.environ.has_key(env_variables[repo]):
                available_projects[repo] = os.environ[env_variables[repo]]
            else:
                if verbose:
                    print ("*** WARNING *** Could not find environment variable ",
                           "%s for project %s, skipping." % (env_variables[repo],repo))


        # Loop through looking for each environment and check its status
        for (name,path) in available_projects.items():
            print "Checking status of %s at path" % name
            print "   %s" % path
            print status(path)
            print "========================================="
    
    elif function == "clone":
        if len(args) == 0:
            repos_list = env_variables.keys()
        else:
            repos_list = args

        print "This will attempt to clone copies of the following Clawpack projects:"
        print repos_list
        ans = raw_input("Ok ? ")
        if ans.lower() not in ['y','yes']:
            print "*** Aborting"
            sys.exit(1)
        for repo in repos_list:
            status = clone(repo,force=force,verbose=verbose)
            if status == 0:
                if verbose:
                    print "Succesfully cloned %s" % repos
            else:
                print >> sys.stderr, "*** Error cloning %s" % repos
                
        
    
        

