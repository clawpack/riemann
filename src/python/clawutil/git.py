#!/usr/bin/env python

r"""Output the git repository status of the available projects

"""

import sys
import os
import subprocess
import tempfile

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

def check_git_status(path):
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


if __name__ == "__main__":

    verbose = False

    # If there are command line arguments, use these for check
    argv = sys.argv
    if len(argv) > 1:
        repo_checks = []
        for argument in argv[1:]:
            if env_variables.has_key(argument):
                repo_checks.append(argument)
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
        print check_git_status(path)
        print "========================================="
