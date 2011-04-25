"""
clawutil/src/git-clone.py

Run this script in a directory named `clawpack` to use git to clone
each project in the Clawpack organization of GitHub.

Usage:
 $ python git-clone.py                  # will attempt to get all projects
 $ python git-clone.py classic amrclaw  # only gets those specified.

It first checks to make sure a directory of each name does not already exist.
"""

import os,sys

def git_clones(repos_list=[]):
    if repos_list == []:
        repos_list = """clawutil doc classic amrclaw geoclaw pyclaw 
                        riemann visclaw
                        petclaw sharpclaw clawpack.github.com""".split()

    print "This will attempt to clone copies of the following Clawpack projects:"
    print repos_list
    ans = raw_input("Ok ? ")
    if ans.lower() not in ['y','yes']:
        print "*** Aborting"
        return

    for repos in repos_list:
        print " "
        if os.path.isdir(repos) or os.path.isfile(repos):
            print "*** %s already exists, not cloning! ***" % repos
        else:
            git_cmd = "git clone git@github.com:clawpack/%s.git" % repos
            print "Executing: %s" % git_cmd
            try:
                os.system(git_cmd)
                print "Succesfully cloned %s" % repos
            except:
                print "*** Error cloning %s" % repos

if __name__=="__main__":
    git_clones(sys.argv[1:])


