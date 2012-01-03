
"""
Generic code for running the fortran version of Clawpack and sending the
results to subdirectory output of the directory from which this is executed.
Execute via
    $ python $CLAW/python/pyclaw/runclaw.py
from a directory that contains a claw.data file and a Clawpack executable.
"""


def runclaw(xclawcmd=None, outdir=None, overwrite=True, restart=False, 
            rundir=None):
    """
    Run the Fortran version of Clawpack using executable xclawcmd, which is
    typically set to 'xclaw', 'xamr', etc.

    If it is not set by the call, get it from the environment variable
    CLAW_EXE.  Default to 'xclaw' if that's not set.
    
    If rundir is None, all *.data is copied from current directory, if a path 
    is given, data files are copied from there instead.
    """

    import os,glob,shutil,time
    verbose = False
    xclawout = None
    xclawerr = None

    if type(overwrite) is str:
        # convert to boolean
        overwrite = (overwrite.lower() in ['true','t'])
    
    if type(restart) is str:
        # convert to boolean
        restart = (restart.lower() in ['true','t'])
    

    if xclawcmd is None:
        # Determine what executable to use from environment variable CLAW_EXE
        # Default to 'xclaw' if it's not set:
        xclawcmd = os.environ.get('CLAW_EXE', 'xclaw')
    

    if outdir is None:
        outdir = '.'
        
    if rundir is None:
        rundir = os.getcwd()
    rundir = os.path.abspath(rundir)
    print "Will take data from ", rundir

    # directory for fort.* files:
    outdir = os.path.abspath(outdir)
    print '== runclaw: Will write output to ',outdir
    
    import pdb; pdb.set_trace()
    
    
    #returncode = clawjob.runxclaw()

    if 1:
        startdir = os.getcwd()
        xdir = os.path.abspath(startdir)
        outdir = os.path.abspath(outdir)
        rundir = os.path.abspath(rundir)
        xclawcmd = os.path.join(xdir,xclawcmd)
        
        try:
            os.chdir(xdir)
        except:
            raise Exception( "==> runclaw: Cannot change to directory xdir = %s" %xdir)
            return 
    
    
        if os.path.isfile(outdir):
            print "==> runclaw: Error: outdir specified is a file"
            return
        
        if (os.path.isdir(outdir) & (not overwrite)):
            # copy the old outdir before possibly overwriting
            tm = time.localtime(os.path.getmtime(outdir))
            year = str(tm[0]).zfill(4)
            month = str(tm[1]).zfill(2)
            day = str(tm[2]).zfill(2)
            hour = str(tm[3]).zfill(2)
            minute = str(tm[4]).zfill(2)
            second = str(tm[5]).zfill(2)
            outdir_backup = outdir + '_%s%s%s-%s%s%s' \
                  % (year,month,day,hour,minute,second)
            if verbose:
                print "==> runclaw: Directory already exists: ",os.path.split(outdir)[1]
                if restart:
                    print "==> runclaw: Copying directory to:      ",os.path.split(outdir_backup)[1]
                else:
                    print "==> runclaw: Moving directory to:      ",os.path.split(outdir_backup)[1]
                time.sleep(1)
            
            try:
                shutil.move(outdir,outdir_backup)
                if restart:
                    shutil.copytree(outdir_backup,outdir)
            except:
                print "==> runclaw: Could not move directory... copy already exists?"
            
            
        if (not os.path.isdir(outdir)):
            try:
                os.mkdir(outdir)
            except:
                print "Cannot make directory ",outdir
                return
    
        try:
            os.chdir(outdir)
        except:
            print '==> runclaw: *** Error in runxclaw: cannot move to outdir = ',\
                  outdir
            raise
            return
    
        fortfiles = glob.glob(os.path.join(outdir,'fort.*'))
        if (overwrite and (not restart)):
            # remove any old versions:
            if verbose:
                print "==> runclaw: Removing all old fort files in ", outdir
            for file in fortfiles:
                os.remove(file)
        elif restart:
            if verbose:
                print "==> runclaw: Restart: leaving original fort files in ", outdir
        else:
            if len(fortfiles) > 1:
                print "==> runclaw: *** Remove fort.* and try again,"
                print "  or use overwrite=True in call to runxclaw"
                print "  e.g., by setting CLAW_OVERWRITE = True in Makefile"
                return
            
        
        try:
            os.chdir(rundir)
        except:
            raise Exception("Cannot change to directory %s" % rundir)
            return 
    
        datafiles = glob.glob('*.data')
        if datafiles == ():
            print "==> runclaw: Warning: no data files found in directory ",rundir
        else:
            if rundir != outdir:
                for file in datafiles:
                    shutil.copy(file,os.path.join(outdir,file))
    
        if xclawout:
            xclawout = open(xclawout,'wb')
        if xclawerr:
            xclawerr = open(xclawerr,'wb')
    
        os.chdir(outdir)
    
        #print "\nIn directory outdir = ",outdir,"\n"
    
        # execute command to run fortran program:
    
        try:
            #print "\nExecuting ",xclawcmd, "  ...  "
            #pclaw = subprocess.Popen(xclawcmd,stdout=xclawout,stderr=xclawerr)
            #print '+++ pclaw started'
                #pclaw.wait()   # wait for code to run
            #returncode = pclaw.returncode
            #print '+++ pclaw done'
            
            returncode = os.system(xclawcmd)
    
            if returncode == 0:
                print "\n==> runclaw: Finished executing\n"
            else:
                print "\n ==> runclaw: *** Runtime error: return code = %s\n " % returncode
        except:
            raise Exception("Could not execute command %s" % xclawcmd)
    
        os.chdir(startdir)

    if returncode != 0:
        print '==> runclaw: *** fortran returncode = ', returncode, '   aborting'
    print '==> runclaw: Done executing %s via pyclaw.runclaw.py' % xclawcmd
    print '==> runclaw: Output is in ', outdir
    

#----------------------------------------------------------

if __name__=='__main__':
    """
    If executed at command line prompt, simply call the function, with
    any argument used as setplot:
    """
    import sys
    args = sys.argv[1:]   # any command line arguments
    runclaw(*args)
