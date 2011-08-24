import os, sys

def send_outdir(outdir="_output"):
    clawdir = os.environ['CLAW'] + '/'
    thisdir = os.getcwd()
    thisdir = thisdir.replace(clawdir,'')
    print "+++ ",thisdir
    if thisdir[0] == '/':
        raise Exception("This directory is not a subdirectory of clawdir = %s" \
                 % clawdir)

    tarfile = thisdir.replace('/','-') + '-' + outdir + '.tar'
    regdir = outdir + '-regression_data'
    if os.path.exists(regdir):
        raise Exception("Directory %s already exists" % regdir)

    os.system("ln -s %s  %s" % (outdir, regdir))
    os.system("tar -cHf %s %s" % (tarfile, regdir))
    os.system("gzip %s" % tarfile)
    tarfile = tarfile + '.gz'

    remote_regdir = "clawpack@homer.u.washington.edu:public_html/regression_data/"

    print "Trying to scp %s to %s" % (tarfile, remote_regdir)

    try:
        os.system("scp %s %s" % (tarfile, remote_regdir))
    except:
        print "*** scp command failed!"

if __name__=="__main__":
    args = sys.argv[1:]
    send_outdir(*args)
