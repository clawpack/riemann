"""
Fetch regression data from an archive on the Clawpack website: 
  http://www.clawpack.org/regression_data/TARFILE
where TARFILE is a tarfile specific to this example (based on path within
apps directory).

To add new data to repository use the save_regression_data.py module.

"""


import os, re, sys
import numpy as np

def fetch_regression_data(target_dir):
    import urllib
    
    clawdir = os.environ['CLAW'] + '/'
    thisdir = os.getcwd()
    thisdir = thisdir.replace(clawdir,'')
    if thisdir[0] == '/':
        raise Exception("This directory is not a subdirectory of clawdir = %s" \
                 % clawdir)

    tarfile = thisdir.replace('/','-') + '-' + target_dir + '.tar.gz'
    regdir = target_dir + '-regression_data'

    url = "http://www.clawpack.org/regression_data"

    print "Trying to retrieve %s \n   from %s" % (tarfile, url)
    
    try:     
        url = os.path.join(url, tarfile)
        urllib.urlretrieve(url, tarfile)
        if os.path.getsize(tarfile) < 500:
            os.system("mv %s tarfile_error.html" % tarfile)
            print "\n*** Error: See tarfile_error.html"
            raise Exception("*** Problem retrieving %s" % tarfile)
    except:
        raise Exception("*** Problem retrieving %s" % tarfile)
    
    try:
        os.system("tar -zxf %s" % tarfile)
        print "Regression data should be in ",regdir
    except:
        raise Exception("*** Problem untarring %s" % tarfile)
    return regdir, tarfile


#---------------------------------------------------

if __name__=="__main__":

    target_dir = sys.argv[1]

    try:
        regression_dir, tarfile = fetch_regression_data(target_dir)
    except:
        raise Exception("Error fetching regression data")
    
    print "\nDownloaded a tar file ", tarfile
    ans = raw_input("  Cleanup by removing? ")
    if ans.lower() in ['y','yes']:
        os.system("rm -rf %s" % tarfile)
        

