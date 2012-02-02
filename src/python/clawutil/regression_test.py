"""
Module to do regression test of output directory against a version
that has been archived.  

Can run at command line from an application directory, e.g.:
  $ python $CLAW/python/pyclaw/regression_test.py output OUTDIR1 OUTDIR2
to compare the data in OUTDIR1 to the regression data in OUTDIR2.

If OUTDIR2 is not specified, will attempt to download comparison data from
an archive on the Clawpack website: 
  http://www.clawpack.org/regression_data/TARFILE
where TARFILE is a tarfile specific to this example (based on path within
apps directory).

To add new data to repository use the save_regression_data.py module.

Need to improve and document better!


"""


import os, re, sys
import numpy as np

def approximateDiff(file_name_1, file_name_2, tolerance):
    """
    Takes the "approximate diff" of two files, allowing for
    differences in real data up to the specified tolerance.

    This isn't as sophisticated as a true diff in terms of
    matching largest common subsections.  It's just meant
    to speed the user's analysis if changing the code results
    in structural or numerical differences.

    Originally written by Jonathan Claridge.
    """
  
    def approximatelyEqual(x,y,tolerance):
        if x==y:
            return True

        try:
            if abs(float(x) - float(y)) < tolerance:
                return True
        except:
            return False


    file_1 = open(file_name_1, 'r')
    file_2 = open(file_name_2, 'r')
    lines_1 = file_1.readlines()
    lines_2 = file_2.readlines()
    difference = False
    any_difference = False
    max_difference = 0.
    max_i = 0
    
    diff_output = []
  
    #==== Check file lengths ====
    #--------------------------------------------------------------------
    # A problem here will indicate that something is structurally wrong.
    #--------------------------------------------------------------------
    if not(len(lines_1) == len(lines_2)):
        diff_output.append("Files are of different length")


    #==== Check line by line ====
    #----------------------------------------------------------------------
    # This is where numerical differences will be highlighted.  Also, if
    # the files are comparable up to a certain point, this will show where
    # they begin to diverge.
    #----------------------------------------------------------------------
    for i in range(min(len(lines_1), len(lines_2))):
        split_1 = lines_1[i].split();
        split_2 = lines_2[i].split();
  
        if len(split_1) == len(split_2):
            #-----------------------------------------------------------
            # If lines have the same number of elements, then check for
            # numerical differences.
            #-----------------------------------------------------------
            for j in range(len(split_1)):
                if not(approximatelyEqual(split_1[j],split_2[j],tolerance)):
                    diff_output.append("  Line " +  str(i+1) + ", element " \
                        + str(j+1) + " differs")
                    diff_output.append("  " + file_name_1.rjust(40) + ": " + split_1[j])
                    diff_output.append("  " + file_name_2.rjust(40) + ": " + split_2[j])
                if not(split_1[j] == split_2[j]):
                    try:
                        x1 = float(split_1[j])
                        x2 = float(split_2[j])
                        max_difference = max(abs(x1-x2), max_difference)
                        max_i = i+1
                    except:
                        max_difference = np.nan
        else:
            #-----------------------------------------------------------
            # If lines have a different number of elements, then print
            # their contents.
            #-----------------------------------------------------------
            diff_output.append("  Line " + str(i+1) + ": number of elements differs")
            diff_output.append("  " + file_name_1.rjust(40) + ": " + lines_1[i])
            diff_output.append("  " + file_name_2.rjust(40) + ": " + lines_2[i])

    return diff_output, max_difference, max_i

    
    
def cmp_dir(outdir, regdir):

    """
    Compare files in two directories, typically:
        outdir = most recent output directory
        regdir = archived regression directory
    Report results both to the terminal and to a file diffs.html in outdir.
    For files that differ, create html files highlighting differences.
    """

    from clawutil import chardiff

    
    files = os.listdir(regdir)
    print "Comparing files in the output directory: ",outdir
    print "          with the regression directory: ", regdir
    
    hfile = open(outdir+'/diffs.html','w')
    hfile.write("""<html>
            <h1>Regression tests</h1>
            Comparing files in the output directory: &nbsp; %s<br>
            with the regression directory: &nbsp; %s<p>
            """ % (outdir,regdir))
                    

        

def compare_plots(plotdir, regression_plotdir, comparison_plotdir="_plots_diff")
                  
    """
    Generate plots in plotdir using the setplot module and then
    create html files showing side-by-side comparison with those in regression_plotdir.
    If the latter directory does not exist, attempt to fetch from
    archived regression data.
    """

    from clawutil import imagediff

    if not os.path.isdir(plotdir):
        print "*** Directory not found: ",plotdir
        return

    if not os.path.isdir(regression_plotdir):
        try:
            regression_plotdir, plotdir_tarfile = fetch_regression_data(plotdir)
        except:
            print "*** Directory not found: ",regression_plotdir
            print "*** Attempt to fetch regression data also failed."
            return

    imagediff.imagediff_dir(plotdir, regression_plotdir, comparison_plotdir)


def fetch_regression_data(outdir):
    import urllib
    
    clawdir = os.environ['CLAW'] + '/'
    thisdir = os.getcwd()
    thisdir = thisdir.replace(clawdir,'')
    if thisdir[0] == '/':
        raise Exception("This directory is not a subdirectory of clawdir = %s" \
                 % clawdir)

    tarfile = thisdir.replace('/','-') + '-' + outdir + '.tar.gz'
    regdir = outdir + '-regression_data'

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

    if sys.argv[1]=="output":
        
        try:
            outdir = sys.argv[2]
        except:
            raise Exception("*** Must specify output directory, and optionally regression directory")
        results_file = "regression_results.txt"
        retrieve = len(sys.argv) == 2
        if not retrieve:
            regression_dir = sys.argv[3]
        else:
            try:
                regression_dir, tarfile = fetch_regression_data(outdir)
            except:
                raise Exception("Error fetching regression data")

        cmp_dir(outdir, regression_dir)
        
    elif sys.argv[1]=="plots":
        try:
            plotdir = sys.argv[2]
        except:
            raise Exception("*** Must specify plots directory, and optionally regression directory")
        comparison_plotdir="_plots_comparison"
        retrieve = len(sys.argv) == 2
        if not retrieve:
            regression_dir = sys.argv[3]
        else:
            try:
                regression_dir, tarfile = fetch_regression_data(plotdir)
            except:
                raise Exception("Error fetching regression data")
        
        outdir = ""  # Assume plots already exist
        compare_plots(outdir, plotdir, regression_dir, comparison_plotdir="_plots_comparison", \
                          setplot="setplot.py", printfigs=False)

    else:
        print "*** Usage:  "
        print "    python regression_test.py output <outdir>"
        print "or  python regression_test.py plots <plotdir>"
        retrieve = False
        
    
    if retrieve:
        print "\nDownloaded a tar file and created directory ", regression_dir
        ans = raw_input("  Cleanup by removing these? ")
        if ans.lower() in ['y','yes']:
            os.system("rm -rf %s %s" % (regression_dir, tarfile))
        

