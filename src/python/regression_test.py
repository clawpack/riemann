import os, re, sys


#-----------------------------------------------------------
#|\""""""""""""""""""""""""|\
#| >    approximateDiff    | >
#|/________________________|/


def approximateDiff(file_name_1, file_name_2, tolerance):
    """
    Takes the "approximate diff" of two files, allowing for
    differences in real data up to the specified tolerance.

    This isn't as sophisticated as a true diff in terms of
    matching largest common subsections.  It's just meant
    to speed the user's analysis if changing the code results
    in structural or numerical differences.
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
        else:
            #-----------------------------------------------------------
            # If lines have a different number of elements, then print
            # their contents.
            #-----------------------------------------------------------
            diff_output.append("  Line " + str(i+1) + ": number of elements differs")
            diff_output.append("  " + file_name_1.rjust(40) + ": " + lines_1[i])
            diff_output.append("  " + file_name_2.rjust(40) + ": " + lines_2[i])

    return diff_output                    
# /|""""""""""""""""""""""""/|
#< |    approximateDiff    < |
# \|________________________\|


    
def compare_directories(outdir, regression_dir, tolerance, results_file):
    
    regdir = regression_dir
    output_object = open(results_file, 'w')

    #==== Locate data files in the regression_data directory ====
    try:
        regression_files = os.listdir(regression_dir)
    except:
        print "*** Problem listing regression_dir: ",regression_dir
        raise 
    data_pattern = re.compile("fort.q\d{4}")
    data_files = []
    
    for file_name in regression_files:
        if data_pattern.match(file_name):
            data_files.append(file_name)

    if data_files == []:
        print "*** No fort.q files found in regression_dir: ",regression_dir


            
    #==== Compare each data file to the same file in _output ====
    difference_detected = False
    failure = False
    for data_file in data_files:
        output_object.write("==> %s: " % data_file)
        try:
            diff_results = approximateDiff(os.path.join(outdir,data_file), \
                           os.path.join(regdir,data_file), tolerance)
            if len(diff_results) > 0:
                output_object.write("\n")
                for line in diff_results:
                    output_object.write(line + "\n")
                    difference_detected = True
            else:
                output_object.write("Identical to tolerance %s\n" % tolerance)
        except:
            #---- If the diff fails, indicate the file with which it crashed ----
            output_object.write("Test failed.\n")
            output_object.write("Could not use approximateDiff\n")
            failure = True
            break
    
    #==== Write an OK message if there were no problems ====
    if not(difference_detected) and not(failure):
        #output_object.write("Files checked:")
        #for file_name in regression_files:
            #output_object.write(" " + file_name)  
        output_object.write("\nNo change detected.\n")
        print "All fort.q files are identical to tolerance ", tolerance
    
    output_object.close()
    

def fetch_regression_data(outdir):
    clawdir = os.environ['CLAW'] + '/'
    thisdir = os.getcwd()
    thisdir = thisdir.replace(clawdir,'')
    print "+++ ",thisdir
    if thisdir[0] == '/':
        raise Exception("This directory is not a subdirectory of clawdir = %s" \
                 % clawdir)

    tarfile = thisdir.replace('/','-') + '-' + outdir + '.tar.gz'
    regdir = outdir + '-regression_data'

    remote_regdir = "clawpack@homer.u.washington.edu:public_html/regression_data/"

    print "Trying to scp %s from %s" % (tarfile, remote_regdir)

    try:
        os.system("scp %s/%s ." % (remote_regdir, tarfile))
    except:
        print "*** scp command failed!"
    
    try:
        os.system("tar -zxf %s" % tarfile)
        print "Regression data should be in ",regdir
    except:
        raise Exception("Problem untarring %s" % tarfile)

if __name__=="__main__":

    tolerance = 1e-6
    outdir = '_output'
    results_file = "regression_results.txt"
    if len(sys.argv)>1:
        regression_dir = sys.argv[1]
    else:
        try:
            fetch_regression_data(outdir)
        except:
            raise Exception("Error fetching regression data")

    compare_directories(outdir, regression_dir, tolerance, results_file)

    print "Regression results are in ",results_file
