#!/usr/bin/env python
__doc__ = r"""
Module to diff two images or directories full of images.

Can by used from the command line via:
   python imagediff.py file1 file2
   python imagediff.py dir1 dir2
Assumes dir1 and dir2 are subdirectories of working directory.

Command line flags include:
    -v, --verbose = Verbose output   
    -h, --help = Display this help
"""

import os, sys, getopt

def imagediff_file(fname1, fname2, verbose):
    
    fname3 = make_imagediff(fname1,fname2,verbose=verbose)
    
    hfile = "imagediff.html"
    html = open(hfile,"w")
    html.write("""
         <html>
         <h2>Comparison of image files</h2>
         Left: %s <br>
         Right: %s
         <p>
         <img src="%s" width=400 border="1"> &nbsp; <img src="%s" width=400 border="1">""" \
      % (fname1,fname2,fname1,fname2))
      
    html.write("""
         <h2>Pixel difference is nonzero where image below is black...</h2>
         <img src="%s" width=400 border="1">""" % fname3)
    html.close()
    print "To view comparision and pixelwise difference, see ", hfile
    
    
         
def make_imagediff(fname1,fname2,fname3='', verbose=False):
    ext1 = os.path.splitext(fname1)[1]
    ext2 = os.path.splitext(fname2)[1]
    if ext1 != ext2:
        print "*** Error: extensions of fname1 and fname2 must agree"
        return
    if ext1 not in ['.png','.gif','.jpg']:
        print "*** Error: image files not recognized"
        return
    
    if fname3 == '':
        fname3 = "_image_diff" + ext1
    
    os.system("convert %s %s -compose subtract -composite -threshold 0.001 -negate %s" \
         % (fname1, fname2, fname3))
    
    if verbose:     
        print "Created pixelwise difference ", fname3
    return fname3
    
    
def imagediff_dir(dir1, dir2, dir3="_image_diff", ext='.png', \
                  overwrite=False, verbose=False):
    
    import filecmp,glob
    
    if dir1[-1] == '/': dir1 = dir1[:-1]
    if dir2[-1] == '/': dir2 = dir2[:-1]
    
    files1 = glob.glob("%s/*%s" % (dir1,ext))    
    files1 = [f.replace(dir1+'/','') for f in files1]
    files2 = glob.glob("%s/*%s" % (dir2,ext))
    files2 = [f.replace(dir2+'/','') for f in files2]
    
    files_both = files1 + files2
    files = []
    for f in files_both:
        if f not in files: files.append(f)
    files.sort()
    
    print "Comparing files in the  directory: ", dir1
    print "               with the directory: ", dir2
    
    if os.path.isdir(dir3):
        if (len(os.listdir(dir3)) > 0)  and (not overwrite):
            ans = raw_input("Ok to overwrite files in %s ?  " % dir3)
            if ans.lower() not in ['y','yes']:
                print "*** Aborting"
                return
    else:
        os.system('mkdir -p %s' % dir3)
    startdir = os.getcwd()

    dir1 = '../' + dir1
    dir2 = '../' + dir2
    
    
    os.chdir(dir3)
            
    hname = '_ImageDiffIndex.html'
    hfile = open(hname,'w')
    hfile.write("""<html>
            <h1>Figure comparison and pixel diff</h1>
            Comparing files in the directory: &nbsp; dir1 = %s<br>
            &nbsp;&nbsp; with the directory: &nbsp; dir2 = %s<br>
            &nbsp;&nbsp; with the image extension: %s<p>
            &nbsp;<p>
            <h2>Files:</h2>
            """ % (dir1,dir2,ext))
            
    f_equal, f_diff, f_other = filecmp.cmpfiles(dir1,dir2,files,False)
    
    for f in files:
        fhtml = os.path.splitext(f)[0] + '.html'  ## Specific to Clawpack _plots
        if not os.path.isfile(os.path.join(dir1,fhtml)): fhtml = f
          
        fname1 = os.path.join(dir1,f)
        fname2 = os.path.join(dir2,f)
        fhtml1 = os.path.join(dir1,fhtml)
        fhtml2 = os.path.join(dir2,fhtml)
        
        if f not in files2:
            hfile.write("""
              <table>
              <tr><td><b>%s</b> only appears in dir1</td></tr>
              <tr>
              <td><a href="%s"><img src="%s" width=350 border="1"></a></td> &nbsp; &nbsp;
              <td><img src="XXmissingXX" width=350 height=250 border="1"></td>
              </tr>
              </table><p>"""  % (f, fhtml1, fname1))
        elif f not in files1:
            hfile.write("""
              <table>
              <tr><td></td><td><b>%s</b> only appears in dir2</td></tr>
              <tr>
              <td><img src="XXmissingXX" width=350 height=250 border="1"></td>&nbsp;&nbsp;
              <td><a href="%s"><img src="%s" width=350 border="1"></a></td> &nbsp; 
              </tr>
              </table><p>"""  % (f, fhtml2, fname2))
            
        elif f in f_equal:
            fname1 = os.path.join(dir1,f)
            hfile.write("""
              <table>
              <tr><td><b>%s</b> are identical in dir1 and dir2</td></tr>
              <tr>
              <td><a href="%s"><img src="%s" width=350 border="1"></a></td> &nbsp;&nbsp; 
              </tr>
              </table><p>"""  % (f, fhtml1, fname1))
            # hfile.write("""
            #               <a href="%s"><img src="%s" width=350 border="1"></a>""" \
            #                 % (fname1,fname1))
        else:
            fname3 = f

            fname3 = make_imagediff(fname1,fname2,fname3,verbose=verbose)
            hfile.write("""
              <table>
              <tr><td><b>%s</b> from dir1</td><td><b>%s</b>  from dir2</td>
                  <td>Pixels that differ between the images</td>  </tr>
              <tr>
              <td><a href="%s"><img src="%s" width=350 border="1"></a></td> &nbsp;&nbsp; 
              <td><a href="%s"><img src="%s" width=350 border="1"></a></td> &nbsp;&nbsp; 
              <td><a href="%s"><img src="%s" width=350 border="1"></a></td>  </tr>
              </table><p>""" \
                % (f,f,fhtml1,fname1,fhtml2,fname2,fname3,fname3))
        
    os.chdir(startdir)
    print "To view diffs, open the file ",os.path.join(dir3,hname)
    

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

if __name__ == "__main__":    
    # Parse input arguments
    argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hv",["help","verbose"])
        except getopt.error, msg:
            raise Usage(msg)

        # Default script parameter values
        verbose = False

        # option processing
        for option, value in opts:
            # Script parameters
            if option in ("-v","--verbose"):
                 verbose = True
            if option in ("-h","--help"):
                raise Usage(help_message)

        # Run diff
        if os.path.isfile(args[0]) and os.path.isfile(args[1]):
            sys.exit(imagediff_file(args[0],args[1],verbose=verbose))
        elif os.path.isdir(args[0]) and os.path.isdir(args[1]):
            sys.exit(imagediff_dir(args[0],args[1],verbose=verbose))
        else:
              raise Usage("Both paths must either be files or directories.")
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        sys.exit(2)
    
    
        
        
        
