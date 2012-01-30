"""
Module to diff two images or directories full of images.
"""
import os, sys

def imagediff_file(fname1, fname2):
    
    fname3 = make_imagediff(fname1,fname2)
    
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
    
    
         
def make_imagediff(fname1,fname2,fname3=''):
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
         
    print "Created pixelwise difference ", fname3
    return fname3
    
def imagediff_dir(dir1, dir2, dir3="imagediff_dir", ext='.png', overwrite=False):
    
    import filecmp,glob
    
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
    # dir1 = os.path.abspath(dir1)
    # dir2 = os.path.abspath(dir2)
    # dir3 = os.path.abspath(dir3)
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
        hfile.write("""<h3>%s <a href="%s">in dir1</a> ... <a href="%s">in dir2</a>  </h3>""" \
                % (f, os.path.join(dir1,fhtml), os.path.join(dir2,fhtml)))
        if f not in files2:
            hfile.write("Only appears in dir1\n")
        elif f not in files1:
            hfile.write("Only appears in dir2\n")
        elif f in f_equal:
            hfile.write("Are identical in dir1 and dir2\n")
        else:
            hfile.write("Images differ in the black pixels in the third imageure...<p>\n")
            fname3 = f
            fname1 = os.path.join(dir1,f)
            fname2 = os.path.join(dir2,f)
            fname3 = make_imagediff(fname1,fname2,fname3)
            hfile.write("""
              <a href="%s"><img src="%s" width=350 border="1"></a> &nbsp; 
              <a href="%s"><img src="%s" width=350 border="1"></a> &nbsp; 
              <a href="%s"><img src="%s" width=350 border="1"></a>""" \
                % (fname1,fname1,fname2,fname2,fname3,fname3))
        
    os.chdir(startdir)
    print "To view diffs, open the file ",os.path.join(dir3,hname)
    
if __name__=="__main__":
    
    args = sys.argv[1:]
    imagediff_dir(*args)
    
    
        
        
        
