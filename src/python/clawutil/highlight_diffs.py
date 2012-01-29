import sys,os

def highlight_diffs_html(fname1, fname2, print_all_lines=True):

    #print_all_lines = True  # print all lines (including those with no changes)

    #fname1 = "_output/fort.q0010"  #"test1.txt"
    #fname2 = "_output_old/fort.q0010"  #test2.txt"
    f1 = open(fname1,'r').readlines()
    f2 = open(fname2,'r').readlines()

    if len(f1) != len(f2):
        print "*** files have different number of lines"
    flen = min(len(f1), len(f2))
    
    table1 = []
    table2 = []
    linenos = []
    changed = []
    
    for i in range(flen):
        line1 = f1[i].replace("\n","")
        line2 = f2[i].replace("\n","")
        badline = False
        if line1==line2:
            line_changed = False
        else:
            line_changed = True
            len_line = max(len(line1),len(line2))
            if (len(line1)<len_line):
                badline = True  # signal break after this line
                line1 = line1.ljust(len_line)  # pad the line
            if (len(line2)<len_line):
                badline = True  # signal break after this line
                line2 = line2.ljust(len_line)  # pad the line

            toggle = []   # keep track of where switch between 
                          # matching and not matching sections
            same = True
            for j in range(len_line):
                if (same and (line1[j] != line2[j])) or \
                      ((not same) and (line1[j] == line2[j])):
                    same = not same
                    toggle.append(j)

            if len(toggle)==0:
                print "*** Error: toggle should be nonempty"
                print "*** Aborting"
                raise Exception()
    
        if print_all_lines and (not line_changed):
            changed.append(False)
            table1.append(line1)
            table2.append(line2)
            linenos.append(i)
        elif line_changed:
            changed.append(True)
            table1line = ""
            table2line = ""
            badline = False
            same = True
            for k in range(len(toggle)+1):
                if k==0:
                    j1 = 0
                    j2 = toggle[0]
                else:
                    j1 = toggle[k-1]
                    if k<len(toggle):
                        j2 = toggle[k]
                    else:
                        j2 = len_line+1
                if same:
                    table1line = table1line + line1[j1:j2]
                    table2line = table2line + line2[j1:j2]
                else:
                    table1line = table1line + \
                           "<span class=yellow>%s</span>" % line1[j1:j2]
                    table2line = table2line + \
                           "<span class=yellow>%s</span>" % line2[j1:j2]
                same = (not same)
            table1.append(table1line)
            table2.append(table2line)
            linenos.append(i)

            # Maybe better not break, might want remainder of file too...
            #if badline: break

    numchanges = sum(changed)

    hname = "diffs_all_lines.html"
    html = open(hname,"w")
    print "Point your browser to: "
    print "  all %s lines with diffs: %s" % (flen, hname)
    html.write("""
        <html>
        <style>
        span.yellow { background-color: yellow; }
        span.red { background-color: #ff8888; }
        span.green { background-color: #88ff88; }
        </style>
        <h2>All lines of files with highlighted diffs (%s lines)</h2>
        <h3>
        See also: <a href="diffs_changed_lines.html">Only changed lines (%s lines)</a>
        </h3>
        <table style="border-collapse: collapse; padding: 50px;">
        <col span="3" style="padding: 50px; background-color: #FFFFFF;
        border: 2px solid #000000;" />\n""" % (flen,numchanges))
    html.write("<td></td><td><b>%s&nbsp;&nbsp;</b></td><td><b>%s&nbsp;&nbsp;</b></td></tr>\n" \
                % (os.path.abspath(fname1),os.path.abspath(fname2))) 
    html.write("<tr><td></td><td>__________________</td><td>__________________</td></tr>\n") 

    for i in range(len(table1)):
        if changed[i]:
            html.write("<tr><td><span class=red>Line %s</span></td>" \
                 % (linenos[i]+1))
        else:
            html.write("<tr><td><span class=green>Line %s</span></td>" \
                 % (linenos[i]+1))
        html.write("<td>%s</td><td>%s</td></tr>\n" % (table1[i],table2[i]))

    html.write("</table>\n")
    #if badline:
        #html.write("<h2>Files disagree after this point --- truncated</h2>\n")
    html.write("</html>\n""")
    html.close()
       
    # Only changed lines:
    hname = "diffs_changed_lines.html"
    html = open(hname,"w")
    print "  only %s lines with changes: %s" % (numchanges, hname)
    html.write("""
        <html>
        <style>
        span.yellow { background-color: yellow; }
        span.red { background-color: #ff8888; }
        span.green { background-color: #88ff88; }
        </style>
        <h2>Displaying only changed lines of files (%s lines)</h2>
        <h3>
        See also: <a href="diffs_all_lines.html">All lines with diffs (%s lines)</a>
        </h3>
        <table style="border-collapse: collapse; padding: 50px;">
        <col span="3" style="padding: 50px; background-color: #FFFFFF;
        border: 2px solid #000000;" />\n""" % (numchanges,flen))
    html.write("<td></td><td><b>%s&nbsp;&nbsp;</b></td><td><b>%s&nbsp;&nbsp;</b></td></tr>\n" \
                % (os.path.abspath(fname1),os.path.abspath(fname2))) 
    html.write("<tr><td></td><td>__________________</td><td>__________________</td></tr>\n") 

    for i in range(len(table1)):
        if changed[i]:
            html.write("<tr><td><span class=red>Line %s</span></td>" \
                 % (linenos[i]+1))
            html.write("<td>%s</td><td>%s</td></tr>\n" % (table1[i],table2[i]))

    html.write("</table>\n")
    #if badline:
        #html.write("<h2>Files disagree after this point --- truncated</h2>\n")
    html.write("</html>\n""")
    html.close()
       
if __name__=="__main__":
    args = sys.argv[1:]
    highlight_diffs_html(*args)
            
            

