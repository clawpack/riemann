import sys
import difflib

def highlight_diffs_html(all_lines=True):

    #all_lines = True  # print all lines (including those with no changes)

    fname1 = "_output/fort.q0010"  #"test1.txt"
    fname2 = "_output_old/fort.q0010"  #test2.txt"
    f1 = open(fname1,'r').readlines()
    f2 = open(fname2,'r').readlines()

    if len(f1) != len(f2):
        print "*** files have different number of lines"
    flen = min(len(f1), len(f2))
    
    S = difflib.SequenceMatcher()
    table1 = []
    table2 = []
    linenos = []
    changed = []
    
    for i in range(flen):
        line1 = f1[i].replace("\n","")
        line2 = f2[i].replace("\n","")
        S.set_seqs(line1, line2)
        opcodes = S.get_opcodes()
        line_equal = ((len(opcodes)==1) and opcodes[0][0]=='equal')
        line_changed = (not line_equal)
    
        if all_lines or line_changed:
            changed.append(line_changed)
            table1line = ""
            table2line = ""
            badline = False
            for opcode in opcodes:
                j1 = opcode[1]
                j2 = opcode[2]
                k1 = opcode[3]
                k2 = opcode[4]
                if (j1 != j1) or (j2 != j2):
                    print "*** Oops, line %s, opcode = %s" % (i, str(opcode))
                    badline = True
                if opcode[0] == "equal":
                    table1line = table1line + line1[j1:j2]
                    table2line = table2line + line2[k1:k2]
                elif opcode[0] == "replace":
                    if j1<j2: table1line = table1line + "<span class=yellow>%s</span>" % line1[j1:j2]
                    if k1<k2: table2line = table2line + "<span class=yellow>%s</span>" % line2[k1:k2]
                    
                #else:
                    #print "*** Oops, line %s, opcode = %s" % (i, str(opcode))
                    #badline = True
            table1.append(table1line)
            table2.append(table2line)
            linenos.append(i)
            if badline: break

    html = open("compare_output.html","w")
    html.write("""
        <html>
        <style>
        span.yellow { background-color: yellow; }
        span.red { background-color: #ff8888; }
        </style>
        <table cellspacing="20" cellpadding="2" border="1">\n""")
    for i in range(len(table1)):
        if changed[i]:
            html.write("<tr><td><span class=red>Line %s</span></td>" % linenos[i])
        else:
            html.write("<tr><td>Line %s</td>" % linenos[i])
        html.write("<td>%s</td><td>%s<td></tr>\n" % (table1[i],table2[i]))

    html.write("</table>\n")
    if badline:
        html.write("<h2>Files disagree after this point --- truncated</h2>\n")
    html.write("</html>\n""")
    html.close()
       
if __name__=="__main__":
    highlight_diffs_html()
            
            

