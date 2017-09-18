def lmp2list(infile):
    fp = open(infile,"r")
    cmd_list = []
    while True:
        line = str(fp.readline())
        if not line:
            break
        line_seq = line.split()
        cmd_list.append(line_seq)
    fp.close()
    return cmd_list
    
def lmp2str(infile):
    cmd_list = lmp2list(infile)
    cmd_str = [' '.join(item) for item in cmd_list]
    return cmd_str

def list2str(cmd_list):
    cmd_str = [' '.join(item) for item in cmd_list]
    return cmd_str    

def str2file(cmd_str,outfile):
    fp = open('%s.lmp' % outfile,'w')    
    for item in cmd_str:
        fp.write(item+'\n')
    fp.close()
        