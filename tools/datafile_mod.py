def insertCharge0(filename):
    g=open(filename+'_full','w+')
    start,end = explore(filename,'Atoms','Bonds')
    linenow = 0
    f=open(filename,'r')
    while True:
        line=f.readline()
        if not line:
            break
        linenow += 1
        if (linenow>start) and (linenow<end):
            insertflag = 1
        else:
            insertflag = 0
        writeCharge(g,line,insertflag)
    f.close()
    g.close()
		
def explore(filename,Key1,Key2):
    f=open(filename,'r')
    linenow = 0
    while True:
        line=f.readline()
        if not line:
            break
        linenow += 1
	if line.startswith(Key1):
            start = linenow
        if line.startswith(Key2):
            end = linenow
    f.close()
    return start,end

			        
def writeCharge(g,line,insertflag):
    if (insertflag):
        line=line.split()
        if len(line)>1:
            line.insert(3,'0')
            line=' '.join(line)+'\n'
        else:
            line='\n'
        g.write(line)
    else:
        g.write(line)
  
def addAtomType(filename,type_num_old,type_num_new):
    mass_start,mass_end = explore(filename,'Masses','Bond Coeffs')
    f=open(filename,'r')
    lines =f.readlines()
    f.close()
    for line in lines:
        if line.split()[-2:] == ['atom','types']:
            index = lines.index(line)
            lines[index] = ' '.join([str(type_num_new),'atom','types'])+'\n'
            break
    masslines = lines[mass_start:mass_end-1]
    massnow = type_num_old + 1
    for j in range(len(masslines)):
        line = masslines[j]
        k = j+mass_start
        if len(line.split())>1:
#            if (line.split()[0] == str(type_num_old)) and (line.split()[0] <str(type_num_new)):
#                lines.insert(linenow,'%d 1.0\n'%massnow)
#                massnow += 1
            if (line.split()[0] == str(type_num_old)):
                insert_start_pos = k+1
    for i in range(type_num_new-type_num_old):
        lines.insert(insert_start_pos+i,'%d 1.0\n'%(massnow+i))
    f.close()
    g=open('data.dreiding','w')
    g.writelines(lines)
    g.close()
    