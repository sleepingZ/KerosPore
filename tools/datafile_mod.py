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

class dataFile:
    def __init__(self,inputfile,outputfile = "data.dreiding"):
        f = open(inputfile,'r')
        lines = f.readlines()
        f.close()
        self.lines = lines
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.basicInfo()
    
    def basicInfo(self):
        lines = self.lines
        for line in lines:
            line_seq = line.split()
            if line_seq[-2:] == ['atom','types']:
                self.nAtomType = int(line_seq[0])
            if line_seq[-2:] == ['bond','types']:
                self.nBondType = int(line_seq[0])
            if line_seq[-2:] == ['angle','types']:
                self.nAngleType = int(line_seq[0])

    def explore(self,Key1,Key2):
        lines = self.lines
        for line in lines:
            if line.startswith(Key1):
                start = lines.index(line) + 1
            if line.startswith(Key2):
                end = lines.index(line) + 1
        return start,end

    def addAtomType(self,type_num_old,type_num_new):
        mass_start,mass_end = self.explore('Masses','Bond Coeffs')
        lines = self.lines
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
                if (line.split()[0] == str(type_num_old)):
                    insert_start_pos = k+1
        for i in range(type_num_new-type_num_old):
            lines.insert(insert_start_pos+i,'%d 1.0\n'%(massnow+i))
    
    def addBondType(self,type_num_old,type_num_new):
        bond_start,bond_end = self.explore('Bond Coeffs','Angle Coeffs')
        lines = self.lines
        for line in lines:
            if line.split()[-2:] == ['bond','types']:
                index = lines.index(line)
                lines[index] = ' '.join([str(type_num_new),'bond','types'])+'\n'
                break
        bondlines = lines[bond_start:bond_end-1]
        bondnow = type_num_old + 1
        for j in range(len(bondlines)):
            line = bondlines[j]
            k = j+bond_start
            if len(line.split())>1:
                if (line.split()[0] == str(type_num_old)):
                    insert_start_pos = k+1
        for i in range(type_num_new-type_num_old):
            lines.insert(insert_start_pos+i,'%d 1.0 1.0\n'%(bondnow+i))

    def addAngleType(self,type_num_old,type_num_new):
        angle_start,angle_end = self.explore('Angle Coeffs','Dihedral Coeffs')
        lines = self.lines
        for line in lines:
            if line.split()[-2:] == ['angle','types']:
                index = lines.index(line)
                lines[index] = ' '.join([str(type_num_new),'angle','types'])+'\n'
                break
        anglelines = lines[angle_start:angle_end-1]
        anglenow = type_num_old + 1
        for j in range(len(anglelines)):
            line = anglelines[j]
            k = j+angle_start
            if len(line.split())>1:
                if (line.split()[0] == str(type_num_old)):
                    insert_start_pos = k+1
        for i in range(type_num_new-type_num_old):
            lines.insert(insert_start_pos+i,'%d 1.0 1.0\n'%(anglenow+i))
    
    def addExtra(self):
        lines = self.lines
        for line in lines:
            if line.strip().split()[-2:] == ['improper','types']:
                index = lines.index(line) + 1
                break
        lines.insert(index,'4 extra special per atom\n')
        lines.insert(index,'6 extra angle per atom\n')
        lines.insert(index,'4 extra bond per atom\n')
            
    def output(self):
        g=open(self.outputfile,'w')
        g.writelines(self.lines)
        g.close()