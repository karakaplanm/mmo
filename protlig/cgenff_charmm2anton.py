#!/usr/bin/env python
###################################################################
# This simple python script takes a cgenff CHARMM stream file generated 
# by paramchemand print out formatted text for corresponding viparr files
# usage: ./cgenff_charmm2anton.py *.str
#
# Please make sure the stream file and viparr ff cgenff_base 
# have the same CGenFF version number
#
# The code are based on older versions from Anton wiki page
# https://wiki.psc.edu/twiki/view/Anton/AddingCustomMoleculesFromCHARMM
# 
# Author: Jing Huang (jing.huang.ff@gmail.com)
# MacKerell lab
# Nov. 2013
###################################################################

import sys, string
from math import fabs

def printfirst(str1):
    print """Add the following text to """+str1+""" :"""
    printlast()

def printlast():
    print """----------------------"""

def printtop(MolTopList):
#   from Anton wiki, made small changes to increase the robustness

    #Molecule Name
    MolName="UNK"
    for line in MolTopList:
        rec = line.split()
        if(rec[0]=="RESI"):
            MolName=rec[1]
    print """   "%s": {"""%(MolName)

    #atoms section
    print """      "atoms": ["""    
    count=0
    Elmts=[0]*32
    for line in MolTopList:
        rec = line.split()
        if(rec[0]=='ATOM'):
            #ATOM H14C HL     0.25
            [AtomName,AtomType,Charge]=rec[1:4]
            el=-1
            if(AtomType[0]=='H'):el=1
            if(AtomType[0]=='C'):el=6
            if(AtomType[0]=='N'):el=7
            if(AtomType[0]=='O'):el=8
            if(AtomType[0]=='P'):el=15
            if(count>0):
                print ", "
            print "         [\"%s\", %d, %s, [\"%s\"]]"%(AtomName,el,Charge,AtomType),
            count=count+1
            Elmts[el]=Elmts[el]+1
    print
    print """      ],"""

    #bonds section
    print """      "bonds": ["""
    count=0
    for line in MolTopList:
        rec = line.split()
        if(rec[0]=='BOND' or rec[0]=='DOUBLE'):
            for i in range(1,len(rec),2):
                if(count>0):
                    print ", "
                print "         [\"%s\", \"%s\"]"%(rec[i],rec[i+1]),
                count=count+1
    print
    print """      ]""",
    
    #impropers if present
    count=0
    for line in MolTopList:
        rec = line.split()
        if(rec[0]=='IMPR' or rec[0]=="IMPH"):
            count=count+1
    if(count>0):
        print ","
        print """      "impropers": ["""
        count=0
        for line in MolTopList:
            rec = line.split()
            if(rec[0]=='IMPR' or rec[0]=="IMPH"):
                for i in range(1,len(rec),4):
                    if(count>0):
                        print ", "
                    print "         [\"%s\", \"%s\", \"%s\", \"%s\"]"%(rec[i],rec[i+1],rec[i+2],rec[i+3]),
                    count=count+1
        print
        print """      ]""",
    print
    print """     }"""


def printparbond(BondsLines):
    for l in BondsLines:
        l1=l.strip()
        f=l1.split()
        memostr=" "
        if f[4]=="!":
            for i in range(5,len(f)):
                memostr = memostr + f[i]+ " "
        print "   {\"type\": [\"%s\", \"%s\"], \"params\": {\"r0\": %s, \"fc\": %s}, \"memo\": \"%s\"},"%(f[0],f[1],f[3],f[2],memostr)
 
def printparangl(AngleLines):
    AnglesHarm=[]
    AnglesUB=[]
    memostr=" "
    for l in AngleLines:
        l1=l.strip()
        if(l1.find("!")>=0):
            l1=l1[0:l1.find("!")]
        f=l1.split()
        if(len(f)==5):
            AnglesHarm.append([f[0],f[1],f[2],f[4],f[3]])
        if(len(f)==7):
            AnglesHarm.append([f[0],f[1],f[2],f[4],f[3]])
            AnglesUB.append([f[0],f[1],f[2],f[6],f[5]])
    if(len(AnglesHarm)>0):
        printfirst("angle_harm")
        for rec in AnglesHarm:
            print "   {\"type\": [\"%s\", \"%s\", \"%s\"], \"params\": {\"theta0\": %s, \"fc\": %s}, \"memo\": \"%s\"},"%(rec[0],rec[1],rec[2],rec[3],rec[4], memostr)
        printlast()
    if(len(AnglesUB)>0):
        printfirst("ureybradley_harm")
        for rec in AnglesUB:
            print "   {\"type\": [\"%s\", \"%s\", \"%s\"], \"params\": {\"r0\": %s, \"fc\": %s}, \"memo\": \"%s\"},"%(rec[0],rec[1],rec[2],rec[3],rec[4], memostr)
        printlast()

def printpardihe(DihedralLines):
    class Dih:
        Phi_eq=None
        C=[0.0]*7
        A=["A1"]*4
        def __init__(self,A):
            self.A=A
            self.C=[0.0]*7
    
    class DihDB:
        DihList=[]
        def GetRecNum(self,A):
            for i in xrange(0,len(self.DihList)):
                if self.DihList[i].A[0]==A[0] and self.DihList[i].A[1]==A[1] and self.DihList[i].A[2]==A[2] and self.DihList[i].A[3]==A[3]:
                    return i
            self.DihList.append(Dih(A))
            return len(self.DihList)-1
        def AddRec(self,A,Kchi,n,delta):
            i=self.GetRecNum(A)
            if self.DihList[i].Phi_eq==None:
                self.DihList[i].Phi_eq=delta
            else:
                if fabs(self.DihList[i].Phi_eq-delta)>0.001:
                    print "Error: Phi_eq is not equal to previous Phi_eq",delta,self.DihList[i].A, self.DihList[i].Phi_eq,self.DihList[i].C
            self.DihList[i].C[0]=self.DihList[i].C[0]+fabs(Kchi) # charmm always have Kchi positive
            self.DihList[i].C[n]=Kchi
    
    m_DihDB=DihDB()
    for l in DihedralLines:
        l1=l.strip()
        if(l1.find("!")>=0):
            l1=l1[0:l1.find("!")]
        f=l1.split()
        # Anton prefer to have Phi_eq to be 0 rather than 180
        a1=float(f[4])
        if f[6]=="180.00":
            f[6]="0.00"
            a1 = -1.0 * a1
        m_DihDB.AddRec([f[0],f[1],f[2],f[3]],a1,int(f[5]),float(f[6]))
        #m_DihDB.AddRec([f[0],f[1],f[2],f[3]],float(f[4]),int(f[5]),float(f[6]))
    memostr=" "
    for i in xrange(0,len(m_DihDB.DihList)):
        print "   {\"type\": [\"%s\", \"%s\", \"%s\", \"%s\"], \"params\": {\"phi0\": %.3f, \"fc0\": %.3f, \"fc1\": %.3f, \"fc2\": %.3f, \"fc3\": %.3f, \"fc4\": %.3f, \"fc5\": %.3f, \"fc6\": %.3f}, \"memo\": \"%s\"},"% \
        (m_DihDB.DihList[i].A[0],m_DihDB.DihList[i].A[1],m_DihDB.DihList[i].A[2],m_DihDB.DihList[i].A[3],\
        m_DihDB.DihList[i].Phi_eq,\
        m_DihDB.DihList[i].C[0],m_DihDB.DihList[i].C[1],m_DihDB.DihList[i].C[2],m_DihDB.DihList[i].C[3],\
        m_DihDB.DihList[i].C[4],m_DihDB.DihList[i].C[5],m_DihDB.DihList[i].C[6],memostr)

def printparimpr(BondsLines):
    for l in BondsLines:
        l1=l.strip()
        f=l1.split()
        memostr=" "
        if f[7]=="!":
            for i in range(8,len(f)):
                memostr = memostr + f[i]+ " "
        print "   {\"type\": [\"%s\", \"%s\", \"%s\", \"%s\"], \"params\": {\"phi0\": %s, \"fc\": %s}, \"memo\": \"%s\"},"%(f[0],f[1],f[2],f[3],f[6],f[4],memostr)
 

if __name__ == "__main__":
    try:
        strfile = open(sys.argv[1],'r')
    except IOError:                     #catch exception
        print ('Files do not exist!\n')

#   remove blank and comment lines
    noblanks = filter(lambda x: len(x.strip())>0, strfile.readlines())
    nocomments = filter(lambda x: x.strip()[0] not in ['*','!'], noblanks)

#   split top and par
    for line in nocomments:
        if (line[0:3]=="END"):
            splitindex=nocomments.index(line)
            break
    toppart = nocomments[0:splitindex+1]
    parpart = nocomments[splitindex+1:]

    printfirst("templates")
    printtop(toppart)
    printlast()

#   split par into different part:
    for line in parpart:
        if (line[0:4]=="BOND"):
            splitind1=parpart.index(line)
        if (line[0:4]=="ANGL"):
            splitind2=parpart.index(line)
        if (line[0:4]=="DIHE"):
            splitind3=parpart.index(line)
        if (line[0:4]=="IMPR"):
            splitind4=parpart.index(line)
#   cgenff stream won't have "CMAP", "NONB", "HBON", "NBFI" terms
        if (line[0:3]=="END"):
            splitind5=parpart.index(line)

    parbond = parpart[splitind1+1:splitind2]
    parangl = parpart[splitind2+1:splitind3]
    pardihe = parpart[splitind3+1:splitind4]
    parimpr = parpart[splitind4+1:splitind5]

    if (len(parbond)>0):
        printfirst("stretch_harm") #printfirst("bonds.Harm")
        printparbond(parbond)
        printlast()
    if (len(parangl)>0):
        printparangl(parangl)
    if (len(pardihe)>0):
        printfirst("dihedral_trig") #printfirst("propers.Proper_Trig")
        printpardihe(pardihe)
        printlast()
    if (len(parimpr)>0):
        printfirst("improper_harm") #printfirst("impropers.Improper_Harm")
        printparimpr(parimpr)
        printlast()

    exit()


