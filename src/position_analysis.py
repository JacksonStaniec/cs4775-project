#homogenous distribution across windows
#iterate across windows, from -200 to 200 by lengths of 20. For each window, find all the words across the different sequences in that window.
#For each window, offset by 1 and count matches of that word. If shifted by 1,2,3,4,5, only count once.
#Calculate expected number of occurrences for word in that window by normalizing total occurrences by multiplying by number of sequences in that window and dividing by total number of sequences
#Do a chi squared test with degrees of freedom equal to number of windows minus one. Observed count minus expected count squared divided by expected count^2
import numpy as np
from scipy.stats import chisquare
import pandas as pd
def do_position_analysis(fastafile, oligomer_length = 6, window_size = 20, bp_to_read = 200):
    seqs_for_analysis = read_fastaFile(filename=fastafile)

    windows = []
    for x in range(0,2*bp_to_read,window_size):
        toAppend = []
        for t in seqs_for_analysis:
            if t[x:x+window_size]!='':
                toAppend.append(t[x:x+window_size])
        windows.append(toAppend)
    foundwordsAll = identifywordsAll_oligomer(windows)
    chisqTable = pd.DataFrame(columns=('motif','pval','counts'))
    for i, word in enumerate(foundwordsAll):
        # observed =[]
        # for window in windows:
        #     observed.append(sum([s.count(word) for s in window]))
        observed = [sum([s.count(word) for s in window]) for window in windows]
        # expected = [sum([s.count(word) for s in seqs_for_analysis]) for window in windows] 
        expected = [len(window)/len(windows) for window in windows] #number of seqs in that window
        expected = [sum(observed)*window/sum(expected) for window in expected]
        #take sum(observed) and spread it according to len(window)
        # assert False
        # do chisq test
        chisq_statistic = chisquare(observed,expected).pvalue
        chisqTable.loc[i]=(word,chisq_statistic,round(sum(observed)))
    return chisqTable
def read_fastaFile(filename):
    with open(filename, "r") as f:
        output = []
        s = ""
        for l in f.readlines():
            if l.strip()[0] == ">":
                # skip the line that begins with ">"
                if s == "": continue
                output.append(s)
                s = ""
            # keep appending to the current sequence when the line doesn't begin
            # with ">"
            else: s += l.strip()
        output.append(s)
        return output
def identifywords_oligomer(seqList,oligomer_length,givetuple=True): #list of seqs in that window
    foundwords =[]
    returned = []
    for seq in seqList:
        for offset in range(len(seq)-oligomer_length):
            foundword = seq[offset:offset+oligomer_length]
            if foundword in foundwords:
                pass
            else:
                foundwords.append(foundword)
                returned.append((foundword,sum([s.count(foundword) for s in seqList]))) #nonoverlapping
        # if givetuple:
            return returned
    return foundwords
def identifywordsAll_oligomer(seqList,oligomer_length,givetuple=True): #list of seqs in that window
    foundwords =[]
    returned = []
    for window in seqList:
        for seq in window:
            for offset in range(len(seq)-oligomer_length):
                foundword = seq[offset:offset+oligomer_length]
                if foundword in foundwords:
                    pass
                else:
                    foundwords.append(foundword)
                    # returned.append((foundword,sum([s.count(foundword) for s in seqList]))) #nonoverlapping
            # if givetuple:
                # return returned
    return foundwords
class positionanalysis:
    def __init__(self,oligomer_length=6,window_size=20,bp_to_read=200):
        self.oligomer_length = oligomer_length
        self.window_size = window_size
        self.bp_to_read = bp_to_read
    @classmethod
    def read_fasta(cls,filename): 
        with open(filename, "r") as f:
            output = []
            s = ""
            for l in f.readlines():
                if l.strip()[0] == ">":
                    # skip the line that begins with ">"
                    if s == "": continue
                    output.append(s)
                    s = ""
                # keep appending to the current sequence when the line doesn't begin
                # with ">"
                else: s += l.strip()
            output.append(s)
            return output
    @classmethod
    def identifywords(cls,seqList,oligomer_length,givetuple=True): #list of seqs in that window
        foundwords =[]
        returned = []
        for seq in seqList:
            for offset in range(len(seq)-oligomer_length):
                foundword = seq[offset:offset+oligomer_length]
                if foundword in foundwords:
                    pass
                else:
                    foundwords.append(foundword)
                    returned.append((foundword,sum([s.count(foundword) for s in seqList]))) #nonoverlapping
            # if givetuple:
                return returned
        return foundwords
    @classmethod
    def identifywordsAll(cls,seqList,oligomer_length,givetuple=True): #list of seqs in that window
        foundwords =[]
        returned = []
        for window in seqList:
            for seq in window:
                for offset in range(len(seq)-oligomer_length):
                    foundword = seq[offset:offset+oligomer_length]
                    if foundword in foundwords:
                        pass
                    else:
                        foundwords.append(foundword)
                        # returned.append((foundword,sum([s.count(foundword) for s in seqList]))) #nonoverlapping
                # if givetuple:
                    # return returned
        return foundwords
    def loadFasta(self,filename):
        self.seqs_for_analysis = positionanalysis.read_fasta(filename)
        self.windows = []
        for x in range(0,2*self.bp_to_read,self.window_size):
            toAppend = []
            for t in self.seqs_for_analysis:
                if t[x:x+self.window_size]!='':
                    toAppend.append(t[x:x+self.window_size])
            self.windows.append(toAppend)
            self.foundwordsAll = positionanalysis.identifywordsAll(self.windows,self.oligomer_length)
        print('fasta load successful')
    def perform_analysis(self):
        self.chisqTable = pd.DataFrame(columns=('motif','pval','counts'))
        for i, word in enumerate(self.foundwordsAll):
            # observed =[]
            # for window in windows:
            #     observed.append(sum([s.count(word) for s in window]))
            observed = [sum([s.count(word) for s in window]) for window in self.windows]
            # expected = [sum([s.count(word) for s in seqs_for_analysis]) for window in windows] 
            expected = [len(window)/len(self.windows) for window in self.windows] #number of seqs in that window
            expected = [sum(observed)*window/sum(expected) for window in expected]
            #take sum(observed) and spread it according to len(window)
            # assert False
            # do chisq test
            chisq_statistic = chisquare(observed,expected).pvalue
            self.chisqTable.loc[i]=(word,chisq_statistic,round(sum(observed)))
        return self.chisqTable.sort_values('pval').reset_index(drop=True)



# words_in_window=[identifywords(window) for window in windows]
# foundwords = [identifywords(window,givetuple=False) for window in windows]
#calculate chisq assuming homogenous



