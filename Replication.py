# ===================================
#Fnding Hidden Messages 
# ===================================

# -----------------------------------
#Pattern Count function

#Input: DNA String(Text), Pattern string (i.e AGG) 
#Output: How a

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count
# -----------------------------------

# -----------------------------------
#Pattern Matching function

#Input: DNA String(Text), Pattern string (i.e AGG) 
#Output: How a

def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    pattern_length = len(Pattern)
    genome_length = le(Genome)
    for char in range(genome_length - pattern_length +1) :
        if Genome[char:char + pattern_length] == Pattern : 
            positions.append(char)
    return positions
# -----------------------------------
#FrequencyMap function

#Goal: Ranging through the text, with a sliding window of length "k". 

#Input: DNA String (Text), k(of k-mer)integer
#Output: Frequency 

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = freq.get(Pattern,0) + 1
    return freq
# -----------------------------------

# -----------------------------------
#FrequentWords function here. 

#Goal: Find the most frequent words/patterns

#Subroutine: FrequencyMap

#Input: A string Text and an integer k 
#Output: Returned list of the most frequent k-mers 


def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m : 
            mer = (key)
            words.append(mer)
    return words
# -----------------------------------

# -----------------------------------
#Reverse function 

#Input: A Text string 
#Output: Reversed Text string 

def Reverse(Pattern):
    rev = ' '
    for char in Pattern : 
        rev = char + rev 
    return rev

# -----------------------------------

# -----------------------------------
#Complement function 

#Input: A DNA string Pattern
#Output: Complementary (not reversed) string, with each nucleotide replaced

def Complement(Pattern):
    comp_lst = ' '
    for char in Pattern :
        if char == 'A':
            comp_lst = comp_lst + 'T'
        elif char == 'T' :
            comp_lst = comp_lst + 'A'
        elif char == 'C':
            comp_lst = comp_lst + 'G'
        elif char == 'G' :
            comp_lst = comp_lst + 'C'
    return comp_lst

# -----------------------------------

# -----------------------------------
#ReverseComplement Function 

#Input: DNA string Pattern
#Output: Reverse complement of Pattern

#Subroutines: Reverse(Pattern), Complement(Pattern)

def ReverseComplement(Pattern):   
    Pattern = Reverse(Pattern) # reverse all letters in a string (get 5'-3' directionality)
    Pattern = Complement(Pattern) #complement each latter in a string
    return Pattern
# -----------------------------------

# -----------------------------------
#PatternMatching function 

#Goal: Finding all occurrences of a pattern in a string 

#Input: Pattern and Genome strings 
#Output: All starting postions of where a Pattern substring appears in Genome

def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    pattern_length = len(Pattern)
    genome_length = len(Genome)
    for char in range(genome_length - pattern_length +1) :
        if Genome[char:char + pattern_length] == Pattern : 
            positions.append(char)
    return positions
# -----------------------------------

# -----------------------------------
#SymbolArray function 

#Input: Genome and symbol strings
#Output: An array of the frequency of symbol in question, per position of Genome

#Subroutine: PatternCount

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array
# -----------------------------------

# -----------------------------------
#FasterSymbolArray

#Only considering two symbols each time the "window" slides 

#Input: Pattern and Genome strings 
#Output: All starting postions of where a Pattern substring appears in Genome

#Subroutine: PatternCount


def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array
# -----------------------------------

# -----------------------------------
#SkewArray Function 

#Goal: To return a list of skew array integer values, per value 

#Skew --> setting Skew[i] = number of occurrences of "G" - number of occurrences of C in the first "i" nucleotides of string Genome 

#Input: Genome string 
#Output: Skew array of Genome as a list 

def SkewArray(Genome):
    Skew = [0]
    for i in range(len(Genome)):
        if Genome[i]=='G':
            Skew.append(Skew[i]+1)
        elif Genome[i] == 'C' :
            Skew.append(Skew[i]-1)
        else :
            Skew.append(Skew[i])
    return Skew
# -----------------------------------

# -----------------------------------
#Minimum Skew function 

#Goal: Finding a position in a genome where the skew diagram attains a minimum 

#Input: Genome String 
#Output: All integers "i'minimizing Skew[i], from values of "i"


#Subroutine: SkewArray 
def MinimumSkew(Genome):
    positions = [] 
    skew_values = SkewArray(Genome)
    min_val = min(skew_values)
    for i in range (len(skew_values)) :
        if skew_values[i]==min_val:
            positions.append(i)
    return positions
# -----------------------------------

# -----------------------------------
#HammingDistance Function 

#Hamming distance --> Total number of mismatches between strings "p" and "q" 

#Mismatch --> If at position "i" in "k-mers" between strings "p" and "q", it is not the same 

#Input: Two strings of equal length 
#Output: Hammingdistance between these strings 

def HammingDistance(p, q):
    dist = 0
    for i in range(len(p)):
        if p[i] != q[i] :
            dist = dist + 1
    return dist
# -----------------------------------
#ApproximatePatternMatching function 

#Goal: Fin all approxiate occurences (at most "d" mismatches) of a pattern in a string

#Input: Pattern and Text strings, d integer
#Output: All starting positions where substring Patterns appear in Text, with at most "d" mismatches 

#Subroutine: HammingDistance
def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] 
    for i in range(len(Text) - len(Pattern)+1):
        sub = Text[i:i+len(Pattern)]
        e = HammingDistance(Pattern, sub)
        if e <= d:
            positions.append(i)
    return positions


# -----------------------------------

# -----------------------------------
#ApproximatePatternCount function 

#Goal: Returning the number of occurrences of Pattern in Text with at most "d" mismatches 

#Input:Pattern and Texxt strings, integer d 
#Output: Number of times Pattern appeas in Text with at most "d" mismatches 

#Subroutine: HammingDistance 

def ApproximatePatternCount(Pattern, Text, d):
    count = 0 
    for i in range(len(Text) - len(Pattern)+1):
        sub = Text[i:i+len(Pattern)]
        e = HammingDistance(Pattern, sub)
        if e <= d:
            count = count +1
    return count

# -----------------------------------
