# Developer:
# Karim Tarek Emam
# ----------------------------------------

# Used Libraries
from Bio import SeqIO
# End of Used Libraries
# Global Variables
seeds = {}
database = []
acceptedseeds = []
allHSP = []
resHSP = []
# End of Global Variables
# Class HSP 
class HSP:
    Mscore = 0
    QuerySeq = ""
    DBSeq = ""
    seed = ""
    SeqID = 0
    QueryStartIndex = 0
    QueryEndIndex = 0
    DataBaseStartIndex = 0
    DataBaseEndIndex = 0

    def CheckMax(self, score):
        # This Function Check if current Score is greater than Mscore and Let it to be new Mscore
        update = False
        if score > self.Mscore:
            self.Mscore = score
            update = True
        return update

    def EndExtend(self, score):
        # This Function if the total Score of extended sequence drops more than the stopValue (5)
        end = False
        StopValue = 5
        diff = self.Mscore - score
        if diff > StopValue:
            end = True
        return end

    def updateHSP(self, seqid,seed,hitq, hitd, siq, eiq, sid, eid):
        # This Function update the saved data for extended sequence if there is New Mscore
        self.SeqID = seqid
        self.seed = seed
        self.QuerySeq = hitq
        self.DBSeq = hitd
        self.QueryStartIndex = siq
        self.QueryEndIndex = eiq
        self.DataBaseStartIndex = sid
        self.DataBaseEndIndex = eid

    def RequiredHSP(self, HSPThreshold):
        # This Function Check that the extended sequence satisfies that is equal or greater than HSPThreshold
        if self.Mscore < HSPThreshold:
            stop = True
        else:
            stop = False
        return stop

    def Printdata(self):
        # This Function Display Result of the Search
        index = str(self.SeqID)
        Score = str(self.Mscore)
        startindexQ = str(self.QueryStartIndex)
        endindexQ = str(self.QueryEndIndex)
        startindexD = str(self.DataBaseStartIndex)
        endindexD = str(self.DataBaseEndIndex)
        print("Sequence " + index)
        print("Hit " + self.seed + ":")
        print("Database Sequence: " + self.DBSeq)
        print("Staring Index : " + startindexD + " , End Index : " + endindexD)
        print("Query Sequence:    " + self.QuerySeq)
        print("Staring Index : " + startindexQ + " , End Index : " + endindexQ)
        print("Score: " + Score)
        print("----------------------------------")
# End of Class HSP 

def Readfile():
    # This Loop gets the Sequences of Fasta file to List database
    for record in SeqIO.parse("Project\sequences.fasta", "fasta"):
        database.append(record.seq)

def ReturnBlosum62Value(c1, c2):
    # Values of each protein in Blosum62
    blosum62 = {
        ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
        ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
        ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
        ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
        ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
        ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
        ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
        ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
        ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
        ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
        ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
        ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
        ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
        ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
        ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
        ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
        ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
        ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
        ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
        ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
        ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
        ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
        ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
        ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
        ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
        ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
        ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
        ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
        ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
        ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
        ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
        ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
        ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
        ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
        ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
        ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
        ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
        ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
        ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
        ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
        ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
        ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
        ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
        ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
        ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
        ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
        ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
        ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
        ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
        ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
        ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
        ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
        ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
        ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
        ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
        ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
        ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
        ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
        ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
        ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
        ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
        ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
        ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
        ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
        ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
        ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
        ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
        ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
        ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
    }
    # To Handle Error of No value found for that Key as for example ('X', 'Q') value is  -1 but we can't find ('Q','X')
    try: 
        return blosum62[(c1, c2)]
    except:
        return blosum62[(c2, c1)]

def LCR(query):
    # This Function to Remove Low compexity Regions
    LQuery = list(query) # Converts the String "Query" to a List which will be easily access each letter
    for i in range(len(LQuery)):
        if i+3 < len(LQuery): # check current index in List plus 3 still in the boundries of the List
            LetterOne = ReturnBlosum62Value(LQuery[i],LQuery[i])
            LetterTwo = ReturnBlosum62Value(LQuery[i+1],LQuery[i+1])
            LetterThree = ReturnBlosum62Value(LQuery[i+2],LQuery[i+2])
            LetterFour = ReturnBlosum62Value(LQuery[i+3],LQuery[i+3])
            # Get the 3 consecutive after the current letter and check if first, third and second, fourh have the same score 
            if LetterOne == LetterThree and LetterTwo == LetterFour:
                LQuery[i] = "X" # if condition satisfy our case then Replace the current Letter with X
    NewQuery = ToStr(LQuery)  #Convert back the List to String and return it
    return NewQuery

def words(query, worldlen):
    # This Function will split the query to some Letters of length of wordlen
    words = [] # create a List to store the generated words
    for i in range(0, len(query)): # Loop through the Query String
        if i + worldlen - 1 < len(query): # Check if current index + worldlen - 1 is within the boundries of Query length
            words.append(query[i:i + worldlen]) # insert the word to words List which contains the current Letter + Letters of wordlen
    return words

def ToStr(s):
    # This Function to convert the List to a String
    str1 = ""
    for element in s: # Loop Through each element in the List and add it to the string
        str1 += element
    return str1

def NeighbourWords(oword):
    # This Function generates the Neighbour of my current word
    ProteinLetters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'R', 'S', 'T', 'V', 'W', 'Y'] # A list contains 20 amino acid that can replace on letter in the word
    mutantword = list(oword) # Convert the string to List
    # The Nested For Loops will loop through each letter in mutantword List and in each iteration will replace the letter with one Letter in ProteinLetters List
    for i in range(len(mutantword)):
        for j in range(len(ProteinLetters)):
            mutantword[i] = ProteinLetters[j]
            mword = ToStr(mutantword)
            score = 0
            for k in range(len(mword)):
                score += ReturnBlosum62Value(mword[k], oword[k]) # Get the whole score for the word after replacement and save it in a dictionary
            seeds.update({(mword, oword): score}) # Save the original and modified word as key for score as value
        mutantword = list(oword) # remove the modification be for enetering the loop for the next Letter to be changed

def wantedseeds(wordthresold):
    # Loop on Every value in the dic and remove (key and value) if value is less than wordthresold
    for k in seeds.copy():
        if seeds[k] < wordthresold:
            del seeds[k]

def findacceptedseeds():
    # This Function Search for the seeds in the Database and store them if found in database
    se = []
    key = list(seeds.keys()) # Create a list of keys of seeds dictionary
    # This Nested loop to loop on every seed and evey sequence in database
    for i in database:
        for j in key:
            if j[0] in i: # Check if seed is a part of this sequence then append it to se List
                se.append(j)
    [acceptedseeds.append(x) for x in se if x not in acceptedseeds]  # Removes duplications in se list by adding to the acceptedseeds list

def getscore(st1, st2):
    # This Function get the score of word compared with the database sequence
    size = len(st1)
    score = 0
    for i in range(size):
        score += ReturnBlosum62Value(st1[i], st2[i])
    return score

def MaxScore(hsp, seqid, seed, thescore, quer, seque, sique, eique, sida, eida):
    # This Function checks the MaxScore if less than the new score so save that score and update saved data to our new data
    # if Mscore greater than the new score check if score doesn't drop greater than stopValue (5)
    stop = False
    if hsp.CheckMax(thescore):
        hsp.updateHSP(seqid, seed, quer, seque, sique, eique, sida, eida)
    else:
        stop = hsp.EndExtend(thescore)
    return stop

def extend(query, wordlen, HSPthreshold):
    # This Function Extends the acceptedseeds to a sequence of score more than or equal to HSPthreshold
    for i in range(len(database)): # Loop over each Sequence in Database
        hsp = HSP() # create an object from class HSP to save data of Sequence that will be extended
        for j in acceptedseeds:
            if j[0] in database[i]:
                # Get the start index and end index of Query and database sequence
                qsi = query.index(j[1])
                dsi = database[i].index(j[0])
                leftd = dsi
                rightd = dsi + wordlen
                leftq = qsi
                rightq = qsi + wordlen
                size = len(query) - 1
                hitq = ""
                hitd = ""
                flag = False
                # Extend the sequence until you can't extend from left or right or score doesn't drop more than StopValue 
                while leftq > 0 and rightq <= size and not flag:
                    leftq = leftq - 1
                    rightq = rightq + 1
                    leftd = leftd - 1
                    rightd = rightd + 1
                    hitq = query[leftq:rightq]
                    hitd = database[i][leftd:rightd]
                    score = getscore(hitq, hitd)
                    flag = MaxScore(hsp, i, j[0], score, hitq, hitd, leftq, rightq, leftd, rightd) # Call this function to check new score and Mscore and check if score drops
                # Extend from left since we can't from the right as query ends from right
                while leftq > 0 and not flag:
                    leftq = leftq - 1
                    leftd = leftd - 1
                    hitq = query[leftq] + hitq
                    hitd = database[i][leftd] + hitd
                    score = getscore(hitq, hitd)
                    flag = MaxScore(hsp, i, j[0], score, hitq, hitd, leftq, rightq, leftd, rightd) # Call this function to check new score and Mscore and check if score drops
                # Decrease rightq and rightd to avoid Error of index not in boundries
                rightq = rightq - 1 
                rightd = rightd - 1
                # Extend from right since we can't from the left as query ends from left
                while rightq < size and not flag:
                    rightq = rightq + 1
                    rightd = rightd + 1
                    hitd = hitd + database[i][rightd]
                    hitq = hitq + query[rightq]
                    score = getscore(hitq, hitd)
                    flag = MaxScore(hsp, i, j[0], score, hitq, hitd, leftq, rightq, leftd, rightd) # Call this function to check new score and Mscore and check if score drops
            # Check Score of extended sequence satisfies the user's HSPthreshold
            if hsp.RequiredHSP(HSPthreshold):
                continue
            else:
                allHSP.append(hsp)

def reworkHSP():
    # This Function removes any duplicates and sort the results in Descending order
    allHSP.sort(key=lambda x: x.Mscore, reverse=True) 
    for i in allHSP:
        if i not in resHSP:
            resHSP.append(i)

def DisplayHSP():
    # This Function Display the result of Search if found 
    if resHSP:
        for i in resHSP:
            i.Printdata()
    else:
        print("Query wasn't found in Database!!!!")


if __name__ == '__main__':
    Readfile()
    Query = input("Enter The Query: ").upper()
    WordLength = int(input("Enter Word Length: "))
    WordThreshold = int(input("Enter Word Threshold: "))
    HSPThreshold = int(input("Enter HSP Threshold: "))
    UpdatedQuery = LCR(Query)
    print(UpdatedQuery)
    aa = words(UpdatedQuery, WordLength)
    for i in aa:
        NeighbourWords(i)
    wantedseeds(WordThreshold)
    findacceptedseeds()
    extend(UpdatedQuery, WordLength, HSPThreshold)
    reworkHSP()
    DisplayHSP()
