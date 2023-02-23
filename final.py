#allow python program to accept system arguments
import sys

if len(sys.argv) < 4:
    print("please provide the input fasta file as well as the alignment flag and the config file for parameters. For more information refer to the README file")
    sys.exit(1)

# Read in the fasta file
with open(sys.argv[1], 'r') as f:
    fasta_lines = f.readlines()

# Extract the sequences from the fasta file and store them in the variables seq1 and seq2 respectively
seq1 = ''
seq2 = ''
current_sequence = ''
for line in fasta_lines:
    line = line.strip()
    if line.startswith('>'):
        if current_sequence == '':
            current_sequence = 'seq1'
        else:
            current_sequence = 'seq2'
    else:
        if current_sequence == 'seq1':
            seq1 += line
        else:
            seq2 += line


class DP_cell:
    def __init__(self, score=0, pointer=None, seq=""):
        self.score = score
        self.pointer = pointer
        self.seq = seq

#for clarity i implemented the two algorithms in two separate functions 

#NEEDLEMAN WUNSCH FOR GLOBAL ALIGNMENT


def needleman_wunsch(s1, s2, match=1, mismatch=-2, gap=-5, gap_ext=-2):
    # Initialize DP table
    dp = [[DP_cell() for j in range(len(s2)+1)] for i in range(len(s1)+1)]
    
    # Initialize first row and column
    for i in range(1, len(s1)+1):
        dp[i][0].score = gap + (i-1) * gap_ext
        dp[i][0].pointer = "up"
        dp[i][0].seq = s1[i-1]
    for j in range(1, len(s2)+1):
        dp[0][j].score = gap + (j-1) * gap_ext
        dp[0][j].pointer = "left"
        dp[0][j].seq = s2[j-1]
    
    # Fill in DP table
    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            match_mismatch = match if s1[i-1] == s2[j-1] else mismatch
            scores = [dp[i-1][j-1].score + match_mismatch,
                      dp[i-1][j].score + gap + gap_ext,
                      dp[i][j-1].score + gap + gap_ext]
            dp[i][j].score = max(scores)
            dp[i][j].pointer = ["diag", "up", "left"][scores.index(dp[i][j].score)]
            dp[i][j].seq = s1[i-1] + s2[j-1]
    
    # Traceback and generate alignments
    align1 = ""
    align2 = ""
    i = len(s1)
    j = len(s2)
    while i > 0 or j > 0:
        if dp[i][j].pointer == "diag":
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            i -= 1
            j -= 1
        elif dp[i][j].pointer == "up":
            align1 = s1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = s2[j-1] + align2
            j -= 1
    
    # Calculate scores and generate report
    score = sum([match if align1[i] == align2[i] else mismatch for i in range(len(align1))])
    gap_count = align1.count("-")
    identity = str(align1.count("|")) + "/" + str(len(align1)) + " (" + str(int(100 * align1.count("|") / len(align1))) + "%)"
    
    report = ""
    report += "Scores:    match = " + str(match) + ", mismatch = " + str(mismatch) + ", h =" + str(gap) + ", g = " + str(gap_ext) + "\n\n"
    report += "Sequence 1 = \"" + s1 + "\", length = " + str(len(s1)) + " characters\n"
    report += "Sequence 2 = \"" + s2 + "\", length = " + str(len(s2)) + " characters\n\n"
    
    pos = 0
    while pos < len(align1):
        string2add = "s1  " + str(pos+1) + "   "
        report += "s1  " + str(pos+1) + "   "
        for i in range(pos, min(pos+60, len(align1))):
            report += align1[i]
        report += "  " + str(min(pos+60, len(align1))) + "\n"
        report += " "*len(string2add)
        for i in range(pos, min(pos+60, len(align1))):
            
            if align1[i] == align2[i]:
                report += "|"
            else:
                report += " "
        report += "\n"
        report += "s2  " + str(pos+1) + "   "
        for i in range(pos, min(pos+60, len(align2))):
            report += align2[i]
        report += "  " + str(min(pos+60, len(align2))) + "\n\n"
        pos += 60

    # Format report and print results
    matches = 0
    mismatches = 0
    opening_gaps = 0
    gap_extensions = 0
    for i in range(len(align1)):
        if align1[i] == align2[i]:
            matches += 1
        elif align1[i] == "-" or align2[i] == "-":
            if i == 0 or align1[i-1] == "-" or align2[i-1] == "-":
                opening_gaps += 1
            else:
                gap_extensions += 1
        else:
            mismatches += 1

    report += "Report:" + "\n\n"
    report += "Global optimal score = " + str(score) + "\n\n"
    report += "Number of:  matches = " + str(matches) + ", mismatches = " + str(mismatches) + ", opening gaps = " + str(opening_gaps) + ", gap extensions = " + str(gap_extensions) + "\n\n"
    report += "Identities = " + str(matches) + "/" + str(len(align1)) + " (" + str(round(matches/len(align1)*100, 2)) + "%), "
    report += "Gaps = " + str(opening_gaps + gap_extensions) + "/" + str(len(align1)) + " (" + str(round((opening_gaps + gap_extensions)/len(align1)*100, 2)) + "%)\n"

    return report



#SMITH WATERMAN FOR LOCAL ALIGNMENT


def smith_waterman(s1, s2, match=1, mismatch=-2, gap=-5, gap_ext=-2):
    n, m = len(s1), len(s2)
    F = [[0] * (m+1) for i in range(n+1)]
    pointer = [[0] * (m+1) for i in range(n+1)]
    
    max_score, max_i, max_j = 0, 0, 0
    
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = F[i-1][j-1] + (match if s1[i-1] == s2[j-1] else mismatch)
            delete = F[i-1][j] + gap
            insert = F[i][j-1] + gap
            F[i][j] = max(0, diag, delete, insert)
            if F[i][j] == diag:
                pointer[i][j] = 1
            elif F[i][j] == delete:
                pointer[i][j] = 2
            elif F[i][j] == insert:
                pointer[i][j] = 3
            else:
                pointer[i][j] = 0
            if F[i][j] > max_score:
                max_score, max_i, max_j = F[i][j], i, j

    align1, align2 = "", ""
    while F[max_i][max_j] != 0:
        if pointer[max_i][max_j] == 1:
            align1 += s1[max_i-1]
            align2 += s2[max_j-1]
            max_i -= 1
            max_j -= 1
        elif pointer[max_i][max_j] == 2:
            align1 += s1[max_i-1]
            align2 += "-"
            max_i -= 1
        elif pointer[max_i][max_j] == 3:
            align1 += "-"
            align2 += s2[max_j-1]
            max_j -= 1

    align1 = align1[::-1]
    align2 = align2[::-1]

    match_count, mismatch_count, gap_open_count, gap_ext_count = 0, 0, 0, 0


    report = "  Alignment:\n\n"
    pos = 0
    while pos < len(align1):
        string2add = "s1  " + str(pos+1) + "   "
        report += string2add
        for i in range(pos, min(pos+60, len(align1))):
            report += align1[i]
        report += "  " + str(min(pos+60, len(align1))) + "\n"
        report += " "*len(string2add)
        for i in range(pos, min(pos+60, len(align1))):
            if align1[i] == align2[i]:
                report += "|"
                match_count += 1
            elif align1[i] == "-" or align2[i] == "-":
                report += " "
                if (i > 0) and (align1[i-1] != "-" and align2[i-1] != "-"):
                    gap_open_count += 1
                else:
                    gap_ext_count += 1
            else:
                report += "."
                mismatch_count += 1
        report += "\n"
        report += "s2  " + str(pos+1) + "   "
        for i in range(pos, min(pos+60, len(align2))):
            report += align2[i]
        report += "  " + str(min(pos+60, len(align2))) + "\n\n"
        pos += 60

    report += "Local score = " + str(max_score) + "\n"
    report += "Number of: matches = " + str(match_count) + ", mismatches = " + str(mismatch_count) + ", opening gaps = " + str(gap_open_count) + ", gap extensions = " + str(gap_ext_count) + "\n\n"
    #report += "Scores:    match = " + str(match) + ", mismatch = " + str(mismatch) + ", h = " + str(gap) + ", g = " + str(gap_ext) + "\n"

    return report


#read the parameters from the param.config file 

with open(sys.argv[3], 'r') as f:
    file_contents = f.read()

# Extract the variables from the file contents
variables = {}
for line in file_contents.split('\n'):
    line = line.strip()
    if line:
        var_name, var_value = line.split('=')
        variables[var_name.strip()] = int(var_value.strip())

# Get the values of the match, mismatch, h, and g variables
#if the value isn't provided the default value is taken for all 
match = variables.get('match', 1)
mismatch = variables.get('mismatch', -2)
h = variables.get('h', -5)
g = variables.get('g', -2)

# Use the variables as needed
print(f'match: {match}')
print(f'mismatch: {mismatch}')
print(f'gap opening penalty: {h}')
print(f'gap extension penalty: {g}')

#check alignment flag and show report statistics 

if sys.argv[2] == '0':
    print(needleman_wunsch(seq1, seq2, match, mismatch, h, g))
elif sys.argv[2] == '1' :
    print(smith_waterman(seq1, seq2, match, mismatch, h, g))
else :
    print ("enter alignment flag : 0 for global alignment 1 for local alignment")


