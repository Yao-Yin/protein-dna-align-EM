inputPath = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\RepeatPeps.lib"
outputPath = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\LINE1RepeatPeps.lib"

wf = open(outputPath, "w")
targetPattern = "LINE/L1"
with open(inputPath, "r") as rf:
    curr_line = rf.readline()
    while curr_line:
        curr_next = rf.readline()
        if curr_line.find(targetPattern) != -1:
            wf.write(curr_line)
            wf.write(curr_next)
        curr_line = rf.readline()

wf.close()