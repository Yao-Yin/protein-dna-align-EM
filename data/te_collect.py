annotation_file = r"D:\\pseudogenome\\hg38_dfam.hits"

# test for get sequences

import requests, os, csv, json
#chr1	DF0000293.4	L1ME1_3end	40.9	3.1e-08	9.5	200	419	-	11674	11484	11694	11480	248956422	26.27
# url = "https://dfam.org/api/families/DF0000293"
#proteinSeqRequest = requests.get(url)
#print(proteinSeqRequest.status_code)
with open(annotation_file, "r") as f:
    for i in range(20):
        print(f.readline())