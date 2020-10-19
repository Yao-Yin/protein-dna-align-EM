# test for get sequences

import requests, sys
 
server = "https://rest.ensembl.org"
ext = "/sequence/id/ENSP00000160740?type=protein"
 
r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
print("hello") 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
 
print(r.text)
